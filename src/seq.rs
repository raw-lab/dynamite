//! Encoded nucleotide sequence and the low-level sequence statistics shared by
//! both gene-calling engines.
//!
//! Bases are packed into small integers so that k-mer indices, GC tests and
//! translation are all cheap and mutually consistent:
//!
//! ```text
//!   A = 0b00   C = 0b01   G = 0b10   T = 0b11
//! ```
//!
//! This is the same 2-bit convention Prodigal uses internally, which means a
//! 6-mer index computed here (`mer_ndx`) is bit-for-bit identical to Prodigal's
//! and the dicodon coding model is directly comparable. `is_gc` reduces to
//! "the two bits differ" (C/G) versus "the two bits match" (A/T).
//!
//! Ambiguous / non-ACGT characters are folded to a concrete base (so that
//! stop/start scanning never panics) but their positions are remembered in an
//! `n` bitmap so that optional run-masking and N-aware translation still work.

/// Numeric code for adenine.
pub const A: u8 = 0;
/// Numeric code for cytosine.
pub const C: u8 = 1;
/// Numeric code for guanine.
pub const G: u8 = 2;
/// Numeric code for thymine.
pub const T: u8 = 3;

/// Map an input byte to a 2-bit base code, returning `(code, is_ambiguous)`.
///
/// `U`/`u` is treated as `T`. IUPAC ambiguity codes are folded to their most
/// representative concrete base (`S`,`B`,`V` -> `G`, everything else -> `A`),
/// matching PHANOTATE's behaviour, and flagged as ambiguous.
#[inline]
pub fn encode_base(b: u8) -> (u8, bool) {
    match b {
        b'A' | b'a' => (A, false),
        b'C' | b'c' => (C, false),
        b'G' | b'g' => (G, false),
        b'T' | b't' | b'U' | b'u' => (T, false),
        b'S' | b's' | b'B' | b'b' | b'V' | b'v' => (G, true),
        _ => (A, true),
    }
}

/// Decode a 2-bit base code back to an upper-case ASCII nucleotide.
#[inline]
pub fn decode_base(code: u8) -> u8 {
    match code & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

/// An encoded contig: forward strand, its reverse complement, and the
/// ambiguity bitmap. Both strands are pre-computed because the engines scan
/// each strand independently many times.
#[derive(Clone)]
pub struct Seq {
    /// Forward-strand base codes (`0..=3`), 5' -> 3'.
    pub fwd: Vec<u8>,
    /// Reverse-complement base codes (`0..=3`), 5' -> 3' of the minus strand.
    pub rev: Vec<u8>,
    /// `true` where the original input base was ambiguous / non-ACGT (forward
    /// coordinates).
    pub n: Vec<bool>,
    /// GC fraction over unambiguous bases.
    pub gc: f64,
}

impl Seq {
    /// Build an encoded sequence from raw ASCII nucleotides.
    pub fn from_ascii(bytes: &[u8]) -> Seq {
        let len = bytes.len();
        let mut fwd = Vec::with_capacity(len);
        let mut n = Vec::with_capacity(len);
        let mut gc_count = 0u64;
        let mut acgt = 0u64;
        for &b in bytes {
            let (code, amb) = encode_base(b);
            if code == C || code == G {
                gc_count += 1;
            }
            if !amb {
                acgt += 1;
            }
            fwd.push(code);
            n.push(amb);
        }
        let rev = revcomp(&fwd);
        let gc = if acgt > 0 {
            gc_count as f64 / acgt as f64
        } else {
            0.5
        };
        Seq { fwd, rev, n, gc }
    }

    /// Length of the contig in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.fwd.len()
    }

    /// Whether the contig is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.fwd.is_empty()
    }

    /// Reconstruct the forward strand as uppercase ASCII (`A/C/G/T`), with
    /// ambiguous positions rendered as `N`. Used to feed sequence to vendored
    /// engines that parse raw ASCII (e.g. the FragGeneScan HMM).
    pub fn to_ascii_upper(&self) -> Vec<u8> {
        self.fwd
            .iter()
            .zip(self.n.iter())
            .map(|(&c, &amb)| {
                if amb {
                    b'N'
                } else {
                    decode_base(c).to_ascii_uppercase()
                }
            })
            .collect()
    }
}

/// Reverse-complement a slice of base codes.
///
/// Complementing in this encoding is simply `3 - code` (A<->T, C<->G); the
/// sequence is then reversed.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| 3 - (c & 0b11)).collect()
}

/// `true` if the base at `i` is G or C.
#[inline]
pub fn is_gc(seq: &[u8], i: usize) -> bool {
    let c = seq[i];
    c == C || c == G
}

/// GC fraction over the inclusive range `[a, b]`.
pub fn gc_content(seq: &[u8], a: usize, b: usize) -> f64 {
    if b < a {
        return 0.0;
    }
    let mut gc = 0.0f64;
    let mut tot = 0.0f64;
    for i in a..=b {
        if is_gc(seq, i) {
            gc += 1.0;
        }
        tot += 1.0;
    }
    if tot == 0.0 {
        0.0
    } else {
        gc / tot
    }
}

/// Index of the `len`-mer starting at `pos`, base-4 little-endian.
///
/// Identical to Prodigal's `mer_ndx`: each base contributes `code << (2*i)`.
/// For `len == 6` the result is in `0..4096`. Positions that would run off the
/// end of the sequence are treated as `A` (code 0).
#[inline]
pub fn mer_ndx(len: usize, seq: &[u8], pos: usize) -> usize {
    let mut ndx = 0usize;
    let slen = seq.len();
    for i in 0..len {
        let p = pos + i;
        let code = if p < slen { seq[p] as usize } else { 0 };
        ndx |= code << (2 * i);
    }
    ndx
}

/// Render a `len`-mer index back to ASCII text (handy for motif reporting).
pub fn mer_text(len: usize, ndx: usize) -> String {
    if len == 0 {
        return "None".to_string();
    }
    let mut s = String::with_capacity(len);
    for i in 0..len {
        let code = ((ndx >> (2 * i)) & 0b11) as u8;
        s.push(decode_base(code) as char);
    }
    s
}

/// Compute the background frequency of every `len`-mer across both strands,
/// matching Prodigal's `calc_mer_bg`. The returned vector has `4^len` entries.
pub fn calc_mer_bg(len: usize, fwd: &[u8], rev: &[u8]) -> Vec<f64> {
    let size = 1usize << (2 * len);
    let mut counts = vec![0u64; size];
    let slen = fwd.len();
    if slen < len {
        return vec![0.0; size];
    }
    let mut glob = 0u64;
    for i in 0..=(slen - len) {
        counts[mer_ndx(len, fwd, i)] += 1;
        counts[mer_ndx(len, rev, i)] += 1;
        glob += 2;
    }
    let denom = glob as f64;
    counts
        .into_iter()
        .map(|c| if denom > 0.0 { c as f64 / denom } else { 0.0 })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_roundtrip() {
        let s = Seq::from_ascii(b"ATGCAANTAG");
        // revcomp of revcomp == original forward
        assert_eq!(revcomp(&s.rev), s.fwd);
    }

    #[test]
    fn gc_test() {
        // GCGC -> all GC
        let s = Seq::from_ascii(b"GCGC");
        assert!((s.gc - 1.0).abs() < 1e-9);
        for i in 0..4 {
            assert!(is_gc(&s.fwd, i));
        }
        // ATAT -> no GC
        let s2 = Seq::from_ascii(b"ATAT");
        assert!(s2.gc.abs() < 1e-9);
    }

    #[test]
    fn mer_index_matches_base4() {
        // ATG = A(0) T(3) G(2) -> 0<<0 | 3<<2 | 2<<4 = 12 + 32 = 44
        let s = Seq::from_ascii(b"ATG");
        assert_eq!(mer_ndx(3, &s.fwd, 0), 0 | (3 << 2) | (2 << 4));
        assert_eq!(&mer_text(3, mer_ndx(3, &s.fwd, 0)), "ATG");
    }
}

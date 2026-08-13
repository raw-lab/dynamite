//! A small, engine-independent open-reading-frame enumerator.
//!
//! Several callers (GeneMark, FragGeneScan, AUGUSTUS, the six-frame translator)
//! need the same primitive: walk all six frames of a contig and emit every
//! maximal stop-to-stop region, together with the best (most upstream) start
//! codon inside it. This module provides that, parameterised by genetic code,
//! minimum length, and whether reads running off the contig ends are allowed
//! (partial genes).

use crate::gencode::GeneticCode;
use crate::seq::Seq;

/// A candidate ORF in 1-based inclusive coordinates (stop codon included for
/// complete genes). `strand` is +1 / -1.
#[derive(Clone, Debug)]
pub struct OrfCand {
    pub begin: usize,
    pub end: usize,
    pub strand: i8,
    /// 0=ATG, 1=GTG, 2=TTG, 3=edge/none.
    pub start_type: u8,
    pub partial_left: bool,
    pub partial_right: bool,
}

/// Enumerate ORFs over both strands.
///
/// * `min_len` — minimum ORF length in nucleotides (stop included).
/// * `require_start` — if true, only ORFs that begin with a start codon are
///   reported (no 5'-partial edge ORFs); if false, edge/partial ORFs are kept.
pub fn find_orfs(
    seq: &Seq,
    code: &GeneticCode,
    min_len: usize,
    require_start: bool,
    allow_partial_ends: bool,
) -> Vec<OrfCand> {
    let mut out = Vec::new();
    let n = seq.len();
    if n < 3 {
        return out;
    }
    // Forward strand uses seq.fwd (2-bit codes); reverse uses seq.rev.
    scan_strand(&seq.fwd, n, 1, code, min_len, require_start, allow_partial_ends, &mut out);
    scan_strand(&seq.rev, n, -1, code, min_len, require_start, allow_partial_ends, &mut out);
    out
}

#[allow(clippy::too_many_arguments)]
fn scan_strand(
    s: &[u8],
    n: usize,
    strand: i8,
    code: &GeneticCode,
    min_len: usize,
    require_start: bool,
    allow_partial_ends: bool,
    out: &mut Vec<OrfCand>,
) {
    // For each of the three frames, track the position just after the previous
    // stop and the first start seen since then.
    for frame in 0..3usize {
        let mut region_start: Option<usize> = if allow_partial_ends { Some(frame) } else { None };
        let mut first_start: Option<usize> = None;
        let mut i = frame;
        while i + 3 <= n {
            let (b0, b1, b2) = (s[i], s[i + 1], s[i + 2]);
            let is_start = code.is_start(b0, b1, b2);
            let is_stop = code.is_stop(b0, b1, b2);
            if is_start && first_start.is_none() {
                first_start = Some(i);
                if region_start.is_none() {
                    region_start = Some(i);
                }
            }
            if is_stop {
                // Close the region [region_start .. i+2].
                if let Some(rs) = region_start {
                    let begin0 = if require_start {
                        first_start
                    } else {
                        first_start.or(Some(rs))
                    };
                    if let Some(b0pos) = begin0 {
                        let stop_end = i + 2; // 0-based inclusive end of stop codon
                        let len = stop_end + 1 - b0pos;
                        if len >= min_len {
                            push_orf(s, n, strand, b0pos, stop_end, first_start.is_none(), false, code, out);
                        }
                    }
                }
                region_start = None;
                first_start = None;
                // Next region begins after this stop.
                if allow_partial_ends {
                    region_start = Some(i + 3);
                }
            }
            i += 3;
        }
        // Trailing region with no terminal stop (3'-partial), if allowed.
        if allow_partial_ends {
            if let Some(rs) = region_start {
                let begin0 = first_start.or(Some(rs));
                if let Some(b0pos) = begin0 {
                    // Extend to the last full codon in this frame.
                    let last = frame + ((n - frame) / 3) * 3;
                    if last >= b0pos + 3 {
                        let stop_end = last - 1;
                        let len = stop_end + 1 - b0pos;
                        if len >= min_len {
                            push_orf(s, n, strand, b0pos, stop_end, first_start.is_none(), true, code, out);
                        }
                    }
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn push_orf(
    s: &[u8],
    n: usize,
    strand: i8,
    begin0: usize,
    end0: usize,
    partial_left: bool,
    partial_right: bool,
    code: &GeneticCode,
    out: &mut Vec<OrfCand>,
) {
    // Determine start type from the first codon (on the strand's coordinates).
    let st = if partial_left {
        3
    } else {
        GeneticCode::start_type(s[begin0], s[begin0 + 1], s[begin0 + 2])
    };
    let _ = code;
    // Convert strand-local 0-based [begin0,end0] to global 1-based [begin,end].
    // Partial flags stay STRAND-relative (partial_left = 5' end, partial_right =
    // 3' end), matching `Gene::fill_sequences`, which reverse-complements the
    // coding strand before applying the start-codon→M and stop-strip rules.
    let (begin, end, pl, pr) = if strand == 1 {
        (begin0 + 1, end0 + 1, partial_left, partial_right)
    } else {
        // reverse strand: position p on rev corresponds to global n-1-p.
        let g_begin = n - 1 - end0; // smaller global coordinate
        let g_end = n - 1 - begin0;
        (g_begin + 1, g_end + 1, partial_left, partial_right)
    };
    out.push(OrfCand {
        begin,
        end,
        strand,
        start_type: st,
        partial_left: pl,
        partial_right: pr,
    });
}

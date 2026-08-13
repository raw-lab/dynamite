//! **FragGeneScan** engine — a thin adapter over the vendored, verbatim
//! [`crate::fgs`] core of FragGeneScanRs (Van der Jeugt et al. 2022, the
//! maintained Rust port of FragGeneScan; Rho et al. 2010).
//!
//! Unlike the other engines this is **not** a from-scratch reimplementation: it
//! drives the real FragGeneScan forward/Viterbi HMM with its original trained
//! models, then maps each predicted gene into DYNAMITE's [`Gene`] using
//! FragGeneScan's own frameshift-aware nucleotide/protein extraction. Output is
//! therefore faithful to FragGeneScanRs, including insertion/deletion handling.
//!
//! `whole_genome = false` runs the short-read / fragment decode (genes may run
//! off either end, sequencing-error frameshifts modelled) — used for `reads`
//! and `long-reads` modes. `whole_genome = true` runs the complete-sequence
//! decode — used when FragGeneScan is forced on assembled contigs.

use crate::fgs::dna::Nuc;
use crate::fgs::hmm;
use crate::fgs::viterbi::viterbi;
use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::seq::Seq;
use std::path::PathBuf;

/// Pull the sequence line out of FragGeneScan's `>header\nSEQUENCE\n` buffer.
fn seq_line(buf: &[u8]) -> String {
    let s = String::from_utf8_lossy(buf);
    s.lines().nth(1).unwrap_or("").to_string()
}

fn gc_of(nt: &str) -> f64 {
    if nt.is_empty() {
        return 0.0;
    }
    let gc = nt
        .bytes()
        .filter(|b| matches!(b, b'g' | b'G' | b'c' | b'C'))
        .count();
    gc as f64 / nt.len() as f64
}

/// Run the FragGeneScan HMM over every contig/read. `whole_genome` selects the
/// complete-sequence (`true`) vs short-read/fragment (`false`) decode.
pub fn call(contigs: &[(String, Seq)], whole_genome: bool) -> (Vec<Gene>, u32) {
    let tt = 11u32; // FragGeneScan's models assume the standard/bacterial code
    let code = GeneticCode::new(tt);
    // A path that cannot exist forces the embedded (`include_bytes!`) models.
    let (global, locals) = match hmm::get_train_from_file(PathBuf::from("/nonexistent/dynamite-fgs-models"), PathBuf::from("complete")) {
        Ok(m) => m,
        Err(_) => return (Vec::new(), tt),
    };
    let dummy = b"s".to_vec();

    let mut genes = Vec::new();
    for (name, s) in contigs {
        if s.is_empty() {
            continue;
        }
        let ascii = s.to_ascii_upper();
        let nseq: Vec<Nuc> = ascii.iter().map(|&b| Nuc::from(b)).collect();
        let head: Vec<u8> = name.bytes().collect();
        let pred = viterbi(&global, &locals, head, nseq, whole_genome);

        for (i, fg) in pred.genes.iter().enumerate() {
            // FragGeneScan's own (insertion/deletion-aware) sequence extraction.
            let mut ntbuf = Vec::new();
            let _ = fg.dna(&mut ntbuf, &dummy, false);
            let mut aabuf = Vec::new();
            let _ = fg.protein(&mut aabuf, &dummy, whole_genome);
            let nt = seq_line(&ntbuf).to_ascii_lowercase();
            let aa = seq_line(&aabuf);
            if aa.is_empty() {
                continue;
            }

            let strand = if fg.forward_strand { 1i8 } else { -1 };
            // Classify the start codon from the coding-strand 5' triplet.
            let start_type = code
                .as_ref()
                .map(|c| {
                    let b = nt.as_bytes();
                    if b.len() >= 3 {
                        let enc = |x: u8| crate::seq::encode_base(x).0;
                        GeneticCode::start_type(enc(b[0]), enc(b[1]), enc(b[2]))
                    } else {
                        let _ = c;
                        3
                    }
                })
                .unwrap_or(3);

            genes.push(Gene {
                contig: name.clone(),
                id: i + 1,
                begin: fg.start,
                end: fg.end,
                strand,
                start_type,
                rbs_motif: "None".into(),
                rbs_spacer: "None".into(),
                partial_left: false,
                partial_right: false,
                score: fg.score,
                cscore: fg.score,
                sscore: 0.0,
                conf: 100.0 / (1.0 + (-fg.score / 20.0).exp()),
                gc: gc_of(&nt),
                trans_table: tt,
                engine: Engine::FragGeneScan,
                aa,
                nt,
            });
        }
    }
    (genes, tt)
}

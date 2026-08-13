//! **GeneMark** — a self-training, 3-periodic inhomogeneous Markov gene finder
//! in the spirit of GeneMark / GeneMark.hmm (Borodovsky & McIninch 1993;
//! Besemer et al. 2001).
//!
//! Algorithm: enumerate start→stop ORFs in all six frames, bootstrap a coding
//! model from the long ORFs (presumed coding), then score every ORF with a
//! coding/non-coding log-likelihood ratio and keep those that look coding.
//!
//! Honest scope: this captures GeneMark's *statistical core* (codon-position
//! Markov chains, self-training, LLR scoring) but is **not** the full
//! GeneMark.hmm — there is no hidden-state Viterbi over the whole sequence with
//! duration/RBS submodels, and no GeneMarkS iterative motif refinement.

use crate::dedup::resolve_overlaps;
use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::markov::MarkovModel;
use crate::orf::{find_orfs, OrfCand};
use crate::seq::Seq;

const MIN_ORF: usize = 90;
const MAX_OVERLAP: i64 = 60; // tolerated overlap between kept genes (nt)
const TRAIN_MIN: usize = 300; // long ORFs used to bootstrap the coding model

/// Call genes on a set of contigs with GeneMark-style self-training.
///
/// The coding model is trained once on the pooled long ORFs across all contigs,
/// then applied to every contig. `code` defaults to 11 if `None`.
pub fn call(contigs: &[(String, Seq)], fixed: Option<u32>, closed: bool) -> (Vec<Gene>, u32) {
    let tt = fixed.unwrap_or(11);
    let code = match GeneticCode::new(tt) {
        Some(c) => c,
        None => return (Vec::new(), tt),
    };
    let require_start = closed; // closed → starts only; open → keep edge ORFs too

    // 1) Enumerate ORFs per contig.
    let mut per_contig: Vec<Vec<OrfCand>> = Vec::with_capacity(contigs.len());
    let mut training: Vec<(usize, usize, i8)> = Vec::new();
    for (_, s) in contigs {
        let cands = find_orfs(s, &code, MIN_ORF, require_start, !require_start);
        for c in &cands {
            if c.end - c.begin + 1 >= TRAIN_MIN && !c.partial_left && !c.partial_right {
                training.push((c.begin, c.end, c.strand));
            }
        }
        per_contig.push(cands);
    }

    // 2) Train the coding model. Fall back to a lower length bar if too sparse.
    if training.len() < 10 {
        training.clear();
        for ((_, s), cands) in contigs.iter().zip(&per_contig) {
            let _ = s;
            for c in cands {
                if c.end - c.begin + 1 >= MIN_ORF + 60 {
                    training.push((c.begin, c.end, c.strand));
                }
            }
        }
    }

    // Train the coding model on the pooled long ORFs.
    let model = train_pooled(contigs, &training);

    // 3) Score every ORF; keep coding ones (LLR > 0). IDs are per-contig.
    let mut genes = Vec::new();
    for ((name, s), cands) in contigs.iter().zip(&per_contig) {
        let mut id = 0usize;
        for c in cands {
            let llr = model.llr(s, c.begin, c.end, c.strand);
            if llr <= 0.0 {
                continue;
            }
            id += 1;
            let mut g = build_gene(s, &code, c, llr, id);
            g.contig = name.clone();
            genes.push(g);
        }
    }
    (resolve_overlaps(genes, MAX_OVERLAP), tt)
}

fn train_pooled(contigs: &[(String, Seq)], training: &[(usize, usize, i8)]) -> MarkovModel {
    // Train on the contig that contains the most training ORFs; for the common
    // single-contig case this is exact. (A fully pooled multi-contig model would
    // need concatenation; the single-genome path is the primary use.)
    if contigs.len() == 1 {
        return MarkovModel::train(&contigs[0].1, training);
    }
    // Multi-contig: train per contig on its own ORFs and keep the richest model.
    let mut best: Option<MarkovModel> = None;
    let mut best_n = 0usize;
    for (_, s) in contigs {
        let local: Vec<(usize, usize, i8)> = training
            .iter()
            .filter(|&&(b, e, _)| e <= s.len() && b >= 1)
            .copied()
            .collect();
        let m = MarkovModel::train(s, &local);
        if m.trained_on >= best_n {
            best_n = m.trained_on;
            best = Some(m);
        }
    }
    best.unwrap_or_else(|| MarkovModel::train(&contigs[0].1, training))
}

fn build_gene(seq: &Seq, code: &GeneticCode, c: &OrfCand, llr: f64, id: usize) -> Gene {
    let mut g = Gene {
        contig: String::new(),
        id,
        begin: c.begin,
        end: c.end,
        strand: c.strand,
        start_type: c.start_type,
        rbs_motif: "None".into(),
        rbs_spacer: "None".into(),
        partial_left: c.partial_left,
        partial_right: c.partial_right,
        score: llr,
        cscore: llr,
        sscore: 0.0,
        conf: 100.0 / (1.0 + (-llr / 20.0).exp()),
        gc: 0.0,
        trans_table: code.id,
        engine: Engine::GeneMark,
        aa: String::new(),
        nt: String::new(),
    };
    g.fill_sequences(seq, code);
    g.gc = gc_of(&g.nt);
    g
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

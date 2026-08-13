//! Greedy overlap resolution for the ORF-enumeration engines.
//!
//! GeneMark, AUGUSTUS and FragGeneScan enumerate every positive-scoring open
//! reading frame in all six frames. Real implementations of those tools run a
//! hidden-state Viterbi decode that yields **one** consistent gene parse and
//! resolves overlaps; DYNAMITE approximates that with a greedy maximum-score
//! selection: take ORFs in descending coding score and keep one only if it does
//! not overlap an already-kept gene (on the same contig) by more than
//! `max_overlap` nucleotides. This collapses nested / duplicate / wrong-strand
//! ORFs while still permitting the short operon-style overlaps seen in compact
//! genomes. It is a heuristic, not a true Viterbi.

use crate::gene::Gene;
use std::cmp::Ordering;
use std::collections::HashMap;

#[inline]
fn overlap(a: &Gene, b: &Gene) -> i64 {
    let lo = a.begin.max(b.begin) as i64;
    let hi = (a.end.min(b.end)) as i64;
    (hi - lo + 1).max(0)
}

/// Resolve overlaps greedily by descending coding score, keeping genomic order
/// and renumbering gene IDs per contig.
pub fn resolve_overlaps(mut genes: Vec<Gene>, max_overlap: i64) -> Vec<Gene> {
    if genes.len() < 2 {
        return genes;
    }
    // Highest coding score first; break ties by longer gene.
    genes.sort_by(|a, b| {
        b.cscore
            .partial_cmp(&a.cscore)
            .unwrap_or(Ordering::Equal)
            .then((b.end - b.begin).cmp(&(a.end - a.begin)))
    });

    let mut kept: Vec<Gene> = Vec::with_capacity(genes.len());
    for g in genes {
        let clashes = kept
            .iter()
            .any(|k| k.contig == g.contig && overlap(k, &g) > max_overlap);
        if !clashes {
            kept.push(g);
        }
    }

    // Restore genomic order and renumber per contig.
    kept.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.begin.cmp(&b.begin)));
    let mut counts: HashMap<String, usize> = HashMap::new();
    for g in &mut kept {
        let c = counts.entry(g.contig.clone()).or_insert(0);
        *c += 1;
        g.id = *c;
    }
    kept
}

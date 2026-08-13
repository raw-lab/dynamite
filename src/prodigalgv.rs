//! **Prodigal-GV** preset — phages, giant viruses (NCLDV / *Nucleocytoviricota*)
//! and crassphages.
//!
//! Wraps the shared Prodigal dynamic-programming core but auto-detects the
//! genetic code among the tables seen across large/▷giant viruses and
//! crassphages (which frequently reassign stop codons — e.g. code 15 reading
//! `TAG`→Q in some crassphages). Detection uses DYNAMITE's standard-code bias
//! so ordinary code-11 genomes are not spuriously reassigned.

use crate::gene::{Engine, Gene};
use crate::prodigal::dp_call;
use crate::seq::Seq;

/// Candidate genetic codes for giant viruses and crassphages.
pub const GV_CANDIDATES: &[u32] = &[11, 1, 4, 15, 25, 6];

/// Run the Prodigal-GV preset. `fixed = Some(table)` forces a code; `None`
/// auto-detects among [`GV_CANDIDATES`]. Treated as metagenome-like so multiple
/// contigs with mixed coding are handled gracefully.
pub fn call(contigs: &[(String, Seq)], fixed: Option<u32>, closed: bool) -> (Vec<Gene>, u32) {
    dp_call(contigs, fixed, GV_CANDIDATES, true, closed, Engine::ProdigalGv)
}

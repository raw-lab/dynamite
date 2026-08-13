//! # DYNAMITE core
//!
//! `dynamite_core` is the engine library behind the DYNAMITE gene caller — an
//! ORF / gene caller for **all domains of life**. It exposes a family of
//! callers that share a common sequence model, genetic-code table, and
//! [`gene::Gene`] output type:
//!
//! * [`phanotate`] — PHANOTATE-style graph / shortest-path caller (phages).
//! * [`prodigal`] — Prodigal dynamic-programming engine; the `prodigal`,
//!   [`prodigalgv`], and [`prodigal_meta`] presets configure it for bacteria /
//!   archaea, phages + giant viruses / crassphages, and metagenomes.
//! * [`genemark`] — self-training 3-periodic Markov caller (GeneMark-style).
//! * [`augustus`] — eukaryotic single-exon caller (AUGUSTUS-style).
//! * [`fraggenescan`] — short-read caller with frameshift handling.
//!
//! Plus helpers: [`sixframe`] (six-frame translation), [`codon`] (codon-bias /
//! density statistics), [`rna`] (heuristic tRNA / rRNA detection), and
//! [`doctor`] (self-test). [`run_records`] ties the callers together, resolving
//! the engine and genetic code from a [`model::Mode`] (or an explicit engine
//! override) and predicting every contig into one gene list.

pub mod augustus;
pub mod banner;
pub mod codon;
pub mod dedup;
pub mod doctor;
pub mod fasta;
pub mod annotate;
pub mod fgs;
pub mod fraggenescan;
pub mod gcframe;
pub mod gencode;
pub mod gene;
pub mod genemark;
pub mod markov;
pub mod model;
pub mod node;
pub mod orf;
pub mod phanotate;
pub mod prodigal;
pub mod prodigal_meta;
pub mod prodigalgv;
pub mod rbs;
pub mod rna;
pub mod seq;
pub mod sixframe;

use anyhow::{anyhow, Result};
use rayon::prelude::*;

use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::model::{config_for, CodeChoice, EngineKind, Mode};
use crate::seq::Seq;

/// Everything needed to run a prediction over a set of records.
#[derive(Clone, Debug)]
pub struct RunConfig {
    /// Organism preset.
    pub mode: Mode,
    /// Explicit engine override (takes precedence over the mode's default).
    pub engine: Option<EngineKind>,
    /// Genetic-code choice (overrides the mode's default when not `Default`).
    pub code: CodeChoice,
    /// Minimum ORF length in nucleotides (PHANOTATE engine). `<= 0` uses the
    /// engine default.
    pub min_len: i64,
    /// Disallow genes that run off the ends of a contig.
    pub closed: bool,
    /// Worker threads (`0` = all available).
    pub threads: usize,
}

impl Default for RunConfig {
    fn default() -> Self {
        RunConfig {
            mode: Mode::Auto,
            engine: None,
            code: CodeChoice::Default,
            min_len: 0,
            closed: false,
            threads: 0,
        }
    }
}

/// The result of a run: every gene called, plus the `(id, length)` of every
/// input contig (needed for GFF/GenBank headers).
pub struct RunOutput {
    pub genes: Vec<Gene>,
    pub contigs: Vec<(String, usize)>,
    /// Which engine ran.
    pub engine: Engine,
    /// The genetic code that was used (resolved, after any auto-detection).
    pub trans_table: u32,
}

fn resolve_fixed(code: &CodeChoice, mode_default: Option<u32>) -> Option<u32> {
    match code {
        CodeChoice::Fixed(id) => Some(*id),
        CodeChoice::Default => mode_default,
        CodeChoice::Auto => None,
    }
}

fn kind_to_engine(k: EngineKind) -> Engine {
    match k {
        EngineKind::Phanotate => Engine::Phanotate,
        EngineKind::Prodigal => Engine::Prodigal,
        EngineKind::ProdigalGv => Engine::ProdigalGv,
        EngineKind::ProdigalMeta => Engine::ProdigalMeta,
        EngineKind::GeneMark => Engine::GeneMark,
        EngineKind::Augustus => Engine::Augustus,
        EngineKind::FragGeneScan => Engine::FragGeneScan,
    }
}

/// Run a prediction over the given FASTA records.
pub fn run_records(records: &[fasta::Record], cfg: &RunConfig) -> Result<RunOutput> {
    let mode_cfg = config_for(cfg.mode);
    let engine = cfg.engine.unwrap_or(mode_cfg.engine);
    let closed = cfg.closed || mode_cfg.closed;
    let min_len = if cfg.min_len <= 0 {
        phanotate::MIN_ORF_LEN
    } else {
        cfg.min_len
    };
    let fixed = resolve_fixed(&cfg.code, mode_cfg.code);

    let contig_list: Vec<(String, usize)> =
        records.iter().map(|r| (r.id.clone(), r.seq.len())).collect();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cfg.threads)
        .build()
        .map_err(|e| anyhow!("failed to build thread pool: {e}"))?;

    // PHANOTATE works directly on ASCII per record; all other engines take
    // 2-bit-encoded contigs.
    if engine == EngineKind::Phanotate {
        let code_id = fixed.unwrap_or(11);
        let code = GeneticCode::new(code_id)
            .ok_or_else(|| anyhow!("unsupported genetic code {code_id}"))?;
        let genes: Vec<Gene> = pool.install(|| {
            records
                .par_iter()
                .flat_map_iter(|r| {
                    let (g, _score) = phanotate::predict(&r.id, &r.seq, &code, min_len, closed);
                    g.into_iter()
                })
                .collect()
        });
        return Ok(RunOutput {
            genes,
            contigs: contig_list,
            engine: Engine::Phanotate,
            trans_table: code_id,
        });
    }

    let contigs: Vec<(String, Seq)> = records
        .iter()
        .map(|r| (r.id.clone(), Seq::from_ascii(&r.seq)))
        .collect();

    let (genes, table) = pool.install(|| match engine {
        EngineKind::Prodigal => prodigal::call(&contigs, fixed, closed),
        EngineKind::ProdigalGv => prodigalgv::call(&contigs, fixed, closed),
        EngineKind::ProdigalMeta => prodigal_meta::call(&contigs, fixed, closed),
        EngineKind::GeneMark => genemark::call(&contigs, fixed, closed),
        EngineKind::Augustus => augustus::call(&contigs, fixed, closed),
        EngineKind::FragGeneScan => {
            let whole = !matches!(cfg.mode, crate::model::Mode::Reads | crate::model::Mode::LongReads);
            fraggenescan::call(&contigs, whole)
        }
        EngineKind::Phanotate => unreachable!(),
    });

    Ok(RunOutput {
        genes,
        contigs: contig_list,
        engine: kind_to_engine(engine),
        trans_table: table,
    })
}

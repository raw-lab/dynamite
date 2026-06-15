//! **Annotate** — a one-shot, Prokka-style annotation bundle.
//!
//! Runs gene calling (any DYNAMITE engine) plus heuristic tRNA/rRNA detection,
//! assigns sequential `locus_tag`s along each contig, and renders the full set
//! of output files Prokka emits — `.gff .gbk .fna .faa .ffn .tsv .txt .tbl .fsa
//! .err .log` — from a single pass.
//!
//! Two files in Prokka's set are produced by NCBI's `tbl2asn` / `table2asn`,
//! which is a separate binary DYNAMITE does not bundle: the `.sqn` (ASN.1
//! Sequin) and the official `.err` discrepancy report. DYNAMITE writes the
//! `.tbl` + `.fsa` that `table2asn` consumes, and a *DYNAMITE* QA report as
//! `.err` (clearly labelled — it is not the NCBI discrepancy report). If
//! `table2asn` is on `PATH`, the CLI runs it to produce the real `.sqn`/`.err`.

use crate::fasta::Record;
use crate::gencode::GeneticCode;
use crate::seq::Seq;
use crate::{rna, run_records, RunConfig};
use anyhow::Result;
use std::fmt::Write as _;

/// Settings for an annotation run.
pub struct AnnotateConfig {
    /// Gene-calling configuration (engine / mode / code / threads).
    pub run: RunConfig,
    /// `locus_tag` prefix, e.g. `DYNAMITE` → `DYNAMITE_00001`.
    pub locus_prefix: String,
    /// Organism name for GenBank / Sequin metadata.
    pub organism: String,
    /// Strain for GenBank / Sequin metadata.
    pub strain: String,
    /// Whether to scan for tRNA / rRNA features.
    pub detect_rna: bool,
    /// Genetic code used for tRNA anticodon → amino-acid mapping.
    pub rna_code: u32,
}

/// Every output file, rendered to a string.
pub struct AnnotationBundle {
    pub gff: String,
    pub gbk: String,
    pub fna: String,
    pub faa: String,
    pub ffn: String,
    pub tsv: String,
    pub txt: String,
    pub tbl: String,
    pub fsa: String,
    pub err: String,
    pub log: String,
    /// Feature counts for the run summary.
    pub n_cds: usize,
    pub n_trna: usize,
    pub n_rrna: usize,
}

/// Unified annotation feature (CDS or RNA).
struct Feat {
    contig: String,
    begin: usize,
    end: usize,
    strand: i8,
    ftype: String, // "CDS" | "tRNA" | "rRNA"
    locus_tag: String,
    product: String,
    nt: String,
    aa: String, // CDS only
    partial_left: bool,
    partial_right: bool,
    trans_table: u32,
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b't',
        b'C' | b'c' => b'g',
        b'G' | b'g' => b'c',
        b'T' | b't' => b'a',
        _ => b'n',
    }
}

/// Lowercase reverse-complement of an ASCII slice.
fn revcomp_ascii(s: &[u8]) -> String {
    s.iter().rev().map(|&b| complement(b) as char).collect()
}

/// 1-based inclusive subsequence on the given strand, lowercase.
fn subseq(ascii: &[u8], begin: usize, end: usize, strand: i8) -> String {
    if begin == 0 || end < begin || end > ascii.len() {
        return String::new();
    }
    let slice = &ascii[begin - 1..end];
    if strand >= 0 {
        slice.iter().map(|&b| b.to_ascii_lowercase() as char).collect()
    } else {
        revcomp_ascii(slice)
    }
}

fn wrap(out: &mut String, seq: &str, width: usize) {
    let b = seq.as_bytes();
    let mut i = 0;
    while i < b.len() {
        let e = (i + width).min(b.len());
        out.push_str(std::str::from_utf8(&b[i..e]).unwrap_or(""));
        out.push('\n');
        i = e;
    }
    if b.is_empty() {
        out.push('\n');
    }
}

/// Map an `rna::RnaFeature` to a (ftype, product) pair.
fn rna_product(kind: &str, note: &str) -> (String, String) {
    if kind.starts_with("tRNA") {
        let p = if note.is_empty() {
            "tRNA".to_string()
        } else {
            note.to_string()
        };
        ("tRNA".to_string(), p)
    } else if kind.starts_with("rRNA") {
        let p = if kind.contains("16S") {
            "16S ribosomal RNA"
        } else if kind.contains("23S") {
            "23S ribosomal RNA"
        } else if kind.contains("5S") {
            "5S ribosomal RNA"
        } else {
            "ribosomal RNA"
        };
        ("rRNA".to_string(), p.to_string())
    } else {
        ("misc_RNA".to_string(), kind.to_string())
    }
}

/// Build the full annotation bundle from input records.
pub fn annotate(records: &[Record], cfg: &AnnotateConfig) -> Result<AnnotationBundle> {
    // 1) CDS calling.
    let out = run_records(records, &cfg.run)?;
    let engine_label = out.engine.label().to_string();
    let cds_table = out.genes.first().map(|g| g.trans_table).unwrap_or(11);

    // contig ascii (uppercased, original order) keyed for sequence emission.
    let contigs: Vec<(String, String, Vec<u8>)> = records
        .iter()
        .map(|r| {
            let up: Vec<u8> = r.seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            (r.id.clone(), r.desc.clone(), up)
        })
        .collect();

    // 2) RNA detection (optional).
    let mut rna_by_contig: Vec<(String, Vec<rna::RnaFeature>)> = Vec::new();
    if cfg.detect_rna {
        let code = GeneticCode::new(cfg.rna_code).unwrap_or_else(|| GeneticCode::new(11).unwrap());
        for (id, _desc, up) in &contigs {
            let s = Seq::from_ascii(up);
            let feats = rna::find_all(&s, &code);
            rna_by_contig.push((id.clone(), feats));
        }
    }

    // 3) Build unified, per-contig feature lists; assign locus_tags in order.
    let mut feats: Vec<Feat> = Vec::new();
    let mut counter: usize = 0;
    let mut n_cds = 0usize;
    let mut n_trna = 0usize;
    let mut n_rrna = 0usize;

    for (cid, _desc, up) in &contigs {
        let mut local: Vec<Feat> = Vec::new();
        // CDS on this contig
        for g in out.genes.iter().filter(|g| &g.contig == cid) {
            local.push(Feat {
                contig: cid.clone(),
                begin: g.begin,
                end: g.end,
                strand: g.strand,
                ftype: "CDS".into(),
                locus_tag: String::new(),
                product: "hypothetical protein".into(),
                nt: g.nt.clone(),
                aa: g.aa.clone(),
                partial_left: g.partial_left,
                partial_right: g.partial_right,
                trans_table: g.trans_table,
            });
        }
        // RNA on this contig
        if let Some((_, rfeats)) = rna_by_contig.iter().find(|(id, _)| id == cid) {
            for rf in rfeats {
                let (ftype, product) = rna_product(&rf.kind, &rf.note);
                local.push(Feat {
                    contig: cid.clone(),
                    begin: rf.begin,
                    end: rf.end,
                    strand: rf.strand,
                    ftype,
                    locus_tag: String::new(),
                    product,
                    nt: subseq(up, rf.begin, rf.end, rf.strand),
                    aa: String::new(),
                    partial_left: false,
                    partial_right: false,
                    trans_table: 0,
                });
            }
        }
        local.sort_by(|a, b| a.begin.cmp(&b.begin).then(a.end.cmp(&b.end)));
        for mut f in local {
            counter += 1;
            f.locus_tag = format!("{}_{:05}", cfg.locus_prefix, counter);
            match f.ftype.as_str() {
                "CDS" => n_cds += 1,
                "tRNA" => n_trna += 1,
                "rRNA" => n_rrna += 1,
                _ => {}
            }
            feats.push(f);
        }
    }

    // ---- renders ----
    let org_full = if cfg.strain.is_empty() {
        cfg.organism.clone()
    } else {
        format!("{} {}", cfg.organism, cfg.strain)
    };

    // .fna — input contigs
    let mut fna = String::new();
    for (id, desc, up) in &contigs {
        if desc.is_empty() {
            let _ = writeln!(fna, ">{id}");
        } else {
            let _ = writeln!(fna, ">{id} {desc}");
        }
        wrap(&mut fna, &String::from_utf8_lossy(up), 70);
    }

    // .fsa — input contigs with Sequin tags (table2asn input)
    let mut fsa = String::new();
    for (id, _desc, up) in &contigs {
        let _ = writeln!(
            fsa,
            ">{id} [organism={}] [strain={}] [gcode={}] [molecule=DNA] [tech=wgs]",
            cfg.organism, cfg.strain, cds_table
        );
        wrap(&mut fsa, &String::from_utf8_lossy(up), 70);
    }

    // .faa — CDS proteins
    let mut faa = String::new();
    for f in feats.iter().filter(|f| f.ftype == "CDS") {
        let _ = writeln!(faa, ">{} {}", f.locus_tag, f.product);
        wrap(&mut faa, &f.aa, 60);
    }

    // .ffn — all transcript nucleotides
    let mut ffn = String::new();
    for f in &feats {
        let _ = writeln!(ffn, ">{} {}", f.locus_tag, f.product);
        wrap(&mut ffn, &f.nt, 60);
    }

    // .gff — master GFF3 + ##FASTA
    let mut gff = String::new();
    gff.push_str("##gff-version 3\n");
    for (id, _d, up) in &contigs {
        let _ = writeln!(gff, "##sequence-region {id} 1 {}", up.len());
    }
    for f in &feats {
        let strand = if f.strand >= 0 { '+' } else { '-' };
        let phase = if f.ftype == "CDS" { "0" } else { "." };
        let _ = writeln!(
            gff,
            "{ctg}\tDYNAMITE\t{ft}\t{b}\t{e}\t.\t{strand}\t{phase}\tID={lt};locus_tag={lt};product={prod}",
            ctg = f.contig,
            ft = f.ftype,
            b = f.begin,
            e = f.end,
            lt = f.locus_tag,
            prod = f.product,
        );
    }
    gff.push_str("##FASTA\n");
    for (id, _d, up) in &contigs {
        let _ = writeln!(gff, ">{id}");
        wrap(&mut gff, &String::from_utf8_lossy(up), 70);
    }

    // .tsv — Prokka-style feature table
    let mut tsv = String::new();
    tsv.push_str("locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct\n");
    for f in &feats {
        let len = f.end.saturating_sub(f.begin) + 1;
        let _ = writeln!(
            tsv,
            "{}\t{}\t{}\t\t\t\t{}",
            f.locus_tag, f.ftype, len, f.product
        );
    }

    // .tbl — NCBI feature table (table2asn input)
    let mut tbl = String::new();
    for (id, _d, _up) in &contigs {
        let _ = writeln!(tbl, ">Feature {id}");
        for f in feats.iter().filter(|f| &f.contig == id) {
            // 5'..3' coordinates in feature orientation; < / > mark partials.
            let (c5, c3) = if f.strand >= 0 {
                (f.begin, f.end)
            } else {
                (f.end, f.begin)
            };
            let p5 = if f.partial_left { "<" } else { "" };
            let p3 = if f.partial_right { ">" } else { "" };
            let _ = writeln!(tbl, "{p5}{c5}\t{p3}{c3}\t{}", f.ftype);
            let _ = writeln!(tbl, "\t\t\tlocus_tag\t{}", f.locus_tag);
            let _ = writeln!(tbl, "\t\t\tproduct\t{}", f.product);
            if f.ftype == "CDS" {
                let _ = writeln!(tbl, "\t\t\tcodon_start\t1");
                let _ = writeln!(tbl, "\t\t\ttransl_table\t{}", f.trans_table);
            }
        }
    }

    // .gbk — GenBank with all features
    let mut gbk = String::new();
    for (id, _d, up) in &contigs {
        let len = up.len();
        let _ = writeln!(
            gbk,
            "LOCUS       {:<16}{} bp    DNA     linear       01-JAN-2025",
            id, len
        );
        let _ = writeln!(
            gbk,
            "DEFINITION  {}, whole genome shotgun sequence.",
            if org_full.is_empty() { "Unknown organism" } else { &org_full }
        );
        let _ = writeln!(gbk, "ACCESSION   {id}");
        let _ = writeln!(gbk, "SOURCE      {}", if cfg.organism.is_empty() { "Unknown" } else { &cfg.organism });
        let _ = writeln!(gbk, "  ORGANISM  {}", if cfg.organism.is_empty() { "Unknown" } else { &cfg.organism });
        let _ = writeln!(gbk, "FEATURES             Location/Qualifiers");
        let _ = writeln!(gbk, "     source          1..{len}");
        if !cfg.organism.is_empty() {
            let _ = writeln!(gbk, "                     /organism=\"{}\"", cfg.organism);
        }
        if !cfg.strain.is_empty() {
            let _ = writeln!(gbk, "                     /strain=\"{}\"", cfg.strain);
        }
        let _ = writeln!(gbk, "                     /mol_type=\"genomic DNA\"");
        for f in feats.iter().filter(|f| &f.contig == id) {
            let loc = if f.strand >= 0 {
                format!("{}..{}", f.begin, f.end)
            } else {
                format!("complement({}..{})", f.begin, f.end)
            };
            let _ = writeln!(gbk, "     {:<16}{}", f.ftype, loc);
            let _ = writeln!(gbk, "                     /locus_tag=\"{}\"", f.locus_tag);
            let _ = writeln!(gbk, "                     /product=\"{}\"", f.product);
            if f.ftype == "CDS" {
                let _ = writeln!(gbk, "                     /codon_start=1");
                let _ = writeln!(gbk, "                     /transl_table={}", f.trans_table);
                // wrapped translation
                let aa = &f.aa;
                let mut first = true;
                let b = aa.as_bytes();
                let mut i = 0;
                while i < b.len() {
                    let e = (i + 58).min(b.len());
                    let chunk = std::str::from_utf8(&b[i..e]).unwrap_or("");
                    if first {
                        let _ = writeln!(gbk, "                     /translation=\"{chunk}");
                        first = false;
                    } else {
                        let _ = writeln!(gbk, "                     {chunk}");
                    }
                    i = e;
                }
                if first {
                    let _ = writeln!(gbk, "                     /translation=\"\"");
                } else {
                    // close the quote on the last written line
                    let trimmed = gbk.trim_end_matches('\n').len();
                    gbk.truncate(trimmed);
                    gbk.push_str("\"\n");
                }
            }
        }
        // ORIGIN
        let _ = writeln!(gbk, "ORIGIN");
        let lower = String::from_utf8_lossy(up).to_ascii_lowercase();
        let bytes = lower.as_bytes();
        let mut pos = 0;
        while pos < bytes.len() {
            let line_end = (pos + 60).min(bytes.len());
            let mut line = String::new();
            let _ = write!(line, "{:>9}", pos + 1);
            let mut j = pos;
            while j < line_end {
                let g_end = (j + 10).min(line_end);
                line.push(' ');
                line.push_str(std::str::from_utf8(&bytes[j..g_end]).unwrap_or(""));
                j = g_end;
            }
            let _ = writeln!(gbk, "{line}");
            pos = line_end;
        }
        let _ = writeln!(gbk, "//");
    }

    // .txt — statistics
    let total_bases: usize = contigs.iter().map(|(_, _, up)| up.len()).sum();
    let mut txt = String::new();
    let _ = writeln!(txt, "organism: {}", if org_full.is_empty() { "Unknown organism" } else { &org_full });
    let _ = writeln!(txt, "engine: {}", engine_label);
    let _ = writeln!(txt, "contigs: {}", contigs.len());
    let _ = writeln!(txt, "bases: {}", total_bases);
    let _ = writeln!(txt, "CDS: {}", n_cds);
    let _ = writeln!(txt, "tRNA: {}", n_trna);
    let _ = writeln!(txt, "rRNA: {}", n_rrna);
    let _ = writeln!(txt, "tmRNA: 0");

    // .err — DYNAMITE QA report (NOT the NCBI discrepancy report)
    let mut err = String::new();
    err.push_str("# DYNAMITE QA report — internal sanity checks only.\n");
    err.push_str("# This is NOT the NCBI discrepancy report; run `table2asn` on the\n");
    err.push_str("# .fsa + .tbl for the official .err / .sqn.\n\n");
    let short: Vec<&Feat> = feats
        .iter()
        .filter(|f| f.ftype == "CDS" && f.aa.len() < 30)
        .collect();
    let partial: Vec<&Feat> = feats
        .iter()
        .filter(|f| f.partial_left || f.partial_right)
        .collect();
    let internal_stop: Vec<&Feat> = feats
        .iter()
        .filter(|f| f.ftype == "CDS" && f.aa.trim_end_matches('*').contains('*'))
        .collect();
    let _ = writeln!(err, "SHORT_CDS (<30 aa): {}", short.len());
    for f in &short {
        let _ = writeln!(err, "  {} {}:{}-{} {} aa", f.locus_tag, f.contig, f.begin, f.end, f.aa.len());
    }
    let _ = writeln!(err, "PARTIAL_FEATURES: {}", partial.len());
    for f in &partial {
        let _ = writeln!(
            err,
            "  {} {}:{}-{} partial={}{}",
            f.locus_tag,
            f.contig,
            f.begin,
            f.end,
            if f.partial_left { "5'" } else { "" },
            if f.partial_right { "3'" } else { "" }
        );
    }
    let _ = writeln!(err, "CDS_WITH_INTERNAL_STOP: {}", internal_stop.len());
    for f in &internal_stop {
        let _ = writeln!(err, "  {} {}:{}-{}", f.locus_tag, f.contig, f.begin, f.end);
    }

    // .log — run record (deterministic; no timestamps)
    let mut log = String::new();
    log.push_str("DYNAMITE annotate\n");
    let _ = writeln!(log, "engine: {}", engine_label);
    let _ = writeln!(log, "locus_tag prefix: {}", cfg.locus_prefix);
    let _ = writeln!(log, "organism: {}", org_full);
    let _ = writeln!(log, "detect_rna: {}", cfg.detect_rna);
    let _ = writeln!(log, "contigs: {}", contigs.len());
    let _ = writeln!(log, "features: CDS={} tRNA={} rRNA={}", n_cds, n_trna, n_rrna);
    log.push_str("outputs: gff gbk fna faa ffn tsv txt tbl fsa err log\n");

    Ok(AnnotationBundle {
        gff,
        gbk,
        fna,
        faa,
        ffn,
        tsv,
        txt,
        tbl,
        fsa,
        err,
        log,
        n_cds,
        n_trna,
        n_rrna,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::Record;
    use crate::model::{CodeChoice, Mode};

    #[test]
    fn bundle_smoke() {
        let seq = b"ATGGCACGTGACATTAAACTGACCGGTAAACGTGCCTTAGATTCAGCTTAA".to_vec();
        let recs = vec![Record {
            id: "c1".into(),
            desc: "test contig".into(),
            seq,
        }];
        let cfg = AnnotateConfig {
            run: RunConfig {
                mode: Mode::Reads, // vendored FGS HMM: no training needed
                engine: None,
                code: CodeChoice::Default,
                min_len: 0,
                closed: false,
                threads: 1,
            },
            locus_prefix: "T".into(),
            organism: "Testus exampli".into(),
            strain: "x1".into(),
            detect_rna: false,
            rna_code: 11,
        };
        let b = annotate(&recs, &cfg).unwrap();
        assert!(b.tsv.starts_with("locus_tag\tftype\tlength_bp"));
        assert!(b.fna.contains(">c1"));
        assert!(b.fsa.contains("[organism=Testus exampli]"));
        assert!(b.gff.contains("##gff-version 3"));
        assert!(b.gff.contains("##FASTA"));
        assert!(b.gbk.contains("LOCUS"));
        assert!(b.gbk.trim_end().ends_with("//"));
        assert!(b.tbl.contains(">Feature c1"));
        assert!(b.txt.contains("contigs: 1"));
    }
}

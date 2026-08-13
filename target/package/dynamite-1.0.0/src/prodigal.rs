//! The Prodigal-style dynamic-programming engine: the DP traceback
//! (`dprog` / `score_connection`), gene extraction (`add_genes`,
//! `tweak_final_starts`), and the high-level self-training / prediction
//! orchestration including DYNAMITE's genetic-code auto-detection.
//!
//! Faithful Rust port of `dprog.c` and parts of `gene.c` (Prodigal /
//! Prodigal-GV, Doug Hyatt et al., GPL-3.0).

use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::node::{
    add_nodes, calc_dicodon_gene, determine_sd_usage, gc_frame_plot, intergenic_mod,
    raw_coding_score, rbs_score, record_gc_bias, record_overlapping_starts, score_nodes,
    sort_nodes, train_starts_sd, Node, Tinf, STOP,
};
use crate::seq::Seq;
use rayon::prelude::*;

const MAX_OPP_OVLP: i64 = 200;
const MAX_NODE_DIST: isize = 500;

/// The training concatenation separator. `TTAATTAATTAA` forces stop codons in
/// all six reading frames, preventing ORFs from spanning contig junctions
/// (exactly as Prodigal's `read_seq_training` does).
pub const TRAIN_GAP: &[u8] = b"TTAATTAATTAA";

/// Confidence from a model score (`calculate_confidence`).
pub fn calculate_confidence(score: f64, start_weight: f64) -> f64 {
    let mut conf;
    if score / start_weight < 41.0 {
        let c = (score / start_weight).exp();
        conf = (c / (c + 1.0)) * 100.0;
    } else {
        conf = 99.99;
    }
    if conf <= 50.0 {
        conf = 50.0;
    }
    conf
}

/// The core dynamic-programming routine (`dprog`). Returns the index of the
/// maximum-scoring terminal node, or `-1`.
pub fn dprog(nodes: &mut [Node], tinf: &Tinf, flag: bool) -> i64 {
    let nn = nodes.len();
    if nn == 0 {
        return -1;
    }
    for n in nodes.iter_mut() {
        n.score = 0.0;
        n.traceb = -1;
        n.tracef = -1;
    }

    for i in 0..nn {
        let i_i = i as isize;
        let mut min: isize = if i_i < MAX_NODE_DIST {
            0
        } else {
            i_i - MAX_NODE_DIST
        };
        let guard = (nodes[i].strand == -1 && nodes[i].typ != STOP)
            || (nodes[i].strand == 1 && nodes[i].typ == STOP);
        if guard && nodes[min as usize].ndx >= nodes[i].stop_val {
            while min >= 0 && nodes[i].ndx != nodes[i].stop_val {
                min -= 1;
            }
        }
        if min < MAX_NODE_DIST {
            min = 0;
        } else {
            min -= MAX_NODE_DIST;
        }
        let mut j = min as usize;
        while j < i {
            score_connection(nodes, j, i, tinf, flag);
            j += 1;
        }
    }

    let mut max_sc = -1.0f64;
    let mut max_ndx: i64 = -1;
    for i in (0..nn).rev() {
        if nodes[i].strand == 1 && nodes[i].typ != STOP {
            continue;
        }
        if nodes[i].strand == -1 && nodes[i].typ == STOP {
            continue;
        }
        if nodes[i].score > max_sc {
            max_sc = nodes[i].score;
            max_ndx = i as i64;
        }
    }
    if max_ndx == -1 {
        return -1;
    }

    // First pass: untangle triple overlaps.
    let mut path = max_ndx;
    while nodes[path as usize].traceb != -1 {
        let nxt = nodes[path as usize].traceb;
        if nodes[path as usize].strand == -1
            && nodes[path as usize].typ == STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].typ == STOP
            && nodes[path as usize].ov_mark != -1
            && nodes[path as usize].ndx > nodes[nxt as usize].ndx
        {
            let ov = nodes[path as usize].ov_mark as usize;
            let tmp = nodes[path as usize].star_ptr[ov];
            let mut i = tmp;
            while nodes[i as usize].ndx != nodes[tmp as usize].stop_val {
                i -= 1;
            }
            nodes[path as usize].traceb = tmp;
            nodes[tmp as usize].traceb = i;
            nodes[i as usize].ov_mark = -1;
            nodes[i as usize].traceb = nxt;
        }
        path = nodes[path as usize].traceb;
    }

    // Second pass: untangle simple overlaps.
    let mut path = max_ndx;
    while nodes[path as usize].traceb != -1 {
        let nxt = nodes[path as usize].traceb;
        if nodes[path as usize].strand == -1
            && nodes[path as usize].typ != STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].typ == STOP
        {
            let mut i = path;
            while nodes[i as usize].ndx != nodes[path as usize].stop_val {
                i -= 1;
            }
            nodes[path as usize].traceb = i;
            nodes[i as usize].traceb = nxt;
        }
        if nodes[path as usize].strand == 1
            && nodes[path as usize].typ == STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].typ == STOP
        {
            let slot = (nodes[path as usize].ndx % 3) as usize;
            let tb = nodes[nxt as usize].star_ptr[slot];
            nodes[path as usize].traceb = tb;
            if tb >= 0 {
                nodes[tb as usize].traceb = nxt;
            }
        }
        if nodes[path as usize].strand == -1
            && nodes[path as usize].typ == STOP
            && nodes[nxt as usize].strand == -1
            && nodes[nxt as usize].typ == STOP
        {
            let slot = (nodes[nxt as usize].ndx % 3) as usize;
            let tb = nodes[path as usize].star_ptr[slot];
            nodes[path as usize].traceb = tb;
            if tb >= 0 {
                nodes[tb as usize].traceb = nxt;
            }
        }
        path = nodes[path as usize].traceb;
        if path < 0 {
            break;
        }
    }

    // Mark forward pointers.
    let mut path = max_ndx;
    while nodes[path as usize].traceb != -1 {
        let tb = nodes[path as usize].traceb;
        nodes[tb as usize].tracef = path;
        path = tb;
    }

    if nodes[max_ndx as usize].traceb == -1 {
        -1
    } else {
        max_ndx
    }
}

/// Score the connection between two nodes (`score_connection`), index form.
fn score_connection(nodes: &mut [Node], p1: usize, p2: usize, tinf: &Tinf, flag: bool) {
    // Copy the scalar fields we need (all Copy).
    let n1_type = nodes[p1].typ;
    let n1_strand = nodes[p1].strand;
    let n1_ndx = nodes[p1].ndx;
    let n1_traceb = nodes[p1].traceb;
    let n1_score = nodes[p1].score;
    let n1_cscore = nodes[p1].cscore;
    let n1_sscore = nodes[p1].sscore;
    let n1_gc = nodes[p1].gc_score;

    let n2_type = nodes[p2].typ;
    let n2_strand = nodes[p2].strand;
    let n2_ndx = nodes[p2].ndx;
    let n2_stop_val = nodes[p2].stop_val;
    let n2_score = nodes[p2].score;
    let n2_cscore = nodes[p2].cscore;
    let n2_sscore = nodes[p2].sscore;
    let n2_gc = nodes[p2].gc_score;

    let mut left = n1_ndx;
    let mut right = n2_ndx;
    let mut ovlp = 0i64;
    let mut maxfr: i64 = -1;
    let mut score = 0.0f64;
    let mut scr_mod = 0.0f64;

    // ---- Invalid connections ----
    if n1_type != STOP && n2_type != STOP && n1_strand == n2_strand {
        return;
    } else if n1_strand == 1 && n1_type != STOP && n2_strand == -1 {
        return;
    } else if n1_strand == -1 && n1_type == STOP && n2_strand == 1 {
        return;
    } else if n1_strand == -1 && n1_type != STOP && n2_strand == 1 && n2_type == STOP {
        return;
    }

    // ---- Edge artifacts ----
    if n1_traceb == -1 && n1_strand == 1 && n1_type == STOP {
        return;
    }
    if n1_traceb == -1 && n1_strand == -1 && n1_type != STOP {
        return;
    }

    let bias = tinf.bias;
    let gc_dot = |g: &[f64; 3]| bias[0] * g[0] + bias[1] * g[1] + bias[2] * g[2];

    // ---- Genes ----
    if n1_strand == n2_strand && n1_strand == 1 && n1_type != STOP && n2_type == STOP {
        if n2_stop_val >= n1_ndx {
            return;
        }
        if n1_ndx % 3 != n2_ndx % 3 {
            return;
        }
        right += 2;
        if !flag {
            scr_mod = gc_dot(&n1_gc);
        } else {
            score = n1_cscore + n1_sscore;
        }
    } else if n1_strand == n2_strand && n1_strand == -1 && n1_type == STOP && n2_type != STOP {
        if nodes[p1].stop_val <= n2_ndx {
            return;
        }
        if n1_ndx % 3 != n2_ndx % 3 {
            return;
        }
        left -= 2;
        if !flag {
            scr_mod = gc_dot(&n2_gc);
        } else {
            score = n2_cscore + n2_sscore;
        }
    }
    // ---- Intergenic ----
    else if n1_strand == 1 && n1_type == STOP && n2_strand == 1 && n2_type != STOP {
        left += 2;
        if left >= right {
            return;
        }
        if flag {
            score = intergenic_mod(&nodes[p1], &nodes[p2], tinf);
        }
    } else if n1_strand == 1 && n1_type == STOP && n2_strand == -1 && n2_type == STOP {
        left += 2;
        right -= 2;
        if left >= right {
            return;
        }
        maxfr = -1;
        let mut maxval = 0.0f64;
        for i in 0..3 {
            let sp = nodes[p2].star_ptr[i];
            if sp == -1 {
                continue;
            }
            let k = sp as usize;
            let n3_stop_val = nodes[k].stop_val;
            let n3_ndx = nodes[k].ndx;
            let cur_ovlp = left - n3_stop_val + 3;
            if cur_ovlp <= 0 || cur_ovlp >= MAX_OPP_OVLP {
                continue;
            }
            if cur_ovlp >= n3_ndx - left {
                continue;
            }
            if n1_traceb == -1 {
                continue;
            }
            if cur_ovlp >= n3_stop_val - nodes[n1_traceb as usize].ndx - 2 {
                continue;
            }
            let val = if flag {
                nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[k], &nodes[p2], tinf)
            } else {
                gc_dot(&nodes[k].gc_score)
            };
            if val > maxval {
                maxfr = i as i64;
                maxval = nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[k], &nodes[p2], tinf);
            }
        }
        if maxfr != -1 {
            let k = nodes[p2].star_ptr[maxfr as usize] as usize;
            if !flag {
                scr_mod = gc_dot(&nodes[k].gc_score);
            } else {
                score = nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[k], &nodes[p2], tinf);
            }
        } else if flag {
            score = intergenic_mod(&nodes[p1], &nodes[p2], tinf);
        }
    } else if n1_strand == -1 && n1_type != STOP && n2_strand == -1 && n2_type == STOP {
        right -= 2;
        if left >= right {
            return;
        }
        if flag {
            score = intergenic_mod(&nodes[p1], &nodes[p2], tinf);
        }
    } else if n1_strand == -1 && n1_type != STOP && n2_strand == 1 && n2_type != STOP {
        if left >= right {
            return;
        }
        if flag {
            score = intergenic_mod(&nodes[p1], &nodes[p2], tinf);
        }
    }
    // ---- Possible operons ----
    else if n1_strand == 1 && n2_strand == 1 && n1_type == STOP && n2_type == STOP {
        if n2_stop_val >= n1_ndx {
            return;
        }
        let slot = (n2_ndx % 3) as usize;
        if nodes[p1].star_ptr[slot] == -1 {
            return;
        }
        let k = nodes[p1].star_ptr[slot] as usize;
        left = nodes[k].ndx;
        right += 2;
        if !flag {
            scr_mod = gc_dot(&nodes[k].gc_score);
        } else {
            score = nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[p1], &nodes[k], tinf);
        }
    } else if n1_strand == -1 && n1_type == STOP && n2_strand == -1 && n2_type == STOP {
        if nodes[p1].stop_val <= n2_ndx {
            return;
        }
        let slot = (n1_ndx % 3) as usize;
        if nodes[p2].star_ptr[slot] == -1 {
            return;
        }
        let k = nodes[p2].star_ptr[slot] as usize;
        left -= 2;
        right = nodes[k].ndx;
        if !flag {
            scr_mod = gc_dot(&nodes[k].gc_score);
        } else {
            score = nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[k], &nodes[p2], tinf);
        }
    }
    // ---- Overlapping opposite-strand 3' ends ----
    else if n1_strand == 1 && n1_type == STOP && n2_strand == -1 && n2_type != STOP {
        if n2_stop_val - 2 >= n1_ndx + 2 {
            return;
        }
        ovlp = (n1_ndx + 2) - (n2_stop_val - 2) + 1;
        if ovlp >= MAX_OPP_OVLP {
            return;
        }
        if (n1_ndx + 2 - n2_stop_val - 2 + 1) >= (n2_ndx - n1_ndx + 3 + 1) {
            return;
        }
        let bnd = if n1_traceb == -1 {
            0
        } else {
            nodes[n1_traceb as usize].ndx
        };
        if (n1_ndx + 2 - n2_stop_val - 2 + 1) >= (n2_stop_val - 3 - bnd + 1) {
            return;
        }
        left = n2_stop_val - 2;
        if !flag {
            scr_mod = gc_dot(&n2_gc);
        } else {
            score = n2_cscore + n2_sscore - 0.15 * tinf.st_wt;
        }
    }

    if !flag {
        score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
    }

    if n1_score + score >= n2_score {
        nodes[p2].score = n1_score + score;
        nodes[p2].traceb = p1 as i64;
        nodes[p2].ov_mark = maxfr;
    }
}

/// Eliminate genes with negative scores (`eliminate_bad_genes`).
fn eliminate_bad_genes(nodes: &mut [Node], dbeg: i64, tinf: &Tinf) {
    if dbeg == -1 {
        return;
    }
    // Walk to the head.
    let mut path = dbeg;
    while nodes[path as usize].traceb != -1 {
        path = nodes[path as usize].traceb;
    }
    let head = path;

    // Apply intergenic modifiers along the chain.
    let mut path = head;
    while nodes[path as usize].tracef != -1 {
        let tf = nodes[path as usize].tracef;
        if nodes[path as usize].strand == 1 && nodes[path as usize].typ == STOP {
            let m = intergenic_mod(&nodes[path as usize], &nodes[tf as usize], tinf);
            nodes[tf as usize].sscore += m;
        }
        if nodes[path as usize].strand == -1 && nodes[path as usize].typ != STOP {
            let m = intergenic_mod(&nodes[path as usize], &nodes[tf as usize], tinf);
            nodes[path as usize].sscore += m;
        }
        path = tf;
    }

    let mut path = head;
    while nodes[path as usize].tracef != -1 {
        let tf = nodes[path as usize].tracef;
        if nodes[path as usize].strand == 1
            && nodes[path as usize].typ != STOP
            && nodes[path as usize].cscore + nodes[path as usize].sscore < 0.0
        {
            nodes[path as usize].elim = true;
            nodes[tf as usize].elim = true;
        }
        if nodes[path as usize].strand == -1
            && nodes[path as usize].typ == STOP
            && nodes[tf as usize].cscore + nodes[tf as usize].sscore < 0.0
        {
            nodes[path as usize].elim = true;
            nodes[tf as usize].elim = true;
        }
        path = tf;
    }
}

/// Minimal gene record used between `add_genes` and `tweak_final_starts`.
#[derive(Clone, Copy)]
struct RawGene {
    begin: i64,
    end: i64,
    start_ndx: i64,
    stop_ndx: i64,
}

/// Extract the gene list from the traceback (`add_genes`).
fn add_genes(nodes: &[Node], dbeg: i64) -> Vec<RawGene> {
    let mut glist = Vec::new();
    if dbeg == -1 {
        return glist;
    }
    let mut path = dbeg;
    while nodes[path as usize].traceb != -1 {
        path = nodes[path as usize].traceb;
    }
    let mut cur = RawGene {
        begin: 0,
        end: 0,
        start_ndx: -1,
        stop_ndx: -1,
    };
    let mut have = false;
    let mut p = path;
    while p != -1 {
        let node = &nodes[p as usize];
        if node.elim {
            p = node.tracef;
            continue;
        }
        if node.strand == 1 && node.typ != STOP {
            cur.begin = node.ndx + 1;
            cur.start_ndx = p;
            have = true;
        }
        if node.strand == -1 && node.typ == STOP {
            cur.begin = node.ndx - 1;
            cur.stop_ndx = p;
            have = true;
        }
        if node.strand == 1 && node.typ == STOP {
            cur.end = node.ndx + 3;
            cur.stop_ndx = p;
            glist.push(cur);
            cur = RawGene {
                begin: 0,
                end: 0,
                start_ndx: -1,
                stop_ndx: -1,
            };
            have = false;
        }
        if node.strand == -1 && node.typ != STOP {
            cur.end = node.ndx + 1;
            cur.start_ndx = p;
            glist.push(cur);
            cur = RawGene {
                begin: 0,
                end: 0,
                start_ndx: -1,
                stop_ndx: -1,
            };
            have = false;
        }
        p = node.tracef;
    }
    let _ = have;
    glist
}

/// Refine close / rare-codon starts (`tweak_final_starts`).
fn tweak_final_starts(genes: &mut [RawGene], nodes: &[Node], tinf: &Tinf) {
    let ng = genes.len();
    let nn = nodes.len() as i64;
    for i in 0..ng {
        let ndx = genes[i].start_ndx;
        if ndx < 0 {
            continue;
        }
        let ndxu = ndx as usize;
        let sc = nodes[ndxu].sscore + nodes[ndxu].cscore;
        let mut igm = 0.0f64;
        if i > 0 {
            let prev_start = genes[i - 1].start_ndx;
            let prev_stop = genes[i - 1].stop_ndx;
            if prev_start >= 0 {
                let ps = nodes[prev_start as usize].strand;
                if nodes[ndxu].strand == 1 && ps == 1 && prev_stop >= 0 {
                    igm = intergenic_mod(&nodes[prev_stop as usize], &nodes[ndxu], tinf);
                } else if nodes[ndxu].strand == 1 && ps == -1 {
                    igm = intergenic_mod(&nodes[prev_start as usize], &nodes[ndxu], tinf);
                }
            }
        }
        if i + 1 < ng {
            let next_start = genes[i + 1].start_ndx;
            let next_stop = genes[i + 1].stop_ndx;
            if next_start >= 0 {
                let ns = nodes[next_start as usize].strand;
                if nodes[ndxu].strand == -1 && ns == 1 {
                    igm = intergenic_mod(&nodes[ndxu], &nodes[next_start as usize], tinf);
                } else if nodes[ndxu].strand == -1 && ns == -1 && next_stop >= 0 {
                    igm = intergenic_mod(&nodes[ndxu], &nodes[next_stop as usize], tinf);
                }
            }
        }

        let mut maxndx = [-1i64; 2];
        let mut maxsc = [0.0f64; 2];
        let mut maxigm = [0.0f64; 2];

        let lo = (ndx - 100).max(0);
        let hi = (ndx + 100).min(nn);
        let mut j = lo;
        while j < hi {
            if j == ndx {
                j += 1;
                continue;
            }
            let ju = j as usize;
            if nodes[ju].typ == STOP || nodes[ju].stop_val != nodes[ndxu].stop_val {
                j += 1;
                continue;
            }
            let mut tigm = 0.0f64;
            let mut skip = false;
            if i > 0 {
                let prev_start = genes[i - 1].start_ndx;
                let prev_stop = genes[i - 1].stop_ndx;
                if prev_start >= 0 {
                    let ps = nodes[prev_start as usize].strand;
                    if nodes[ju].strand == 1 && ps == 1 && prev_stop >= 0 {
                        if nodes[prev_stop as usize].ndx - nodes[ju].ndx > crate::node::MAX_SAM_OVLP {
                            skip = true;
                        } else {
                            tigm = intergenic_mod(&nodes[prev_stop as usize], &nodes[ju], tinf);
                        }
                    } else if nodes[ju].strand == 1 && ps == -1 {
                        if nodes[prev_start as usize].ndx - nodes[ju].ndx >= 0 {
                            skip = true;
                        } else {
                            tigm = intergenic_mod(&nodes[prev_start as usize], &nodes[ju], tinf);
                        }
                    }
                }
            }
            if !skip && i + 1 < ng {
                let next_start = genes[i + 1].start_ndx;
                let next_stop = genes[i + 1].stop_ndx;
                if next_start >= 0 {
                    let ns = nodes[next_start as usize].strand;
                    if nodes[ju].strand == -1 && ns == 1 {
                        if nodes[ju].ndx - nodes[next_start as usize].ndx >= 0 {
                            skip = true;
                        } else {
                            tigm = intergenic_mod(&nodes[ju], &nodes[next_start as usize], tinf);
                        }
                    } else if nodes[ju].strand == -1 && ns == -1 && next_stop >= 0 {
                        if nodes[ju].ndx - nodes[next_stop as usize].ndx > crate::node::MAX_SAM_OVLP {
                            skip = true;
                        } else {
                            tigm = intergenic_mod(&nodes[ju], &nodes[next_stop as usize], tinf);
                        }
                    }
                }
            }
            if skip {
                j += 1;
                continue;
            }

            let cs = nodes[ju].cscore + nodes[ju].sscore;
            if maxndx[0] == -1 {
                maxndx[0] = j;
                maxsc[0] = cs;
                maxigm[0] = tigm;
            } else if cs + tigm > maxsc[0] {
                maxndx[1] = maxndx[0];
                maxsc[1] = maxsc[0];
                maxigm[1] = maxigm[0];
                maxndx[0] = j;
                maxsc[0] = cs;
                maxigm[0] = tigm;
            } else if maxndx[1] == -1 || cs + tigm > maxsc[1] {
                maxndx[1] = j;
                maxsc[1] = cs;
                maxigm[1] = tigm;
            }
            j += 1;
        }

        for k in 0..2 {
            let mndx = maxndx[k];
            if mndx == -1 {
                continue;
            }
            let m = mndx as usize;
            let far = (nodes[m].ndx - nodes[ndxu].ndx).abs();
            if nodes[m].tscore < nodes[ndxu].tscore
                && maxsc[k] - nodes[m].tscore >= sc - nodes[ndxu].tscore + tinf.st_wt
                && nodes[m].rscore > nodes[ndxu].rscore
                && nodes[m].uscore > nodes[ndxu].uscore
                && nodes[m].cscore > nodes[ndxu].cscore
                && far > 15
            {
                maxsc[k] += nodes[ndxu].tscore - nodes[m].tscore;
            } else if far <= 15
                && nodes[m].rscore + nodes[m].tscore > nodes[ndxu].rscore + nodes[ndxu].tscore
                && !nodes[ndxu].edge
                && !nodes[m].edge
            {
                if nodes[ndxu].cscore > nodes[m].cscore {
                    maxsc[k] += nodes[ndxu].cscore - nodes[m].cscore;
                }
                if nodes[ndxu].uscore > nodes[m].uscore {
                    maxsc[k] += nodes[ndxu].uscore - nodes[m].uscore;
                }
                if igm > maxigm[k] {
                    maxsc[k] += igm - maxigm[k];
                }
            } else {
                maxsc[k] = -1000.0;
            }
        }

        let mut mndx: i64 = -1;
        for k in 0..2 {
            if maxndx[k] == -1 {
                continue;
            }
            if mndx == -1 && maxsc[k] + maxigm[k] > sc + igm {
                mndx = k as i64;
            } else if mndx >= 0
                && maxsc[k] + maxigm[k] > maxsc[mndx as usize] + maxigm[mndx as usize]
            {
                mndx = k as i64;
            }
        }
        if mndx != -1 {
            let chosen = maxndx[mndx as usize];
            let c = chosen as usize;
            if nodes[c].strand == 1 {
                genes[i].start_ndx = chosen;
                genes[i].begin = nodes[c].ndx + 1;
            } else {
                genes[i].start_ndx = chosen;
                genes[i].end = nodes[c].ndx + 1;
            }
        }
    }
}

/// Train a model on a (possibly concatenated) sequence under one genetic code.
pub fn train(seq: &Seq, trans_table: u32, closed: bool) -> Tinf {
    let code = GeneticCode::new(trans_table).expect("validated table");
    let mut tinf = Tinf::new(trans_table, seq.gc);

    let mut nodes = add_nodes(seq, &code, closed);
    sort_nodes(&mut nodes);
    if nodes.is_empty() {
        return tinf;
    }
    let gc = gc_frame_plot(seq);
    record_gc_bias(&gc, &mut nodes, &mut tinf);
    record_overlapping_starts(&mut nodes, &tinf, false);
    let ipath = dprog(&mut nodes, &tinf, false);
    calc_dicodon_gene(&mut tinf, seq, &nodes, ipath);
    raw_coding_score(seq, &mut nodes, &tinf);
    rbs_score(seq, &mut nodes, &tinf);
    train_starts_sd(seq, &mut nodes, &mut tinf);
    determine_sd_usage(&mut tinf);
    tinf
}

/// Predict genes on one contig with a trained model. Returns the genes and the
/// summed coding score (used by code auto-detection).
pub fn predict(seq: &Seq, contig: &str, tinf: &Tinf, closed: bool, is_meta: bool) -> (Vec<Gene>, f64) {
    let code = tinf.code();
    let mut nodes = add_nodes(seq, &code, closed);
    sort_nodes(&mut nodes);
    if nodes.is_empty() {
        return (Vec::new(), 0.0);
    }
    score_nodes(seq, &mut nodes, tinf, closed, is_meta);
    record_overlapping_starts(&mut nodes, tinf, true);
    let ipath = dprog(&mut nodes, tinf, true);
    eliminate_bad_genes(&mut nodes, ipath, tinf);
    let mut raw = add_genes(&nodes, ipath);
    tweak_final_starts(&mut raw, &nodes, tinf);

    let mut coding_total = 0.0f64;
    let mut genes = Vec::with_capacity(raw.len());
    for (k, g) in raw.iter().enumerate() {
        if g.start_ndx < 0 || g.stop_ndx < 0 {
            continue;
        }
        let sidx = g.start_ndx as usize;
        let stopidx = g.stop_ndx as usize;
        let strand = nodes[sidx].strand;
        let (begin, end) = if g.begin <= g.end {
            (g.begin, g.end)
        } else {
            (g.end, g.begin)
        };
        let begin = begin.max(1);
        let end = end.min(seq.len() as i64);
        if end <= begin {
            continue;
        }

        let cscore = nodes[sidx].cscore;
        let sscore = nodes[sidx].sscore;
        coding_total += cscore;
        let total = cscore + sscore;
        let conf = calculate_confidence(total, tinf.st_wt);

        let start_type = match nodes[sidx].typ {
            0 => 0u8,
            1 => 1,
            2 => 2,
            _ => 3,
        };
        let partial_left = begin == 1 && (nodes[if strand == 1 { sidx } else { stopidx }].edge);
        let partial_right = end == seq.len() as i64
            && (nodes[if strand == 1 { stopidx } else { sidx }].edge);

        let mut gene = Gene {
            contig: contig.to_string(),
            id: k + 1,
            begin: begin as usize,
            end: end as usize,
            strand: if strand == 1 { 1 } else { -1 },
            start_type,
            rbs_motif: "None".to_string(),
            rbs_spacer: "None".to_string(),
            partial_left,
            partial_right,
            score: total,
            cscore,
            sscore,
            conf,
            gc: nodes[sidx].gc_cont,
            trans_table: tinf.trans_table,
            engine: Engine::Dynamic,
            aa: String::new(),
            nt: String::new(),
        };
        gene.fill_sequences(seq, &code);
        genes.push(gene);
    }
    // Renumber sequentially.
    for (k, g) in genes.iter_mut().enumerate() {
        g.id = k + 1;
    }
    (genes, coding_total)
}

/// Build a single training sequence from many contigs, separated by
/// stop-forcing gaps (`TTAATTAATTAA`).
pub fn build_training_seq(contigs: &[(String, Seq)]) -> Seq {
    if contigs.len() == 1 {
        return contigs[0].1.clone();
    }
    let mut bytes: Vec<u8> = Vec::new();
    for (i, (_, s)) in contigs.iter().enumerate() {
        if i > 0 {
            bytes.extend_from_slice(TRAIN_GAP);
        }
        for &c in &s.fwd {
            bytes.push(crate::seq::decode_base(c));
        }
    }
    Seq::from_ascii(&bytes)
}

/// Choose the genetic code by self-training under each candidate.
///
/// Naively picking the highest coding score is biased: codes with fewer stop
/// codons (4, 15, 25 — two stops each) let ORFs read through what would be a
/// stop under the standard code (11 / 1 — three stops), which inflates apparent
/// coding density even on ordinary genomes. To avoid spurious reassignments we
/// treat the standard code as a baseline and only adopt an alternative code
/// when it beats the baseline by a clear margin. Genuine alternative-code
/// genomes clear this bar easily, because using the wrong (standard) code
/// fragments many real genes at read-through stops and sharply lowers the
/// score; mere read-through inflation does not.
pub fn autodetect_code(train_seq: &Seq, candidates: &[u32], closed: bool, is_meta: bool) -> (u32, Tinf) {
    /// An alternative code must exceed the standard code's score by this factor.
    const ALT_MARGIN: f64 = 1.10;

    let mut results: Vec<(u32, Tinf, f64)> = Vec::new();
    for &tt in candidates {
        if GeneticCode::new(tt).is_none() {
            continue;
        }
        let tinf = train(train_seq, tt, closed);
        let (_g, score) = predict(train_seq, "__train__", &tinf, closed, is_meta);
        results.push((tt, tinf, score));
    }
    if results.is_empty() {
        let tinf = train(train_seq, 11, closed);
        return (11, tinf);
    }

    // Baseline = the standard code if available (prefer 11, then 1).
    let baseline_idx = results
        .iter()
        .position(|(tt, _, _)| *tt == 11)
        .or_else(|| results.iter().position(|(tt, _, _)| *tt == 1));

    // Highest-scoring candidate overall.
    let mut best_idx = 0usize;
    for (i, (_, _, s)) in results.iter().enumerate() {
        if *s > results[best_idx].2 {
            best_idx = i;
        }
    }

    let chosen = match baseline_idx {
        None => best_idx, // no standard code among candidates; take the best
        Some(bi) => {
            if best_idx == bi {
                bi
            } else {
                let base_score = results[bi].2;
                let best_score = results[best_idx].2;
                // Require a clear margin (handle non-positive baselines safely).
                let threshold = if base_score > 0.0 {
                    base_score * ALT_MARGIN
                } else {
                    base_score // any improvement over a non-positive baseline counts
                };
                if best_score > threshold {
                    best_idx
                } else {
                    bi
                }
            }
        }
    };

    let (tt, tinf, _) = results.swap_remove(chosen);
    (tt, tinf)
}

// ---------------------------------------------------------------------------
// Preset entry points
// ---------------------------------------------------------------------------

/// Run the DP engine over a set of contigs with an explicit genetic-code policy
/// and engine tag.
///
/// * `fixed = Some(table)` forces a genetic code; `None` auto-detects among
///   `candidates` (with DYNAMITE's standard-code bias).
/// * Per-contig prediction is parallelised with rayon.
///
/// Returns `(genes, resolved_table)`. This is the shared core behind the
/// `prodigal`, `prodigalgv`, and `prodigal-meta` presets.
pub fn dp_call(
    contigs: &[(String, Seq)],
    fixed: Option<u32>,
    candidates: &[u32],
    is_meta: bool,
    closed: bool,
    engine: Engine,
) -> (Vec<Gene>, u32) {
    if contigs.is_empty() {
        return (Vec::new(), fixed.unwrap_or(11));
    }
    let train_seq = build_training_seq(contigs);
    let (code_id, tinf) = match fixed {
        Some(id) => (id, train(&train_seq, id, closed)),
        None => autodetect_code(&train_seq, candidates, closed, is_meta),
    };
    let genes: Vec<Gene> = contigs
        .par_iter()
        .flat_map_iter(|(id, s)| {
            let (mut g, _score) = predict(s, id, &tinf, closed, is_meta);
            for x in &mut g {
                x.engine = engine;
            }
            g.into_iter()
        })
        .collect();
    (genes, code_id)
}

/// **Prodigal** preset — bacteria and archaea: a single self-trained genome
/// under the standard bacterial/archaeal code 11 (override with `fixed`).
pub fn call(contigs: &[(String, Seq)], fixed: Option<u32>, closed: bool) -> (Vec<Gene>, u32) {
    dp_call(contigs, Some(fixed.unwrap_or(11)), &[11], false, closed, Engine::Prodigal)
}

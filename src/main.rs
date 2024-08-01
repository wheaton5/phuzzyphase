#[macro_use]

extern crate clap;
extern crate bio;
extern crate hashbrown;
extern crate rand;
extern crate rayon;
extern crate threadpool;
extern crate statrs;
extern crate petgraph;
extern crate sanitize_filename;
extern crate log;
extern crate env_logger;
extern crate pprof;

use log::{debug, error, log_enabled, info, trace, Level};

use std::sync::mpsc::channel;
use std::thread::Thread;


use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::process::Command;
use std::{thread, time};
use std::str;
use std::cmp::Ordering;
use threadpool::ThreadPool;
use std::sync::{Arc, Barrier};

use statrs::distribution::{Binomial};
use statrs::distribution::Discrete;
use bio::alignment::pairwise::banded;
use bio::io::fasta;
use bio::utils::TextSlice;
use petgraph::unionfind::UnionFind;

use statrs::distribution::Multinomial;

use failure::{Error, ResultExt};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::Format;
use rust_htslib::bcf::{self, Read as BcfRead};
use std::convert::TryInto;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use crate::statrs::distribution::DiscreteCDF;

use hashbrown::{HashMap, HashSet};
use std::collections::BinaryHeap;

use clap::App;

const K: usize = 6; // kmer match length
const W: usize = 20; // Window size for creating the band
const MATCH: i32 = 2; // Match score
const MISMATCH: i32 = -4; // Mismatch score
const GAP_OPEN: i32 = -4; // Gap open score
const GAP_EXTEND: i32 = -2; // Gap extend score

fn main() {
    env_logger::init();
    let guard = pprof::ProfilerGuard::new(100).unwrap();

    let result = _main();
    if let Ok(report) = guard.report().build() {
        let file = File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };
    if let Err(v) = result {
        let fail = v.as_fail();
        println!("Phasstphase error. v{}.", env!("CARGO_PKG_VERSION"));
        println!("Error: {}", fail);

        for cause in fail.iter_causes() {
            println!("Info: caused by {}", cause);
        }

        println!("\n{}", v.backtrace());

        println!("If you think this is bug in Phasstphase, please file a bug at https://github.com/wheaton5/phasstphase, and include the information above and the command-line you used.");
        std::process::exit(1)
    }
}

fn _main() -> Result<(), Error> {
    println!("Welcome to phasst phase!");
    let params = load_params();
    //if params.long_read_bam == None && params.linked_read_bam == None && params.hic_bam == None {
    //    eprintln!("Must supply at least one bam");
    //    std::process::exit(1);
    //}
    if Path::new(&params.output).is_dir() {
        eprintln!("restarting from partial output in {}", params.output);
    } else {
        fs::create_dir(params.output.to_string())?;
    }

    let fai = params.fasta.to_string() + ".fai";
    let fa_index_iter = fasta::Index::from_file(&fai)
        .expect(&format!("error opening fasta index: {}", fai))
        .sequences();
    let mut chroms: Vec<String> = Vec::new();
    
    let mut starts: Vec<usize> = Vec::new();
    let mut ends: Vec<usize> = Vec::new();
    if let Some(chrom) = params.chrom {
        chroms.push(chrom);
        starts.push(params.start.expect("could not unwrap region start"));
        ends.push(params.end.expect("could not unwrap region end"));
    } else {
        for chrom in fa_index_iter {
            chroms.push(chrom.name.to_string());
            starts.push(0);
            ends.push(chrom.len as usize);
        }
    }
   
    //let vcf_reader = bcf::IndexedReader::from_path(params.vcf.to_string())?;
    let mut chunks: Vec<ThreadData> = Vec::new();
    for (i, chrom) in chroms.iter().enumerate() {
        //if chrom.chars().count() > 5 { continue; }
        let data = ThreadData {
            index: i,
            long_read_bam: match &params.long_read_bam {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            hic_bam: match &params.hic_bam {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            fasta: params.fasta.to_string(),
            chrom: chrom.to_string(),
            vcf: params.vcf.to_string(),
            min_mapq: params.min_mapq,
            min_base_qual: params.min_base_qual,
            start: starts[i],
            end: ends[i],
            window: params.window,
            output: params.output.to_string(),
            vcf_out: format!("{}/chrom_{}.vcf.gz", params.output, sanitize_filename::sanitize(chrom)),
            phased_vcf_out: format!("{}/phased_chrom_{}.vcf.gz", params.output, sanitize_filename::sanitize(chrom)),
            vcf_out_done: format!("{}/chrom_{}.vcf.done", params.output, sanitize_filename::sanitize(chrom)),
            phased_vcf_done: format!("{}/phased_chrom_{}.vcf.done", params.output, sanitize_filename::sanitize(chrom)),
            phasing_window: params.phasing_window,
            seed: params.seed,
            ploidy: params.ploidy,
            hic_phasing_posterior_threshold: params.hic_phasing_posterior_threshold,
            long_switch_threshold: params.long_switch_threshold,
        };
        chunks.push(data);
    }
    println!("{} chromosome chunks and {} threads",chunks.len(), params.threads);
    //let pool = rayon::ThreadPoolBuilder::new()
    //    .num_threads(params.threads)
    //    .build()
    //    .unwrap();
    //let results: Vec<_> = pool.install(|| {
    //    chunks
    //        .par_iter()
    //        .map(|rec_chunk| phase_chunk(&rec_chunk))
    //        .collect()
    //});
    //for data in chunks {
    //   pool.spawn(move || phase_chunk(&data).expect("thread failed"));
    //}
    let pool = ThreadPool::new(params.threads);
    let (tx, rx) = channel();
    let final_output = format!("{}/phasstphase.vcf.gz",params.output);
    let mut cmd: Vec<String> = Vec::new();
    cmd.push("concat".to_string());
    cmd.push("-o".to_string());
    cmd.push(final_output.to_string());
    let jobs = chunks.len();
    for data in chunks {
        cmd.push(format!("{}/phased_chrom_{}.vcf.gz", data.output, data.chrom));
        let tx = tx.clone();
        pool.execute(move|| {
            phase_chunk(&data).expect("thread failed");
            tx.send(1).expect("wait for me");
        }
        );
    }
    println!("before barrier");
    assert_eq!(rx.iter().take(jobs).fold(0,|a,b|a+b),jobs);
    println!("after barrier");
    println!("merging final vcf");
    let result = Command::new("bcftools")
        .args(&cmd)
        .status()
        .expect("bcftools failed us");
    let result = Command::new("tabix")
        .args(&["-p", "vcf", &final_output])
        .status()
        .expect("couldn't tabix index final vcf");
    
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct PhaseBlock {
    id: usize,
    start_index: usize,
    start_position: usize,
    end_index: usize,
    end_position: usize,
}

struct MetaPhaseBlock {
    id: usize,
    phase_blocks: Vec<PhaseBlock>,
}

fn maximization(
    molecules: &Vec<Vec<Allele>>,
    posteriors: &Vec<Vec<f32>>,
    cluster_centers: &mut Vec<Vec<f32>>,
    min_index: &mut usize,
    max_index: &mut usize,
    molecule_support: &mut Vec<Vec<f32>>,
) -> f32 {
    let mut updates: HashMap<usize, Vec<(f32, f32)>> = HashMap::new(); // variant index to vec across
    let mut variant_molecule_count: HashMap<usize, usize> = HashMap::new();
    *min_index = std::usize::MAX;
    *max_index = 0;
    // haplotype clusters to a tuple of (numerator, denominator)
    for molecule_index in 0..molecules.len() {
        // if molecule does not support any haplotype over another, dont use it in maximization
        let mut different = false;
        for haplotype in 0..cluster_centers.len() {
            //if (posteriors[molecule_index][haplotype] - posteriors[molecule_index][0]).abs() < 0.01 { // but why would this happen?
            if posteriors[molecule_index][haplotype] == posteriors[molecule_index][0] {
                different = true;
            }
        }
        //if !different {
            //println!("debug, molecule does not support any haplotype over another");
        //    continue;
        //}
        for allele in &molecules[molecule_index] {
            let variant_index = allele.index;
            let alt = allele.allele;
            let count = variant_molecule_count.entry(variant_index).or_insert(0);
            *count += 1;
            for haplotype in 0..cluster_centers.len() {
                let posterior = posteriors[molecule_index][haplotype];
                let numerators_denominators = updates
                    .entry(variant_index)
                    .or_insert(vec![(0.0, 0.0); cluster_centers.len()]);
                if alt {
                    numerators_denominators[haplotype].0 += posterior;
                    numerators_denominators[haplotype].1 += posterior;
                } else {
                    numerators_denominators[haplotype].1 += posterior;
                }
            }
        }
    }
    let mut updated_variants = 0;
    let mut total_change = 0.0;
    for (variant_index, haplotypes) in updates.iter() {
        if variant_molecule_count.get(variant_index).unwrap() > &3 { // TODO dont hard code, make parameters
            updated_variants += 1;
            //TODO Dont hard code stuff
            *min_index = (*min_index).min(*variant_index);
            *max_index = (*max_index).max(*variant_index);
            for haplotype in 0..haplotypes.len() {
                let numerators_denominators = haplotypes[haplotype];
                let allele_fraction = numerators_denominators.0 / numerators_denominators.1;
                let cluster_value = allele_fraction
                    .max(0.001)
                    .min(0.999);
                total_change +=
                    (cluster_value - cluster_centers[haplotype][*variant_index]).abs();
                
                cluster_centers[haplotype][*variant_index] = cluster_value;
    
                molecule_support[haplotype][*variant_index] = numerators_denominators.1;

            }
        }
    }
    total_change
}

// molecules is vec of molecules each a vec of alleles
// cluster centers is vec of haplotype clusters  by variant loci
// posteriors (return value) is molecules by clusters to posterior probability
fn expectation(
    molecules: &Vec<Vec<Allele>>,
    cluster_centers: &Vec<Vec<f32>>,
    posteriors: &mut Vec<Vec<f32>>,
) -> (bool, f32, f32) {
    let mut delta = 0.0;
    //let mut posteriors: Vec<Vec<f32>> = Vec::new();
    let mut any_different = false;
    let mut log_likelihood: f32 = 0.0;
    for (moldex, molecule) in molecules.iter().enumerate() {
        let mut log_probs: Vec<f32> = Vec::new(); // for each haplotype
        for haplotype in 0..cluster_centers.len() {
            let mut log_prob = 0.0; // log(0) = probability 1
            for allele in molecule.iter() {
                let lp;
                if allele.allele {
                    // alt allele
                    lp = cluster_centers[haplotype][allele.index].ln();
                    log_prob += cluster_centers[haplotype][allele.index].ln(); // adding in log space, multiplying in probability space
                } else {
                    lp = (1.0 - cluster_centers[haplotype][allele.index]).ln();
                    log_prob += (1.0 - cluster_centers[haplotype][allele.index]).ln();
                }
            }
            
            log_probs.push(log_prob);
        }
        let bayes_log_denom = log_sum_exp(&log_probs);
        log_likelihood += bayes_log_denom;
        let mut mol_posteriors: Vec<f32> = Vec::new();

        for log_prob in log_probs {
            mol_posteriors.push((log_prob - bayes_log_denom).exp());
        }

        for haplotype in 0..cluster_centers.len() {
            if mol_posteriors[haplotype] != mol_posteriors[0] {
                any_different = true;
            }
            //error!("previous posterior moldex {} hap {} {} to {} for a difference of {}", 
            //    moldex, haplotype, posteriors[moldex][haplotype], mol_posteriors[haplotype], 
            //    (posteriors[moldex][haplotype] - mol_posteriors[haplotype]).abs());
            delta += (posteriors[moldex][haplotype] - mol_posteriors[haplotype]).abs();
            posteriors[moldex][haplotype] = mol_posteriors[haplotype];
        }

        //posteriors.push(mol_posteriors);

    }
    (!any_different, log_likelihood, delta)
}

fn log_sum_exp(p: &Vec<f32>) -> f32 {
    let max_p: f32 = p.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let sum_rst: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn log_sum_exp64(p: &Vec<f64>) -> f64 {
    let max_p: f64 = p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let sum_rst: f64 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn phase_chunk(data: &ThreadData) -> Result<(), Error> {
    println!("checking for file {}", data.phased_vcf_done);
    if Path::new(&data.phased_vcf_done).exists() {
        println!("phasing complete for chrom {}. Delete .done or output directory if you want to rerun", data.chrom);
        return Ok(());
    } else { println!("could not find file, continuing to phasing");}
    println!("checking for file {}", data.vcf_out_done);
    if !Path::new(&data.vcf_out_done).exists() {
        println!("could not find file, gathing variant assignments");
        get_all_variant_assignments(data).expect("we error here at get all variant assignments");
    } else { println!("found file, continuing from previous partial run");}
    let mut vcf_reader = bcf::IndexedReader::from_path(format!("{}", data.vcf_out.to_string()))
        .expect("could not open indexed vcf reader on output vcf");
    let chrom = vcf_reader
        .header()
        .name2rid(data.chrom.as_bytes())
        .expect("can't get chrom rid, make sure vcf and bam and fasta contigs match!");
    let vcf_info = inspect_vcf(&mut vcf_reader, &data);
    let (mut cluster_centers, mut molecule_support) =  init_cluster_centers(vcf_info.num_variants, &data);
    let mut window_start: usize = data.start;
    let mut window_end: usize = window_start + data.phasing_window;
    //let mut position_to_index: HashMap<usize, usize> = HashMap::new();
    //let mut position_so_far: usize = 0;
    let mut phase_blocks: Vec<PhaseBlock> = Vec::new();
    let mut putative_phase_blocks: Vec<PhaseBlock> = Vec::new();
    let mut phase_block_start: usize = 0;
    let mut last_attempted_index: usize = 0;
    let mut in_phaseblock = false;
    let mut last_window_start: Option<usize> = None;
    let seed = [data.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    'outer: while window_start < data.end
        && window_start < vcf_info.final_position as usize
    {
        if let Some(last_start) = last_window_start {
            if last_start >= window_start {
                error!("last window start is >= current start {} >= {}", last_start, window_start);
            }
        }
        last_window_start = Some(window_start);
        let mut cluster_center_delta: f32 = 10.0;
        vcf_reader
            .fetch(chrom, window_start as u64, Some(window_end as u64))
            .expect("some actual error");
        let (molecules, first_var_index, last_var_index) = get_read_molecules(&mut vcf_reader, &vcf_info, READ_TYPE::HIFI);
        let mut iteration = 0;
        let mut min_index: usize = 0;
        let mut max_index: usize = 0;
        let mut last_cluster_center_delta = cluster_center_delta;
        let mut last_posterior_delta = 0.0;
        info!("next window {}-{}", window_start, window_end);
        let mut posteriors: Vec<Vec<f32>> = Vec::new();
        for moldex in molecules.iter() {
            let mut post: Vec<f32> = Vec::new();
            for hap in cluster_centers.iter() {
                post.push(0.0);
            }
            posteriors.push(post);
        }
        while cluster_center_delta > 0.01 {
            
            let (breaking_point, _log_likelihood, posterior_delta) = expectation(&molecules, &cluster_centers, &mut posteriors);

            if in_phaseblock && breaking_point {
                in_phaseblock = false;
                while cluster_centers[0][last_attempted_index] == cluster_centers[1][last_attempted_index] && 
                    vcf_info.variant_positions[last_attempted_index] > window_start {
                    last_attempted_index -= 1;
                }

                // WE CHANGED THIS AFTER PREVIOUS TESTING
                if vcf_info.variant_positions[last_attempted_index] < window_start {
                    last_attempted_index += 1;
                }
                // END
                
                putative_phase_blocks.push(PhaseBlock{
                    id: putative_phase_blocks.len(),
                    start_index: phase_block_start,
                    start_position: vcf_info.variant_positions[phase_block_start],
                    end_index: last_attempted_index,
                    end_position: vcf_info.variant_positions[last_attempted_index]
                });
                
                // WE CHANGED THIS AFTER PREVIOUS TESTING
                phase_block_start = last_attempted_index + 1;
                if phase_block_start >= vcf_info.variant_positions.len() {
                    break 'outer;
                }
                window_start = vcf_info.variant_positions[phase_block_start];
                // END
                /** WE COMMENTED THIS OUT AFTER PREVIOUS TESTING
                phase_block_start = last_attempted_index;
                last_window_start = Some(window_start);
                if let Some(previous_window_start) = last_window_start {
                    while previous_window_start >= window_start {
                        phase_block_start += 1;
                        
                        if phase_block_start >= vcf_info.variant_positions.len() {
                            break 'outer;
                        }
                        window_start = vcf_info.variant_positions[phase_block_start];
                    }
                } else {
                    phase_block_start += 1;
                    if phase_block_start >= vcf_info.variant_positions.len() {
                        break 'outer;
                    } // we are at the end of the region or chromosome, we are done
                    window_start = vcf_info.variant_positions[phase_block_start];
                }
                **/
                eprintln!("in phaseblock but hit breakpoint, reseting window start to {}", window_start);
                window_end = window_start + data.phasing_window;
                for haplotype in 0..cluster_centers.len() {
                    cluster_centers[haplotype][phase_block_start] =
                        rng.gen::<f32>().min(0.98).max(0.02);
                }
                in_phaseblock = false;
                continue 'outer;
            } else if !in_phaseblock && breaking_point {
                
                eprintln!("not in a phaseblock, but hit breaking point. resetting window");
                for haplotype in 0..cluster_centers.len() {
                    cluster_centers[haplotype][phase_block_start] = 0.5;
                }
                if first_var_index == usize::MAX { // there were no variants in this window
                    // we need to find a new starting point
                    eprintln!("\tthere were no variants, scanning forward for window start");
                    for (index, position) in vcf_info.variant_positions.iter().enumerate() {
                        phase_block_start = index;
                        if *position > window_end {
                            break;
                        }
                    }
                    //eprintln!("unable to start phaseblock, resetting");
                } else {
                    if last_var_index >= vcf_info.variant_positions.len() {
                        eprintln!("\tthere were variants and last_var_index was {} at {}, next is {} at {} and current window is {}-{}",
                        last_var_index, vcf_info.variant_positions[last_var_index], 
                        last_var_index + 1, vcf_info.variant_positions[last_var_index + 1], window_start, window_end);
                    }
                    
                    //eprintln!("unphased variants {}-{} positions {}-{}", first_var_index, last_var_index,
                    //vcf_info.variant_positions[first_var_index], vcf_info.variant_positions[last_var_index]);
                    for (index, position) in vcf_info.variant_positions.iter().enumerate() {
                        phase_block_start = index;
                        if *position > window_start {
                            break;
                        }
                    }
                    //phase_block_start = last_var_index + 1;
                }
                
                if phase_block_start >= vcf_info.variant_positions.len() {
                    break 'outer;
                } // we are at the end of the region or chromosome
                last_window_start = Some(window_start);
                window_start = vcf_info.variant_positions[phase_block_start];
                eprintln!("reseting window start to {}", window_start);
                window_end = window_start + data.phasing_window;
                let seed = [data.seed; 32];
                let mut rng: StdRng = SeedableRng::from_seed(seed);
                for haplotype in 0..cluster_centers.len() {
                    cluster_centers[haplotype][phase_block_start] =
                        rng.gen::<f32>().min(0.98).max(0.02);
                }
                continue 'outer;
            }
            cluster_center_delta = maximization(
                &molecules,
                &posteriors,
                &mut cluster_centers,
                &mut min_index,
                &mut max_index,
                &mut molecule_support
            );
            last_cluster_center_delta = cluster_center_delta;
            last_posterior_delta = posterior_delta;
            if max_index != 0 {
                last_attempted_index = max_index;
                if !in_phaseblock {
                    phase_block_start = min_index;
                }
                in_phaseblock = true;
            }

            iteration += 1;
        }
        
        window_start += data.phasing_window / 4;
        window_start = window_start.min(vcf_info.final_position as usize);
        window_end += data.phasing_window / 4;
        window_end = window_end.min(vcf_info.final_position as usize);
    } // end 'outer loop
    println!("DONE phasing long reads! thread {} chrom {}", data.index, data.chrom);
    if in_phaseblock {
        putative_phase_blocks.push(PhaseBlock{
            id: putative_phase_blocks.len(),
            start_index: phase_block_start,
            start_position: vcf_info.variant_positions[phase_block_start],
            end_index: last_attempted_index,
            end_position: vcf_info.variant_positions[last_attempted_index]
        });
    }

    let mut new_phaseblock_id: usize = 0;
    for putative_phase_block in putative_phase_blocks {
        let mut cut_blocks = test_long_switch(putative_phase_block.start_index, putative_phase_block.end_index, &mut cluster_centers, &vcf_info, &mut vcf_reader, &data);
        for mut phase_block in cut_blocks {
            phase_block.id = new_phaseblock_id;
            new_phaseblock_id += 1;
            phase_blocks.push(phase_block);
        }
    }
    
    debug!("DONE long switch test! thread {} chrom {}", data.index, data.chrom);
    let mut total_gap_length = 0;
    for (id, phase_block) in phase_blocks.iter().enumerate() {
        let mut gap = phase_block.start_position;
        if id > 0 {
            gap = phase_block.start_position - phase_blocks[id - 1].end_position;
        }
        total_gap_length += gap;
        info!(
            "phase block {} from {}-{}, {}-{} length {} with gap from last of {}",
            id,
            phase_block.start_position,
            phase_block.end_position,
            phase_block.start_index,
            phase_block.end_index,
            phase_block.end_position - phase_block.start_position,
            gap
        );
    }
    if phase_blocks.len() > 0 {
        info!("and final gap of {}", data.end - phase_blocks[phase_blocks.len() - 1].end_position);
    }
    
    info!("with total gap length of {}", total_gap_length);
    // get phaseblock N50... 
    let mut sizes: Vec<usize> = Vec::new();
    let mut total: usize = 0;
    for phase_block in phase_blocks.iter() {
        let length = phase_block.end_position - phase_block.start_position;
        sizes.push(length);
        total += length;
    }
    sizes.sort_by(|a, b| b.cmp(a));
    let mut so_far = 0;
    for size in sizes {
        so_far += size;
        if so_far > total/2 {
            println!("N50 phase blocks from long reads for chrom {} is {}", data.chrom, size);
            break;
        }
    }

    let allele_phase_block_id: HashMap<usize, usize> = phase_phaseblocks(data, &mut cluster_centers, &phase_blocks, &vcf_info); 
    println!("DONE hic phasing! thread {} chrom {}", data.index, data.chrom);
    
    output_phased_vcf(data, cluster_centers, phase_blocks, &vcf_info, &molecule_support, &allele_phase_block_id);
    let result = Command::new("bcftools")
        .args(&["index", "-f", &data.vcf_out])
        .status()
        .expect("bcftools failed us");
    fs::File::create(data.phased_vcf_done.to_string()).expect("cant create .done file. are the permissions wrong?");
    println!("thread {} chrom {}finished", data.index, data.chrom);
    Ok(())
}

fn test_long_switch(start_index: usize, end_index: usize, 
    cluster_centers: &mut Vec<Vec<f32>>, vcf_info: &VCF_info, 
    vcf_reader: &mut bcf::IndexedReader, data: &ThreadData) -> Vec<PhaseBlock> {
    let mut to_return: Vec<PhaseBlock> = Vec::new();
    //if data.ploidy > 2 {
    //to_return.push(PhaseBlock{
    //    start_index: start_index,
    //    start_position: vcf_info.variant_positions[start_index],
    //    end_index: end_index + 1,
    //    end_position: vcf_info.variant_positions[end_index],
    //    id: 0,
    //});
    //return to_return; // currently not doing test_long_switch for polyploid
    //}



    let pairings = pairings(data.ploidy);
    let log_prior = (1.0/(pairings.len() as f32)).ln();
    let chrom = vcf_reader
        .header()
        .name2rid(data.chrom.as_bytes())
        .expect("can't get chrom rid, make sure vcf and bam and fasta contigs match!");
    let mut phase_block_start = start_index;
    let mut cluster_center_copies: Vec<Vec<Vec<f32>>> = Vec::new(); // pairings by cluster center copies
    for pairing in pairings.iter() {
        cluster_center_copies.push(swap_full(&cluster_centers, &pairing));
    }
    let one_million: usize = 1000000;
    let mut which_million: usize = vcf_info.variant_positions[start_index] / one_million;
    // TODO REALLY THIS TIME URGENT -- if we cut a phaseblock, we need to only test going forward, not read back into old phaseblock
    let mut min_position: u64 = vcf_info.variant_positions[start_index] as u64;
    for breakpoint in start_index..end_index {
        
        let mut log_likelihoods: Vec<f32> = Vec::new();
        let position = vcf_info.variant_positions[breakpoint];
        /*
        if position / one_million > which_million {
            println!("{}", position);
            println!("{} > {} == {}", position / one_million, which_million, position / one_million > which_million);
            which_million = position / one_million;
            println!(" which million is now {}", which_million);
        } 
        */
        let vcf_fetch_start = (position as u64).checked_sub(15000).unwrap_or(0).max(min_position); // so position - 15000 but not going negative or lower than min_position
        let vcf_fetch_end = Some(position as u64 + 15000); 
        vcf_reader
            .fetch(chrom, vcf_fetch_start, vcf_fetch_end) // TODO dont hard code things
            .expect("could not fetch in vcf");
        let (molecules, first_var_index, last_var_index) = get_read_molecules(vcf_reader, &vcf_info, READ_TYPE::HIFI);
        for (index, pairing) in pairings.iter().enumerate() {
            // TODO just deleted this, but might need it back
            //swap(cluster_centers, breakpoint, &pairing, 50);
            swap_copied(&mut cluster_center_copies[index], breakpoint, &pairing);
            let mut posteriors: Vec<Vec<f32>> = Vec::new();
            for moldex in molecules.iter() {
                let mut post: Vec<f32> = Vec::new();
                for hap in cluster_centers.iter() {
                    post.push(0.0);
                }
                posteriors.push(post);
            }
            let (_break, log_likelihood, post_delta) = expectation(&molecules, &cluster_center_copies[index], &mut posteriors);
            log_likelihoods.push(log_likelihood + log_prior);
            // TODO just deleted this, if we screwed up we might need it back
            //swap(cluster_centers, breakpoint, &pairing, 50); // reversing the swap
        }
        let log_bayes_denom = log_sum_exp(&log_likelihoods);
        let log_posterior = log_likelihoods[0] - log_bayes_denom;
        let posterior = log_posterior.exp();
        if posterior < data.long_switch_threshold {
            // end phase block and add to to_return
            println!("breaking phaseblock at {}", position);
            to_return.push(PhaseBlock {
                id: to_return.len(),
                start_index: phase_block_start,
                start_position: vcf_info.variant_positions[phase_block_start],
                end_index: breakpoint+1,
                end_position: vcf_info.variant_positions[breakpoint],
            });
            if breakpoint < end_index {
                min_position = vcf_info.variant_positions[breakpoint + 1] as u64;
            }
            phase_block_start = breakpoint + 1;
            let start_position = vcf_info.variant_positions[start_index];
            let end_position = vcf_info.variant_positions[end_index];
            let position = vcf_info.variant_positions[breakpoint];
            //vcf_reader
            //    .fetch(chrom, (position as u64).checked_sub(15000).unwrap_or(0), Some(position as u64))
            //    .expect("could not fetch in vcf");
            //let (molecules, first_var_index, last_var_index)  = get_read_molecules(vcf_reader, &vcf_info, READ_TYPE::HIFI);
            trace!("HIT POTENTIAL LONG SWITCH ERROR. phase block from indexes {}-{}, positions {}-{}, posterior {}, breakpoint {} position {} with {} molecules", 
                start_index, end_index, start_position, end_position, posterior, breakpoint, position, molecules.len());
            let mut posteriors: Vec<Vec<f32>> = Vec::new();
            for moldex in molecules.iter() {
                let mut post: Vec<f32> = Vec::new();
                for hap in cluster_centers.iter() {
                    post.push(0.0);
                }
                posteriors.push(post);
            }
            let (_break, log_likelihood, post_delta) = expectation(&molecules, &cluster_centers, &mut posteriors);
            //eprintln!("mol posteriors {:?}", posteriors);

            //trace!("hap1 {:?}", &cluster_centers[0][breakpoint..(breakpoint+10)]);
            //trace!("hap2 {:?}", &cluster_centers[1][breakpoint..(breakpoint+10)]);
        
        }
        /*** else {
            let start_position = vcf_info.variant_positions[start_index];
            let end_position = vcf_info.variant_positions[end_index];
            let position = vcf_info.variant_positions[breakpoint];
            vcf_reader
                .fetch(chrom, (position as u64).checked_sub(15000).unwrap_or(0), Some(position as u64))
                .expect("could not fetch in vcf");
                let (molecules, first_var_index, last_var_index)  = get_read_molecules(vcf_reader, &vcf_info, READ_TYPE::HIFI);
            let (_break, posteriors, log_likelihood) = expectation(&molecules, &cluster_centers);
            eprintln!("yay we did it right??? posterior {}, molecules {} log likelihoods {:?}", posterior, molecules.len(), log_likelihoods);
        }  ***/
    }
    if to_return.len() == 0 {
        if vcf_info.variant_positions.len() > 0 {
            to_return.push(PhaseBlock{
                id: to_return.len(),
                start_index: start_index,
                start_position: vcf_info.variant_positions[start_index],
                end_index: end_index+1,
                end_position: vcf_info.variant_positions[end_index]});
        }
    } else if to_return[to_return.len()-1].end_index != end_index {
        to_return.push(PhaseBlock {
            id: to_return.len(),
            start_index: phase_block_start,
            start_position: vcf_info.variant_positions[phase_block_start],
            end_index: end_index+1,
            end_position: vcf_info.variant_positions[end_index],
        });
    }
    to_return
}

fn swap_copied(cluster_centers_copy: &mut Vec<Vec<f32>>, breakpoint: usize, pairing: &Vec<(usize, usize)>) {
    let mut touched = [false;32];
    for (hap1, hap2) in pairing {
        if touched[*hap1] {continue;}
        touched[*hap1] = true;
        touched[*hap2] = true;
        let tmp = cluster_centers_copy[*hap1][breakpoint];
        cluster_centers_copy[*hap1][breakpoint] = cluster_centers_copy[*hap2][breakpoint];
        cluster_centers_copy[*hap2][breakpoint] = tmp;
    }
}

fn swap_full(cluster_centers: &Vec<Vec<f32>>, pairing: &Vec<(usize, usize)>) -> Vec<Vec<f32>> {
    let mut cluster_centers_copy = cluster_centers.clone();
    let mut touched = [false;32];
    for (hap1, hap2) in pairing {
        if touched[*hap1]  { continue; }
        touched[*hap1] = true;
        touched[*hap2] = true;
        for i in 0..cluster_centers[0].len() {
            let tmp = cluster_centers[*hap1][i];
            cluster_centers_copy[*hap1][i] = cluster_centers[*hap2][i];
            cluster_centers_copy[*hap2][i] = tmp;            
        }
    }
    cluster_centers_copy
}

fn swap(cluster_centers: &mut Vec<Vec<f32>>, breakpoint: usize, pairing: &Vec<(usize, usize)>, length: usize) {
    let mut touched = [false;32];
    for (hap1, hap2) in pairing {
        if touched[*hap1]  { continue; }
        touched[*hap1] = true;
        touched[*hap2] = true;
        for locus in (breakpoint+1)..(breakpoint+1+length) {
            if locus < cluster_centers[0].len() {
                let tmp = cluster_centers[*hap1][locus];
                cluster_centers[*hap1][locus] = cluster_centers[*hap2][locus];
                cluster_centers[*hap2][locus] = tmp;
            }
        }
    }
}


fn phase_phaseblocks(data: &ThreadData, cluster_centers: &mut Vec<Vec<f32>>, 
    phase_blocks: &Vec<PhaseBlock>, vcf_info: &VCF_info) -> HashMap<usize, usize> {
    let mut vcf_reader = bcf::IndexedReader::from_path(format!("{}", data.vcf_out.to_string()))
        .expect("could not open indexed vcf reader on output vcf");
    let chrom = vcf_reader
        .header()
        .name2rid(data.chrom.as_bytes())
        .expect("cant get chrom rid");
    let mut allele_phase_block_id: HashMap<usize, usize> = HashMap::new();
    
    match &data.hic_bam {
        Some(hic) => {}, // if has hic continue as normal
        None => {
            for (id, pb) in phase_blocks.iter().enumerate() {
                for allele_id in pb.start_index..(pb.end_index + 1) { 
                    allele_phase_block_id.insert(allele_id, pb.start_position); 
                }    
            }
            return allele_phase_block_id;
        }
    }
    let vcf_info: VCF_info = inspect_vcf(&mut vcf_reader, &data);
    let mut phase_block_ids: HashMap<usize,usize> = HashMap::new();
    println!("entering hic phasing with {} phase blocks", phase_blocks.len());
    for (id, phase_block) in phase_blocks.iter().enumerate() {
        //eprintln!("phase block {} from {}-{}",id, phase_block.start_index, phase_block.end_index);
        for i in phase_block.start_index..(phase_block.end_index+1) {
            //eprintln!("\tinserting phase block id {} for vardex {}", id, i);
            phase_block_ids.insert(i,id);
        }
    }
    vcf_reader
        .fetch(chrom, data.start as u64, Some(data.end as u64))
        .expect("some actual error");
    let (hic_reads, _, _) = get_read_molecules(&mut vcf_reader, &vcf_info, READ_TYPE::HIC);
    println!("{} hic reads hitting > 1 variant", hic_reads.len());

    // Okay what do I really want? Maybe a map from (pb1_id, pb2_id) to a map from (allele1, allele2) to counts array????
    // OR map pb1_id to map from pb2_id to map from (allele1, allele2) to counts array???
    // benefits of 1: fewer hashmaps
    // benefits of 2: can query all phaseblocks for a given phaseblock
    // okay, wrap this in a struct
    let mut phaseblock_allele_pair_counts: PhaseBlockAllelePairCounts = PhaseBlockAllelePairCounts::new();

    // Eventually I will have a max heap of posterior probability phasing and a pair of metaphaseblocks
    // I will pop the heap and merge these two phaseblocks if they have high enough probability
    // but then I need to remove from the heap any element containing either of these metaphaseblocks
    // or I can keep a set of ids that have already been merged?? and ignore them if the heap pops with 
    // one already in the merged set? 

    // let mut all_counts: HashMap<(usize, usize), HashMap<(u8,u8), usize>> = HashMap::new();
    // ok that is a map from (phase_block_id1, phase_block_id2) to a map from (pb1_hap, pb2_hap) to counts
    let mut allele_pair_counts: HashMap<(usize, usize), [u64; 4]> = HashMap::new(); // allele_index1, allele_index2 => [(alt,alt),(alt,ref),(ref,alt),(ref,ref)]
    // and a map from (allele_index1, allele_index2) to  counts array for [1/1, 1/0, 0/1, 0/0]
    for hic_read in hic_reads {
        for i in 0..hic_read.len() {
            for j in (i+1)..hic_read.len() {
                let allele1 = hic_read[i];
                let allele2 = hic_read[j];
                if !phase_block_ids.contains_key(&allele1.index) || !phase_block_ids.contains_key(&allele2.index) { continue; }
                let phase_block1 = phase_block_ids.get(&allele1.index).expect(&format!("if you are reading this, i screwed up, allele index {} not in a phase block", allele1.index));
                let phase_block2 = phase_block_ids.get(&allele2.index).expect("why didnt the previous one fail first?");
                if phase_block1 != phase_block2 {
                    let minpb = phase_block1.min(phase_block2);
                    let maxpb = phase_block1.max(phase_block2);
                    let pb_allele_pair_counts: &mut AllelePairCounts = phaseblock_allele_pair_counts.counts.entry((*minpb, *maxpb)).or_insert(AllelePairCounts::new());
                    if allele1.index < allele2.index {
                        let counts = allele_pair_counts.entry((allele1.index, allele2.index)).or_insert([0;4]);
                        let pbcounts: &mut [u64; 4] = pb_allele_pair_counts.counts.entry((allele1.index, allele2.index)).or_insert([0;4]);
                        if allele1.allele && allele2.allele { // alt and alt
                            counts[0] += 1;
                            pbcounts[0] += 1;
                        } else if allele1.allele && !allele2.allele { // alt and ref
                            counts[1] += 1;
                            pbcounts[1] += 1;
                        } else if !allele1.allele && allele2.allele { // ref and alt
                            counts[2] += 1;
                            pbcounts[2] += 1;
                        } else { // ref and ref
                            counts[3] += 1;
                            pbcounts[3] += 1;
                        }
                    } else {
                        let counts = allele_pair_counts.entry((allele2.index, allele1.index)).or_insert([0;4]);
                        let pbcounts: &mut [u64; 4] = pb_allele_pair_counts.counts.entry((allele2.index, allele1.index)).or_insert([0;4]);
                        if allele2.allele && allele1.allele { // alt and alt
                            counts[0] += 1;
                            pbcounts[0] += 1;
                        } else if allele2.allele && !allele1.allele { // alt and ref
                            counts[1] += 1;
                            pbcounts[1] += 1;
                        } else if !allele2.allele && allele1.allele { // ref and alt
                            counts[2] += 1;
                            pbcounts[2] += 1;
                        } else { // ref and ref
                            counts[3] += 1;
                            pbcounts[3] += 1;
                        }
                    }
                }
            }
        }
    }
    println!("DONE building allele pair counts (size {})! thread {} chrom {}", allele_pair_counts.len(), data.index, data.chrom);
    
    let mut phase_block_pair_phasing_log_likelihoods: HashMap<(usize, usize), HashMap<usize, f64>> = HashMap::new();
    // this is a map from (phase_block_id1, phase_block_id2) to a map from a pairing index to a log likelihood
    // okay now i have allele_pair_counts which will contribute log likelihoods to phaseblock pairs
    let all_possible_pairings = pairings(data.ploidy); // get all pairings
    let log_phasing_prior = (1.0/(all_possible_pairings.len() as f64)).ln();
    let error = 0.2; // TODO do not hard code
    // each pairing implies a multinomial distribution on each pair of alleles
    for ((allele1_index, allele2_index), counts) in allele_pair_counts.iter() {
        let mut total_counts: u64 = 0;
        for count in counts.iter() {
            total_counts += *count as u64;
        }
        let phase_block1 = phase_block_ids.get(&allele1_index).expect("if you are reading this, i screwed up");
        let phase_block2 = phase_block_ids.get(&allele2_index).expect("why didnt the previous one fail first?");
        let min = phase_block1.min(phase_block2);
        let max = phase_block1.max(phase_block2);
        let phase_block_log_likelihoods = phase_block_pair_phasing_log_likelihoods.entry((*min, *max)).or_insert(HashMap::new());
        for (pairing_index, haplotype_pairs) in all_possible_pairings.iter().enumerate() {
            let log_likelihood = phase_block_log_likelihoods.entry(pairing_index).or_insert(log_phasing_prior);
            let mut pair_probabilities: [f64;4] = [error; 4];
            let mut total = 0.0; // for normalization to sum to 1
            for (phase_block1_hap, phase_block2_hap) in haplotype_pairs.iter() {
                let phase_block1_allele_frac = cluster_centers[*phase_block1_hap][*allele1_index] as f64;
                let phase_block2_allele_frac = cluster_centers[*phase_block2_hap][*allele2_index] as f64;
                pair_probabilities[0] += phase_block1_allele_frac * phase_block2_allele_frac;
                pair_probabilities[1] += phase_block1_allele_frac * (1.0 - phase_block2_allele_frac);
                pair_probabilities[2] += (1.0 - phase_block1_allele_frac) * phase_block2_allele_frac;
                pair_probabilities[3] += (1.0 - phase_block1_allele_frac) * (1.0 - phase_block2_allele_frac);
            }
            for psuedocount in pair_probabilities.iter() { total += psuedocount; }
            for index in 0..4 { pair_probabilities[index] /= total; } // normalize to sum to 1.0
            let multinomial_distribution = Multinomial::new(&pair_probabilities, total_counts).unwrap();
            let ln_pmf = multinomial_distribution.ln_pmf(counts);
            *log_likelihood += ln_pmf;
            println!("allele {}, allele {}, counts {:?} pairing {} probs {:?} adding ln_pmf {} which is prob {} for a total log likelihood of {}",
                allele1_index, allele2_index, counts, pairing_index, pair_probabilities, ln_pmf, ln_pmf.exp(), log_likelihood);
        }
    }
    println!("DONE creating phase block pair log likelihoods (size {})! thread {} chrom {}", phase_block_pair_phasing_log_likelihoods.len(), data.index, data.chrom);
    

    // How do we want to go about this?
    // 1. Try to phase largest phase blocks together, merge if possible, then continue
    // 2. Or find phase blocks that are big enough to be phased and try to phase the ones that are next to one another, merge, then continue
    // 3. Or find all pairwise posterior probabilities of phasing... merge highest posterior probability phasing phaseblocks, then compute 
    // phasing of the newly merged metaphaseblock to other phaseblocks.. continuing to merge highest posterior probability phasing blocks 
    // until the highest posterior is below some threshold.
    // okay lets do #3
    let mut meta_phaseblocks: Vec<MetaPhaseBlock> = Vec::new();
    let mut meta_phaseblock_map: HashMap<usize, MetaPhaseBlock> = HashMap::new();
    let mut first_phaseblock_start = 0;
    for (index, pb) in phase_blocks.iter().enumerate() {
        // TODO changes that im not sure about
        //if index == 0 {
        //    first_phaseblock_start = pb.start_position;
        //}
        meta_phaseblocks.push(MetaPhaseBlock { id: pb.start_position, phase_blocks: vec![*pb] });// pb.id, phase_blocks: vec![*pb] });
        meta_phaseblock_map.insert(index, MetaPhaseBlock { id: pb.start_position, phase_blocks: vec![*pb]});//pb.id, phase_blocks: vec![*pb] });
        for allele_id in pb.start_index..(pb.end_index+1) { 
            allele_phase_block_id.insert(allele_id, pb.start_position);//first_phaseblock_start); 
        }
    }
    let mut pairwise_posteriors_heap = 
        pairwise_phasing(&allele_pair_counts, &phase_block_ids, cluster_centers, data);

    while pairwise_posteriors_heap.len() != 0 {
        let phase_block_pair_phasing = pairwise_posteriors_heap.pop().unwrap();
        println!("in pairwise_posteriors_heap loop, {}",pairwise_posteriors_heap.len());
        if phase_block_pair_phasing.log_posterior_phasing.exp() <= data.hic_phasing_posterior_threshold {
            break;
        }
    
        // check if either ids are already merged and if so, continue
        if !meta_phaseblock_map.contains_key(&phase_block_pair_phasing.id1) || 
            !meta_phaseblock_map.contains_key(&phase_block_pair_phasing.id2) {
            continue;
        }
        let meta_block1 = meta_phaseblock_map.get(&phase_block_pair_phasing.id1).unwrap();
        let meta_block2 = meta_phaseblock_map.get(&phase_block_pair_phasing.id2).unwrap();
        let mut pbmin1 = usize::MAX;
        let mut pbmax1 = 0;
        let mut pbmin2 = usize::MAX;
        let mut pbmax2 = 0;
        for pb in &meta_block1.phase_blocks {
            pbmin1 = pbmin1.min(pb.start_position);
            pbmax1 = pbmax1.max(pb.end_position);
        }
        for pb in &meta_block2.phase_blocks {
            pbmin2 = pbmin2.min(pb.start_position);
            pbmax2 = pbmax2.max(pb.end_position);
        }
        let pbsize1 = pbmax1 - pbmin1;
        let pbsize2 = pbmax2 - pbmin2;
        let span = pbmax1.max(pbmax2) - pbmin1.min(pbmin2);
        println!("merging phaseblocks {} and {} for total span of {} with probability {}", meta_block1.id, meta_block2.id, span, phase_block_pair_phasing.log_posterior_phasing.exp());
        if meta_block1.phase_blocks.len() == 1 && meta_block2.phase_blocks.len() == 1 {
            println!("\tcounts {:?}", phase_block_pair_phasing.counts);
        }
        // merge the phaseblocks
        let new_meta_phaseblock = merge_phaseblocks(meta_block1, 
            meta_block2, 
            &phase_block_pair_phasing.phase_pairing, 
            cluster_centers, &mut allele_phase_block_id);
        
        // remove both from meta_phaseblock_map
        meta_phaseblock_map.remove(&phase_block_pair_phasing.id1);
        meta_phaseblock_map.remove(&phase_block_pair_phasing.id2);
        
        // compute phasing posteriors of newly merged phaseblock to all existing phaseblocks
        let all_new_phasings: Vec<MetaBlockPairPhasing> = get_new_phasings(&meta_phaseblock_map, 
            &new_meta_phaseblock, &cluster_centers, &data,
            &allele_pair_counts);
        // and add the new pairing phasings to the heap
        for meta_block_phasing in all_new_phasings {
            pairwise_posteriors_heap.push(meta_block_phasing);
        }
        meta_phaseblock_map.insert(new_meta_phaseblock.id, new_meta_phaseblock);
        
    }
    
    
    
    println!("done with hic");
    return allele_phase_block_id;
}

fn get_new_phasings(meta_phaseblock_map: &HashMap<usize, MetaPhaseBlock>, 
    new_meta_phaseblock: &MetaPhaseBlock,  
    cluster_centers: &Vec<Vec<f32>>, 
    data: &ThreadData, allele_pair_counts: &HashMap<(usize, usize), [u64;4]>) -> Vec<MetaBlockPairPhasing> {
    let mut to_return: Vec<MetaBlockPairPhasing> = Vec::new();
    for (id, meta_phase_block) in meta_phaseblock_map {
        let (pairing, posterior) = phase_two_phaseblocks(data, allele_pair_counts, 
        cluster_centers, new_meta_phaseblock, meta_phase_block);
        let meta_block_pair_phasing = MetaBlockPairPhasing{
            log_posterior_phasing: posterior.ln(),
            id1: meta_phase_block.id,
            id2: new_meta_phaseblock.id,
            counts: vec![],
            phase_pairing: pairing,
        };
        to_return.push(meta_block_pair_phasing);
    }
    return to_return;
}

fn pairwise_phasing(allele_pair_counts: &HashMap<(usize, usize), [u64; 4]>, 
    phase_block_ids: &HashMap<usize, usize>,
    cluster_centers: &Vec<Vec<f32>>, data: &ThreadData) -> BinaryHeap<MetaBlockPairPhasing> {
    let mut phase_block_pair_phasing_log_likelihoods: HashMap<(usize, usize), Vec<f64>> = HashMap::new();
    // this is a map from (phase_block_id1, phase_block_id2) to a map from a pairing index to a log likelihood
    
    let mut pair_posterior_heap: BinaryHeap<MetaBlockPairPhasing> = BinaryHeap::new();
    let all_possible_pairings = pairings(data.ploidy); // get all pairings
    let log_phasing_prior = (1.0/(all_possible_pairings.len() as f64)).ln();
    let error = 0.2; // TODO do not hard code
    let mut phase_block_pair_counts: HashMap<(usize, usize), Vec<u64>> = HashMap::new();
    // map from pb1,pb2 to Vec<usize> which is the counts of 0,0 pairing then 0,1 pairing, then 1,0 then 1,1
    // notes to self and josh: allele_pair_counts is a map from (allele1, allele2) to an array of
    // counts for [(alt,alt), (alt,ref), (ref,alt), (ref,ref)]
    // each pairing implies a multinomial distribution on each pair of alleles
    for ((allele1_index, allele2_index), counts) in allele_pair_counts.iter() {
        let mut total_counts: u64 = 0;
        for count in counts.iter() {
            total_counts += *count as u64;
        }
        let phase_block1 = phase_block_ids.get(&allele1_index).expect("if you are reading this, i screwed up");
        let phase_block2 = phase_block_ids.get(&allele2_index).expect("why didnt the previous one fail first?");
        
        let min = phase_block1.min(phase_block2);
        let max = phase_block1.max(phase_block2);
        let phase_block_counts = phase_block_pair_counts.entry((*min, *max)).or_insert(vec![0; 2_usize.pow(data.ploidy as u32)]);
        
        if (cluster_centers[0][*allele1_index] > 0.9 && cluster_centers[0][*allele2_index] > 0.9 &&
            cluster_centers[1][*allele1_index] < 0.1 && cluster_centers[1][*allele2_index] < 0.1) || 
            (cluster_centers[0][*allele1_index] < 0.1 && cluster_centers[0][*allele2_index] < 0.1 &&
            cluster_centers[1][*allele1_index] > 0.9 && cluster_centers[1][*allele2_index] > 0.9) {
                // we are in cis
                phase_block_counts[0] += counts[0];
                phase_block_counts[1] += counts[1];
                phase_block_counts[2] += counts[2];
                phase_block_counts[3] += counts[3];
        } else if (cluster_centers[0][*allele1_index] > 0.9 && cluster_centers[0][*allele2_index] < 0.1 &&
            cluster_centers[1][*allele1_index] < 0.1 && cluster_centers[1][*allele2_index] > 0.9) || 
            (cluster_centers[0][*allele1_index] < 0.1 && cluster_centers[0][*allele2_index] > 0.9 &&
            cluster_centers[1][*allele1_index] > 0.9 && cluster_centers[1][*allele2_index] < 0.1) {
                // we are in trans
                phase_block_counts[0] += counts[1];
                phase_block_counts[3] += counts[2];
                phase_block_counts[1] += counts[0];
                phase_block_counts[2] += counts[3];
            }

        let phase_block_log_likelihoods = 
            phase_block_pair_phasing_log_likelihoods.entry((*min, *max)).or_insert(Vec::new());

        for (pairing_index, haplotype_pairs) in all_possible_pairings.iter().enumerate() {
            if phase_block_log_likelihoods.len() == pairing_index {
                phase_block_log_likelihoods.push(log_phasing_prior);
            }
            //
            let mut pair_probabilities: [f64;4] = [error; 4];
            let mut total = 0.0; // for normalization to sum to 1
            for (phase_block1_hap, phase_block2_hap) in haplotype_pairs.iter() {
                let phase_block1_allele_frac = cluster_centers[*phase_block1_hap][*allele1_index] as f64;
                let phase_block2_allele_frac = cluster_centers[*phase_block2_hap][*allele2_index] as f64;
                pair_probabilities[0] += phase_block1_allele_frac * phase_block2_allele_frac;
                pair_probabilities[1] += phase_block1_allele_frac * (1.0 - phase_block2_allele_frac);
                pair_probabilities[2] += (1.0 - phase_block1_allele_frac) * phase_block2_allele_frac;
                pair_probabilities[3] += (1.0 - phase_block1_allele_frac) * (1.0 - phase_block2_allele_frac);
            }
            for psuedocount in pair_probabilities.iter() { total += psuedocount; }
            for index in 0..4 { pair_probabilities[index] /= total; } // normalize to sum to 1.0
            let multinomial_distribution = Multinomial::new(&pair_probabilities, total_counts).unwrap();
            phase_block_log_likelihoods[pairing_index] += multinomial_distribution.ln_pmf(counts);
        }
    }
    for ((phaseblock1, phaseblock2), pairing_log_likelihoods) in phase_block_pair_phasing_log_likelihoods {
        
        let mut max_log_posterior: f64 = f64::MIN;
        let mut max_pairing = 0;
        let log_denominator: f64 = log_sum_exp64(&pairing_log_likelihoods);
        println!("pb {} and {} pairing_log_likelihoods {:?}, log denom {}", 
            phaseblock1, phaseblock2, pairing_log_likelihoods, log_denominator);

        for (i, log_likelihood) in pairing_log_likelihoods.iter().enumerate() {
            let log_post = *log_likelihood - log_denominator;
            if log_post > max_log_posterior {
                max_log_posterior = log_post;
                max_pairing = i;
            }
        }
        println!("\t log_post {} posterior {}", max_log_posterior, max_log_posterior.exp());
        let mut pairing: Vec<(usize, usize)> = Vec::new();
        for (h1, h2) in &all_possible_pairings[max_pairing] { 
            pairing.push((*h1, *h2)); 
        }
        let mut copyme: Vec<u64> = Vec::new();
        for count in phase_block_pair_counts.get(&(phaseblock1, phaseblock2)).unwrap_or(&vec![0;0]).iter() {
            copyme.push(*count);
        }
        if max_log_posterior.exp() < data.hic_phasing_posterior_threshold {
            continue;
        }
        pair_posterior_heap.push(MetaBlockPairPhasing{
            log_posterior_phasing: max_log_posterior,
            id1: phaseblock1,
            id2: phaseblock2,
            counts: copyme,
            phase_pairing: pairing,
        });
    }
    println!("pair posterior heap hic phasing size {}",pair_posterior_heap.len());

    return pair_posterior_heap;
}

struct MetaBlockPairPhasing {
    log_posterior_phasing: f64,
    id1: usize,
    id2: usize,
    counts: Vec<u64>,
    phase_pairing: Vec<(usize, usize)>,
}

impl Ord for MetaBlockPairPhasing {
    fn cmp(&self, other: &Self) -> Ordering {
        self.log_posterior_phasing.partial_cmp(&other.log_posterior_phasing).unwrap()
    }
}

impl PartialOrd for MetaBlockPairPhasing {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.log_posterior_phasing.partial_cmp(&other.log_posterior_phasing)
    }
}

impl PartialEq for MetaBlockPairPhasing {
    fn eq(&self, other: &Self) -> bool {
        self.log_posterior_phasing == other.log_posterior_phasing
    }
}

impl Eq for MetaBlockPairPhasing { 
}

struct PhaseBlockAllelePairCounts {
    counts: HashMap<(usize, usize), AllelePairCounts>, // phaseblock_id1, phaseblock_id2 -> allelepaircounts
}

impl PhaseBlockAllelePairCounts {
    fn new() -> PhaseBlockAllelePairCounts {
        PhaseBlockAllelePairCounts { counts: HashMap::new(), }
    }
}

struct AllelePairCounts {
    counts: HashMap<(usize, usize), [u64; 4]>, // allele1_id, allele2_id -> [alt/alt, alt/ref, ref/alt, ref/ref] counts
}

impl AllelePairCounts {
    fn new() -> AllelePairCounts {
        AllelePairCounts { counts: HashMap::new(), }
    }
}

fn merge_phaseblocks(phaseblock1: &MetaPhaseBlock, 
    phaseblock2: &MetaPhaseBlock, 
    best_pairing: &Vec<(usize, usize)>, 
    cluster_centers: &mut Vec<Vec<f32>>,  
    allele_phase_block_id: &mut HashMap<usize, usize>) -> MetaPhaseBlock {
    let mut to_return: MetaPhaseBlock = MetaPhaseBlock { 
        id: 0, phase_blocks: Vec::new() // set id later
    };
    let mut first_phaseblock_start = 0;
    for (index, pb) in phaseblock1.phase_blocks.iter().enumerate() {
        if index == 0 {
            first_phaseblock_start = pb.start_position;
            to_return.id = first_phaseblock_start;
        }
        to_return.phase_blocks.push(*pb);
        for allele_id in pb.start_index..(pb.end_index+1) {
            allele_phase_block_id.insert(allele_id, first_phaseblock_start);
        }
    }
    for pb in phaseblock2.phase_blocks.iter() {
        to_return.phase_blocks.push(*pb);
        for allele_id in pb.start_index..(pb.end_index+1) {
            allele_phase_block_id.insert(allele_id, to_return.id);
        }
    }
    let mut tmp_cluster_center: Vec<Vec<Vec<f32>>> = Vec::new(); // phaseblock, haplotype, index
    for (pb_index, pb) in phaseblock2.phase_blocks.iter().enumerate() {
        tmp_cluster_center.push(Vec::new());
        for haplotype in 0..cluster_centers.len() { 
            tmp_cluster_center[pb_index].push(Vec::new());
            for index in pb.start_index..(pb.end_index+1) {
                tmp_cluster_center[pb_index][haplotype].push(0.0);
            }
        }
    }
    
    
    for (pb1_hap, pb2_hap) in best_pairing {
        for (pb_index, pb) in phaseblock2.phase_blocks.iter().enumerate() {
            for index in pb.start_index..(pb.end_index +1) {
                tmp_cluster_center[pb_index][*pb1_hap][index - pb.start_index] = cluster_centers[*pb2_hap][index];
            }
        }
    }
    for (pb_index, pb) in phaseblock2.phase_blocks.iter().enumerate() {
        for haplotype in 0..cluster_centers.len() {
            for index in pb.start_index..(pb.end_index + 1) {
                cluster_centers[haplotype][index] = tmp_cluster_center[pb_index][haplotype][index - pb.start_index];
            } // TODO double check end index is inclusive
        }
    }
    return to_return;
}

fn phase_two_phaseblocks(data: &ThreadData, allele_pair_counts: &HashMap<(usize, usize), [u64;4]>, 
    cluster_centers: &Vec<Vec<f32>>, 
    phaseblock1: &MetaPhaseBlock, phaseblock2: &MetaPhaseBlock) -> (Vec<(usize, usize)>, f64) {
    let pairings = pairings(data.ploidy);
    let log_phasing_prior = (1.0/pairings.len() as f64).ln();
    let mut pairing_log_likelihoods: Vec<f64> = Vec::new();
    for _ in &pairings { pairing_log_likelihoods.push(log_phasing_prior); }
    let mut pairing_posteriors: Vec<f64> = Vec::new();
    let error = 0.2; //TODO dont hard code
    // get allele pairs
    for pb1 in &phaseblock1.phase_blocks {
        for allele1 in pb1.start_index..pb1.end_index {
            for pb2 in &phaseblock2.phase_blocks {
                for allele2 in pb2.start_index..pb2.end_index {
                    let min = allele1.min(allele2);
                    let max = allele1.max(allele2);
                    if let Some(counts) = allele_pair_counts.get(&(min, max)) { // is there a better way of doing this?
                        // maybe if i had allele_pair_counts for each phase_block?? then I could iterate instead of doing
                        // this n^2 business. Oh well, lets code it and see if its slow, optimize only if slow.
                        let mut total_counts = 0;
                        for count in counts.iter() { total_counts += *count as u64; }
                        for (pairing_index, haplotype_pairs) in pairings.iter().enumerate() {
                            let mut pair_probabilities: [f64;4] = [error; 4];
                            let mut total = 0.0;
                            for (pb1_hap, pb2_hap) in haplotype_pairs {
                                let pb1_allele_frac = cluster_centers[*pb1_hap][allele1] as f64;
                                let pb2_allele_frac = cluster_centers[*pb2_hap][allele2] as f64;
                                pair_probabilities[0] += pb1_allele_frac * pb2_allele_frac; // both alt
                                pair_probabilities[1] += pb1_allele_frac * (1.0 - pb2_allele_frac); // alt/ref
                                pair_probabilities[2] += (1.0 - pb1_allele_frac) * pb2_allele_frac; // ref/alt
                                pair_probabilities[3] += (1.0 - pb1_allele_frac) * (1.0 - pb2_allele_frac); // ref/ref
                            }
                            for prob in pair_probabilities.iter() {
                                total += prob;
                            }
                            for i in 0..pair_probabilities.len() {
                                pair_probabilities[i] /= total; // probabilities must sum to 1
                            }
                            let multinomial = Multinomial::new(&pair_probabilities, total_counts).unwrap();
                            pairing_log_likelihoods[pairing_index] += multinomial.ln_pmf(counts);
                        }
                    }
                }
            }
        }
    } 
    let mut max_posterior = 0.0;
    let mut max_pairing = 0;
    let log_denominator = log_sum_exp64(&pairing_log_likelihoods);
    for (i, log_likelihood) in pairing_log_likelihoods.iter().enumerate() {
        let post = (*log_likelihood - log_denominator).exp();
        pairing_posteriors.push(post);
        if post > max_posterior {
            max_posterior = post;
            max_pairing = i;
        }
    }
    let mut pairing: Vec<(usize, usize)> = Vec::new();
    for (h1,h2) in &pairings[max_pairing] { pairing.push((*h1, *h2)); }
    return (pairing, max_posterior);
}




fn pairings(size: usize) -> Vec<Vec<(usize, usize)>> {
    let pairings_so_far: Vec<(usize, usize)> = Vec::new();
    let left: Vec<usize> = (0..size).collect();
    let right: Vec<usize> = (0..size).collect();
    return generate_pairings(pairings_so_far, left, right);
}


// generates all bipartite graph pairings
fn generate_pairings(pairings_so_far: Vec<(usize, usize)>, left_remaining: Vec<usize>, right_remaining: Vec<usize>) -> Vec<Vec<(usize, usize)>> {
    assert!(left_remaining.len() == right_remaining.len());
    let mut to_return: Vec<Vec<(usize, usize)>> = Vec::new();
    if left_remaining.len() == 0 {
        to_return.push(pairings_so_far);
        return to_return;
    }
    let left = &left_remaining[0];
    for right in right_remaining.iter() {
        let mut so_far = pairings_so_far.to_vec();
        so_far.push((*left, *right));
        let mut new_left: Vec<usize> = Vec::new();
        for l in left_remaining.iter() {
            if l != left {
                new_left.push(*l);
            }
        }
        let mut new_right: Vec<usize> = Vec::new();
        for r in right_remaining.iter() {
            if r != right {
                new_right.push(*r);
            }
        }
        let new_pairings = generate_pairings(so_far, new_left, new_right);
        for pairing in new_pairings {
            to_return.push(pairing);
        }
    }
    return to_return;
}

fn _binomial_test(cis: usize, trans: usize) -> f64 {
    let min = cis.min(trans) as u64;
    let n = Binomial::new(0.5, (cis+trans) as u64).unwrap();
    let p_value = n.cdf(min) * 2.0;
    p_value
}

fn output_phased_vcf(
    data: &ThreadData,
    cluster_centers: Vec<Vec<f32>>,
    phase_blocks: Vec<PhaseBlock>,
    vcf_info: &VCF_info,
    molecule_support: &Vec<Vec<f32>>,
    allele_phase_block_id: &HashMap<usize, usize>,
) {
    let mut vcf_reader = bcf::IndexedReader::from_path(format!("{}", data.vcf_out.to_string()))
        .expect("could not open indexed vcf reader on output vcf");
    let header_view = vcf_reader.header();
    let mut new_header = bcf::header::Header::from_template(header_view);
    new_header.push_record(br#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="phase set id">"#);
    let cluster_center_string = format!(r#"##FORMAT=<ID=CC,Number={},Type=Float,Description="phasstphase cluster center values">"#, cluster_centers.len());
    let molecule_support_string = format!(r#"##FORMAT=<ID=MS,Number={},Type=Float,Description="weighted molecule support for this haplotype at this location">"#, cluster_centers.len());
    new_header.push_record(&cluster_center_string.as_bytes());
    new_header.push_record(&molecule_support_string.as_bytes());
    let mut vcf_writer = bcf::Writer::from_path(
        data.phased_vcf_out.to_string(),
        &new_header,
        false,
        Format::Vcf,
    ).expect("could not open vcf writer");
    let mut index: usize = 0;
    let mut index_to_phase_block: HashMap<usize, usize> = HashMap::new();
    for (id, phase_block) in phase_blocks.iter().enumerate() {
        for i in phase_block.start_index..(phase_block.end_index + 1) {
            index_to_phase_block.insert(i, id);
        }
    }
    for (i, r) in vcf_reader.records().enumerate() {
        let mut rec = r.expect("could not unwrap vcf record");
        let mut new_rec = vcf_writer.empty_record();
        copy_vcf_record(&mut new_rec, &rec);
        //println!("{}",index);
        match allele_phase_block_id.get(&index) {
            Some(id) => {
                new_rec.push_format_integer(b"PS", &[*id as i32]).expect("you did it again, pushing your problems down to future you");
            },
            None => {eprintln!("no PS for index {}",index); },
        }
        let mut cluster_center: Vec<f32> = Vec::new();
        for h in 0..cluster_centers.len() {
            cluster_center.push(cluster_centers[h][i]);
        }
        new_rec.push_format_float(b"CC", &cluster_center).expect("blame past self");
        let mut mol_support: Vec<f32> = Vec::new();
        for h in 0..cluster_centers.len() {
            mol_support.push(molecule_support[h][i]);
        }
        new_rec.push_format_float(b"MS", &mol_support).expect("failed push format");
        //let phase_block_id = index_to_phase_block.get(&index).expect("i had it coming");
        let genotypes = infer_genotype(&cluster_centers, index, &vcf_info);
        new_rec.push_genotypes(&genotypes)
            .expect("i did expect this error");
        vcf_writer.write(&new_rec).expect("could not write record");
        index += 1;
    }
}

fn infer_genotype(cluster_centers: &Vec<Vec<f32>>, index: usize, vcf_info: &VCF_info) -> Vec<GenotypeAllele> {
    let mut genotypes: Vec<GenotypeAllele> = Vec::new();
    let mut phased = true;
    let mut all_point5 = true;
    let lower = 0.15;
    let upper = 0.85;
    for haplotype in 0..cluster_centers.len() {
        if cluster_centers[haplotype][index] != 0.5 { all_point5 = false;}
        if cluster_centers[haplotype][index] > upper {
        } else if cluster_centers[haplotype][index] < lower {
        } else { phased = false; }
    }
    if all_point5 {
        // we need the genotypes
        let sample = &vcf_info.genotypes[index];
        for g in sample {
            genotypes.push(*g);
        }
        return genotypes;
    }
    if !phased {
        let sample = &vcf_info.genotypes[index];
        for g in sample {
            genotypes.push(*g);
        }
        return genotypes;
    }
    for haplotype in 0..cluster_centers.len() {
        if cluster_centers[haplotype][index] > upper {
                genotypes.push(GenotypeAllele::Phased(1));
        } else if cluster_centers[haplotype][index] < lower {
                genotypes.push(GenotypeAllele::Phased(0));
        }    
    }
    genotypes
}



enum READ_TYPE {
    HIFI, HIC,
}

fn get_read_molecules(vcf: &mut bcf::IndexedReader, vcf_info: &VCF_info, read_type: READ_TYPE) -> (Vec<Vec<Allele>>, usize, usize) {
    let mut molecules: HashMap<String, Vec<Allele>> = HashMap::new();
    let ref_tag;
    let alt_tag;
    match read_type {
        READ_TYPE::HIFI => {
            ref_tag = b"RM"; 
            alt_tag = b"AM";
        },
        READ_TYPE::HIC => {
            ref_tag = b"RH";
            alt_tag = b"AH";
        }
    }
    let mut last_var_index = 0;
    let mut first_var_index = usize::MAX;
    for rec in vcf.records() {
        let rec = rec.expect("couldnt unwrap record");
        let pos = rec.pos();
        let var_index = vcf_info
            .position_to_index
            .get(&(pos as usize))
            .expect("please don't do this to me");
        last_var_index = *var_index;
        if first_var_index == usize::MAX { first_var_index = *var_index; }

        match rec.format(ref_tag).string() { 
            Ok(rec_format) => {
                for rec in rec_format.iter() {
                    let ref_mols = std::str::from_utf8(rec).expect(":").to_string();
                    for read in ref_mols.split(";") {
                        let mol = molecules.entry(read.to_string()).or_insert(Vec::new());
                        mol.push(Allele {
                            index: *var_index,
                            allele: false,
                        });
                    }
                }
            }
            Err(_) => (),
        }
        match rec.format(alt_tag).string() {
            Ok(rec_format) => {
                for rec in rec_format.iter() {
                    let alt_mols = std::str::from_utf8(rec).expect(":").to_string();
                    for read in alt_mols.split(";") {
                        let mol = molecules.entry(read.to_string()).or_insert(Vec::new());
                        mol.push(Allele {
                            index: *var_index,
                            allele: true,
                        });
                    }
                }
            }
            Err(_) => (),
        }
    }
    let mut to_return: Vec<Vec<Allele>> = Vec::new();
    for (index, (_read_name, alleles)) in molecules.iter().enumerate() {
        if alleles.len() < 3 {
            continue;
        }
        let mut mol: Vec<Allele> = Vec::new();
        for allele in alleles {
            mol.push(*allele);
        }
        to_return.push(mol);
    }
    (to_return, first_var_index, last_var_index)
}

#[derive(Clone, Copy)]
struct Allele {
    index: usize,
    allele: bool, // alt is true, ref is false
}

struct VCF_info {
    num_variants: usize,
    final_position: i64,
    variant_positions: Vec<usize>,
    position_to_index: HashMap<usize, usize>,
    genotypes: Vec<Vec<GenotypeAllele>>,
}

fn inspect_vcf(vcf: &mut bcf::IndexedReader, data: &ThreadData) -> VCF_info {
    let chrom = vcf
        .header()
        .name2rid(data.chrom.as_bytes())
        .expect("cant get chrom rid");
    let mut genotypes: Vec<Vec<GenotypeAllele>> = Vec::new();
    vcf.fetch(chrom, 0, None).expect("could not fetch in vcf");
    let mut position_to_index: HashMap<usize, usize> = HashMap::new();
    let mut num = 0;
    let mut last_pos = 0;
    let mut variant_positions: Vec<usize> = Vec::new();
    for r in vcf.records() {
        let rec = r.expect("could not unwrap vcf record");
        last_pos = rec.pos();
        let gts = rec.genotypes().expect("couldnt unwrap genotypes");
        genotypes.push(gts.get(0).to_vec());
        variant_positions.push(rec.pos() as usize);
        position_to_index.insert(rec.pos() as usize, num);
        num += 1;
    }
    println!("vcf has {} variants for chrom {}", num, chrom);
    VCF_info {
        num_variants: num,
        final_position: last_pos,
        variant_positions: variant_positions,
        position_to_index: position_to_index,
        genotypes: genotypes,
    }
}

fn init_cluster_centers(num: usize, data: &ThreadData) -> (Vec<Vec<f32>>, Vec<Vec<f32>>) {
    let seed = [data.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let mut cluster_centers: Vec<Vec<f32>> = Vec::new();
    let mut molecule_supports: Vec<Vec<f32>> = Vec::new();
    for k in 0..data.ploidy {
        let mut cluster_center: Vec<f32> = Vec::new();
        let mut molecule_support: Vec<f32> = Vec::new();
        for v in 0..num {
            cluster_center.push(0.5);
            molecule_support.push(0.0);
        }
        if cluster_center.len() > 0 {
            cluster_center[0] = rng.gen::<f32>().min(0.98).max(0.02);
            cluster_centers.push(cluster_center);
            molecule_supports.push(molecule_support);
        }
    }
    (cluster_centers, molecule_supports)
}

fn get_all_variant_assignments(data: &ThreadData) -> Result<(), Error> {
    eprintln!("get all variant assignemnts");
    let mut long_read_bam_reader = match &data.long_read_bam {
        Some(x) => Some(bam::IndexedReader::from_path(x).expect("could not open bam, maybe no index?")),
        None => None,
    };
    let mut hic_bam_reader = match &data.hic_bam {
        Some(x) => Some(bam::IndexedReader::from_path(x).expect("could not open bam, maybe no index2?")),
        None => None,
    };
    let mut fasta = fasta::IndexedReader::from_file(&data.fasta).expect("cannot open fasta file");

    let mut vcf_reader = bcf::IndexedReader::from_path(data.vcf.to_string()).expect("could not load index for vcf... looking for .csi file");
    let header_view = vcf_reader.header();
    let mut new_header = bcf::header::Header::from_template(header_view);
    new_header.push_record(br#"##fileformat=VCFv4.2"#);
    new_header.push_record(br#"##FORMAT=<ID=AM,Number=1,Type=String,Description="alt molecules long reads">"#);
    new_header.push_record(br#"##FORMAT=<ID=RM,Number=1,Type=String,Description="ref molecules long reads">"#);
    new_header.push_record(br#"##FORMAT=<ID=AH,Number=1,Type=String,Description="alt molecules hic">"#);
    new_header.push_record(br#"##FORMAT=<ID=RH,Number=1,Type=String,Description="ref molecules hic">"#);

    {
        // creating my own scope to close later to close vcf writer
        eprintln!("writing vcf {} using {} as header template",data.vcf_out.to_string(), data.vcf.to_string());
        let mut vcf_writer =
            bcf::Writer::from_path(data.vcf_out.to_string(), &new_header, false, Format::Vcf).expect("cant open vcf writer");
        let chrom = vcf_reader.header().name2rid(data.chrom.as_bytes()).expect("could not read vcf header");
        println!("fetching {}:{}-{}",chrom, data.start, data.end);
        match vcf_reader.fetch(chrom, data.start as u64, Some(data.end as u64)) {
            Ok(_) => {
                let mut total = 0;
                let mut hets = 0;
                for (i, _rec) in vcf_reader.records().enumerate() {
                    total += 1;
                    let rec = _rec.expect("cant unwrap vcf record");
                    let pos = rec.pos();
                    let alleles = rec.alleles();
                    let mut new_rec = vcf_writer.empty_record();
                    copy_vcf_record(&mut new_rec, &rec);
                    if alleles.len() != 3 {
                        continue; // ignore multi allelic sites
                    }
                    let reference = std::str::from_utf8(alleles[0]).expect("this really shouldnt fail");
                    let alternative = std::str::from_utf8(alleles[1]). expect("this really shouldnt fail2");
                    let rec_chrom = rec.rid().expect("could not unwrap vcf record id");

                    let genotypes = rec.genotypes().expect("cant get genotypes");
                    let genotype = genotypes.get(0); // assume only 1 and get the first one
                    if is_heterozygous(genotype) {
                        hets += 1;
                        get_variant_assignments(
                            &data.chrom,
                            pos as usize,
                            reference.to_string(),
                            alternative.to_string(),
                            data.min_mapq,
                            data.min_base_qual,
                            &mut long_read_bam_reader,
                            &mut hic_bam_reader,
                            &mut fasta,
                            data.window,
                            data.start,
                            data.end,
                            &mut vcf_writer,
                            &mut new_rec,
                        );
                    }
                }
                println!(
                    "done, saw {} records of which {} were hets in chrom {}",
                    total, hets, data.chrom
                );
            },
            Err(_e) => (),
        }; // skip to chromosome for this thread
        
    } //creating my own scope to close vcf
    eprintln!("finished scope which closes the file");

    let result = Command::new("bcftools")
        .args(&["index", "-f", &data.vcf_out])
        .status()
        .expect("bcftools failed us");

    

    fs::File::create(data.vcf_out_done.to_string()).expect("cant create .done file. are the permissions wrong?");
    Ok(())
}

fn copy_vcf_record(new_rec: &mut bcf::record::Record, rec: &bcf::record::Record) {
    new_rec.set_rid(rec.rid());
    new_rec.set_pos(rec.pos());
    new_rec.set_id(&rec.id()).expect("set id failed");
    for filter in rec.filters() {
        new_rec.push_filter(&filter).expect("push filter failed");
    }
    let mut alleles = rec.alleles();
    alleles.pop();
    new_rec
        .set_alleles(&alleles)
        .expect("could not write alleles to new record???");
    new_rec.set_qual(rec.qual());
    let header = rec.header();
    for header_record in header.header_records() {
        match header_record {
            bcf::header::HeaderRecord::Filter { key, values } => {}
            bcf::header::HeaderRecord::Info { key, values } => {
                let mut format = FORMAT {
                    Id: "blah".to_string(),
                    Type: FORMAT_TYPE::Integer,
                };
                for (x, y) in values {
                    match x.as_str() {
                        "ID" => format.Id = y,
                        "Type" => match y.as_str() {
                            "Integer" => format.Type = FORMAT_TYPE::Integer,
                            "Float" => format.Type = FORMAT_TYPE::Float,
                            "String" => format.Type = FORMAT_TYPE::String,
                            "Char" => format.Type = FORMAT_TYPE::Char,
                            &_ => (),
                        },
                        &_ => (),
                    }
                }
                match format.Type {
                    FORMAT_TYPE::Integer => match rec.info(&format.Id.as_bytes()).integer() {
                        Ok(rec_format) => {
                            for thingy in rec_format.iter() {
                                new_rec
                                    .push_info_integer(&format.Id.as_bytes(), thingy)
                                    .expect("fail1");
                            }
                        }
                        Err(_) => (),
                    },
                    FORMAT_TYPE::Float => match rec.info(&format.Id.as_bytes()).float() {
                        Ok(rec_format) => {
                            for thingy in rec_format.iter() {
                                new_rec
                                    .push_info_float(&format.Id.as_bytes(), thingy)
                                    .expect("fail1");
                            }
                        }
                        Err(_) => (),
                    },
                    FORMAT_TYPE::String => match rec.info(&format.Id.as_bytes()).string() {
                        Ok(rec_format) => {
                            new_rec
                                .push_info_string(
                                    &format.Id.as_bytes(),
                                    &rec_format.expect("blerg"),
                                )
                                .expect("fail1");
                        }
                        Err(_) => (),
                    },
                    FORMAT_TYPE::Char => match rec.info(&format.Id.as_bytes()).string() {
                        Ok(rec_format) => {
                            new_rec
                                .push_info_string(
                                    &format.Id.as_bytes(),
                                    &rec_format.expect("blerg2"),
                                )
                                .expect("fail1");
                        }
                        Err(_) => (),
                    },
                }
            }
            bcf::header::HeaderRecord::Format { key, values } => {
                let mut format = FORMAT {
                    Id: "blah".to_string(),
                    Type: FORMAT_TYPE::Integer,
                };
                for (x, y) in values {
                    match x.as_str() {
                        "ID" => format.Id = y,
                        "Type" => match y.as_str() {
                            "Integer" => format.Type = FORMAT_TYPE::Integer,
                            "Float" => format.Type = FORMAT_TYPE::Float,
                            "String" => format.Type = FORMAT_TYPE::String,
                            "Char" => format.Type = FORMAT_TYPE::Char,
                            &_ => (),
                        },
                        &_ => (),
                    }
                }
                match format.Type {
                    FORMAT_TYPE::Integer => match rec.format(&format.Id.as_bytes()).integer() {
                        Ok(rec_format) => {
                            for thingy in rec_format.iter() {
                                new_rec
                                    .push_format_integer(&format.Id.as_bytes(), thingy)
                                    .expect("noooooooooo");
                            }
                        }
                        Err(_) => (),
                    },
                    FORMAT_TYPE::Float => match rec.format(&format.Id.as_bytes()).float() {
                        Ok(rec_format) => {
                            for thingy in rec_format.iter() {
                                new_rec
                                    .push_format_float(&format.Id.as_bytes(), thingy)
                                    .expect("fail1");
                            }
                        }
                        Err(_) => (),
                    },
                    FORMAT_TYPE::String => {
                        if format.Id == "GT".to_string() {
                            let sample_count = header.sample_count();
                            for i in 0..sample_count {
                                let gt = rec.genotypes().expect("lkjlkj").get(i as usize);
                                new_rec.push_genotypes(&gt).expect("nope");
                            }
                        } else {
                            match rec.format(&format.Id.as_bytes()).string() {
                                Ok(rec_format) => {
                                    new_rec
                                        .push_format_string(&format.Id.as_bytes(), &rec_format)
                                        .expect("fail1");
                                }
                                Err(_) => (),
                            }
                        }
                    }
                    FORMAT_TYPE::Char => match rec.format(&format.Id.as_bytes()).string() {
                        Ok(rec_format) => {
                            new_rec
                                .push_format_string(&format.Id.as_bytes(), &rec_format)
                                .expect("fail1");
                        }
                        Err(_) => (),
                    },
                }
            }
            bcf::header::HeaderRecord::Contig { key, values } => {}
            bcf::header::HeaderRecord::Structured { key, values } => {}
            bcf::header::HeaderRecord::Generic { key, value } => {}
        }
    }
}

struct FORMAT {
    Id: String,
    Type: FORMAT_TYPE,
}

#[derive(Debug)]
enum FORMAT_TYPE {
    Integer,
    Float,
    String,
    Char,
}

fn get_variant_assignments<'a>(
    chrom: &String,
    pos: usize,
    ref_allele: String,
    alt_allele: String,
    min_mapq: u8,
    min_base_qual: u8,
    long_reads: &mut Option<bam::IndexedReader>,
    hic_reads: &mut Option<bam::IndexedReader>,
    fasta: &mut fasta::IndexedReader<std::fs::File>,
    window: usize,
    start: usize,
    end: usize,
    vcf_writer: &mut bcf::Writer,
    vcf_record: &mut bcf::record::Record,
) {
    //println!("pos {} window {} start {}",pos, window, start);
    //if (pos + window) > end {
    if pos > end {
        return;
    }
    match long_reads {
        Some(bam) => {
            let (ref_string, alt_string) = get_read_assignments(chrom.to_string(), pos, bam, 
                fasta, &ref_allele, &alt_allele, min_base_qual, min_mapq, window);
            vcf_record
                .push_format_string(b"AM", &[alt_string.as_bytes()])
                .expect("blarg");
            vcf_record
                .push_format_string(b"RM", &[ref_string.as_bytes()])
                .expect("gggg");
        }
        None => (),
    }
    match hic_reads {
        Some(bam) => {
            let (ref_string, alt_string) = get_read_assignments(chrom.to_string(), pos, bam, 
                fasta, &ref_allele, &alt_allele, min_base_qual, min_mapq, window);
            vcf_record
                .push_format_string(b"AH", &[alt_string.as_bytes()])
                .expect("blarg2");
            vcf_record
                .push_format_string(b"RH", &[ref_string.as_bytes()])
                .expect("gggg2");
        }
        None => (),
    }
    vcf_writer.write(vcf_record).expect("nope");
}

fn get_read_assignments(
    chrom: String,
    pos: usize,
    bam: &mut bam::IndexedReader,
    fasta: &mut fasta::IndexedReader<std::fs::File>,
    ref_allele: &String,
    alt_allele: &String,
    min_base_qual: u8,
    min_mapq: u8,
    window: usize,
) -> (String, String) {
    let tid = bam
        .header()
        .tid(chrom.as_bytes())
        .expect(&format!("cannot find chrom tid {}", chrom));
    bam.fetch((tid, pos as u32, (pos + 1) as u32))
        .expect("blah"); // skip to region of bam of this variant position
    let ref_start = (pos.checked_sub(window)).unwrap_or(0) as u64;
    let ref_end = (pos.checked_add(window).unwrap_or(pos).checked_add(ref_allele.len())).unwrap_or(pos) as u64;
    let padding = 15;
    let ref_start_padded = ref_start.checked_sub(padding).unwrap_or(0);// TODO dont hard code
    let ref_end_padded = ref_end.checked_add(padding).unwrap_or(ref_end); // TODO dont hard code
    //eprintln!("padded {}-{}",ref_start_padded,ref_end_padded);
    
    fasta
        .fetch(&chrom, ref_start_padded, ref_end_padded)
        .expect("fasta fetch failed");
    let mut ref_sequence: Vec<u8> = Vec::new();
    fasta
        .read(&mut ref_sequence)
        .expect("failed to read fasta sequence");
    //if pos == 69505 { println!("yes, this one"); }
    let mut alt_sequence: Vec<u8> = Vec::new();
    for i in 0..(window + padding as usize) as usize {
        alt_sequence.push(ref_sequence[i]);
    }
    for base in alt_allele.as_bytes() {
        alt_sequence.push(*base);
    }
    for i in (window as usize + ref_allele.len() + padding as usize)..ref_sequence.len() {
        alt_sequence.push(ref_sequence[i]);
    }
    let mut read_names_ref: Vec<String> = Vec::new();
    let mut read_names_alt: Vec<String> = Vec::new();
    
    let mut ambiguous_count: f32 = 0.0;
    let mut total: f32 = 0.0;
    for _rec in bam.records() {
        let rec = _rec.expect("cannot read bam record");
        if rec.mapq() < min_mapq {
            continue;
        }
        if rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        let mut read_start: Option<usize> = None; // I am lazy and for some reason dont know how to do things, so this is my bad solution
        let mut read_end: usize = rec.seq_len();
        let mut min_bq = 93;
        let qual = rec.qual();
        for pos_pair in rec.aligned_pairs() {
            if (pos_pair[1] as u64) >= ref_start && (pos_pair[1] as u64) < ref_end {
                if pos_pair[1] as usize >= pos
                    && pos_pair[1] as usize <= pos + ref_allele.len().max(alt_allele.len())
                {
                    min_bq = min_bq.min(qual[pos_pair[0] as usize]);
                }
                read_end = pos_pair[0] as usize;
                if read_start == None {
                    read_start = Some(pos_pair[0] as usize); // assign this value only the first time in this loop that outter if statement is true
                                                             // getting the position in the read that corresponds to the reference region of pos-window
                }
            }
        }

        if min_bq < min_base_qual {
            continue;
        }
        if read_start == None {
            //println!(
            //    "what happened, read start {:?} read end {}",
            //    read_start, read_end
            //);
            continue;
        }
        let read_start = read_start.expect("why read start is none");
        let seq = rec.seq().as_bytes()[read_start..read_end].to_vec();
        let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
        let mut aligner = banded::Aligner::new(GAP_OPEN, GAP_EXTEND, score, K, W);
        let ref_alignment = aligner.semiglobal(&seq, &ref_sequence);
        let alt_alignment = aligner.semiglobal(&seq, &alt_sequence);
        total += 1.0;
        if ref_alignment.score > alt_alignment.score {
            read_names_ref.push(std::str::from_utf8(rec.qname()).expect("wtff").replace(":","_").replace(";","-").to_string());
        } else if alt_alignment.score > ref_alignment.score {
            read_names_alt.push(std::str::from_utf8(rec.qname()).expect("wtf").replace(":","_").replace(";","-").to_string());
        } else {
            ambiguous_count += 1.0;
        }
        
        //if pos == 69505 { println!("scores {} and {}",ref_alignment.score, alt_alignment.score);
        //println!("read_names_ref {}, read_names_alt {}", read_names_ref.len(), read_names_alt.len()); }
    }
    if ambiguous_count / total > 0.15 { // TODO dont hard code stuff
        read_names_ref = Vec::new();
        read_names_alt = Vec::new();
    }
    let concat_ref = read_names_ref.join(";");
    let concat_alt = read_names_alt.join(";");
    (concat_ref, concat_alt)
}

fn pretty_print(aln: &bio::alignment::Alignment, x: &Vec<u8>, y: &Vec<u8>) {
    let mut s1: Vec<u8> = Vec::new();
    let mut s2: Vec<u8> = Vec::new();
    let mut s3: Vec<u8> = Vec::new();
    let mut x_index: usize = 0;
    let mut y_index: usize = 0;
    if aln.xstart > aln.ystart {
        for i in 0..(aln.xstart-aln.ystart) {
            s2.push(b' ');
            s1.push(x[x_index]);
            x_index += 1;
            s3.push(b' ');
        }
    } else if aln.ystart > aln.xstart {
        for i in 0..(aln.ystart-aln.xstart) {
            s2.push(b' ');
            s3.push(y[y_index]);
            y_index += 1;
            s1.push(b' ');
        }
    }

    for op in &aln.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                s1.push(x[x_index]);
                x_index += 1;
                s2.push(b'|');
                s3.push(y[y_index]);
                y_index += 1;
            },
            bio::alignment::AlignmentOperation::Subst => {
                s1.push(x[x_index]);
                x_index += 1;
                s3.push(y[y_index]);
                y_index += 1;
                s2.push(b'x');
            },
            bio::alignment::AlignmentOperation::Del => {
                s1.push(b'-');
                s2.push(b'-');
                s3.push(y[y_index]);
                y_index += 1;
            },
            bio::alignment::AlignmentOperation::Ins => {
                s1.push(x[x_index]);
                x_index += 1;
                s2.push(b'-');
                s3.push(b'-');
            },
            _ => {},

        }
    }
    if x_index < x.len() {
        for xi in x_index..x.len() {
            s1.push(x[xi]);
        }
    }
    if y_index < y.len() {
        for yi in y_index..y.len() {
            s3.push(y[yi]);
        }
    }
    println!("{}",str::from_utf8(&s1).expect("couldnt unwrap"));
    println!("{}",str::from_utf8(&s2).expect("couldnt unwrap2"));
    println!("{}",str::from_utf8(&s3).expect("couldnt unwrap3"));
}

fn is_heterozygous(gt: bcf::record::Genotype) -> bool {
    if gt[0] == bcf::record::GenotypeAllele::Unphased(0)
        && gt[1] == bcf::record::GenotypeAllele::Unphased(1)
    {
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(1)
        && gt[1] == bcf::record::GenotypeAllele::Unphased(0)
    {
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(0)
        && gt[1] == bcf::record::GenotypeAllele::Phased(1)
    {
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(1)
        && gt[1] == bcf::record::GenotypeAllele::Phased(0)
    {
        return true;
    }
    return false;
}

struct ThreadData {
    index: usize,
    long_read_bam: Option<String>,
    hic_bam: Option<String>,
    fasta: String,
    chrom: String,
    start: usize,
    end: usize,
    vcf: String,
    min_mapq: u8,
    min_base_qual: u8,
    window: usize,
    output: String,
    vcf_out: String,
    phased_vcf_out: String,
    vcf_out_done: String,
    phased_vcf_done: String,
    phasing_window: usize,
    seed: u8,
    ploidy: usize,
    hic_phasing_posterior_threshold: f64,
    long_switch_threshold: f32,
}

#[derive(Clone)]
struct Params {
    hic_bam: Option<String>,
    long_read_bam: Option<String>,
    output: String,
    seed: u8,
    threads: usize,
    fasta: String,
    ploidy: usize,
    vcf: String,
    min_mapq: u8,
    min_base_qual: u8,
    window: usize,
    phasing_window: usize,
    hic_phasing_posterior_threshold: f64,
    long_switch_threshold: f32,
    chrom: Option<String>,
    start: Option<usize>,
    end: Option<usize>,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();

    let output = params.value_of("output").unwrap();
    let hic_bam = match params.value_of("hic_bam") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    let long_read_bam = match params.value_of("long_read_bam") {
        Some(x) => Some(x.to_string()),
        None => None,
    };

    let threads = params.value_of("threads").unwrap();
    let threads = threads.to_string().parse::<usize>().unwrap();

    let seed = params.value_of("seed").unwrap();
    let seed = seed.to_string().parse::<u8>().unwrap();

    let min_mapq = params.value_of("min_mapq").unwrap();
    let min_mapq = min_mapq.to_string().parse::<u8>().unwrap();

    let min_base_qual = params.value_of("min_base_qual").unwrap();
    let min_base_qual = min_base_qual.to_string().parse::<u8>().unwrap();

    let hic_phasing_posterior_threshold = params.value_of("hic_phasing_posterior_threshold").unwrap();
    let hic_phasing_posterior_threshold = hic_phasing_posterior_threshold.to_string().parse::<f64>().unwrap();


    let fasta = params.value_of("fasta").unwrap();

    let ploidy = params.value_of("ploidy").unwrap();
    let ploidy = ploidy.to_string().parse::<usize>().unwrap();
    assert!(ploidy != 1, "what are you doing trying to phase a haploid genome, ploidy can't be 1 here");

    let allele_alignment_window = params.value_of("allele_alignment_window").unwrap();
    let allele_alignment_window = allele_alignment_window
        .to_string()
        .parse::<usize>()
        .unwrap();

    let phasing_window = params.value_of("phasing_window").unwrap();
    let phasing_window = phasing_window.to_string().parse::<usize>().unwrap();

    let long_switch_threshold = params.value_of("long_switch_threshold").unwrap();
    let long_switch_threshold = long_switch_threshold.to_string().parse::<f32>().unwrap();

    let vcf = params.value_of("vcf").unwrap();
    let region = params.value_of("region");
    let chrom;
    let start;
    let end;
    if let Some(region) = region {
        let toks = region.split(":").collect::<Vec<&str>>();
        chrom = Some(toks[0].to_string());
        let tmp = toks[1].to_string();
        let toks2 = tmp.split("-").collect::<Vec<&str>>();
        let tmp = toks2[0].to_string();
        start = Some(tmp.parse::<usize>().expect("couldn't parse region"));
        let tmp = toks2[1].to_string();
        end = Some(tmp.parse::<usize>().expect("couldn't parse region"));
    } else {
        chrom = None;
        start = None;
        end = None;
    }



    Params {
        output: output.to_string(),
        long_read_bam: long_read_bam,
        hic_bam: hic_bam,
        threads: threads,
        seed: seed,
        fasta: fasta.to_string(),
        ploidy: ploidy,
        vcf: vcf.to_string(),
        min_mapq: min_mapq,
        min_base_qual: min_base_qual,
        window: allele_alignment_window,
        phasing_window: phasing_window,
        hic_phasing_posterior_threshold: hic_phasing_posterior_threshold,
        long_switch_threshold: long_switch_threshold,
        chrom: chrom,
        start: start,
        end: end,
    }
}

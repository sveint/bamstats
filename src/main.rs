use std::collections::BTreeMap;

extern crate rust_htslib;
extern crate bio;

use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::io;
use std::env;
use std::string;
use bio::io::fasta;

const FASTA_SEQS: &'static [ &'static str ] =  &[
    "1 dna:chromosome chromosome:GRCh37:1:1:249250621:1",
    "2 dna:chromosome chromosome:GRCh37:2:1:243199373:1",
    "3 dna:chromosome chromosome:GRCh37:3:1:198022430:1",
    "4 dna:chromosome chromosome:GRCh37:4:1:191154276:1",
    "5 dna:chromosome chromosome:GRCh37:5:1:180915260:1",
    "6 dna:chromosome chromosome:GRCh37:6:1:171115067:1",
    "7 dna:chromosome chromosome:GRCh37:7:1:159138663:1",
    "8 dna:chromosome chromosome:GRCh37:8:1:146364022:1",
    "9 dna:chromosome chromosome:GRCh37:9:1:141213431:1",
    "10 dna:chromosome chromosome:GRCh37:10:1:135534747:1",
    "11 dna:chromosome chromosome:GRCh37:11:1:135006516:1",
    "12 dna:chromosome chromosome:GRCh37:12:1:133851895:1",
    "13 dna:chromosome chromosome:GRCh37:13:1:115169878:1",
    "14 dna:chromosome chromosome:GRCh37:14:1:107349540:1",
    "15 dna:chromosome chromosome:GRCh37:15:1:102531392:1",
    "16 dna:chromosome chromosome:GRCh37:16:1:90354753:1",
    "17 dna:chromosome chromosome:GRCh37:17:1:81195210:1",
    "18 dna:chromosome chromosome:GRCh37:18:1:78077248:1",
    "19 dna:chromosome chromosome:GRCh37:19:1:59128983:1",
    "20 dna:chromosome chromosome:GRCh37:20:1:63025520:1",
    "21 dna:chromosome chromosome:GRCh37:21:1:48129895:1",
    "22 dna:chromosome chromosome:GRCh37:22:1:51304566:1",
    "X dna:chromosome chromosome:GRCh37:X:1:155270560:1",
    "Y dna:chromosome chromosome:GRCh37:Y:2649521:59034049:1"
];

fn main() {
    // Load BAM file from 1st argument
    let mut bamfile: String = string::String::new();
    if let Some(arg1) = env::args().nth(1) {
        bamfile = arg1;
    }
    let bam = bam::Reader::new(&bamfile).ok().expect("Error opening bam.");

    let mut ratio_counter: BTreeMap<u32, u32> = BTreeMap::new();

    let mut idx: u64 = 0;
    let mut processed_chr = 9999;

    // pileup over all covered sites
    for p in bam.pileup() {
        let pileup = p.ok().expect("Error reading BAM file.");

        if pileup.depth() < 20 {
            continue;
        }
        //println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());

        let mut alignment_counter: BTreeMap<u8, u32> = BTreeMap::new();
        for alignment in pileup.alignments() {
            let record = alignment.record();
            if record.is_unmapped() || record.is_duplicate() {
                continue;
            }
            let seq = record.seq();
            let base = seq[alignment.qpos()];

            match alignment_counter.get_mut(&base) {
                Some(value) => { *value += 1; continue; }
                None => {}
            }
            alignment_counter.insert(base, 1);
        }

        let mut sum: u32 = 0;
        let mut max: u32 = 0;
        for (base, &number) in alignment_counter.iter() {
            sum += number;
            if number > max {
                max = number;
            }
        }
        if sum > 0 {
            idx += 1;
            if processed_chr != pileup.tid() {
                processed_chr = pileup.tid();
                println!("Processing chromosome: {}. Processed positions so far: {}", processed_chr + 1, idx);
            }
            // Count minority as pct of total
            let ratio = (sum - max) * 100 / sum;
            match ratio_counter.get_mut(&ratio) {
                Some(value) => { *value += 1; continue; }
                None => {}
            }
            ratio_counter.insert(ratio, 1);
        }

    }

    for (ratio, &number) in ratio_counter.iter() {
        println!("Ratio {}: {}", ratio, number);
    }
}

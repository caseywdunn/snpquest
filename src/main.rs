use clap::Parser;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufRead;
use std::io::{Read, Write};
use std::path::Path;
use std::path::PathBuf;
use std::vec;

type Kmer = u64;
type MapType = rustc_hash::FxHashMap<Kmer, u64>;

// A nested datatype to hold count of snp variation:
// - snp site, kmer with the last two bits masked to 0.
// - a vector of length n where n is the number of samples
// - an array of length 4, corresponding to A, C, G, T variants at the snp
// - a value that is the count of each variant

/// index `00` is count of `A`
/// index `01` is count of `C`
/// index `10` is count of `G`
/// index `11` is count of `T`
type VariantCount = [u16; 4];

// Will be length of samples
type SampleVariantCounts = Vec<VariantCount>;

// Key is snp site, kmer with the last two bits masked to 0
type Variants = FxHashMap<Kmer, SampleVariantCounts>;

struct Sample {
    name: String,
    path: String,
    k: usize,
    kmers: Vec<Kmer>,
}

#[derive(serde::Serialize, Deserialize)]
struct SnpSet {
    sample_names: Vec<String>,
    k: usize,
    min_count: usize,
    ploidy: usize,
    snps: Vec<Kmer>,
}

/// A collection of kmer counting and analysis tools
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// minimum k-mer count
    #[arg(long, default_value_t = 2)]
    min_count: usize,

    /// ploidy, used to remove snps whose number of variants are more than expected
    /// in one or more samples
    #[arg(long, default_value_t = 2)]
    ploidy: usize,

    /// discard fraction of reads with highest k-mer counts
    #[arg(long, default_value_t = 0.05)]
    discard_fraction: f64,

    /// Name of the sample, could be species and sample ID eg Nanomia-bijuga-YPMIZ035039.
    /// Will be used as the prefix for output files
    #[arg(short, long, default_value_t = String::from("sample") )]
    sample: String,

    /// Name of file with previously identified snp sites. Will discover and write
    /// the snps to a file if none is provided
    #[arg(long, default_value_t = String::from("") )]
    snp_file: String,

    /// kmer prefix for thinning kmers
    #[arg(short, long, default_value_t = String::from("") )]
    prefix: String,

    /// Directory for output files, defaults to current directory
    #[arg(short, long, default_value_t = String::from("./") )]
    outdir: String,

    /// Input files. If no files are specified, data will be read
    /// from stdin. These are from jellyfish dump -c, eg Physalia_111156.kmers .
    /// Sample names are parsed from the file names
    #[arg()]
    input: Vec<String>,

    /// Verbosity
    #[arg(long, default_value_t = 0)]
    verbosity: usize,
}

/// Convert a kmer string to a u64
/// * `00` represents `A`
/// * `01` represents `C`
/// * `10` represents `G`
/// * `11` represents `T`
///
fn string_to_kmer(s: &str) -> Kmer {
    let mut kmer: Kmer = 0;
    for c in s.chars() {
        kmer <<= 2;
        match c {
            'A' => kmer |= 0,
            'C' => kmer |= 1,
            'G' => kmer |= 2,
            'T' => kmer |= 3,
            _ => panic!("Invalid character in kmer string: {}", c),
        }
    }
    kmer
}

fn discover_snp_sites(samples: &[Sample], k: usize, ploidy: usize) -> SnpSet {
    println!("Finding snps... ");
    // Create a Kmer bitmask for all but the last two bits
    let mut kmer_mask: Kmer = 0;
    for _ in 0..(k - 1) {
        kmer_mask |= 0b11;
        kmer_mask <<= 2;
    }

    let mut variants: Variants = Variants::default();

    for i in 0..samples.len() {
        let sample = &samples[i];
        for kmer in sample.kmers.iter() {
            // Mask the last two bits of the kmer
            let snp_site = kmer & kmer_mask;
            // Get the last two bits of the kmer
            let variant = kmer & 0b11;
            // If the snp site is not in the map, add it
            variants
                .entry(snp_site)
                .or_insert_with(|| vec![VariantCount::default(); samples.len()]);
            // Increment the count for the variant
            let sample_variants = variants.get_mut(&snp_site).unwrap();
            sample_variants[i][variant as usize] += 1;
        }
    }

    println!(
        "  Number of unique kmers ingested across all samples: {}",
        variants.len()
    );

    let mut snps: Vec<Kmer> = vec![];

    // Loop through the variants. If there are more than one across all samples and no more than
    // the ploidy in any one sample, identify the most common variant at the site and add
    // it to snps

    'kmer: for (snp_site, sample_variants) in variants.iter() {
        let mut variants_all_samples: VariantCount = [0; 4];
        // Count the number of variants at the site
        for sample in sample_variants.iter() {
            // if more than ploidy entries are greater than 0, break
            if sample.iter().filter(|&&x| x > 0).count() > ploidy {
                continue 'kmer;
            }
            // Add the sample's variants to the total
            for i in 0..4 {
                variants_all_samples[i] += sample[i];
            }
        }

        if variants_all_samples.iter().filter(|&&x| x > 0).count() > 1 {
            // Find the most common variant at the site
            // Keep the smallest variant if there is a tie
            let mut max_variant: usize = 0;
            for i in 1..4 {
                if variants_all_samples[i] > variants_all_samples[max_variant] {
                    max_variant = i;
                }
            }
            // Add the snp site to the list of snps
            snps.push(snp_site | max_variant as Kmer);
        }
    }

    snps.sort();

    SnpSet {
        sample_names: samples.iter().map(|s| s.name.clone()).collect(),
        k,
        min_count: 2,
        ploidy,
        snps,
    }
}

fn main() {
    let start_run = std::time::Instant::now();

    // Ingest command line arguments
    let args = Args::parse();
    let outdir = if args.outdir.is_empty() {
        "."
    } else {
        &args.outdir
    };

    // Print the program name and version
    println!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    // Print the arguments
    println!("{:?}", args);

    // Parse the outdir path and sample, create directories if necessary
    let path = PathBuf::from(&args.outdir);

    // Create the output directory if it does not exist
    let directory = format!("{}/", path.to_str().unwrap());
    std::fs::create_dir_all(&directory).unwrap();

    let mut k: usize = 0;

    let start = std::time::Instant::now();
    println!("Ingesting kmers...");
    std::io::stdout().flush().unwrap();

    // Create a vector of samples
    let mut samples: Vec<Sample> = vec![];

    for file_name in args.input.iter() {
        println!(" Reading {}...", file_name);
        // get the sample name from the file name by removing the path and extension, if any
        let sample_name = Path::new(file_name)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        let mut line_n: u64 = 0;
        // Open the file for buffered reading
        let file_path = Path::new(&file_name);
        let file = std::fs::File::open(file_path).unwrap();
        let reader = std::io::BufReader::new(file);

        let mut kmer_counts = MapType::default();

        // Iterate over the lines of the file
        for line in reader.lines() {
            line_n += 1;
            let line = line.unwrap();
            let mut parts = line.split_whitespace();
            let kmer_string = parts.next().unwrap();
            let count: u64 = parts.next().unwrap().parse().unwrap();

            // continue if the kmer count is less than the minimum count
            if count < args.min_count as u64 {
                continue;
            }

            // continue if the kmer does not start with the prefix
            if !args.prefix.is_empty() && !kmer_string.starts_with(&args.prefix) {
                continue;
            }

            let kmer = string_to_kmer(kmer_string);

            kmer_counts.insert(kmer, count);
            if k == 0 {
                k = kmer_string.len();
            } else if k != kmer_string.len() {
                panic!("kmer length mismatch: {} != {}", k, kmer_string.len());
            }
        }
        println!("  Number of processed kmers: {}", line_n);
        println!("  Number of ingested kmers: {}", kmer_counts.len());
        println!(
            "  Total count of ingested kmers: {}",
            kmer_counts.values().sum::<u64>()
        );
        println!("  kmer length: {}", k);

        // remove the most common kmers
        let mut counts: Vec<_> = kmer_counts.iter().collect();
        counts.sort_by(|a, b| b.1.cmp(a.1));
        let n_to_discard = (args.discard_fraction * counts.len() as f64).round() as usize;
        println!("  Discarding {} kmers with highest counts", n_to_discard);
        let keys_to_remove: Vec<Kmer> = counts
            .iter()
            .take(n_to_discard)
            .map(|(kmer, _)| **kmer)
            .collect();
        for kmer in keys_to_remove {
            kmer_counts.remove(&kmer);
        }
        println!("  Number of kmers after discarding: {}", kmer_counts.len());

        let mut kmers: Vec<Kmer> = kmer_counts.keys().copied().collect();
        kmers.sort();

        samples.push(Sample {
            name: sample_name.to_string(),
            path: file_name.to_string(),
            k,
            kmers,
        });
    }

    // Verify that all samples have the same kmer length
    let k = samples[0].k;
    for sample in samples.iter() {
        if sample.k != k {
            panic!("kmer length mismatch across samples: {} != {}", k, sample.k);
        }
    }

    println!("Ingesting samples done, time: {:?}", start.elapsed());

    let mut snp_set = SnpSet {
        sample_names: vec![],
        k: 0,
        min_count: 0,
        ploidy: 0,
        snps: vec![],
    };

    // If a snp file is provided, read it
    if !args.snp_file.is_empty() {
        println!("Reading snp file: {}", args.snp_file);
        let mut file = File::open(args.snp_file).expect("Failed to open file");
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer).expect("Failed to read file");

        snp_set = bincode::deserialize(&buffer).expect("Failed to deserialize");
    } else {
        // Discover the snp sites
        let start = std::time::Instant::now();
        println!("Discovering snp sites...");
        snp_set = discover_snp_sites(&samples, k, args.ploidy);
        println!("Discovering snp sites done, time: {:?}", start.elapsed());

        // Write the snp set to a file
        print!("Writing snp file: {}.snps... ", args.sample);
        let mut file_path = PathBuf::from(&outdir);
        file_path.push(format!("{}.snps", args.sample));
        let mut file = File::create(file_path).expect("Failed to create file");
        let encoded: Vec<u8> = bincode::serialize(&snp_set).expect("Failed to serialize");
        file.write_all(&encoded).expect("Failed to write file");
        println!("done");
    }

    println!("Number of snps: {}", snp_set.snps.len());

    println!("Total run time: {:?}", start_run.elapsed());
}

#[cfg(test)]
mod tests {
    use super::*;

    // These functions are used in the tests
    fn create_test_samples() -> Vec<Sample> {
        let sample1 = Sample {
            name: "sample1".to_string(),
            path: "/path/to/sample1".to_string(),
            k: 3,
            kmers: vec![
                // Two variants at one site
                // One shared with sample 2
                0b001111, 0b001101, // Singleton, should be excluded
                0b101110, // Exceed ploidy, should be excluded
                0b111100, 0b111101, 0b111110,
            ],
        };

        let sample2 = Sample {
            name: "sample2".to_string(),
            path: "/path/to/sample2".to_string(),
            k: 3,
            kmers: vec![
                // Two variants at one site, unique to sample2
                // Smaller should be retained
                0b001001, 0b001011, // Should be included due to variation in sample1
                0b001111, // Should be excluded due to ploidy in sample1
                0b111111,
            ],
        };

        vec![sample1, sample2]
    }

    // Tests
    #[test]
    fn test_test() {
        // Testing using the node indices from the HashMap
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn test_discover_snp_sites() {
        let samples = create_test_samples();
        let k = 3;
        let ploidy = 2;
        let snp_set = discover_snp_sites(&samples, k, ploidy);
        let snps = snp_set.snps;
        assert_eq!(snps.len(), 2);

        // Should be in order
        assert_eq!(snps[0], 0b001001);
        assert_eq!(snps[1], 0b001111);
    }
}

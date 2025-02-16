use std::fs::File;
use std::io::BufRead;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::path::Path;
use std::vec;

use clap::Parser;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Deserialize;

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

// Bit packed encoding of one or more nucleotides at a given site
// It is a serialized base plane
// 0000 No variant observed
// 0001 A
// 0010 C
// 0100 G
// 1000 T
// These can be summed, so for example 0011 is A and C observed
type Locus = u8;

type Genotype = Vec<Locus>;

#[derive(Clone)]
struct Sample {
    name: String,
    path: String,
    k: usize,
    kmers: Vec<Kmer>,
    genotype: Genotype,
}

#[derive(serde::Serialize, Deserialize, Clone)]
struct SnpSet {
    sample_names: Vec<String>,
    k: usize,
    min_count: usize,
    ploidy: usize,
    prefix: String,
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

    /// sample frequency minimum. The snp site must be present in at least this fraction of samples
    #[arg(long, default_value_t = 1.0)]
    freq_min: f64,

    /// Name of the run, could be clade eg Physalia.
    /// Will be used as the prefix for output files
    #[arg(short, long, default_value_t = String::from("run") )]
    run_name: String,

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

fn discover_snp_sites(
    samples: &[Sample],
    k: usize,
    ploidy: usize,
    freq_min: f64,
    prefix: String,
) -> SnpSet {
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

    // Loop through the variants. Criteria for inclusion are:
    // - More than one across all samples
    // - No more than the ploidy in any one sample
    // - At least one variant in at least `freq_min` fraction of samples
    // If met, identify the most common variant at the site and add it to snps

    let n_samples_min = (freq_min * samples.len() as f64).round() as usize;
    println!(
        "  Minimum number of samples for a site to be considered a snp: {}",
        n_samples_min
    );

    let mut n_rejected_exceed_ploidy = 0;
    let mut n_rejected_below_freq = 0;

    'kmer: for (snp_site, sample_variants) in variants.iter() {
        let mut variants_all_samples: VariantCount = [0; 4];
        // Count the number of variants at the site
        let mut n_samples_with_site = 0;
        for sample in sample_variants.iter() {
            if sample.iter().copied().sum::<u16>() > 0 {
                n_samples_with_site += 1;
            }
            // if more than ploidy entries are greater than 0, break
            if sample.iter().filter(|&&x| x > 0).count() > ploidy {
                n_rejected_exceed_ploidy += 1;
                continue 'kmer;
            }
            // Add the sample's variants to the total
            for i in 0..4 {
                variants_all_samples[i] += sample[i];
            }
        }

        if n_samples_with_site < n_samples_min {
            n_rejected_below_freq += 1;
            continue;
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
            snps.push(snp_site | max_variant as Kmer);
        }
    }

    println!("  Number of snps: {}", snps.len());
    println!(
        "  Number of snps rejected due to exceeding ploidy: {}",
        n_rejected_exceed_ploidy
    );
    println!(
        "  Number of snps rejected due to insufficient sampling across samples: {}",
        n_rejected_below_freq
    );

    snps.sort();

    SnpSet {
        sample_names: samples.iter().map(|s| s.name.clone()).collect(),
        k,
        min_count: 2,
        ploidy,
        prefix: prefix,
        snps,
    }
}

fn kmer_to_string(kmer: Kmer, k: usize) -> String {
    let mut s = String::new();
    for i in 0..k {
        let base = (kmer >> (2 * (k - i - 1))) & 0b11;
        s.push(match base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => panic!("Invalid base in kmer: {}", base),
        });
    }
    s
}

fn snp_caller(snp_set: &SnpSet, sample: &mut Sample) {
    let mask: Kmer = !0u64 & !0b11;
    let kmer_set: FxHashSet<Kmer> = sample.kmers.iter().copied().collect();
    sample.genotype.clear();

    for &snp in &snp_set.snps {
        let site = snp & mask;

        // Check for variants by constructing kmers for each possible base
        let mut locus: Locus = 0;

        for variant in 0..4 {
            let candidate = site | variant as Kmer;
            if kmer_set.contains(&candidate) {
                locus |= match variant {
                    0 => 0b0001, // A
                    1 => 0b0010, // C
                    2 => 0b0100, // G
                    3 => 0b1000, // T
                    _ => 0,
                };
            }
        }

        sample.genotype.push(locus);
    }
}

fn locus_to_string(locus: &Locus) -> String {
    let mut s = String::new();
    
    if locus & 0b0001 > 0 {
        s.push('A');
    }

    if locus & 0b0010 > 0 {
        s.push('C');
    }

    if locus & 0b0100 > 0 {
        s.push('G');
    }

    if locus & 0b1000 > 0 {
        s.push('T');
    }
    s
}

fn write_vcf_from_sample(sample: &Sample, snp_set: &SnpSet, outdir: &str) {
    let mut file_path = PathBuf::from(&outdir);
    file_path.push(format!("{}.vcf", sample.name));

    // Open the file for writing
    let mut file = File::create(file_path).expect("Failed to create file");

    // Write the header
    writeln!(file, "##fileformat=VCFv4.2").unwrap();
    writeln!(file, "##source=snpquest").unwrap();
    writeln!(file, "##contig=<ID=1>").unwrap(); // Add contig definition for chromosome 1
    writeln!(file, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").unwrap();
    writeln!(file, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">").unwrap();
    writeln!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample.name).unwrap();

    // Write the records
    for (i, &snp) in snp_set.snps.iter().enumerate() {
        let kmer_string = kmer_to_string(snp, snp_set.k);

        // The last character is the major variant, which we will use as the ref
        let reference = kmer_string.chars().last().unwrap().to_string();

        // Construct the alternate alleles
        // The reference allele is always at index 0, and the alternates follow in a consistent order
        let alleles = vec!["A".to_string(), "C".to_string(), "G".to_string(), "T".to_string()];
        let ref_index = alleles.iter().position(|x| x == &reference).unwrap();

        // Remove the reference allele from the alternates
        let mut alts = alleles.clone();
        alts.remove(ref_index);

        let alts_string = alts.join(",");

        // The ID is the kmer_string with the last character set to X
        let id = kmer_string[0..snp_set.k - 1].to_string() + "X";

        // String with the observed variants in this sample in this snp
        let locus = sample.genotype[i];
        let variants: Vec<char> = locus_to_string(&locus).chars().collect();

        // Map variants to their corresponding allele indices
        let variants_numeric: Vec<String> = variants
            .iter()
            .map(|v| {
                alleles
                    .iter()
                    .position(|x| x == &v.to_string())
                    .unwrap_or_else(|| panic!("Invalid variant: {}", v))
                    .to_string()
            })
            .collect();

        // If there are more than two variants, panic
        if variants_numeric.len() > 2 {
            panic!("More than two variants at site {}: {}", i, variants_numeric.len());
        }

        // Construct the genotype string
        let genotype_text = match variants_numeric.len() {
            0 => "./.:255,255,255".to_string(), // Not called
            1 => {
                if variants_numeric[0] == "0" {
                    format!("{}/{}:10,10,100", variants_numeric[0], variants_numeric[0]) // Homozygous major
                } else {
                    format!("{}/{}:100,100,10", variants_numeric[0], variants_numeric[0]) // Homozygous minor
                }
                
            }
            2 => format!("{}/{}:10,10,100", variants_numeric[0], variants_numeric[1]), // Heterozygote
            _ => unreachable!(),
        };

        // Other fields
        let position = i + 1;
        let chromosome = "1";

        // Write the record
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t30\tPASS\t.\tGT:PL\t{}",
            chromosome, position, id, reference, alts_string, genotype_text
        ).unwrap();
    }
}

fn write_pseudogenome (snp_set: &SnpSet, outdir: &str, run_name: String) {
    let mut file_path = PathBuf::from(&outdir);
    file_path.push(format!("{}.pseudogenome.fasta", run_name));

    // Open the file for writing
    let mut file = File::create(file_path).expect("Failed to create file");

    // Write the header
    writeln!(file, ">1").unwrap();

    // Write the sequence
    for snp in &snp_set.snps {
        let kmer_string = kmer_to_string(*snp, snp_set.k);
        let major_variant = kmer_string.chars().last().unwrap();
        write!(file, "{}", major_variant).unwrap();
    }
}

fn main() {
    let start_run = std::time::Instant::now();

    let args = Args::parse();
    let outdir = if args.outdir.is_empty() {
        "."
    } else {
        &args.outdir
    };

    println!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    println!("{:?}", args);

    let path = PathBuf::from(&args.outdir);
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
            let (kmer_string, count_str) = line.split_once(' ').unwrap();
            let count: u64 = count_str.parse().unwrap();

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
            genotype: vec![],
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
        prefix: String::from(""),
        snps: vec![],
    };

    // If a snp file is provided, read it
    if !args.snp_file.is_empty() {
        println!("Reading snp file: {}", args.snp_file);
        let mut file = File::open(args.snp_file.clone()).expect("Failed to open file");
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer).expect("Failed to read file");

        snp_set = bincode::deserialize(&buffer).expect("Failed to deserialize");
    } else {
        // Discover the snp sites
        let start = std::time::Instant::now();
        println!("Discovering snp sites...");
        snp_set = discover_snp_sites(&samples, k, args.ploidy, args.freq_min, args.prefix);
        println!("Discovering snp sites done, time: {:?}", start.elapsed());


    }
    println!("Number of snps: {}", snp_set.snps.len());

    // Loop over the samples and call the genotype
    println!("Calling genotypes...");
    for sample in samples.iter_mut() {
        snp_caller(&snp_set, sample);
        let n_called = sample.genotype.iter().filter(|&&x| x > 0).count();
        println!(
            "  Sample {}: {} sites called with {} kmers, {} variants",
            sample.name,
            n_called,
            sample.kmers.len(),
            sample.genotype.len()
        );
    }
    println!("Calling genotypes done.");

    // Write the vcf files
    print!("Writing vcf files... ");
    for sample in samples.iter() {
        write_vcf_from_sample(sample, &snp_set, outdir);
    }
    println!("done");

    // Write the pseudogenome
    print!("Writing pseudogenome... ");
    write_pseudogenome(&snp_set, outdir, args.run_name.clone());
    println!("done");

    if args.snp_file.is_empty() {
        // Write the snp set to a file
        print!("Writing snp file: {}.snps... ", args.run_name);
        let mut file_path = PathBuf::from(&outdir);
        file_path.push(format!("{}.snps", args.run_name));
        let mut file = File::create(file_path).expect("Failed to create file");
        let encoded: Vec<u8> = bincode::serialize(&snp_set).expect("Failed to serialize");
        file.write_all(&encoded).expect("Failed to write file");
        println!("done");
    }

    println!("Total run time: {:?}", start_run.elapsed());
}

#[cfg(test)]
mod tests {
    use super::*;

    // Functions used in the tests
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
        let snp_set = discover_snp_sites(&samples, k, ploidy, 0.0, String::from(""));
        let snps = snp_set.snps;
        assert_eq!(snps.len(), 2);

        // Should be in order
        assert_eq!(snps[0], 0b001001);
        assert_eq!(snps[1], 0b001111);
    }
}

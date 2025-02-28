use std::fs::File;
use std::io::BufRead;
use std::io::{Read, Write};
use std::path::Path;
use std::path::PathBuf;
use std::vec;

use clap::Parser;
use rustc_hash::FxHashMap;
use serde::Deserialize;

type Kmer = u64;
type MapType = rustc_hash::FxHashMap<Kmer, u16>;

// Notes
// - To adjust for nucleotide biases in the canonical variant, we will take the reverse
// complement of the variant site of odd variant kmers

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

#[derive(Clone, Debug, Default)]
struct Locus {
    a_count: u16,
    c_count: u16,
    g_count: u16,
    t_count: u16,
}


fn locus_to_base(locus: &Locus) -> String {
    // Convert Locus to IUPAC ambiguity code based on which bases are present
    let mut code = 0u8;
    
    if locus.a_count > 0 { code |= 0b0001; }
    if locus.c_count > 0 { code |= 0b0010; }
    if locus.g_count > 0 { code |= 0b0100; }
    if locus.t_count > 0 { code |= 0b1000; }
    
    match code {
        0b0000 => "N".to_string(),
        0b0001 => "A".to_string(),
        0b0010 => "C".to_string(),
        0b0100 => "G".to_string(),
        0b1000 => "T".to_string(),
        0b0011 => "M".to_string(), // A or C
        0b0101 => "R".to_string(), // A or G
        0b1001 => "W".to_string(), // A or T
        0b0110 => "S".to_string(), // C or G
        0b1010 => "Y".to_string(), // C or T
        0b1100 => "K".to_string(), // G or T
        _ => panic!("Unexpected combination of bases: {}", code),
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'Y' => 'R', // C or T
            'R' => 'Y', // A or G
            'S' => 'S', // C or G
            'W' => 'W', // A or T
            'K' => 'M', // G or T
            'M' => 'K', // A or C
            'B' => 'V', // C or G or T
            'V' => 'B', // A or C or G
            'D' => 'H', // A or G or T
            'H' => 'D', // A or C or T
            'N' => 'N', // A or C or G or T
            _ => panic!("Invalid nucleotide: {}", c),
        })
        .collect()
}

type Genotype = Vec<Locus>;

#[derive(Clone)]
struct Sample {
    name: String,
    path: String,
    k: usize,
    kmers: MapType,
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

struct NucleotideCounts {
    a: u64,
    c: u64,
    g: u64,
    t: u64,
}

impl NucleotideCounts {
    fn add(&mut self, other: &NucleotideCounts) {
        self.a += other.a;
        self.c += other.c;
        self.g += other.g;
        self.t += other.t;
    }

    fn get_pi(&self) -> (f64, f64, f64, f64, u64) {
        let total = self.a + self.c + self.g + self.t;
        if total == 0 {
            return (0.0, 0.0, 0.0, 0.0, 0);
        }
        let a_freq = self.a as f64 / total as f64;
        let c_freq = self.c as f64 / total as f64;
        let g_freq = self.g as f64 / total as f64;
        let t_freq = self.t as f64 / total as f64;
        (a_freq, c_freq, g_freq, t_freq, total)
    }

    fn ingest_kmer_string(&mut self, kmer_string: &str, count: u16) {
        let count = count as u64; // Convert to u64 for internal counting
        for c in kmer_string.chars() {
            match c {
                'A' => self.a += count,
                'C' => self.c += count,
                'G' => self.g += count,
                'T' => self.t += count,
                _ => panic!("Invalid character in kmer string: {}", c),
            }
        }
    }
}

/// SNP analyses from raw sequence reads
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// minimum k-mer count
    #[arg(long, default_value_t = 2)]
    min_count: u16,

    /// ploidy, used to remove snps whose number of variants are more than expected
    /// in one or more samples
    #[arg(long, default_value_t = 2)]
    ploidy: usize,

    /// discard fraction of reads with highest k-mer counts
    #[arg(long, default_value_t = 0.05)]
    discard_fraction: f64,

    /// minimum entropy. kmers with entropy less than this will be discarded
    #[arg(long, default_value_t = 1.5)]
    entropy_min: f64,

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

fn shannon_entropy(kmer: &str) -> f64 {
    let mut counts = [0u32; 4]; // For A, C, G, T

    for base in kmer.chars() {
        match base {
            'A' => counts[0] += 1,
            'C' => counts[1] += 1,
            'G' => counts[2] += 1,
            'T' => counts[3] += 1,
            _ => {} // optionally handle invalid bases if needed
        }
    }

    let k = kmer.len() as f64;
    let mut entropy = 0.0;

    for &count in &counts {
        if count > 0 {
            let p = count as f64 / k;
            entropy -= p * p.log2();
        }
    }

    entropy
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
    for i in 0..(k - 1) {
        kmer_mask <<= 2;
        kmer_mask |= 0b11;
    }
    kmer_mask <<= 2; // Shift once more to mask last two bits

    let mut variants: Variants = Variants::default();

    for i in 0..samples.len() {
        let sample = &samples[i];
        for (kmer, count) in sample.kmers.iter() {
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
            sample_variants[i][variant as usize] += count;
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
        prefix,
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
    let mask: Kmer = !0b11;
    sample.genotype.clear();

    for &snp in &snp_set.snps {
        let site = snp & mask;

        // Check for variants by constructing kmers for each possible base
        let mut locus = Locus::default();

        for variant in 0..4 {
            let candidate = site | variant as Kmer;
            if let Some(count) = sample.kmers.get(&candidate) {
                match variant {
                    0 => locus.a_count = *count, // A
                    1 => locus.c_count = *count, // C
                    2 => locus.g_count = *count, // G
                    3 => locus.t_count = *count, // T
                    _ => {}
                };
            }
        }

        sample.genotype.push(locus);
    }
}

fn locus_to_string(locus: &Locus) -> String {
    let mut s = String::new();

    if locus.a_count > 0 {
        s.push('A');
    }

    if locus.c_count > 0 {
        s.push('C');
    }

    if locus.g_count > 0 {
        s.push('G');
    }

    if locus.t_count > 0 {
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
    writeln!(
        file,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )
    .unwrap();
    writeln!(
        file,
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">"
    )
    .unwrap();
    writeln!(
        file,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        sample.name
    )
    .unwrap();

    // Write the records
    for (i, &snp) in snp_set.snps.iter().enumerate() {
        let kmer_string = kmer_to_string(snp, snp_set.k);

        // The last character is the major variant, which we will use as the ref
        let mut reference = kmer_string.chars().last().unwrap().to_string();

        // If the kmer is odd, take the reverse complement
        if i % 2 == 1 {
            let ref_char = match reference.chars().next().unwrap() {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => panic!("Invalid character in kmer string: {}", reference),
            };
            reference = ref_char.to_string();
        }

        // The ID is the kmer_string with the last character set to X
        let id = kmer_string[0..snp_set.k - 1].to_string() + "X";

        // String with the observed variants in this sample in this snp
        let locus = &sample.genotype[i];
        let mut variants: Vec<char> = locus_to_string(locus).chars().collect();

        // Take reverse complement of the variants if the index is odd
        if i % 2 == 1 {
            for variant in &mut variants {
                *variant = match *variant {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => panic!("Invalid character in kmer string: {}", variant),
                };
            }
        }

        // If there are more than one variants, sort them according to the following criteria:
        // - If one is the same as the reference, it is first
        // - All other variants are sorted in alphabetical order
        variants.sort_by(|a, b| {
            if a == &reference.chars().next().unwrap() {
                std::cmp::Ordering::Less
            } else if b == &reference.chars().next().unwrap() {
                std::cmp::Ordering::Greater
            } else {
                a.cmp(b)
            }
        });

        // If there are more than two variants, panic
        if variants.len() > 2 {
            panic!(
                "More than two variants at site, ploidy > 2 not yet supported {}: {}",
                i,
                variants.len()
            );
        }

        // Construct the alt and genotype string
        let genotype_string: String;
        let alts_string: String;

        // There is a single variant, and it is the same as the reference
        // homozygous major
        if variants.len() == 1 && variants[0] == reference.chars().next().unwrap() {
            alts_string = ".".to_string();
            genotype_string = "0/0:10".to_string();
        }
        // There is a single variant, and it is different from the reference
        // homozygous minor
        else if variants.len() == 1 && variants[0] != reference.chars().next().unwrap() {
            alts_string = variants[0].to_string();
            genotype_string = "1/1:10".to_string();
        }
        // There are two variants, one is the same as the reference
        // heterozygous major/minor
        else if variants.len() == 2 && variants[0] == reference.chars().next().unwrap() {
            alts_string = variants[1].to_string();
            genotype_string = "0/1:10,10,100".to_string();
        }
        // There are two variants, neither is the same as the reference
        // heterozygous minor/minor
        else if variants.len() == 2 && variants[0] != reference.chars().next().unwrap() {
            alts_string = format!("{},{}", variants[0], variants[1]);
            genotype_string = "1/2:100,10,10".to_string();
        }
        // Unsupported case
        else {
            panic!(
                "Unexpected case!\n  kmer_string:{}\n  ref:{}\n  variants:{}",
                kmer_string,
                reference,
                variants.iter().collect::<String>()
            );
        }

        // Other fields
        let position = i + 1;
        let chromosome = "1";

        // Write the record
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t30\tPASS\t.\tGT:PL\t{}",
            chromosome, position, id, reference, alts_string, genotype_string
        )
        .unwrap();
    }
}

fn write_fasta_from_sample(sample: &Sample, snp_set: &SnpSet, outdir: &str) {
    let mut file_path = PathBuf::from(&outdir);
    file_path.push(format!("{}.fasta", sample.name));

    // Open the file for writing
    let mut file = File::create(file_path).expect("Failed to create file");

    // Write the header with sample name
    writeln!(file, ">{}", sample.name).unwrap();
    
    // Write the sequence
    // Reverse complement the sequence if the index is odd
    // If homozygous, write the resolved base
    // If heterozygous, write the ambiguous base
    for (i, &snp) in snp_set.snps.iter().enumerate() {
        let locus = &sample.genotype[i];
        let mut base = locus_to_base(locus);
        // Adjust the base if the index is odd
        if i % 2 == 1 {
            base = reverse_complement(&base);
        }
        write!(file, "{}", base).unwrap();
    }
    writeln!(file).unwrap();

}

fn write_pseudogenome(snp_set: &SnpSet, outdir: &str, run_name: String) {
    let mut file_path = PathBuf::from(&outdir);
    file_path.push(format!("{}.pseudogenome.fasta", run_name));

    // Open the file for writing
    let mut file = File::create(file_path).expect("Failed to create file");

    // Write the header
    writeln!(file, ">1").unwrap();

    // Write the sequence
    for (i, snp) in snp_set.snps.iter().enumerate() {
        let kmer_string = kmer_to_string(*snp, snp_set.k);
        let major_variant = kmer_string.chars().last().unwrap();
        // Adjust the major variant if the index is odd
        let major_variant = if i % 2 == 1 {
            match major_variant {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => panic!("Invalid character in kmer string: {}", major_variant),
            }
        } else {
            major_variant
        };
        write!(file, "{}", major_variant).unwrap();
    }
}

fn snp_major_nucleotide_count(snp_set: &SnpSet, adjusted: bool) -> NucleotideCounts {
    let mut counts = NucleotideCounts {
        a: 0,
        c: 0,
        g: 0,
        t: 0,
    };

    for (i, snp) in snp_set.snps.iter().enumerate() {
        let kmer_string = kmer_to_string(*snp, snp_set.k);
        let major_variant = kmer_string.chars().last().unwrap();

        // Adjust the major variant if the index is odd
        let major_variant = if adjusted && i % 2 == 1 {
            match major_variant {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => panic!("Invalid character in kmer string: {}", major_variant),
            }
        } else {
            major_variant
        };
        match major_variant {
            'A' => counts.a += 1,
            'C' => counts.c += 1,
            'G' => counts.g += 1,
            'T' => counts.t += 1,
            _ => panic!("Invalid character in kmer string: {}", major_variant),
        }
    }

    counts
}

fn process_kmer_file(file_name: &str, prefix: &str, discard_fraction: f64, entropy_min: f64, min_count: u16) -> (Sample, NucleotideCounts, usize) {
    println!(" Reading {}...", file_name);
    // Get the sample name from the file name by removing the path and extension
    let sample_name = Path::new(file_name)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");

    let mut line_n: u64 = 0;
    // Open the file for buffered reading
    let file_path = Path::new(file_name);
    let file = std::fs::File::open(file_path).unwrap();
    let reader = std::io::BufReader::new(file);

    let mut kmer_counts = MapType::default();

    let mut nuc_counts_sample = NucleotideCounts {
        a: 0,
        c: 0,
        g: 0,
        t: 0,
    };

    let mut k: usize = 0;

    let mut entropies: Vec<f64> = vec![]; 
    
    // Iterate over the lines of the file
    for line in reader.lines() {
        line_n += 1;
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("Error reading line {}: {}", line_n, e);
                continue;
            }
        };

        // Attempt to split the line into two fields
        let (kmer_string, count_str) = match line.split_once(' ') {
            Some((kmer, count)) => (kmer, count),
            None => {
                eprintln!("Invalid line format at line {}: {}", line_n, line);
                continue;
            }
        };

        // Attempt to parse the count
        let count: u16 = match count_str.parse() {
            Ok(c) => c,
            Err(e) => {
                eprintln!(
                    "Invalid count at line {}: {} (error: {})",
                    line_n, count_str, e
                );
                continue;
            }
        };

        // continue if the count is less than the minimum
        if count < min_count {
            continue;
        }

        // continue if the kmer does not start with the prefix
        if !prefix.is_empty() && !kmer_string.starts_with(prefix) {
            continue;
        }

        let entropy = shannon_entropy(kmer_string);
        entropies.push(entropy);

        // continue if the entropy is less than the minimum
        if entropy < entropy_min {
            continue;
        }

        let kmer = string_to_kmer(kmer_string);

        kmer_counts.insert(kmer, count);
        if k == 0 {
            k = kmer_string.len();
        } else if k != kmer_string.len() {
            panic!("kmer length mismatch: {} != {}", k, kmer_string.len());
        }

        nuc_counts_sample.ingest_kmer_string(kmer_string, count);
    }

    // Print statistics for this sample
    let (a_freq, c_freq, g_freq, t_freq, total) = nuc_counts_sample.get_pi();
    println!("  Number of processed kmers: {}", line_n);
    println!("  Number of ingested nucleotides: {}", total);
    println!("  Number of ingested kmers: {}", kmer_counts.len());
    println!(
        "  Total count of ingested kmers: {}",
        kmer_counts.values().map(|&c| c as u64).sum::<u64>()
    );
    println!(
        "  Sample {} kmer pi: {:.3} A, {:.3} C, {:.3} G, {:.3} T",
        sample_name, a_freq, c_freq, g_freq, t_freq
    );
    println!(
        "  Sample {} expected genome pi: {:.3} A, {:.3} C, {:.3} G, {:.3} T",
        sample_name,
        (a_freq + t_freq) / 2.0,
        (c_freq + g_freq) / 2.0,
        (c_freq + g_freq) / 2.0,
        (a_freq + t_freq) / 2.0
    );
    println!("  kmer length: {}", k);
    println!("  kmer entropy:");
    println!("    Mean kmer entropy: {:.3}", entropies.iter().sum::<f64>() / entropies.len() as f64);
    println!("    Median kmer entropy: {:.3}", {
        let mut entropies = entropies.clone();
        entropies.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mid = entropies.len() / 2;
        if entropies.len() % 2 == 0 {
            (entropies[mid - 1] + entropies[mid]) / 2.0
        } else {
            entropies[mid]
        }
    });
    println!("    Fraction < 1.0: {:.3}", entropies.iter().filter(|&&x| x < 1.0).count() as f64 / entropies.len() as f64);
    println!("    Fraction < 1.5: {:.3}", entropies.iter().filter(|&&x| x < 1.5).count() as f64 / entropies.len() as f64);

    // Remove the most common kmers
    let kmer_counts = filter_high_abundance_kmers(kmer_counts, discard_fraction);
    
    // Create the sample
    let sample = Sample {
        name: sample_name.to_string(),
        path: file_name.to_string(),
        k,
        kmers: kmer_counts,
        genotype: vec![],
    };

    (sample, nuc_counts_sample, k)
}

fn filter_high_abundance_kmers(mut kmer_counts: MapType, discard_fraction: f64) -> MapType {
    let mut counts: Vec<_> = kmer_counts.iter().collect();
    counts.sort_by(|a, b| b.1.cmp(a.1));
    let n_to_discard = (discard_fraction * counts.len() as f64).round() as usize;
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
    kmer_counts
}

fn print_nucleotide_stats(nuc_counts: &NucleotideCounts, label: &str) {
    let (a_freq, c_freq, g_freq, t_freq, total) = nuc_counts.get_pi();
    println!(
        "  {} kmer pi: {:.3} A, {:.3} C, {:.3} G, {:.3} T",
        label, a_freq, c_freq, g_freq, t_freq
    );
    println!(
        "  {} expected genome pi: {:.3} A, {:.3} C, {:.3} G, {:.3} T",
        label,
        (a_freq + t_freq) / 2.0,
        (c_freq + g_freq) / 2.0,
        (c_freq + g_freq) / 2.0,
        (a_freq + t_freq) / 2.0
    );
    println!("  Total number of ingested nucleotides: {}", total);
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

    // Create output directory if it doesn't exist
    let path = PathBuf::from(&args.outdir);
    let directory = format!("{}/", path.to_str().unwrap());
    std::fs::create_dir_all(&directory).unwrap();

    let start = std::time::Instant::now();
    println!("Ingesting kmers...");
    std::io::stdout().flush().unwrap();

    // Create a vector of samples
    let mut samples: Vec<Sample> = vec![];
    let mut nuc_counts_all = NucleotideCounts { a: 0, c: 0, g: 0, t: 0 };
    let mut k: usize = 0;

    // Process each input file
    for file_name in args.input.iter() {
        let (sample, nuc_counts_sample, sample_k) = process_kmer_file(file_name, &args.prefix, args.discard_fraction, args.entropy_min, args.min_count);
        
        if k == 0 {
            k = sample_k;
        }
        
        nuc_counts_all.add(&nuc_counts_sample);
        samples.push(sample);
    }

    // Verify that all samples have the same kmer length
    let k = samples[0].k;
    for sample in samples.iter() {
        if sample.k != k {
            panic!("kmer length mismatch across samples: {} != {}", k, sample.k);
        }
    }

    let (a_freq, c_freq, g_freq, t_freq, total) = nuc_counts_all.get_pi();
    print_nucleotide_stats(&nuc_counts_all, "All");
    println!("Ingesting samples done, time: {:?}", start.elapsed());

    // Load or discover SNP sites
    let snp_set = if !args.snp_file.is_empty() {
        println!("Reading snp file: {}", args.snp_file);
        let mut file = File::open(&args.snp_file).expect("Failed to open file");
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer).expect("Failed to read file");
        bincode::deserialize(&buffer).expect("Failed to deserialize")
    } else {
        // Discover the snp sites
        let start = std::time::Instant::now();
        println!("Discovering snp sites...");
        let set = discover_snp_sites(&samples, k, args.ploidy, args.freq_min, args.prefix.clone());
        println!("Discovering snp sites done, time: {:?}", start.elapsed());
        set
    };
    println!("Number of snps: {}", snp_set.snps.len());

    // Loop over the samples and call the genotype
    println!("Calling genotypes...");
    for sample in samples.iter_mut() {
        snp_caller(&snp_set, sample);
        let n_called = sample.genotype.iter().filter(|x| {
            x.a_count > 0 || x.c_count > 0 || x.g_count > 0 || x.t_count > 0
        }).count();
        println!(
            "  Sample {}: {} sites called with {} kmers, {} variants",
            sample.name,
            n_called,
            sample.kmers.len(),
            sample.genotype.len()
        );
    }
    println!("Calling genotypes done.");

    // get pi
    let counts = snp_major_nucleotide_count(&snp_set, false);
    let (a_freq, c_freq, g_freq, t_freq, _total) = counts.get_pi();
    println!(
        "Unadjusted major variant pi: {:.3} A, {:.3} C, {:.3} G, {:.3} T",
        a_freq, c_freq, g_freq, t_freq
    );

    let counts = snp_major_nucleotide_count(&snp_set, true);
    let (a_freq, c_freq, g_freq, t_freq, _total) = counts.get_pi();
    println!(
        "Adjusted major variant pi: {:.3} A, {:.3} C, {:.3} G, {:.3} T",
        a_freq, c_freq, g_freq, t_freq
    );

    // Write the vcf files
    print!("Writing vcf and fasta files... ");
    for sample in samples.iter() {
        write_vcf_from_sample(sample, &snp_set, outdir);
        write_fasta_from_sample(sample, &snp_set, outdir);
    }
    println!("done");

    // Write the pseudogenome
    print!("Writing pseudogenome... ");
    write_pseudogenome(&snp_set, outdir, args.run_name.clone());
    println!("done");

    if args.snp_file.is_empty() {
        // Write the snp set to a file
        print!("Writing snp file: {}.snps... ", args.run_name);
        let mut file_path = PathBuf::from(outdir);
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
        let mut kmers1 = FxHashMap::default();
        kmers1.insert(0b001111, 1); // Two variants at one site, one shared with sample 2
        kmers1.insert(0b001101, 1);
        kmers1.insert(0b101110, 1); // Exceed ploidy, should be excluded
        kmers1.insert(0b111100, 1);
        kmers1.insert(0b111101, 1);
        kmers1.insert(0b111110, 1);

        let sample1 = Sample {
            name: "sample1".to_string(),
            path: "/path/to/sample1".to_string(),
            k: 3,
            kmers: kmers1,
            genotype: vec![],
        };

        let mut kmers2 = FxHashMap::default();
        kmers2.insert(0b001001, 1); // Two variants at one site, unique to sample2
        kmers2.insert(0b001011, 1);
        kmers2.insert(0b001111, 1); // Should be included due to variation in sample1
        kmers2.insert(0b111111, 1); // Should be excluded due to ploidy in sample1

        let sample2 = Sample {
            name: "sample2".to_string(),
            path: "/path/to/sample2".to_string(),
            k: 3,
            kmers: kmers2,
            genotype: vec![],
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

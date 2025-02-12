use clap::Parser;
use std::io::BufRead;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
use std::vec;
use rustc_hash::FxHashMap;

type Kmer = u64;
type MapType = rustc_hash::FxHashMap<Kmer, u64>;


struct Sample {
    name: String,
    path: String,
    k: usize,
    kmers: Vec<Kmer>,
}

/// A collection of kmer counting and analysis tools
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// minimum k-mer count
    #[arg(long, default_value_t = 2)]
    min_count: usize,

    /// discard fraction of reads with highest k-mer counts
    #[arg(long, default_value_t = 0.05)]
    discard_fraction: f64,

    /// Name of the sample, could be species and sample ID eg Nanomia-bijuga-YPMIZ035039. 
    /// Will be used as the prefix for output files
    #[arg(short, long, default_value_t = String::from("sample") )]
    sample: String,

    /// kmer prefix for thinning kmers
    #[arg(short, long, default_value_t = String::from("") )]
    prefix: String,

    /// Directory for output files, defaults to current directory
    #[arg(short, long, default_value_t = String::from("./") )]
    outdir: String,

    /// Input files, fastq. If no files are specified, data will be read
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

fn main() {
    let start_run = std::time::Instant::now();

    // Ingest command line arguments
    let args = Args::parse();

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
            let kmer = string_to_kmer(kmer_string);
            let count: u64 = parts.next().unwrap().parse().unwrap();
            kmer_counts.insert(kmer, count);
            if k == 0 {
                k = kmer_string.len();
            }
            else if k != kmer_string.len() {
                panic!("kmer length mismatch: {} != {}", k, kmer_string.len());
            }
        }
        println!("  Number of processed kmers: {}", line_n);
        println!("  Number of ingested kmers: {}", kmer_counts.len());
        println!("  Total count of ingested kmers: {}", kmer_counts.values().sum::<u64>());
        println!("  kmer length: {}", k);

        // remove the most common kmers
        let mut counts: Vec<_> = kmer_counts.iter().collect();
        counts.sort_by(|a, b| b.1.cmp(a.1));
        let n_to_discard = (args.discard_fraction * counts.len() as f64).round() as usize;
        println!("  Discarding {} kmers with highest counts", n_to_discard);
        let keys_to_remove: Vec<Kmer> = counts.iter().take(n_to_discard).map(|(kmer, _)| **kmer).collect();
        for kmer in keys_to_remove {
            kmer_counts.remove(&kmer);
        }
        println!("  Number of kmers after discarding: {}", kmer_counts.len());

        let mut kmers: Vec<Kmer> = kmer_counts.keys().copied().collect();
        kmers.sort();

        samples.push(Sample {
            name: sample_name.to_string(),
            path: file_name.to_string(),
            k: k,
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

    println!("Ingesting samples done.");

    print!("Combining kmers across samples... ");
    let mut kmers_all: Vec<Kmer> = vec![];
    for sample in samples.iter() {
        kmers_all.extend(sample.kmers.iter());
    }
    kmers_all.sort();
    kmers_all.dedup();
    println!("done.");

    print!("Finding snps... ");
    // Create a Kmer bitmask for all but the last two bits
    let mut kmer_mask: Kmer = 0;
    for _ in 0..(k - 1) {
        kmer_mask |= 0b11;
        kmer_mask <<= 2;
    }

    let mut snps: Vec<Kmer> = vec![];

    // Iterate over the kmers and find those that are variable in the sample set
    for i in 0..kmers_all.len() {
        let mut is_variable = false;
        if i > 0 {
            let kmer = kmers_all[i];
            let kmer_prev = kmers_all[i - 1];
            if (kmer & kmer_mask) == (kmer_prev & kmer_mask) {
                is_variable = true;
            }
        }

        if i < kmers_all.len() - 1 {
            let kmer = kmers_all[i];
            let kmer_next = kmers_all[i + 1];
            if (kmer & kmer_mask) == (kmer_next & kmer_mask) {
                is_variable = true;
            }
        }

        if is_variable {
            snps.push(kmers_all[i]);
        }
    }
    println!("done.");

    println!("Number of snps: {}", snps.len());




    println!("Total run time: {:?}", start_run.elapsed());
}

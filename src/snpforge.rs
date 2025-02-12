use clap::Parser;
use std::io::BufRead;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

use crate::kmer::Chunk;
use crate::kmer::KmerCounts;

mod kmer;

const N_READS_PER_BATCH: u64 = 1000;

/// A collection of kmer counting and analysis tools
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// k-mer length
    #[arg(short, default_value_t = 15)]
    k: usize,

    /// Maximum value for histogram
    #[arg(long, default_value_t = 10000)]
    histo_max: u64,

    /// Number of chunks to divide the data into
    #[arg(short, default_value_t = 100)]
    n: usize,

    /// Maximum number of reads to process
    #[arg(short = 'm', long, default_value_t = 0)]
    max_reads: u64,

    /// Name of the sample, could be species and sample ID eg Nanomia-bijuga-YPMIZ035039. 
    /// Will be used as the prefix for output files
    #[arg(short, long, default_value_t = String::from("sample") )]
    sample: String,

    /// Directory for output files, defaults to current directory
    #[arg(short, long, default_value_t = String::from("./") )]
    outdir: String,

    /// Input files, fastq. If no files are specified, data will be read
    /// from stdin. This can be used to uncompress a gz file and send them
    /// to sharkmer.
    #[arg()]
    input: Option<Vec<String>>,

    /// Verbosity
    #[arg(long, default_value_t = 0)]
    verbosity: usize,
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

    let k = args.k;

    // Check that the arguments are valid
    assert!(
        k < 32,
        "k must be less than 16"
    );
    assert!(k > 0, "k must be greater than 0");
    assert!(k % 2 == 1, "k must be odd");
    assert!(args.histo_max > 0, "histo_max must be greater than 0");
    assert!(args.n > 0, "n must be greater than 0");

    // Ingest the fastq data
    let start = std::time::Instant::now();
    print!("Ingesting reads...");
    std::io::stdout().flush().unwrap();

    let mut chunks: Vec<kmer::Chunk> = Vec::new();
    for _ in 0..args.n {
        chunks.push(Chunk::new(&k));
    }
    let mut chunk_index: usize = 0;

    let mut reads: Vec<kmer::Read> = Vec::new();
    let mut n_reads_read: u64 = 0;
    let mut n_bases_read: u64 = 0;

    match &args.input {
        Some(input_files) => {
            // read from one or more files
            'processing_files: for file_name in input_files.iter() {
                let mut line_n: u64 = 0;
                // Open the file for buffered reading
                let file_path = Path::new(&file_name);
                let file = std::fs::File::open(file_path).unwrap();
                let reader = std::io::BufReader::new(file);
                // Iterate over the lines of the file
                for line in reader.lines() {
                    line_n += 1;
                    if line_n % 4 == 2 {
                        // This is a sequence line
                        let line = line.unwrap();
                        n_bases_read += line.len() as u64;
                        let new_reads = kmer::seq_to_reads(&line);
                        reads.extend(new_reads);
                        n_reads_read += 1;

                        // If we have read enough reads, ingest them into current chunk
                        if n_reads_read % N_READS_PER_BATCH == 0 {
                            chunks[chunk_index].ingest_reads(&reads);
                            chunk_index += 1;
                            if chunk_index == args.n {
                                chunk_index = 0;
                            }
                            reads.clear();
                        }
                    }
                    if args.max_reads > 0 && n_reads_read >= args.max_reads {
                        break 'processing_files;
                    }
                }
            }
        }
        None => {
            // read from stdin
            let stdin = std::io::stdin();
            // Lock stdin for exclusive access
            let handle = stdin.lock();

            let mut line_n: u64 = 0;

            // Create a buffer for reading lines
            for line in handle.lines() {
                line_n += 1;
                if line_n % 4 == 2 {
                    // This is a sequence line
                    let line = line.unwrap();
                    n_bases_read += line.len() as u64;
                    let new_reads = kmer::seq_to_reads(&line);
                    reads.extend(new_reads);
                    n_reads_read += 1;

                    // If we have read enough reads, ingest them into current chunk
                    if n_reads_read % N_READS_PER_BATCH == 0 {
                        chunks[chunk_index].ingest_reads(&reads);
                        chunk_index += 1;
                        if chunk_index == args.n {
                            chunk_index = 0;
                        }
                        reads.clear();
                    }
                }
                if args.max_reads > 0 && n_reads_read >= args.max_reads {
                    break;
                }
            }
        }
    }

    // Ingest any remaining reads
    chunks[chunk_index].ingest_reads(&reads);
    reads.clear();

    println!(" done");

    let mut n_reads_ingested: u64 = 0;
    let mut n_bases_ingested: u64 = 0;
    let mut n_kmers_ingested: u64 = 0;
    for chunk in chunks.iter() {
        n_reads_ingested += chunk.get_n_reads();
        n_bases_ingested += chunk.get_n_bases();
        n_kmers_ingested += chunk.get_n_kmers();
    }

    // Print some stats
    println!("  Read {} reads", n_reads_read);
    println!("  Read {} bases", n_bases_read);
    println!("  Ingested {} subreads", n_reads_ingested);
    println!("  Ingested {} bases", n_bases_ingested);
    println!("  Ingested {} kmers", n_kmers_ingested);
    println!("  Time to ingest reads: {:?}", start.elapsed());

    // Create the histograms
    print!("Consolidating chunks and creating histograms...");
    let start = std::time::Instant::now();
    std::io::stdout().flush().unwrap();

    let mut kmer_counts: KmerCounts = KmerCounts::new(&k);
    let mut histos: Vec<kmer::Histogram> = Vec::with_capacity(args.n);

    // Iterate over the chunks
    for chunk in chunks.iter() {
        kmer_counts.extend(chunk.get_kmer_counts());

        let histo = kmer::Histogram::from_kmer_counts(&kmer_counts, &args.histo_max);

        histos.push(histo);
    }
    println!(" done, time: {:?}", start.elapsed());

    let n_hashed_kmers: u64 = kmer_counts.get_n_kmers();
    println!(
        "  {} unique kmers with a total count of {} were found",
        kmer_counts.get_n_unique_kmers(),
        n_hashed_kmers
    );

    if n_hashed_kmers != n_kmers_ingested {
        panic!(
            "The total count of hashed kmers ({}) does not equal the number of ingested kmers ({})",
            n_hashed_kmers, n_kmers_ingested,
        );
    }

    // Write the histograms to a tab delimited file, with the first column being the count
    // Skip the first row, which is the count of 0. Do not include a header
    print!("Writing histograms to file...");
    std::io::stdout().flush().unwrap();
    let mut file = std::fs::File::create(format!("{}{}.histo", directory, args.sample)).unwrap();
    for i in 1..args.histo_max as usize + 2 {
        let mut line = format!("{}", i);
        for histo in histos.iter() {
            let histo_vec = kmer::Histogram::get_vector(histo);
            line = format!("{}\t{}", line, histo_vec[i]);
        }
        line = format!("{}\n", line);
        file.write_all(line.as_bytes()).unwrap();
    }

    println!(" done");

    // Write the final histogram to a file, ready for GenomeScope2 etc...
    print!("Writing final histogram to file...");
    std::io::stdout().flush().unwrap();

    let last_histo = &histos[histos.len() - 1];
    let last_histo_vec = kmer::Histogram::get_vector(last_histo);

    let mut file =
        std::fs::File::create(format!("{}{}.final.histo", directory, args.sample)).unwrap();
    for i in 1..args.histo_max as usize + 2 {
        let mut line = format!("{}", i);

        line = format!("{}\t{}", line, last_histo_vec[i]);

        line = format!("{}\n", line);
        file.write_all(line.as_bytes()).unwrap();
    }

    println!(" done");

    let n_singleton_kmers = last_histo_vec[1];
    let n_unique_kmers_histo: u64 = last_histo.get_n_unique_kmers();
    let n_kmers_histo: u64 = last_histo.get_n_kmers();

    println!("  {} unique kmers in histogram", n_unique_kmers_histo);
    println!("  {} kmers in histogram", n_kmers_histo);

    if n_kmers_histo != n_kmers_ingested {
        panic!(
            "The total count of kmers in the histogram ({}) does not equal the total expected count of kmers ({})",
            n_kmers_histo,
            n_kmers_ingested,
        );
    }

    if n_unique_kmers_histo != kmer_counts.get_n_unique_kmers() {
        panic!(
            "The total count of unique kmers in the histogram ({}) does not equal the total count of hashed kmers ({})",
            n_unique_kmers_histo,
            kmer_counts.get_n_unique_kmers(),
        );
    }

    print!("Writing stats to file...");
    std::io::stdout().flush().unwrap();
    let mut file_stats =
        std::fs::File::create(format!("{}{}.stats", directory, args.sample)).unwrap();
    let mut line = format!("arguments\t{:?}\n", args);
    line = format!("{}kmer_length\t{}\n", line, args.k);
    line = format!("{}n_reads_read\t{}\n", line, n_reads_read);
    line = format!("{}n_bases_read\t{}\n", line, n_bases_read);
    line = format!("{}n_subreads_ingested\t{}\n", line, reads.len());
    line = format!("{}n_bases_ingested\t{}\n", line, n_bases_ingested);
    line = format!("{}n_kmers\t{}\n", line, n_kmers_ingested);
    line = format!(
        "{}n_multi_kmers\t{}\n",
        line,
        n_kmers_ingested - n_singleton_kmers
    );
    line = format!("{}n_singleton_kmers\t{}\n", line, n_singleton_kmers);

    file_stats.write_all(line.as_bytes()).unwrap();
    println!(" done");

    println!("Total run time: {:?}", start_run.elapsed());
}

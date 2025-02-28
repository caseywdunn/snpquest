# snpquest

In progress. Use at your own risk.

Reference free snp calling from raw reads.

## Overview



## Installation

First, [install the rust build tools](https://www.rust-lang.org/tools/install).

Then clone this repository and build the executable:

    git clone https://github.com/caseywdunn/snpquest.git  # Or download and expand the zip file
    cd snpquest
    cargo build --release

The binary will be in `target/release/snpquest`. You can move it to a directory in your path, or run it from the current directory.

    ./target/release/snpquest --help

## Usage

### Preparing files

For each sample, you will:

- Decompress the fastq.gz files
- Optionally subset the reads, for example take only the first 100M read pairs
- Trim them to remove adapters and low quality regions
- Count the kmers

Below is an example workflow for a single sample.

```bash
    # Decompress the fastq files and retain only the first 400M lines,
    # which (since there are 4 lines per read) is 100M read pairs
    mkdir reads
    zcat raw_reads/sample01_R1.fastq.gz | head -n 400000000 > reads/sample01_R1.fastq
    zcat raw_reads/sample01_R2.fastq.gz | head -n 400000000 > reads/sample01_R2.fastq

    # Trim the reads
    echo ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" > TruSeq3-PE.fa
    trimmomatic PE \
        -threads 4 \
        reads/sample01_R1.fastq \
        reads/sample01_R2.fastq \
        reads/sample01_trimmed_paired_R1.fastq \
        reads/sample01_trimmed_unpaired_R1.fastq \
        reads/sample01_trimmed_paired_R2.fastq \
        reads/sample01_trimmed_unpaired_R2.fastq \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:15 TRAILING:15 MINLEN:50
    
    # Count the kmers using kmer length 31
    mkdir kmers
    jellyfish count -C -m 31 -s 1000000000 -t 4 \
          -o kmers/sample01.jf \
          reads/sample01_trimmed_paired_R1.fastq reads/sample01_trimmed_paired_R1.fastq

    # Dump the kmers to a text file
    jellyfish dump -c -L 2 -o kmers/sample01.kmers kmers/sample01.jf
```

Run this workflow for each sample.

### Run snpquest

Once you have prepared the files for each sample, you can run snpquest.

```bash
    # Run snpquest
    ../target/release/snpquest \
      --outdir output_snp/ \
      --freq-min 1.0 \
      --run-name Physalia_k{params.thisk} \
      kmers/*.kmers
```

This will create a directory `output_snp/` with vnf and fasta files for each sample.

### PCA

You can run a PCA on the vnf files as follows:

```bash

```

### Build a phylogeny

You can build a phylogeny as follows:

```bash
    cd output_snp
    cat *.fasta > all.fasta
    iqtree2 -s all.fasta -B 1000 -nt AUTO -st DNA
```

## Development

Common development tasks:

    cargo test
    cargo clippy
    cargo clippy --fix
    cargo fmt
    cargo build # Debug
    cargo build --release
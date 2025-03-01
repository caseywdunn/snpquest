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

### PCA and Admixture

You can run a PCA on the vnf files as follows:

```bash
ml PLINK
ml BCFtools
ml ADMIXTURE

# Retain sites where one allele is the reference allele, and compress the files
for file in *.vcf; do bcftools view -i 'GT~"0"' -Oz -o "${file}.gz" "$file"; done

for file in *.vcf.gz; do tabix -p vcf "$file"; done  

NAME=pilot

bcftools merge -oz -o $NAME.merged.vcf.gz *.vcf.gz  > $NAME.bcftools.log

plink2 --vcf $NAME.merged.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out $NAME.plink --bad-freqs --max-alleles 2 > $NAME.plink.log

num_clusters=5
admixture --cv $NAME.plink.bed $num_clusters > $NAME.admixture.log
```

Then, to view the PCA in R:

```r
library(ggplot2)
library(ggrepel)

pca <- read.table("pilot.plink.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca)-2)))

ggplot(pca, aes(x=PC1, y=PC2, label=FID)) +
  geom_point(size=3) +
  geom_text_repel(size=3, max.overlaps=20, box.padding=0.5, point.padding=0.5) +
  theme_minimal() +
  labs(title="PCA Plot", x="PC1", y="PC2")
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
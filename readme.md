# snpquest

Reference free snp calling from raw reads.

## Installation

First, [install the rust build tools](https://www.rust-lang.org/tools/install).

Then clone this repository and build the executable:

    git clone https://github.com/caseywdunn/snpquest.git  # Or download and expand the zip file
    cd snpquest
    cargo build --release

## Development

Common development tasks:

    cargo test
    cargo clippy
    cargo clippy --fix
    cargo fmt
    cargo build # Debug
    cargo build --release
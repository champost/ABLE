## Scripts

This folder contains `python` scripts which define the demographic models M1, M2 and M6 models for use with [`dadi`](https://doi.org/10.1371/journal.pgen.1000695) and to convert `fasta` files into the `pseudo_MS` format.

### Python code for a [`dadi`](https://doi.org/10.1371/journal.pgen.1000695) analysis

The following files can be imported into the main python script which performs the `dadi` inference.

- `config_M1_dadi.py`
- `config_M2_dadi.py`
- `config_M6_dadi.py`


### Python code for converting `fasta` files into the `pseudo_MS` format

- `block_cutter_fasta.py`

Usage : `python block_cutter_fasta.py 2000 1800 1 *.fa`

Here, we are converting a set of aligned `fasta` files and reading in "2kb" blocks sequentially with an overall tolerance of 200bp of missing information. The `1` signifies that we would like to keep the accepted amount of sequence (which would be at least 1800bp) as a single block for every consecutive stretch of 2kb.

The `pseudo_MS` files contained in the `data` folder have been generated with the following command lines.

- `Orangutan_500bp_blocks.txt` : `python block_cutter_fasta.py 2000 1800 4 *.fa`
- `Orangutan_2kb_blocks.txt` : `python block_cutter_fasta.py 2000 1800 1 *.fa`



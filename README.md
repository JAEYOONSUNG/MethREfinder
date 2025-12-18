# MethREfinder

MethREfinder is a computational tool designed to identify methylation-sensitive restriction enzyme recognition sites within a given DNA sequence. This package leverages the REBASE database to cross-reference known restriction enzymes and determine their sensitivity to DNA methylation modifications. Researchers can use MethREfinder to assess the compatibility of restriction enzymes with methylated DNA regions, aiding in the selection of appropriate enzymes for molecular cloning, epigenetic studies, and DNA methylation analysis.

## Key Features
- **Identification of restriction enzyme recognition sites within user-provided DNA sequences.**
- **Detection of methylation-sensitive recognition sites based on REBASE methylation data.**
- **Support for common methylation types, including CpG and adenine methylation.**
- **Streamlined and efficient sequence analysis workflow.**
- **Informed enzyme selection for experiments involving methylated DNA.**

## Anaylsis flow
![MethREfinder](https://github.com/user-attachments/assets/7e0dd79a-ccb2-440a-9d7e-3ebfc1e45b00)


## Installation
### Requirements

The MethREfinder is supported for macOS, Linux and Windows machines, which can provide an environment for using R.
It requires R version >=4.2.1 for release, and R version >=4.3 for devel.

```r
# Install from GitHub (example, adjust based on your repo)
install.packages("devtools")
devtools::install_github("JAEYOONSUNG/MethREfinder")
```

## Getting Started
#### Parameter Descriptions

- **target_seq**: A DNA sequence provided as a character string. The sequence can include standard nucleotides (`A, T, G, C`) as well as ambiguous bases following the [IUPAC nucleotide code](https://www.bioinformatics.org/sms/iupac.html).
- **mod_type**: The type of DNA modification. Currently, the following three major methylation types are supported:
  - `4mC` (N4-methylcytosine)
  - `5mC` (5-methylcytosine)
  - `6mA` (N6-methyladenine)
- **mod_position**: The position of the modified base within **target_seq** (1-based indexing).
- **window_sizes**: A vector of integers representing window lengths. These windows are generated based on the position of the modification (`mod_position`) and are designed to simulate the recognition sequence lengths of typical **Type II restriction enzymes**.

```r
library(MethREfinder)

# Fetch REBASE methylation sensitivity data
methylation_table <- fetch_rebase_methylation()

# Extract recognition sequence windows from a pattern
windows <- extract_all_windows_RecSeq(
  sequence = "CGAANNNNNNNTARC",
  mod_position = 4,
  mod_type = "6mA",
  direction = "forward",
  window_sizes = c(4, 5, 6)
)

# Match target sequence against REBASE methylation-sensitive enzymes
matched_enzymes <- match_enzyme_sequences(
  mod_type = "6mA",
  target_seq = "CGAANNNNNNNTARC"
)
```
## Quick start
```r
run_MethREfinder(target_seq = "CGAANNNNNNNTARC", mod_type = "6mA", window_sizes = c(4,5,6), direction = "forward", mod_position = 4)
run_MethREfinder(target_seq = "CCWGG", mod_type = "5mC", window_sizes = c(4,5,6), direction = "forward", mod_position = 2)
```

### Advantages
- **Facilitates the selection of restriction enzymes for DNA samples with potential methylation modifications.**
- **Reduces trial-and-error in enzyme selection for molecular cloning and epigenetics.**
- **Integrates REBASE methylation sensitivity data directly into your sequence analysis pipeline.**

## License
**MIT License**

# MethREfinder

MethREfinder is a computational tool designed to identify methylation-sensitive restriction enzyme recognition sites within a given DNA sequence. This package leverages the REBASE database to cross-reference known restriction enzymes and determine their sensitivity to DNA methylation modifications. Researchers can use MethREfinder to assess the compatibility of restriction enzymes with methylated DNA regions, aiding in the selection of appropriate enzymes for molecular cloning, epigenetic studies, and DNA methylation analysis.

## Key Features
- **Identification of restriction enzyme recognition sites within user-provided DNA sequences.**
- **Detection of methylation-sensitive recognition sites based on REBASE methylation data.**
- **Support for common methylation types, including CpG and adenine methylation.**
- **Streamlined and efficient sequence analysis workflow.**
- **Informed enzyme selection for experiments involving methylated DNA.**

## Installation

```r
# Install from GitHub (example, adjust based on your repo)
install.packages("devtools")
devtools::install_github("JAEYOONSUNG/MethREfinder")
```

## Getting Started

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
```


## Dependencies
- **R (>= 4.0.0)**
- **xlsx**
- **httr**

## Advantages
- **Facilitates the selection of restriction enzymes for DNA samples with potential methylation modifications.**
- **Reduces trial-and-error in enzyme selection for molecular cloning and epigenetics.**
- **Integrates REBASE methylation sensitivity data directly into your sequence analysis pipeline.**

## License
**MIT License**

# NGSreads2microsats
Design microsats from NGS reads. 
The pipeline is decribed in [ngs2microsat.md](ngs2microsat.md)


## Summary

This repository present a method to design microsatellites primers from any species using raw NGS reads even at low coverage.
It can also be used to detect microsatellites in any genome. It is centered around [MSATCOMMANDER](https://code.google.com/archive/p/msatcommander/downloads), adding a few convenient cleaning and selecting steps as well as showing that it can be ran easily from raw reads by creating a fast assembly.

It works in 3 steps:

1. Assembly from raw reads
2. [MSATCOMMANDER](https://code.google.com/archive/p/msatcommander/downloads) identification of microsats.
3. Selection of unique and non-overlapping microsatellites

## Dependencies


python 2.7 with the following packages:
* Biopython
* pandas
* argparse

[MSATCOMMANDER](https://code.google.com/archive/p/msatcommander/downloads)

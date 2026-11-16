# scMut Benchmark

A reproducible benchmark of variant calling tools for single-cell RNA sequencing (scRNA-seq) data, covering both 10x Genomics and Smart-seq2 platforms.

This repository provides the analysis code, tool configurations, benchmarking procedures, simulation methods, and consensus-based variant evaluation used in the study.

## Documentation

The analysis is organized into four Bookdown chapters:

1. **Variant Calling** (`01-variant-calling.Rmd`)  
   Commands, parameters, and configurations for variant calling tools across different sequencing platforms and analysis levels.

2. **Benchmark Evaluation** (`02-benchmark-evaluation.Rmd`)  
   Variant filtering, reference set construction, and performance evaluation using RTG Tools, including precision, recall, and F1 score.

3. **Simulations** (`03-simulations.Rmd`)  
   Simulation procedures for assessing variant detection performance under different sequencing depths and variant allele frequencies.

4. **Consensus Analysis** (`04-consensus.Rmd`)  
   Construction and evaluation of consensus variant call sets using different combinations of variant calling tools.
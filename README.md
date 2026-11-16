# scMut benchmark

A four-chapter Bookdown with four calling/filtering workflows and two auxiliary workflows.

- `run_10x_pooled.py`
- `run_10x_celltype.py`
- `run_10x_singlecell.py`
- `run_smartseq.py`
- `run_simulation.py`
- `run_consensus.py`

Edit `config/workflow.json`, inspect a sample run with `--dry-run`, then run the selected callers. Each stage uses the configured conda environment. Workflow commands and original analysis code appear in the Rmd chapters.

Build with `Rscript render.R` from this directory and open `_book/index.html`. Benchmark analysis is not executed when rendering. Figures are the supplied PDFs; PNG files are display copies.

The original Rmd files and the two supplied Smart-seq2 simulation scripts are the method sources. See `metadata/displayed_code.tsv`, `metadata/workflow_sources.tsv` and `metadata/parameter_conflicts.tsv`.

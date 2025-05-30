# muscat.Replication
Results replication of "muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data" article

# Code organisation

## Benchmarking

Code in `R/` and `_targets.R` for the benchmark

### Visualise the workflow

``` r
tar_visnetwork()
```

### Launch benchmark

``` r
tar_make()
```

## Benchmark results analysis

Run scripts in `benchmarkAnalysisR/`

## LPS downstream

Run script in `LPSR/`

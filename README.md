# CSV (Updating)
This project contains the code of Cdf Smoothing via Virtual points (CSV).

Kindly make sure any libraries required for [LIPP](https://github.com/Jiacheng-WU/lipp), [SALI](https://github.com/cds-ruc/SALI) and [ALEX](https://github.com/microsoft/ALEX) are installed.

First make sure any download the neccessary datasets:
```
./data_download.sh

```

Prepare the datasets
```
make
./data_prep

```

Run read-only benchmarks using

```
./benchmark_lipp.sh
./benchmark_sali.sh
./benchmark_alex.sh

```

Run read-write benchmarks using

```
./benchmark_lipp_insert.sh
./benchmark_sali_insert.sh
./benchmark_alex_insert.sh

```

## Acknowledgements

- The benchmarking framework was built on top of the work performed by [Bachfischer](https://github.com/Bachfischer/LogarithmicErrorRegression).

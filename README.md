# Metareport for backupgenes and GraphPrioritizer results

This respotory is designed to load data results, including gold standards and kernel metrics, from specific locations and then process them based on the input argument (`ld` or `rp`). This data dependes from results obtained in [BackupGenes](https://github.com/federedef/BackupGenes.git) and [GraphPrioritizer](https://github.com/federedef/GraphPrioritizer.git) repositories so path have to be adjusted correctly.

## Description

The `daemon.sh` script performs two main tasks depending on the input argument:

1. **Mode `ld`**: Loads result files from a shared results directory into a local directory `data_results`.
2. **Mode `rp`**: Generates a report based on the data loaded in mode `ld` and saves it in a subdirectory with the current date inside `meta_reports`.
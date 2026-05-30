# DA-FAP skeletal muscle degeneration atlas reproducibility repository

This repository accompanies the manuscript:

**Perturbation-guided transcriptomic atlases prioritize degeneration-associated fibro-adipogenic progenitor programs in skeletal muscle**

It contains the analysis scripts, configuration files, curated dataset metadata, compact processed metadata, processed result tables, supplementary tables, and release manifests used for the DA-FAP analysis. The project is based entirely on public transcriptomic datasets.

## Repository contents

- `scripts/`: modular Python analysis scripts for dataset discovery, download/standardization, single-cell integration, FAP subclustering, trajectory, ligand-receptor scoring, GRN-regression perturbation, spatial support, bulk validation, robustness analysis, and manuscript/figure assembly.
- `config/`: workflow configuration.
- `environment/`: Python/R session and package version records.
- `data/metadata/`: dataset inventory and download manifests.
- `data/processed/`: compact processed metadata tables used for reproducibility checks.
- `results/tables/`: analysis tables and supplementary tables.
- `results/figures/`: placeholder directory; full figure files are provided in the manuscript submission package and will be included in the Zenodo archive.
- `release_packages/`: compact release archive for scripts plus large-output manifests.
- `release_packages/large_data_manifest/`: manifest and SHA-256 checksums for large processed objects that are not stored directly in GitHub because of file-size limits.

## Large processed objects

The full processed `.h5ad` objects and large data archive are approximately 2.9 GB compressed. These files exceed ordinary GitHub repository limits and are therefore tracked here by file manifest and SHA-256 checksums. They will be archived on Zenodo upon manuscript acceptance. All large objects can also be regenerated from the public source datasets using the scripts and configuration provided here.

## Reproduction outline

1. Create the software environment from `environment/environment.yml`.
2. Review the dataset inventory in `data/metadata/dataset_metadata.csv`.
3. Run scripts sequentially from `scripts/00_setup/` through the analysis modules.
4. Regenerate figures and tables from `results/` and `scripts/12_figure_generation/` or the revision scripts in `scripts/18_revision/`.

The analysis is computational and hypothesis-generating. GRN perturbation results should be interpreted as candidate regulator prioritization rather than experimental knockout evidence.

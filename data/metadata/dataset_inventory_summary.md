# Phase 1 Dataset Inventory Summary

Generated: 2026-05-29T13:17:32

## Counts

- Total candidate records: 27
- Inclusion decisions: {'include': 12, 'maybe': 12, 'exclude': 3}
- Organisms: {'human': 15, 'mouse': 11, 'rabbit': 1}
- Data types: {'scRNA-seq / snRNA-seq': 1, 'snRNA-seq': 2, 'bulk RNA-seq': 7, 'scRNA-seq': 13, 'spatial transcriptomics': 4}

## Phase Gate

- Included scRNA-seq/snRNA-seq candidates: 5
- Included bulk validation candidates: 6
- Included spatial transcriptomics candidates: 1

Minimum discovery gate status: PASS

Spatial status: available for supportive analysis

## Included Datasets

- CELLxGENE_E-MTAB-13874_human_ageing: scRNA-seq / snRNA-seq (Homo sapiens); Primary human aging atlas with processed h5ad objects and FAP/stromal subsets; best anchor for DA-FAP discovery.
- GSE167186_snRNA: snRNA-seq (Homo sapiens); Human aging/sarcopenia snRNA-seq with usable processed and raw data; directly supports FAP state analysis.
- GSE167186_bulk: bulk RNA-seq (Homo sapiens); Useful paired bulk validation and deconvolution test, but not independent from the GSE167186 single-nucleus study.
- GSE143704: scRNA-seq (Homo sapiens); Human skeletal muscle reference atlas; useful for annotation and baseline FAP/stromal comparison, but not a primary degeneration cohort.
- GSE130646: scRNA-seq (Homo sapiens); Small but valuable human FAP reference and bulk deconvolution signature source.
- GSE164471: bulk RNA-seq (Homo sapiens); Strong independent human bulk RNA-seq aging validation cohort.
- GSE111017: bulk RNA-seq (Homo sapiens); Large human sarcopenia bulk RNA-seq validation resource; phenotype harmonization needed across subseries.
- GSE25941: bulk RNA-seq (Homo sapiens); Independent human aging validation; microarray but well-described sample groups.
- GSE38718: bulk RNA-seq (Homo sapiens); Independent human aging validation cohort; microarray, smaller sample size.
- GSE225766: spatial transcriptomics (Mus musculus); Best currently identified spatial transcriptomics support for fibro-adipogenic degenerative niches; mouse model, so use for mechanistic support.
- GSE138826: scRNA-seq (Mus musculus); Strong mouse regeneration reference with FAP dynamics; useful for trajectory and communication support.
- GSE282285: bulk RNA-seq (Mus musculus); Directly relevant sorted FAP/MuSC bulk RNA-seq support for aging-associated FAP expansion and IL6/SPP1 signaling; mouse mechanistic support.

## Caution

Human spatial transcriptomics specifically for skeletal muscle aging/degeneration was not confirmed in this Phase 1 pass. Current spatial support is mouse dystrophic/regeneration centered. Treat spatial conclusions as cross-species mechanistic support unless a human spatial dataset is added later.

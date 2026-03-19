#### <p align="center">BRAIN-MAGNET: A novel functional genomics atlas coupled with convolutional neural networks facilitates clinical interpretation of disease relevant variants in non-coding regulatory elements</p>

<div align=center><img src="https://github.com/user-attachments/assets/0c847ee6-a48a-43a6-85d8-cec5ed7bf896" width="40%"></div>

Code and resources from BRAIN-MAGNET. For more information check out our [paper](https://doi.org/10.1101/2024.04.13.24305761).

## Content

* [Quick start](https://github.com/ruizhideng/BRAIN-MAGNET/tree/main?tab=readme-ov-file#1-querying-specific-variants)

* [Application of this model](https://github.com/ruizhideng/BRAIN-MAGNET/tree/main?tab=readme-ov-file#2-visualize-the-dataset-in-ucsc)

* [Citation](https://github.com/ruizhideng/BRAIN-MAGNET/tree/main?tab=readme-ov-file#3-citation)


## 1. Quick start

### (1) Querying specific variants

BRAIN-MAGNET predictions for all possible SNPs from NSC NCREs (~100 million), you can easily score your interested variants from our pre-scored data.

#### Install the latest tabix:

```
conda install bioconda::htslib
````

#### Access a specific region from the remote file without downloading the entire large dataset:

````
wget https://huggingface.co/datasets/RuizhiDeng/BRAIN-MAGNET/resolve/main/BRAIN_MAGNET_scores_hg38.txt.bgz.tbi
tabix https://huggingface.co/datasets/RuizhiDeng/BRAIN-MAGNET/resolve/main/BRAIN_MAGNET_scores_hg38.txt.bgz chr1:191050-191060
````

The output contains the following columns:

| chrom | pos | NCRE | Category | TargetGene | cb score | percentile_all score | percentile_each score |

Below is an example of the expected output format:

```
chr17	2594176	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0142849991522962	91.7301440397596	98.1387478849408
chr17	2594177	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0207236691756407	97.3736678634194	98.9847715736041
chr17	2594178	chr17:2593888-2594479	Category_5	PAFAH1B1	0.016062167105265	93.9489634151448	98.6463620981388
chr17	2594179	chr17:2593888-2594479	Category_5	PAFAH1B1	0.00186096147954231	25.9731169853305	48.7309644670051
chr17	2594180	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0269353554851841	99.1836781785742	99.492385786802
chr17	2594181	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0358390458795475	99.8680421037927	100
chr17	2594182	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0272240349091589	99.2280490245149	99.8307952622673
chr17	2594183	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0253191021981183	98.885096409902	99.1539763113367
chr17	2594184	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0271620684134541	99.2186880527201	99.6615905245347
chr17	2594185	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0182490776109626	95.8971137701526	98.8155668358714
chr17	2594186	chr17:2593888-2594479	Category_5	PAFAH1B1	0.0262525200308301	99.0681267062369	99.3231810490694
```

#### If you need to run multiple queries, consider downloading the files locally for faster access.

```
wget https://huggingface.co/datasets/RuizhiDeng/BRAIN-MAGNET/resolve/main/BRAIN_MAGNET_scores_hg38.txt.bgz
wget https://huggingface.co/datasets/RuizhiDeng/BRAIN-MAGNET/resolve/main/BRAIN_MAGNET_scores_hg38.txt.bgz.tbi
```
Once downloaded, you can proceed to score the data as needed:

````
tabix BRAIN_MAGNET_scores_hg38.txt.bgz chr1:191050-191060
````

#### Tabix also offers the -R option, allowing you to score multiple regions efficiently using a BED file.

regions.bed is like below:

```
chr17  2594176  2594151
chr17  2594152  2594157
```

Then use tabix to request the scores from the specific regions:

```
tabix BRAIN_MAGNET_scores_hg38.txt.bgz -R regions.bed
```

#### (2) Visualize the dataset in UCSC

The UCSC tracks of data are available [here](https://genome.ucsc.edu/s/BarakatLab/BrainMagnet_NSC_ESC_cb_scoreshg38).

#### (3) The code generates the figures of the paper

Check the codes of figures: `analysis/NSC_ChIP-STARR-seq.R`

#### (4) How to prioritize the variants from your GS data?

We provide some general recommendations:

(a) Population allele frequency (AF): Verify that your variants are rare by checking their AF in population databases such as gnomAD, GEL or similar resources. Prioritize variants with low (depending on disease frequency) or absent frequency in the general population.

(b) Regulatory region annotation: Determine whether your variants are located within our NCREs atlas. Pay special attention to variants that fall within high category NCREs.

(c) Functional impact prediction: Assess the potential regulatory impact of the variants. Prioritize those with high cb scores or disruption of our predicted crucial motif sites (you could trace the motif sites from our provided [UCSC tracks](https://genome.ucsc.edu/s/BarakatLab/BrainMagnet_NSC_ESC_cb_scoreshg38)).

(d) Gene-Phenotype relevance: Evaluate whether the predicted target genes of the NCREs explain the patient's clinical phenotype.

### 2. Application of this model

Play with the model: `examples/training.ipynb`

Scan motifs: `examples/scan_motifs.ipynb`

Generate cb scores for each nucleotide: `examples/cb_scores.ipynb`

## 3. Prepare training data from an activity table (user-friendly)

You can start from **your own** `Enhancer_activity.txt`-like file:

- **Input format**: first 3 columns are BED (`chrom`, `start`, `end`), followed by **any number of activity columns** (one per cell/assay).
- **Outputs**:
  - `Enhancer.fa` (forward only)
  - `train_set/Sequences_{Train,Valid,Test}.fa`
  - `train_set/Sequences_activity_{Train,Valid,Test}.txt` (same activity columns as your input, or a subset)
  - optional reverse-complement augmentation appended to each split (headers like `>START-END Reversed:`)

### Install

From the repo root:

```bash
python -m pip install .
```

To train the CNN model, install training extras:

```bash
python -m pip install .[train]
```

If you want `prepare_data` to fetch sequences from a genome FASTA (uses `pyfaidx`):

```bash
python -m pip install .[prepare-genome]
```

### Run

```bash
brain-magnet prepare_data \
  --enhancer-activity /path/to/Enhancer_activity.txt \
  --sequences-fasta /path/to/Enhancer.fa \
  --out-dir /path/to/output_folder
```

If you cannot (or don’t want to) install the package, you can run it directly from the repo:

```bash
python3 scripts/prepare_enhancer_data.py \
  --enhancer-activity /path/to/Enhancer_activity.txt \
  --sequences-fasta /path/to/Enhancer.fa \
  --out-dir /path/to/output_folder
```

If you prefer to fetch sequences from a genome FASTA directly, you can pass `--genome-fasta`,
but that requires an extra Python dependency (`pyfaidx`) in your environment.

Optional:

- choose a subset of activity columns to write:

```bash
brain-magnet prepare_data \
  --enhancer-activity /path/to/Enhancer_activity.txt \
  --genome-fasta /path/to/hg38.fa \
  --out-dir /path/to/output_folder \
  --target-cols "NSC_log2_enrichment,ESC_log2_enrichment"
```

The outputs will be written to:

- `output_folder/Enhancer.fa`
- `output_folder/train_set/Sequences_{Train,Valid,Test}.fa`
- `output_folder/train_set/Sequences_activity_{Train,Valid,Test}.txt`

## 4. Training notebook

The notebook `examples/training.ipynb` is parameterized so other users can run it by editing one config cell:

- `TRAIN_SET_DIR`: path to your prepared `train_set/`
- `TARGET_COLUMN`: activity column to predict
- `OUTPUT_DIR`: where checkpoints/plots/predictions are written

## 5. BRAIN-MAGNET helper utilities

For convenience, we include a minimal set of sequence I/O and encoding helpers under:

- `src/brain_magnet/brain_magnet_helper/IOHelper.py`
- `src/brain_magnet/brain_magnet_helper/SequenceHelper.py`

These functions are adapted from [DeepSTARR](https://github.com/bernardo-de-almeida/DeepSTARR/tree/main/DeepSTARR/Neural_Network_DNA_Demo/helper).

## 6. Training CLI (wraps the notebook)

You can run the same training pipeline as `examples/training.ipynb` from the command line:

```bash
brain-magnet train \
  --train-set-dir /path/to/train_set \
  --target-column NSC_log2_enrichment \
  --output-dir /path/to/output_training/NSC_log2_enrichment
```

It will write:

- `output_dir/models/checkpoint_<target-column>.pth`
- `output_dir/metrics.json` and `output_dir/history.json`
- `output_dir/preds_targets/<target-column>/{preds,targets}_<target-column>_{Train,Valid,Test}.npy`

## 7. Citation
```
@article{deng2024dna,
  title={BRAIN-MAGNET: A novel functional genomics atlas coupled with convolutional neural networks facilitates clinical interpretation of disease relevant variants in non-coding regulatory elements},
  author={R. Deng, E. Perenthaler, A. Nikoncuk, S. Yousefi, K. Lanko, R. Schot, M. Maresca, E. Medico-Salsench, L. E. Sanderson, M. J. Parker, W. F.J. van Ijcken, J. Park, M. Sturm, T. B. Haack, Genomics England Research Consortium, G. V. 
  Roshchupkin, E. Mulugeta, T. S. Barakat},
  journal={medRxiv},
  pages={2024.04.13.24305761},
  year={2024},
  doi={10.1101/2024.04.13.24305761}
}
```
rd_APP: please download data at https://figshare.com/s/c577b8b70d2c7c2e8faa

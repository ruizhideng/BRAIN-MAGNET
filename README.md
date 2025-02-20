## BRAIN-MAGNET: A novel functional genomics atlas coupled with convolutional neural networks facilitates clinical interpretation of disease relevant variants in non-coding regulatory elements.

<div align=center><img src="https://github.com/user-attachments/assets/0c847ee6-a48a-43a6-85d8-cec5ed7bf896" width="40%"></div>

For more information check out our [paper](https://doi.org/10.1101/2024.04.13.24305761), and the UCSC tracks of data are available [here](https://genome.ucsc.edu/s/BarakatLab/BrainMagnet_NSC_ESC_cb_scoreshg38).

### 1. Querying specific variants

BRAIN-MAGNET predictions for all possible SNPs from NSC NCREs (~1 billion)

Install the latest tabix:

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
### 2. Train BRAIN-MAGNET on your data


rd_APP: please download data at https://figshare.com/s/c577b8b70d2c7c2e8faa

## BRAIN-MAGNET: A novel functional genomics atlas coupled with convolutional neural networks facilitates clinical interpretation of disease relevant variants in non-coding regulatory elements.

<div align=center><img src="https://github.com/user-attachments/assets/0c847ee6-a48a-43a6-85d8-cec5ed7bf896" width="40%"></div>

BRAIN-MAGNET predictions for all possible SNPs from NSC NCREs (~1 billion)

For more information check out our paper: https://doi.org/10.1101/2024.04.13.24305761

The UCSC tracks of data: https://genome.ucsc.edu/s/BarakatLab/BrainMagnet_NSC_ESC_cb_scoreshg38

### Querying specific variants

Install the latest tabix:

In your current conda environment (might be slow):

```
conda install -c bioconda -c conda-forge htslib=1.18
````

Query a specific region, from the remote file:

````
tabix https://xxx/BRAIN_MAGNET_scores_hg38.txt.bgz chr1:191050-191060
````

The output has the following columns:

| chrom | pos | NCRE | Category | TargetGene | cb score | percentile_all score | percentile_each score |

and your output like below:

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

If you want to do many queries you might want to first download the files locally

```
wget https://xxx/BRAIN_MAGNET_scores_hg38.txt.gz
wget https://xxx/BRAIN_MAGNET_scores_hg38.txt.gz.tbi
```
and then score:

````
tabix BRAIN_MAGNET_scores_hg38.txt.bgz chr1:191050-191060
````

rd_APP: please download data at https://figshare.com/s/c577b8b70d2c7c2e8faa

# plasmid-gc

## Reproducing Enterobacterales population analysis

### Bioinformatics pipeline

The bioinformatics pipeline for this component was run via Nextflow in `main.nf` and `nextflow.config`. It processes a sample of bacterial genomes as follows:

1. Calculate the length and GC content of each contig
2. Run MOB-typer and PlasmidFinder on all plasmids
3. Cluster the plasmids with pling
4. Run mlst and mashtree on all chromosomes

In total, the our sample contained 5,136 contigs, of which 1,454 chromosomes and 3,682 were plasmids. These were located in the `data/assemblies` directory. Also within the `assemblies` directory was the sample sheet `sample.csv`. This lists whether each contig was a chromosome or a plasmid, and was used for feeding the contigs into different processes e.g. all the plasmids into pling and all the chromosomes into mlst. Also necessary is the `scripts` directory, which contains `contig_statistics.py`.

> Versus the original [study](https://doi.org/10.7554/eLife.85302), 4 genomes were dropped because the chromosomes appeared to be fragmented (`389ef33b-8677-40b6-b160-3102d31d3e2b`, `7ecf1e2e-3d06-4974-9de1-52e485d6d84e`, `cd636592-f493-4f64-9c0e-b6356dac56aa`, and `e525c3c4-6baa-42b8-a7c4-bd821511dea5`).

Your file structure should look like this:

```
plasmid-gc/
├── nextflow.config                       # Nextflow config
├── main.nf                               # Nextflow script
├── scripts/
│   ├── contig_statistics.py
│   └── model.R
├── data/                                 # Input data
│   ├── sample.csv
│   └── assemblies/
│       ├── contig1.fasta
│       ├── contig2.fasta
│       └── ...
└── results/                              # Outputs are written in here
```

Before running, ensure that the `nextflow.config` file has the correct paths and resource requirements for your circumstances. To execute the pipeline, make sure you're in `plasmig-gc`, then run

```
nextflow run main.nf
```

### Statistical modelling

The results from the pipeline can then be processed using `scripts/model.R`. 

## Reproducing *Rhizobium leguminosarum* analysis

### Bioinformatics

You first need to download the data from NCBI with

```
datasets download genome accession GCF_000009265.1 --include gbff --filename refseq_gca9265.zip
unzip refseq_gca9265.zip
```

You also need to download [Supplementary File "pnas.2009094117.sd02.xlsx" from Wheatley et. al (2020)](https://pmc.ncbi.nlm.nih.gov/articles/instance/7519234/bin/pnas.2009094117.sd02.xlsx). Place them both in the same directory, alongside the scripts in `scripts/rhizobium_gc.py` and `scripts/rhizobium_model.R`. 

Your file structure should look like this:

```
rhizobium_analysis/
├── ncbi_datasets/
├── pnas.2009094117.sd02.xlsx                     
├── rhizobium_gc.py 
└── rhizobium_model.R
```

You can then run
```
python3 rhizobium_gc.py
```

This will write the file `rhizobium_gc.tsv`.



# Taxon Abundance

## Outline
This pipeline is used to estimate the relative abundances of reads originating from specific taxa in a sequenced library. The taxonomic level of interest can be specified using the `kraken2` labels:

- `O`: Order
- `F`: Family
- `G`: Genus
- `S`: Species
- `S1`: One level below species

The [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy) is used to define the taxa.

## Usage

```
nextflow run BCCDC-PHL/taxon-abundance \
  --fastq_input <fastq_input_dir> \
  --outdir <output_dir>
```

Alternatively, a 'samplesheet.csv' file may be provided with fields ID, R1, R2:

```
ID,R1,R2
sample-01,/path/to/sample-01_R1.fastq.gz,/path/to/sample-01_R2.fastq.gz
sample-02,/path/to/sample-02_R1.fastq.gz,/path/to/sample-02_R2.fastq.gz
...
```

```
nextflow run BCCDC-PHL/taxon-abundance \
  --samplesheet_input </path/to/samplesheet.csv> \
  --outdir </path/to/outdir> 
```

Reads can be binned by taxonomic group, and extracted to separate output files using the `--extract_reads` flag.
When using this flag, a threshold is applied on the percentage of reads assigned to the taxonomic group, below which
reads are not extracted. The default threshold is 1%. It can be modified using the `--extract_reads_threshold` flag.

For example, to extract reads for any taxonomic group present at 0.1% or more:


```
nextflow run BCCDC-PHL/taxon-abundance \
  --fastq_input <fastq_input_dir> \
  --extract_reads \
  --extract_reads_threshold 0.1 \
  --outdir </path/to/outdir> 
```

## Outputs

An output directory will be created for each sample. Within those directories,

```
<sample_id>
├── <sample_id>_fastp.csv
├── <sample_id>_fastp.json
├── <sample_id>_kraken2.txt
├── <sample_id>_S_bracken_abundances.csv
└── <sample_id>_S_top_5.csv
```

If the `--extract_reads` flag is used, a directory named `extracted_reads_by_taxid` is created, with sub-directories named using the NCBI taxonomy ID that those reads were assigned to.

```
<sample_id>
├── extracted_reads_by_taxid
│   ├── 0
│   │   ├── <sample_id>-taxid-0_R1.fastq.gz
│   │   └── <sample_id>-taxid-0_R2.fastq.gz
│   └── 1768
│       ├── <sample_id>-taxid-1768_R1.fastq.gz
│       └── <sample_id>-taxid-1768_R2.fastq.gz
├── <sample_id>_fastp.csv
├── <sample_id>_fastp.json
├── <sample_id>_kraken2.txt
├── <sample_id>_S_bracken_abundances.csv
└── <sample_id>_S_top_5.csv
```

The main outputs are:

1. `<sample_id>_<taxonomic_level>_top_5.csv`

```csv
sample_id,taxonomy_level,abundance_1_name,abundance_1_ncbi_taxonomy_id,abundance_1_num_assigned_reads,abundance_1_fraction_total_reads,abundance_2_name,abundance_2_ncbi_taxonomy_id,abundance_2_num_assigned_reads,abundance_2_fraction_total_reads,abundance_3_name,abundance_3_ncbi_taxonomy_id,abundance_3_num_assigned_reads,abundance_3_fraction_total_reads,abundance_4_name,abundance_4_ncbi_taxonomy_id,abundance_4_num_assigned_reads,abundance_4_fraction_total_reads,abundance_5_name,abundance_5_ncbi_taxonomy_id,abundance_5_num_assigned_reads,abundance_5_fraction_total_reads
DRR161190,S,Mycobacterium novum,2492438,2014826,0.96731,Mycolicibacter sinensis,875328,22690,0.01089,Mycolicibacter terrae,1788,7829,0.00376,Mycolicibacter minnesotensis,1118379,2839,0.00136,Mycolicibacter hiberniae,29314,2236,0.00107
DRR161192,S,Mycobacterium gallinarum,39689,1859062,0.99582,Mycolicibacterium rhodesiae,36814,1630,0.00087,Mycolicibacterium gadium,1794,1620,0.00087,Mycobacterium paragordonae,1389713,1255,0.00067,Mycobacterium conspicuum,44010,379,0.00020
DRR161197,S,Mycobacterium simiae,1784,2068577,0.99913,Mycobacterium cookii,1775,211,0.00010,Mycobacterium mantenii,560555,185,0.00009,Mycobacterium heidelbergense,53376,178,0.00009,Mycolicibacterium fallax,1793,146,0.00007
DRR161199,S,Mycobacterium cookii,1775,2347528,0.99881,Mycobacterium conspicuum,44010,249,0.00011,Mycobacterium tuberculosis,1773,157,0.00007,Mycobacterium kubicae,120959,157,0.00007,Klebsiella michiganensis,1134687,153,0.00007
```

2. `<sample_id>_fastp.csv`

```csv
sample_id,total_reads_before_filtering,total_reads_after_filtering,total_bases_before_filtering,total_bases_after_filtering,q20_bases_before_filtering,q20_bases_after_filtering,q30_bases_before_filtering,q30_bases_after_filtering,adapter_trimmed_reads,adapter_trimmed_bases
DRR161190,4221278,4196834,422127800,413313356,386783970,380201636,359261515,353924424,282820,6371132
DRR161192,3774078,3756172,377407800,369681342,345579006,339764987,321136725,316409787,248262,5936806
DRR161197,4175044,4158344,417504400,408656988,385097573,378269116,359615124,353971246,315786,7178202
DRR161199,4748928,4728502,474892800,463574412,436236542,427545934,406140193,398995239,389594,9277186
```

If the `--versioned_outdir` is used, then a sub-directory will be created below each sample, named with the pipeline name and minor version:

```
<sample_id>
└── taxon-abundance-v0.1-output
    ├── <sample_id>_fastp.csv
    ├── <sample_id>_fastp.json
    ├── <sample_id>_kraken2.txt
    ├── <sample_id>_S_bracken_abundances.csv
    └── <sample_id>_S_top_5.csv
```

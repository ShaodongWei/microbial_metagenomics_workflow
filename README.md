mamba has to be installed also, due to newer version of snakemake 
use '--use-conda --conda-frontend conda', if mamba is problematic to run, e.g. parameter not recognized, then choose to use conda. so far I cannot make mamba work. 

use channel_priority: strict inside yaml

using mamba will lead to --no-default-packages not recognized error. 
This is because Mamba 2.0.0, released last week, no longer supports the --no-default-packages option that Snakemake passes when creating an environment.

we should use snakemake binning --cores 40 --use-conda --conda-frontend conda --resources gpu=2, so that binning can have GPU

# Overview
This pipeline is designed to preprocess microbial metagenomics/isolates genomics sequencing data by applying a series of data manipulation including quality control, host contamination removal, assembling of reads, taxonomy annotation, and functional profiling. 

# Getting Started
## Prerequisites
### Before using the pipeline, ensure you have the following installed:
Snakemake: [Installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) 

conda: [installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) 


# How to Use the Pipeline
## 1. Clone the Repository
First, clone the repository containing the Snakemake pipeline to your local machine:

```bash
git clone https://github.com/ShaodongWei/microbial_metagenomics_workflow.git

```
## 2. Set up the configuration file 

```
input_directory: "/path/to/fastq" # files has to be paired fastq fields, with format of file_1.fastq or file_2.fastq
output_directory: "/path/to/output"
host_genome: "/path/to/host_genomic.fasta"
threads: 40
assembler: "megahit"  # Options: spades or megahit
assembler_mode: 'meta' # if use spades please choose 'meta' or 'isolate'; if use assembler megahit, it is only for metagenomics. 
run_binning: true # use 'true' or 'false'. The supported binner is SemiBin2. 
metaphlan_db: '/path/to/metaphlan_db' # directory to donwload metaphlan4 database. 
functional_profiling: true # 'true' (using HUMAnN) or 'false'
chocophlan_db: '/data/chocophlan_db' # chocophlan database for humann
```

## 3. Run the pipeline using snakemake
```
conda config --set channel_priority flexible # Set channel priority to be flexible for conda
```
### Run the entire pipeline 
```
snakemake --cores threads_number --use-conda --conda-frontend conda # Conda will install dependencies automatically. All steps will be executed sequentially. 
```

### Run a specific step 
```
snakemake --list # show all steps

snakemake step_name --cores threads_number --use-conda --conda-frontend conda # Run a specific step 

snakemake --cores threads_number --use-conda --conda-frontend conda --dryrun # using dryrun to check what steps will be executed
```
## 4. Steps in the workflow 
### Step 1, quality control
This step is to quality control your single raw fastq files, including quality filtering, adapter trimming, etc. Details can be found [here](https://github.com/OpenGene/fastp)
```
snakemake fastp_trim --cores threads_number --use-conda --conda-frontend conda
```
### Step 2, remove host DNA contamination
This step is to remove the potential host DNA contamination using bowtie2. 
```
snakemake remove_host --cores threads_number --use-conda --conda-frontend conda
```
### Step 3, assembling of short reads into contigs
This step is to assemble (single sample) the short reads using either spades (meta or isolate mode) or megahit.  You have to specify the assembler in the config.yaml file. 
```
snakemake spades_assembly --cores threads_number --use-conda --conda-frontend conda
or
snakemake megahit_assembly --cores threads_number --use-conda --conda-frontend conda
```
### Step 4, binning
This step is to group assembled contigs into bins (species) in each sample using SemiBin2. You have to specify if to run this optional step. 
```
snakemake binning --cores threads_number --use-conda --conda-frontend conda
```
### Step 5, taxonomy annotation
This step is to annotate the taxonomic composition of the samples using metaphlan4, it returns both taxonomy and abundances. You need to specify the metaphlan database location in the configure file. 
```
snakemake taxonomy --cores threads_number --use-conda
```
### Step 6, functional profiling 
This step uses HUMAnN3 to annotate the metabolic function and the pathways involved in the sample. You have to specify if to run this optional step. 
```
snakemake functional_profiling --cores threads_number --use-conda
```
## 5. Troubleshooting

### 5.1 Directory locked
This normally happens when you run snakemake multiple times, you can unlock it by deleting the lock files 
```
rm .snakemake/locks
```
### 5.2 Run a specific step, without dependency steps
For example, if you want to run mapping specifically, without running all previous steps, you can use touch to update the timestamp of the .done file for the previous step which is quality_control rule to let Snakemake think it’s up-to-date and previous step has been done. 
```
touch output/.fastp_trim.done
```


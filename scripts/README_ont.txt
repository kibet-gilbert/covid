Analysis of ONT data:

Requirements:
1. Fast5 files
2. Flowcell type (e.g FLO-MIN106)
3. Sequencing Ligation kit (e.g SQK-LSK109)
4. Barcoding kit (e.g EXP-NBD196 EXP-NBD104 EXP-NBD114 SQK-NBD110-24 SQK-NBD110-96)

Software:
1. Basecalling (e.g guppy-basecaller)
2. Demultiplexing (e.g guppy-barcoder)
3. Downstream analysis pipeline: includes variant-calling pipelines, Genome assembly and ...

Harware:
For basecalling/demultiplexing a PC/HPC with GPU (specific types) or just use CPUs (slower)

Running scripts:
The scripts described here were optimized for ilri HPC:

Basecalling: basecallgpu.5.0.11.hpc.sbatch
	To run it "cd" into the directory containing the "fast5/" directory 
	$ sbatch path-to-script/basecallgpu.5.0.11.hpc.sbatch

variant calling: covid_run-2.2.nanopore.sbatch
	Prerequisites: Install nextflow and singularity. You can opt to run the pipeline with other package managers e.g conda, Docker, Podman and others supported by nextflow. Here we use singularity.
	Setting up: Run `nextflow run nf-core/viralrecon -r 2.2 -profile test,singularity` - downloads and installs all singularity images .
	To analyse your data: cd into the directory containing sub-directories: fast5, basecall, demultiplex.
	$ sbatch -w <node> path-to-script/covid_run-2.2.nanopore.sbatch

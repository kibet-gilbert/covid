# Summary of nf-core/viralrecon Data Analysis.
============================================

## 1. To analyse some datasets
---------------------------
1. Open the script path-to/covid/scripts/covid_run.sbatch and set the variables: FRead_suffix and RRead_suffix to your desired forward read and reverse read suffix. The defaults are FRead_suffix="L001_R1_001.fastq.gz" and RRead_suffix="L001_R2_001.fastq.gz" 
2. "cd" to path-to/covid/data/<your_data_dir> and execute the sbatch command below: $ sbatch -w compute05 ../../scripts/covid_run.sbatch

N/B: compute05 should work okay i.e 48CPUs and 1.5TB disk space.

The nf-core/viralrecon, a nextflow pipeline has been set up to run with -profile singularity option. For more on options that can be set see PARAMETER configuration section.

Help message: Execute the following line for the usage:

$ nextflow run nf-core/viralrecon -r 1.1.0 --help -profile singularity

This workflow is made up of 64 process(The equivalent of functions in bash/python): some can be deactivated as need. See Processes section for more.

## 2. Organization of the covid analysis directory:
------------------------------------------------

The covid directory is structured into four directories: data,scripts,viralrecon and work.(see more below and in the organization section)

2.1. data
---------
The data dir has: core_data, covid<date>, test_data*. The core_data includes all necessary data files like reference genome, Artic primer data, gff files.... The covid<data> dirs contains our sequences and the analysis results from them (more below). 

2.1.1. covid<date>
------------------
PS: Naming will change to <date>covid in the next runs.
This dir has all symbolic links to fastq.gz sequencing data. It has results from the sequence analysis as well in diferent directories/files:
Files
-----
2.1.1.1 <inputfile>.mutations.tsv	: symplified vcf as a tab-separeted file of all mutations in all samples
2.1.1.2 <inputfile>_aa.mutations.tsv	:Symplified amino acid mutaion file
2.1.1.3 <inputfile>_analysis.*		: Any other custom result file excel/pdf
Directories
-----------
alignment	: Musltiple Sequence Alignment (MUSCLE) results dir
mutaions	: Symplified VCF files in a tab-separeted format
nextclade	: nextclade.js analysis results	
phylogeny	: phylogeny (iqtree) analysis results
pangolin	: Pangolin analysis results
results		: Variant calling results - Contains more dirs (see next)
2.1.1.4 results
---------------
assembly	: De-novo assembly results - don't have any results when skipped
multiqc		: multiqc results in .html and data.yaml format
pipeline_info	: All pipeline execution reports in .txt and html format
preprocess	: Pre-processing results: fastQC and fastp (QC & trimming)
variants	: All variant call results. 
	ivar	- All variant call files - Not annotated
	ivar/snpeff - All variant call files annotated by sneff 
	ivar/quast - consensus sequence summary stats by quast
	bam - bam idex files sorted.bam & sorted.bai
	bam/samtoos_stats - Samtool stats
	bam/picard_metrics - picard stats on dublicates in the reads
	bam/mosdepth - mosdepth stats on amplicon and genome coverage
	bam/mosdepth/amplicon/plots - log Plots on how amplicons/primers output looks (pdf) and an overall heatmap (pdf). All have accompanying data in TSV format
	bam/mosdepth/genome/plots -log plots and data for log genome coverage vs position/loci on the genome.

2.2. scripts
------------
The scripts dir has the analysis sbatch scripts and auxilliary scripts for pre-analysis and post-analysis processing: covid_run_<unique-feature>.sbatch

2.3. viralrecon
---------------
The viralrecon dir has been git and has the source code of the pipeline... The idea is to eventually run the pipeline directly from the code when needed (not possible as of now).

2.4. work
---------
The work dir has the conda environment (work/conda/<env>) which can be copied to the working dir (/var/scratch/${USER}/work/conda/<env>)

Organization.
============
.(covid)
├── data
│   ├── core_data
│   ├── covid_01-04-2021
│   │   ├── alignment
│   │   │   ├── covid_01-04-2021_con.aln
│   │   │   ├── covid_01-04-2021_con.fasta
│   │   │   └── covid_01-04-2021_headers
│   │   ├── covid_01-04-2021_aa.mutations.tsv
│   │   ├── covid_01-04-2021_analysis.mutations.pdf
│   │   ├── covid_01-04-2021_analysis.mutations.xlsx
│   │   ├── covid_01-04-2021.mutations.tsv
│   │   ├── covid_01-04-2021_pangolin.tsv
│   │   ├── mutations
│   │   │   ├── COVC21058_S21.snpEff.vcf.tsv
│   │   │   |
│   │   │   |
│   │   │   |
│   │   │   ├── COVC23453_S20.snpEff.vcf.tsv
│   │   │   ├── MoH-Cov-6_S6.snpEff.vcf.tsv
│   │   │   └── Undetermined_S0.snpEff.vcf.tsv
│   │   ├── nextclade
│   │   │   ├── covid_01-04-2021_con.json
│   │   │   ├── covid_01-04-2021_con.tsv
│   │   │   ├── nextstrain_.svg
│   │   │   ├── nextstrain_tree.nexus
│   │   │   └── nextstrain_tree.nwk
│   │   ├── phylogeny
│   │   │   ├── covid_01-04-2021_con.bionj
│   │   │   ├── covid_01-04-2021_con.ckp.gz
│   │   │   ├── covid_01-04-2021_con.contree
│   │   │   ├── covid_01-04-2021_con.iqtree
│   │   │   ├── covid_01-04-2021_con.mldist
│   │   │   ├── covid_01-04-2021_con.model.gz
│   │   │   ├── covid_01-04-2021_con.splits.nex
│   │   │   ├── covid_01-04-2021_con.treefile
│   │   │   └── covid_01-04-2021_con.ufboot
│   │   ├── results
│   │   │   ├── assembly
│   │   │   ├── multiqc
│   │   │   ├── pipeline_info
│   │   │   ├── preprocess
│   │   │   └── variants
│   │   │       ├── bam
│   │   │       │   ├── log
│   │   │       │   ├── mosdepth
│   │   │       │   │   ├── amplicon
│   │   │       │   │   │   └── plots
│   │   │       │   │   └── genome
│   │   │       │   │       └── plots
│   │   │       │   ├── picard_metrics
│   │   │       │   └── samtools_stats
│   │   │       └── ivar
│   │   │           ├── bcftools_stats
│   │   │           ├── consensus
│   │   │           │   └── base_qc
│   │   │           ├── log
│   │   │           ├── quast
│   │   │           │   └── AF0.75
│   │   │           │       ├── aligned_stats
│   │   │           │       ├── basic_stats
│   │   │           │       ├── contigs_reports
│   │   │           │       │   └── minimap_output
│   │   │           │       ├── genome_stats
│   │   │           │       └── icarus_viewers
│   │   │           └── snpeff
│   │   ├── samplesheet.csv
│   │   └── slurm-716659.out
│   ├── covid_11-02-2021
│   ├── covid-run3_09-04-2021
│   ├── general_variants
│   ├── test_data
│   ├── test_data00
│   └── test_data01
├── README
├── scripts
│   ├── aa_codes
│   ├── covid_env.yml
│   ├── Covid_exporatory.ipynb
│   ├── covid_run_denovo.sbatch
│   ├── covid_run_kranken.sbatch
│   ├── covid_run_nomqc.sbatch
│   ├── covid_run.sbatch
│   ├── covid_run_test.sbatch
│   ├── Notes.txt
│   ├── pangolin
│   ├── process_files_bak.sh
│   └── process_files.sh
├── viralrecon
│   ├── assets
│   ├── bin
│   ├── CHANGELOG.md
│   ├── CITATIONS.md
│   ├── CODE_OF_CONDUCT.md
│   ├── conf
│   ├── Dockerfile
│   ├── docs
│   ├── environment.yml
│   ├── LICENSE
│   ├── main.nf
│   ├── nextflow.config
│   └── README.md
└── work
    └── conda

## Processes:
=========
1. GUNZIP_FASTA
2. GUNZIP_GFF
3. CHECK_SAMPLESHEET
4. *SRA_FASTQ_FTP
5. *SRA_FASTQ_DUMP
6. CAT_FASTQ
7. FASTQC
8. FASTP
9. BOWTIE2_INDEX
10. MAKE_SNPEFF_DB
11. BOWTIE2
12. SORT_BAM
13. IVAR_TRIM
14. PICARD_MARKDUPLICATES
15. PICARD_METRICS
16. MOSDEPTH_GENOME
17. MOSDEPTH_AMPLICON
18. MOSDEPTH_AMPLICON_PLOT
19. SAMTOOLS_MPILEUP
20. *VARSCAN2
21. *VARSCAN2_CONSENSUS
22. *VARSCAN2_SNPEFF
23. *VARSCAN2_QUAST
24. IVAR_VARIANTS
25. IVAR_CONSENSUS
26. IVAR_SNPEFF
27. IVAR_QUAST
28. *BCFTOOLS_VARIANTS
29. *BCFTOOLS_CONSENSUS
30. *BCFTOOLS_SNPEFF
31. *BCFTOOLS_QUAST
32. *BCFTOOLS_ISEC
33. **MAKE_BLAST_DB
34. **SPADES
35. **SPADES_BLAST
36. **SPADES_ABACAS
37. **SPADES_PLASMIDID
38. **SPADES_QUAST
39. **SPADES_VG
40. **SPADES_SNPEFF
41. **METASPADES
42. **METASPADES_BLAST
43. **METASPADES_ABACAS
44. **METASPADES_PLASMIDID
45. **METASPADES_QUAST
46. **METASPADES_VG
47. **METASPADES_SNPEFF
48. **UNICYCLER
49. **UNICYCLER_BLAST
50. **UNICYCLER_ABACAS
51. **UNICYCLER_PLASMIDID
52. **UNICYCLER_QUAST
53. **UNICYCLER_VG
54. **UNICYCLER_SNPEFF
55. **MINIA
56. **MINIA_BLAST
57. **MINIA_ABACAS
58. **MINIA_PLASMIDID
59. **MINIA_QUAST
60. **MINIA_VG
61. **MINIA_SNPEFF
62. get_software_versions
63. MULTIQC
64. output_documentation

PARAMETER configurations:
  // Options: Generic
  input = './samplesheet.csv'
  protocol = 'amplicon'
  amplicon_fasta = ${DATADIR}/GCA_009858895.3_ASM985889v3_genomic.200409.fna.gz
  amplicon_bed = ${DATADIR}/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz

  // Options: SRA download
  save_sra_fastq = false
  skip_sra = false

  // Options: Reference genomes
  genome = false
  save_reference = false

  // Options: Read Trimming
  cut_mean_quality = 20
  qualified_quality_phred = 20
  unqualified_percent_limit = 10
  min_trim_length = 50
  skip_adapter_trimming = false
  skip_amplicon_trimming = false
  save_trimmed = false

  // Options: Variant calling
  callers = 'varscan2,ivar,bcftools'
  min_mapped_reads = 1000
  ivar_trim_noprimer = false
  ivar_trim_min_len = 20
  ivar_trim_min_qual = 20
  ivar_trim_window_width = 4
  filter_dups = false
  filter_unmapped = false
  mpileup_depth = 0
  min_base_qual = 20
  min_coverage = 10
  min_allele_freq = 0.25
  max_allele_freq = 0.75
  varscan2_strand_filter = true
  amplicon_left_suffix = '_LEFT'
  amplicon_right_suffix = '_RIGHT'
  save_align_intermeds = false
  save_mpileup = false
  skip_markduplicates = false
  skip_picard_metrics = false
  skip_mosdepth = false
  skip_snpeff = false
  skip_variants_quast = false
  skip_variants = false

  // Options: QC
  skip_fastqc = false
  skip_multiqc = false

  // Boilerplate options
  outdir = './results'
  publish_dir_mode = 'copy'
  name = false
  multiqc_config = false
  email = false
  email_on_fail = false
  max_multiqc_email_size = 25.MB
  plaintext_email = false
  monochrome_logs = false
  help = false
  tracedir = "${params.outdir}/pipeline_info"
  custom_config_version = 'master'
  custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
  hostnames = false
  config_profile_description = false
  config_profile_contact = false
  config_profile_url = false

  // Options: Kraken2
  kraken2_db = 'https://zenodo.org/record/3738199/files/kraken2_human.tar.gz'
  kraken2_db_name = 'human'
  kraken2_use_ftp = false
  save_kraken2_fastq = false
  skip_kraken2 = false

  // Options: De novo assembly
  assemblers = 'spades,metaspades,unicycler,minia'
  minia_kmer = 31
  skip_blast = false
  skip_abacas = false
  skip_plasmidid = false
  skip_vg = false
  skip_assembly_quast = false
  skip_assembly = false


  // Defaults only, expecting to be overwritten
  max_memory = 40.GB
  max_cpus = 30
  max_time = 240.h

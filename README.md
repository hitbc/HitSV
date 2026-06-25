# HitSV: Maximizing discovery of structural variants across sequencing technologies

## Overview

HitSV (Original name: gcSV) is a comprehensive structural variant detection method that supports multiple sequencing data types and analysis modes, enabling accurate detection and genotyping of structural variants in genomes.

Key Features:
- Supports structural variant detection and genotyping using pure long-read sequencing (LRS) data
- Supports structural variant detection and genotyping using pure short-read sequencing (SRS) data
- Supports hybrid structural variant detection by combining low-depth LRS data with SRS data
- Precisely reconstructs inserted sequences through local sequence assembly
- Supports structural variant detection from ultra-low depth (less than 5X) third-generation sequencing data
- Supports multiple sequencing platforms, including ASM (T2T), ONT, PacBio, Illumina, BGI-T7, etc.
- Provides strategies for complex structural variant (CSV) analysis of population samples

The structural variant (SV) detection results of HitSV tool on HG002/3/4/5/6/7 and the 1000 Genomes Project datasets are stored at https://github.com/hitbc/HitSV_call_results

Table of Contents:

* [Introduction](#introduction)
* [Installation](#installation)
  * [1. Directly Obtain and Deploy Precompiled Static Binary](#1-directly-obtain-and-deploy-precompiled-static-binary)
  * [2. Compile HitSV from Source Code](#2-compile-gcsv-from-source-code)
* [Usage](#usage)
  * [1. Detect and Genotype Structural Variants Using Pure LRS Data](#1-detect-and-genotype-structural-variants-using-pure-lrs-data)
    * [Whole Genome SV Detection](#whole-genome-sv-detection)
    * [Chromosome 1 SV Detection](#chromosome-1-sv-detection)
    * [Specific Region SV Detection (chr2: 11,500,000–18,500,000)](#specific-region-sv-detection-chr2-1150000018500000)
  * [2. Detect and Genotype Structural Variants Using Pure SRS Data](#2-detect-and-genotype-structural-variants-using-pure-srs-data)
    * [Whole Genome Detection](#whole-genome-detection)
    * [Specific Region SV Detection (chrX: 11,500,000–18,500,000)](#specific-region-sv-detection-chrx-1150000018500000)
  * [3. Hybrid SV Detection Combining Low-depth LRS Data and SRS Data](#3-hybrid-sv-detection-combining-low-depth-lrs-data-and-srs-data)
  * [4. Multi-process Processing](#4-multi-process-processing)
    * [Divide the Entire Genome into Different Regions for Simultaneous Variant Detection](#divide-the-entire-genome-into-different-regions-for-simultaneous-variant-detection)
    * [Merge Variant Detection Result Files into a Complete File](#merge-variant-detection-result-files-into-a-complete-file)
  * [5. LRS error mode preset for local-assembly](#5-lrs-error-mode-preset-for-local-assembly)
  * [6. Contig-to-Reference Alignment Parameters (Applicable to Both LRS and SRS)](#6-contig-to-reference-alignment-parameters-applicable-to-both-lrs-and-srs)
  * [7. Output Parameters](#7-output-parameters)
  * [8. Local Assembly for Specified Intervals](#8-local-assembly-for-specified-intervals)
    * [Generate VCF Format Results](#generate-vcf-format-results)
    * [Generate BED Format Results](#generate-bed-format-results)
  * [9. Structural Variant Analysis for Population Samples](#9-structural-variant-analysis-for-population-samples)
    * [9.1 Population Sample Complex Structural Variant (CSV) Clustering Analysis](#91-population-sample-complex-structural-variant-csv-clustering-analysis)
    * [9.2 Population Sample Complex Structural Variant Region - MEI-TR Array Detection](#92-population-sample-complex-structural-variant-region-mei-tr-array-detection)
    * [9.3 Population Sample Complex Structural Variant Region - Nested Variant Detection](#93-population-sample-complex-structural-variant-region-nested-variant-detection)
* [Demo](#demo)
* [Changelog](#Changelog)
* [License](#license)


## Introduction

- HitSV is a sequence alignment-based structural variant detection tool (taking BAM files as input).
- Detects and genotypes structural variants using pure LRS (long-read sequencing) data or assembly result data (ASM)
- Detects and genotypes structural variants using pure SRS (short-read sequencing) data
- Performs hybrid structural variant detection and genotyping by combining (low-depth) LRS data with SRS data
- Precisely reconstructs inserted sequences through local sequence assembly
- Supports structural variant detection from ultra-low depth (less than 5X) third-generation sequencing data

HitSV supports data from common second-generation sequencing platforms (Illumina, BGI-T7, etc.) as well as data from common third-generation sequencing platforms, including ONT, PacBio sequencing data, or global assembly result data (such as ASM, T2T).
HitSV demonstrates high detection rates and genotyping accuracy when using high-quality third-generation sequencing data (such as ONT-Q26, PacBio HiFi, etc.); at the same time, HitSV can also analyze third-generation sequencing data with higher error rates (error rates of approximately 3%-5%).

## Installation
HitSV provides precompiled static linked binary files. Use the following commands for quick deployment and execution:

### 1. Directly Obtain and Deploy Precompiled Static Binary

```bash
get latest release at "https://github.com/hitbc/HitSV/releases" 
./HitSV --help
```

### 2. Compile HitSV from Source Code

```bash
git clone https://github.com/hitbc/HitSV/
cd ./Release
make clean
make all -j 12
./HitSV --help
```

## Usage

HitSV supports multiple data inputs for SV detection, including pure LRS data, pure SRS data, or a combination of (low-depth) LRS data and SRS data.

In addition, you can quickly obtain SV detection results for specific regions by setting the local region to be detected. HitSV does not provide built-in multi-threaded parallelization; users can divide entire chromosomes into multiple different regions and run multiple HitSV instances in parallel to achieve multi-process parallelization.

For different types of third-generation sequencing data, HitSV provides different preset parameter settings to optimize detection and genotyping accuracy. "ONT_Q20" and "HIFI" can be used for data with error rates below 1%, and "ERR_PRONE" can be used for data with error rates between 1% and 6%.

### 1. Detect and Genotype Structural Variants Using Pure LRS Data

Introduce LRS data or assembly result data (ASM) for detection by setting the -l parameter.

#### Whole Genome SV Detection

```bash
HitSV call -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### Chromosome 1 SV Detection

Parameters -S and -E represent the start and end chromosome IDs, respectively. Note that chromosome IDs are 0-based.

```bash
HitSV call -S 0 -E 0 -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### Specific Region SV Detection (chr2: 11,500,000–18,500,000)

```bash
HitSV call -S 1 -E 1 -s 11500000 -F 18500000 -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

Note: Chromosome IDs and positions are 0-based. By default, whole-genome SV detection is performed.

### 2. Detect and Genotype Structural Variants Using Pure SRS Data

Before analysis, we recommend using `HitSV srs_fa_stat` to pre-calculate local repeat complexity information for the reference genome. Pre-built results for common human reference genomes grch38 and hg37d5 are available at [GRCh38.stat.txt.gz](https://github.com/hitbc/HitSV/blob/main/demo_HitSV/GRCh38.stat.txt.gz) and [hs37d5.stat.txt.gz](https://github.com/hitbc/HitSV/blob/main/demo_HitSV/hs37d5.stat.txt.gz), which need to be decompressed before use.

When conducting large-scale sequence research or analyzing complete genomes, use `HitSV srs_trans_reads` to preprocess SRS data to reduce I/O requirements and accelerate analysis; this command does not affect variant detection results but significantly improves analysis speed.

#### Whole Genome Detection

```bash
HitSV srs_fa_stat ref.fa > ref.stat.txt
HitSV srs_trans_reads ref.fa sample.SRS.bam TL.bam 
samtools sort --output-fmt=BAM -o TL.sort.bam TL.bam
samtools index TL.sort.bam
HitSV call -n sample.SRS.bam -L TL.sort.bam -r ref.fa -I ref.stat.txt -o output.vcf 2> /dev/null
```

#### Specific Region SV Detection (chrX: 11,500,000–18,500,000)

```bash
HitSV call -S 22 -E 22 -s 11500000 -F 18500000 -n sample.SRS.bam -I ref.stat.txt -r ref.fa -o output.vcf 2> /dev/null
```

### 3. Hybrid SV Detection Combining Low-depth LRS Data and SRS Data

Provide both LRS and SRS datasets in a single HitSV-call to leverage their complementary advantages:

```bash
HitSV call -l sample.LRS.bam -n sample.SRS.bam -L TL.sort.bam -r ref.fa -I ref.stat.txt -o output.vcf 2> /dev/null
```

### 4. Multi-process Processing

HitSV does not directly support multi-threaded variant detection. Multi-process parallelization can be achieved by dividing chromosomes into different regions and running multiple HitSV instances in parallel.

The script HitSV_multy_process.sh provides an example showing how to divide the entire genome into 24 different regions and run 24 HitSV instances in parallel to achieve multi-process parallel variant detection.

```bash
bash HitSV_multy_process.sh sample.LRS.bam ./WORK_DIR ONT_Q26 ref.fa
```

Using the following multi-process processing method:

#### Divide the Entire Genome into Different Regions for Simultaneous Variant Detection

```bash
HitSV call -S 0 -E 0 -l sample.LRS.bam -r ref.fa -o chr1.vcf 2> /dev/null
.....
HitSV call -S 22 -E 22 -H -l sample.LRS.bam -r ref.fa -o chrX.vcf 2> /dev/null
HitSV call -S 23 -E 23 -H -l sample.LRS.bam -r ref.fa -o chrY.vcf 2> /dev/null
```

#### Merge Variant Detection Result Files into a Complete File

```bash
cat chr1.vcf > output.vcf
.....
cat chrX.vcf >> output.vcf
cat chrY.vcf >> output.vcf
```

This method is applicable to all analysis types (pure LRS, pure SRS, hybrid).

### 5. LRS error mode preset for local-assembly

- `--LRS_preset ONT_Q20` or `-p ONT_Q20`: Used to preset the type of input third-generation sequencing data. The default value is high-quality ONT dataset (ONT_Q20), and other available preset parameters include "ASM", "HIFI", and "ERR_PRONE". "ASM" is used for assembled result data, "ONT_Q20" and "HIFI" can be used for data with error rates below 1%, and "ERR_PRONE" can be used for data with error rates between 1% and 10%.

### 6. Contig-to-Reference Alignment Parameters (Applicable to Both LRS and SRS)

--Contig_aln_preset asm5 or -P asm5: Specifies the preset alignment parameters for contig-to-reference realignment via ksw2. This setting controls alignment sensitivity according to the expected sequence divergence between the sample and the reference genome, and is used for both long-read (LRS) and short-read (SRS) data analysis. Available presets are:

asm5 (default): uses -A1 -B19 -O39,81 -E3,1. Optimized for species with low sequence divergence, such as human.

asm10: uses -A1 -B9 -O16,41 -E2,1. Suitable for most species, including soybean, mouse, and monkey.

asm20: uses -A1 -B4 -O6,26 -E2,1. Recommended for species with high sequence divergence, such as zebrafish or maize.
The default value is asm5.

### 7. Output Parameters

- `-B, --output_format` : Specifies the format of the SV callset. Accepted values are VCF (default) and PURE_STR. The VCF option produces a standard VCF v4.2 file with HitSV‑specific annotations in the INFO field. The PURE_STR option writes a simplified tab‑delimited table containing only the core SV coordinates and types, which may be more convenient for custom downstream parsing.

- `-f, --random_phasing` : Enables random phasing of heterozygous SVs. When set to 1 (default), unphased heterozygous genotypes (0/1) are randomly resolved to either 0|1 or 1|0, while homozygous alternative genotypes (1/1) are output as 1|1. All heterozygous variants residing on the same locally assembled contig are phased consistently. Set to 0 to output unphased genotypes without random assignment.

- `-m, --MIN_SV_len` : Minimum length (in bp) for a structural variant to be reported. Variants shorter than this threshold are discarded.

- `-v, --Small_var` : Controls the output of small variants (SNPs and INDELs) flanking detected SVs. When enabled (1, default), small variants that can be unambiguously phased to the same local haplotype as the SV are included in the output. This option is automatically ignored (treated as 0) in ERR_PRONE mode, where the high base‑calling error rate makes reliable small‑variant calling infeasible. Set to 0 to suppress small‑variant output in all modes.

- `-o, --output` : Path to the output file. If not specified, results are written to standard output (stdout).

- `-H, --no-header` : Suppress the VCF header (meta‑information lines and column header line) in the output. This flag is useful when concatenating results from multiple runs or when the header is not required for subsequent analysis.

Note on ERR_PRONE mode: Because the -v option relies on accurate base‑level alignments, it is automatically disabled when using -p ERR_PRONE.

### 8. Local Assembly for Specified Intervals

When specifying intervals with the following parameters (-b --FC_BED), HitSV generates local assembly sequences for specific intervals. This function only applies to ASM and LRS datasets, not to SRS datasets, and the length of each specified interval cannot exceed 50,000 bp.

#### Generate VCF Format Results

The generated results are stored in VCF format, with contigs stored in the INFO field:

```bash
HitSV call -b region.bed -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### Generate BED Format Results

The generated local contig results are stored in BED format for analysis of completed assembly data:

```bash
HitSV call -B BED -b region.bed -l sample.ASM.bam -r ref.fa -o output.bed 2> /dev/null
```

### 9. Structural Variant Analysis for Population Samples

This function is used to analyze the structural variant structure of population samples within specific chromosomal intervals. It is particularly designed to address the detection of complex structural variants (CSV). HitSV directly analyzes the local haplotype sequence information (local contigs) of population samples instead of being based on VCF, avoiding the loss of some details when generating VCF, thereby better detecting and analyzing the complex structural variant composition of population samples. This function only applies to ASM and LRS datasets, not to SRS datasets, and the length of each specified interval cannot exceed 50,000 bp.

#### 9.1 Population Sample Complex Structural Variant (CSV) Clustering Analysis

General steps:

1. Perform single-sample variant detection for all samples;
2. Based on the single-sample variant detection results, identify genomic intervals of interest;
3. Use the (-b --FC_BED) parameter to perform local assembly on each sample in each genomic interval of interest to generate local contigs;
4. Annotate the structure of all local contigs using the TRF algorithm and repeat Masker algorithm;
5. Based on the clustering algorithm for similar contigs, cluster similar contigs together and analyze the complex structural variant composition of population samples.

Refer to "cohort_csv_analysis.md" for the complete workflow.

#### 9.2 Population Sample Complex Structural Variant Region - MEI-TR Array Detection

1. Re-call local contig sequences for each sample based on BED files;
2. Contig sequence RM and TRF annotation;
3. SV-TR array detection;

Refer to "py/mei_tr_array_single_region_detection.py" for the complete workflow.

Based on workflow 9.1, obtain the variant annotation results from repeat Masker: [region].fa.out, then execute:

```bash
python3 ./py/mei_tr_array_single_region_detection.py region.fa.out
```

#### 9.3 Population Sample Complex Structural Variant Region - Nested Variant Detection

1. Re-call local contig sequences for each sample based on BED files;
2. Contig sequence RM and TRF annotation;
3. Nested variant detection;

Obtain nested variant detection results based on the result file "nestedStructures.log" from workflow 9.1.

## Demo

Demo data provides a fragment (1:869000-870000) from GIAB's HG002 real dataset, including demo Pacbio HiFi dataset (average 10x depth) and demo Illumina data (average 60x depth).

Due to the small size of the demo data, the insertion fragment size distribution of the Illumina dataset cannot be effectively calculated. Therefore, an additional file `SRS_HG002_stat.json` is provided to indicate the depth and insertion fragment size distribution statistics of the Illumina dataset (-T SRS_HG002_stat.json). In general usage scenarios (such as WGS), there is no need to provide this file or parameter.

The parameters `-S 0 -E 0 -s 0 -F 1000000` limit the actual region for variant detection to 1:1-1000000.

```bash
cd ./demo_HitSV/
# Obtain reference analysis files
zcat hs37d5.stat.txt.gz > hs37d5.stat.txt
# LRS SV detection
HitSV call -S 0 -E 0 -s 0 -F 1000000 -r hs37d5_1_0_1000000.fa -o lrs_demo.vcf -p HIFI -l LRS_HG002_1_869000_870000_10X_demo.bam
# SRS SV detection
HitSV call -S 0 -E 0 -s 0 -F 1000000 -T SRS_HG002_stat.json -r hs37d5_1_0_1000000.fa -o srs_demo.vcf -I hs37d5.stat.txt -n  SRS_HG002_1_869000_870000_60X_demo.bam
# Hybrid SV detection
HitSV call -S 0 -E 0 -s 0 -F 1000000 -T SRS_HG002_stat.json -r  hs37d5_1_0_1000000.fa -o hybrid_demo.vcf -I hs37d5.stat.txt -n SRS_HG002_1_869000_870000_60X_demo.bam -l LRS_HG002_1_869000_870000_10X_demo.bam -p HIFI
```

## Changelog

### v2.0.2 2026-04-09

Fixed bug: A bug was identified, fixed, and the relevant update has been pushed to GitHub.
When performing variant calling on next-generation sequencing (NGS) data, a local repeat rate file R for the reference genome is required as input. Previously, if this file was not provided, the program would still run, but errors might occur under certain circumstances. Now, when processing NGS data, the program will check for this file and terminate execution if the file is not found, thus avoiding potential issues.

### v2.0.4 2026-06-24

Implementation updates for diverse datasets. To support datasets with broader read-quality profiles and reference-to-sample divergence, we added several implementation options to HitSV. First, an optional contig-polishing step was added after local assembly to improve candidate contigs reconstructed from lower-accuracy long reads, such as ONT reads generated by older chemistries or ultralong-read protocols. Second, candidate-contig realignment was extended to support minimap2 assembly-alignment presets corresponding to different sequence-divergence levels. The ASM5 setting was used as the default for human datasets, whereas ASM10 or ASM20 was used when evaluating non-human species with larger divergence between the analyzed assembly and the reference genome. These options were used only to improve alignment robustness and allele reconstruction across datasets; the same repeat-aware signature-recognition framework was retained across human and non-human benchmarks.

## License

This project is protected by the [GNU GPL v3](LICENSE) license.

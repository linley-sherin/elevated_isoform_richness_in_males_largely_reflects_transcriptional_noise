# Elevated isoform richness in males largely reflects transcriptional noise

Code for manuscript https://doi.org/10.64898/2026.08.10.744030.

## Directory structure

```
01_build_longread_transcriptome/   PacBio Iso-Seq workflow, orthogonal short-read data, QC
02_build_shortread_transcriptome/  Short-read transcriptome built with StringTie
03_data_analysis/                  analyses and figures
tools/                             third-party tools (not tracked)
data/                              inputs and outputs, one folder per section (not tracked)
```

Scripts are numbered in the order they are run, and each writes to a matching folder in `data/`
(e.g. `01_build_longread_transcriptome/09_illumina_star.sh` → `data/01_build_longread_transcriptome/09_star/`).
Helper R code called by the scripts is in each section's `lib/` folder. Figures are numbered after the script that makes them.

### 01_build_longread_transcriptome

| Steps | |
|---|---|
| 01–06 | Iso-Seq: skera, lima, refine, fix bam sample names, cluster, map and collapse |
| 07–10 | Illumina short-read support: cutadapt, Trimmomatic, STAR, kallisto |
| 11–12 | SQANTI3 QC and filter |
| 13–16 | Post-SQANTI QC: ORF prediction (TD2), redundancy removal, long-read kallisto, edgeR expression filter |

### 02_build_shortread_transcriptome

| Steps | |
|---|---|
| 01 | StringTie assembly from the STAR alignments |
| 02 | ORF prediction (TD2) |
| 03 | Transcript filtering |
| 04 | kallisto quantification |
| 05 | edgeR expression filter |

### 03_data_analysis

| Steps | |
|---|---|
| 01–03 | SQANTI results, ORF status, build mastertable of isoform metrics |
| 04 | Differential gene and isoform expression (edgeR) |
| 05–07 | Overview plots, saturation curves, data correlation |
| 08 | Null model of splicing stochasticity |
| 09–11 | Alternative splicing events (SUPPA2) |
| 12–13 | Isoform metrics by sex, isoform switching |

## Usage

Run scripts from the project folder, e.g.

```bash
bash 01_build_longread_transcriptome/07_illumina_cutadapt.sh
```

Each script's header lists what it reads and writes. Shell scripts skip outputs that already exist.

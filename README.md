# PKD1 Variant Calling Pipeline

Analysis pipeline for Labcorp’s LR-PCR assay to detect PKD1 variants.

## Workflow
This workflow uses HiFi reads directly from the Revio instrument (UBAM format). PCR duplicates are marked with `pbmarkdup` and then aligned to “human_g1k_v37_decoy.fasta” using `pbmm2`. Small variants are called with `DeepVariant` and variant calls are indexed with bcftools. Structural variants are called with `sawfish`; however, only small variants are used for validation. Both small and structural variants are jointly phased using `HiPhase`. Phased variant calls are filtered using `bcftools` to only include “PASS” variants (i.e., no low-quality or RefCall variants) and also filtered to only exon regions using “PKD1.exons.bed”. Finally, average sequencing depth over exon regions is calculated using `mosdepth` and the “PKD1.exons.bed” file. 

![](Flowchart.png)

# Install Dependencies
## Tools
| Tool                                                         | Version|
|--------------------------------------------------------------|--------|
| [pbmarkdup](https://github.com/PacificBiosciences/pbmarkdup) | 1.1.0  |
| [pbmm2](https://github.com/PacificBiosciences/pbmm2)         | 1.16.0 |
| [deepvariant](https://github.com/google/deepvariant)         | 1.8.0  |
| [bcftools](https://github.com/samtools/bcftools)             | 1.21   |
| [sawfish](https://github.com/PacificBiosciences/sawfish)     | 0.12.9 |
| [hiphase](https://github.com/PacificBiosciences/hiphase)     | 1.4.5  |
| [mosdepth](https://github.com/brentp/mosdepth)               | 0.3.3  |
| [hap.py](https://github.com/Illumina/hap.py)                 | v0.3.8 |

## Pull singularity containers
```
singularity pull docker://google/deepvariant:1.8.0
singularity pull bcftools_v1.21.sif docker://quay.io/biocontainers/bcftools@sha256:e32f39e2c692c4a6ec8ba4a48b1f9b5b6c01616aba46affdaca0b943014e794e
singularity pull docker://pkrusche/hap.py:v0.3.8
```

## Create conda environment with PacBio tools
```
conda env create -f pkd1.yaml
```

`pkd1.yaml`:
```
name: pkd1
channels:
  - bioconda
dependencies:
  - pbmarkdup=1.1.0
  - pbmm2=1.16.0
  - sawfish=0.12.9
  - hiphase=1.4.5
  - mosdepth=0.3.3
```

```
conda activate pkd1
```

# Code

## Mark Duplicates
```
pbmarkdup -j 16 --log-level INFO \
	--log-file 01_SAMPLE.mkdup.log \
	../00_hifi_reads/SAMPLE.hifi_reads.bam \
	01_SAMPLE.mkdup.bam
```

## Align HiFi Reads
```
pbmm2 align -j 16 --preset HIFI --sort \
	$REFS/human_g1k_v37_decoy.fasta \
	01_SAMPLE.mkdup.bam \
	02_SAMPLE.pbmm2.bam
```

## Call Small Variants
```
singularity run -B $REFS -B $PWD \
	deepvariant_1.8.0.sif \
	/opt/deepvariant/bin/run_deepvariant \
	--model_type PACBIO \
	--ref $REFS/human_g1k_v37_decoy.fasta \
	--reads $PWD/02_SAMPLE.pbmm2.bam \
	--output_vcf $PWD/03_SAMPLE.dv.vcf.gz \
	--output_gvcf $PWD/03_SAMPLE.dv.g.vcf.gz \
	--num_shards 16
```

## Index Small Variants
```
singularity run -B $PWD bcftools_v1.21.sif bcftools index 03_SAMPLE.dv.vcf.gz
```

## Call Structural Variants
```
sawfish discover --threads 16 \
	--ref $REFS/human_g1k_v37_decoy.fasta \
	--bam 02_SAMPLE.pbmm2.bam \
	--output-dir 04_SAMPLE.sawfish_discover

sawfish joint-call --threads 16 \
	--sample 04_SAMPLE.sawfish_discover \
	--output-dir 04_SAMPLE.sawfish_call
```

## Phase Variants
```
hiphase \
	--threads 16 \
	--reference $REFS/human_g1k_v37_decoy.fasta \
	--bam 02_SAMPLE.pbmm2.bam \
	--output-bam 05_SAMPLE.haplotagged.bam \
	--vcf 03_SAMPLE.dv.vcf.gz \
	--output-vcf 05_SAMPLE.dv.phased.vcf.gz \
	--vcf 04_SAMPLE.sawfish_call/genotyped.sv.vcf.gz \
	--output-vcf 05_SAMPLE.sawfish.phased.vcf.gz \
	--stats-file 05_hiphase.stats.csv \
	--blocks-file 05_hiphase.blocks.tsv \
	--summary-file 05_hiphase.summary.tsv
```

## Filter PASS Variants within Region-of-Interest\
```
singularity run -B $PWD -B $REFS bcftools_v1.21.sif \
	bcftools view -R $REFS/PKD1.exons.bed -f PASS \
	05_SAMPLE.dv.phased.vcf.gz \
	-o 06_SAMPLE.dv.phased.pass.vcf.gz -Oz

singularity run -B $PWD -B $REFS bcftools_v1.21.sif \
	bcftools view -R $REFS/PKD1.exons.bed -f PASS \
	05_SAMPLE.sawfish.phased.vcf.gz \
	-o 06_SAMPLE.sawfish.phased.pass.vcf.gz -Oz
```

## Count Total Variants
```
singularity run -B $PWD bcftools_v1.21.sif \
	bcftools view -v snps 06_SAMPLE.dv.phased.pass.vcf.gz | \
	grep -vc '^#' > 06_SAMPLE.dv.phased.pass.snvs.count

singularity run -B $PWD bcftools_v1.21.sif \
	bcftools view -v indels 06_SAMPLE.dv.phased.pass.vcf.gz | \
	grep -vc '^#' > 06_SAMPLE.dv.phased.pass.indels.count

singularity run -B $PWD bcftools_v1.21.sif \
	bcftools view 06_SAMPLE.sawfish.phased.pass.vcf.gz | \
	grep -vc '^#' > 06_SAMPLE.sawfish.phased.pass.count
```

## Calculate Coverage over Region-of-Interest
```
mosdepth -t 8 -b $REFS/PKD1.exons.bed \
	06_SAMPLE 05_SAMPLE.haplotagged.bam

grep 16_region 06_SAMPLE.mosdepth.summary.txt | \
	cut -f 4 > 06_SAMPLE.mosdepth.region.cov
```

## Compare Variants to Expected/Truth
```
singularity exec -B $REFS -B $PWD \
	hap.py_v0.3.8.sif /opt/hap.py/bin/hap.py \
	$PWD/00_expected_variants.norm.vcf \
	$PWD/06_SAMPLE.dv.phased.pass.vcf.gz \
	-f $REFS/PKD1.exons.bed \
	-r $REFS/human_g1k_v37_decoy.fasta \
	-o 07_SAMPLE.dv.phased.pass.happy --set-gt hom
```

## Consolidate Results
```
grep PASS 07_SAMPLE.dv.phased.pass.happy.summary.csv | sed 's/,/\t/g' > 07_SAMPLE.dv.phased.pass.happy.summary.tsv

paste 06_SAMPLE.dv.phased.pass.snvs.count \
	06_SAMPLE.dv.phased.pass.indels.count \
	06_SAMPLE.sawfish.phased.pass.count \
	06_SAMPLE.mosdepth.region.cov \
	07_SAMPLE.dv.phased.pass.happy.summary.tsv \
	> 08_SAMPLE.summary.tsv
```

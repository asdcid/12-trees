# 12-trees

## Introduction

## Data
We have eight branches: A,B,C,D,E,F,G,H, each branch has 3 bio replicates, each bio replicates has 2 technical replicates.

## Method
### 1 Map reads to genome reference using ngm

```
#R1, R2 are the illumina paired-end reads
ID=$sampleID
RGLB=$(zcat --force $R1 | head -n 1 | cut -d: -f3,4 --output-delimiter=.)
ngm \
    --color \
    -t $threads \
    -r $genome_fasta \
    -p \
    -1 $R1 \
    -2 $R2 \
    --rg-id $ID \
    --rg-sm $ID \
    --rg-lb $RGLB \
    --rg-pl illumina \
    -o $outputSam

samtools sort -@ $threads -O bam -o $outputBam $outputSam
```

--------------
### 2 Mark duplicates
### 2.1 Merge different technical replicates according to sample (ABCDEFGH): A, B

```
$GATK MergeSamFiles\ 
 -I $outputBam_techRep_1.sort.bam \
 -I $outputBam_tecRep_2.sort.bam \
 -O $merge.bam
```

### 2.2 Sort BAM files

```
$GATK SortSam \ 
    -I $merge.bam \ 
    -O $sample.sort.bam \ 
    -SO coordinate
``` 

### 2.3 Make duplicateds

```
$GATK MarkDuplicates \ 
    -I $sample.sort.bam \ 
    -O $sample.sort.dedup.bam \
    -M $metrics_sample.txt
```

### 2.4 Build index

```
$GATK BuildBamIndex \ 
    -I $sample.sort.dedup.bam
```
--------------
### NOTE
The next step should be "Base Quality Score Recalibration" (BQSR), but it requires ***knownSites***. According to GATK suggests, I run the first call variants without BSQR, then use the variants result to call the recalibrated data, and using the recalibrated data to call variants again. 

--------------

### 3. HaplotypeCaller
### 3.1 Create dict for genome reference

```
$GATK CreateSequenceDictionary \
    -R $genome_fasta \
    -O ${genome_fasta/.fasta/.dict}
```
--------------

### 3.2 Call variants per-sample with HaplotypeCaller

```
#heterozygosity is different in different species
$GATK HaplotypeCaller \
  -R $genome_fasta \
  --emit-ref-confidence GVCF \
  -I $sample.sort.dedup.bam \
  -stand-call-conf 50 \
  --heterozygosity 0.015 \
  -O $sample.g.vcf.gz
```

### 3.3 Consolidate GVCFs (joint genotype)
This step may not necessary for 8 branches somatic variation call if we merge the 3 bio replicates together previously. 

### 3.3.1 Combine all sample's `gvcf` files
```
for sample in $samples
do
    sample_gvcfs=${sample_gvcfs}"-V ${sample_gvcfs_outputDir}/${sample}.g.vcf.gz"
done

$GATK CombineGVCFs \
      -R $genome_fasta \
      ${sample_gvcfs} \
      -O $combined.g.vcf.gz
```

### 3.3.2 Using GenotypeGVCFs to call joint genotype

```
gatk GenotypeGVCFs \
    -R $genome_fasta \
    --heterozygosity 0.015 \
    -V $combined.g.vcf.gz \
    -O $genotypeGVCFs.vcf.gz
```

--------------

### 4 Filter Variants (Hard filters)


>**Parametersfor hard-filters**
>
>**1. QualByDepth (QD)** 
>This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-reference samples.
>
>**2. FisherStrand (FS)**
>Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.
>
>**3. StrandOddsRatio (SOR)**
>The StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. Higher values indicate more strand bias.
>
>**4. RMSMappingQuality (MQ)**
>This is the Root Mean Square of the mapping quality of the reads across all samples.
>The score of *best mapping* for BWA is 60.  
>
>**5. MappingQualityRankSumTest (MQRankSum)**
>This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele). Note that the mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.
>
>**6. ReadPosRankSumTest (ReadPosRankSum)**
>This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.
>
>**7. ExcessHet** Optional, related annotation: inbreedingCoeff, https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_annotator_ExcessHet.php
>Describes the heterozygosity of the called samples, giving a probability of excess heterozygosity being observed
>Phred-scaled p-value for exact test of excess heterozygosity.
This annotation estimates the probability of the called samples exhibiting excess heterozygosity with respect to the null hypothesis that the samples are unrelated. The higher the score, the higher the chance that the variant is a technical artifact or that there is consanguinuity among the samples. In contrast to Inbreeding Coefficient, there is no minimal number of samples for this annotation. If samples are known to be related, a pedigree file can be provided so that the calculation is only performed on founders and offspring are excluded.


***The filter should be more strict for BQSR, and can be normal for the last round (after get the final BQSR result)***

### 4.1a SNP hard filters

```
# Extract SNP from call set
$GATK SelectVariants \
    -select-type SNP \
    -V $genotypeGVCFs.vcf.gz \
    -O $genotypeGVCFs.snp.vcf.gz

# Plot the figure, determine parameters for filtering SNPs

# Apply the filter to the SNP call set
$GATK VariantFiltration \
    -V $genotypeGVCFs.snp.vcf.gz \
    --filter-expression "DP > 480.0 || QD < 2.0 || MQ < 40.0 || FS > 50.0 || ExcessHet > 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O $genotypeGVCFs.snp.filter.vcf.gz

# Plot the figure, double check the filter performance

```

### 4.1b Indel hard filters

```
# Extract Indel from call set
$GATK SelectVariants \
    -select-type INDEL \
    -V $genotypeGVCFs.vcf.gz \
    -O $genotypeGVCFs.indel.vcf.gz

# Plot the figure, determine parameters for filtering SNPs

# Apply the filter to the SNP call set
$GATK VariantFiltration \
    -V $genotypeGVCFs.indel.vcf.gz \
    --filter-expression "DP > 480.0 || QD < 2.0 || FS > 50.0 || ExcessHet > 40.0 || SOR > 3.0 || ReadPosRankSum < -20.0" \
    --filter-name "Filter" \
    -O $genotypeGVCFs.indel.filter.vcf.gz

# Plot the figure, double check the filter performance

```

### 4.2 Exclude filter sites
Since VariantFilteration only mark the read as "Filter" (the name is defined by yourself) and PASS (the name is defined by GATK) in fields[6], get the "Pass reads"

```
$GATK SelectVariants \
    --exclude-filtered \
    -V $genotypeGVCFs.XXX.filter.vcf.gz \
    -O $genotypeGVCFs.xxx.removeFilter.vcf.gz
```

### 4.3 Merge SNP and Indel

```
$GATK MergeVcfs \
    -I $genotypeGVCFs.snp.removeFilter.vcf.gz \
    -I $genotypeGVCFs.indel.removeFilter.vcf.gz \
    -O $genotypeGVCFs.filter.vcf.gz
```
--------------
**Biallelic sites**
Only include the biallelic site, specific locus in a genome that contains two observed alleles (-m2 -M2: remove multiallelic site, a specific locus in a genome that contains three or more observed alleles).

In GATK:
>True multiallelic sites are not observed very frequently unless you look at very large cohorts, so they are often taken as a sign of a noisy region where artifacts are likely.

--------------
### 4.4 Remove multiallelic sites

```
#get the biallelic sites
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    -o $genotypeGVCFs.filter.biallelic.vcf \
    $genotypeGVCFs.filter.vcf.gz


bgzip $genotypeGVCFs.filter.biallelic.vcf

$GATK IndexFeatureFile \
    -F $genotypeGVCFs.filter.biallelic.vcf.gz

```

--------------

### 5 Evaluating the quality of  variant callset

GATK introduction
 https://gatkforums.broadinstitute.org/gatk/discussion/6308/evaluating-the-quality-of-a-variant-callset

>**1. Number of Indels & SNPs**
>It counts only biallelic sites and filters out multiallelic sites
>Compare the number of Indels and SNPs between different sample.
>
>**2. Indel Ratio**
>The indel ratio is determined to be the total number of insertions divided by the total number of deletions; this tool does not include filtered variants in this calculation. Usually, the indel ratio is around 1, as insertions occur typically as frequently as deletions. However, in rare variant studies, indel ratio should be around 0.2-0.5. For example, an indel ratio of ~0.95 indicates that these variants are not likely to have a bias affecting their insertion/deletion ratio.
>
>**3. TiTv Ratio**: This metric is the ratio of **t**rans**i**tion (Ti, A<->G or C<->T) to **t**rans**v**ersion (Tv, A<->C，A<->T，G<->C和G<->T) SNPs. 
>
>If the distribution of transition and transversion mutations were random (i.e. without any biological influence) we would expect a ratio of **0.5**. This is simply due to the fact that there are twice as many possible transversion mutations than there are transitions. 
>
>However, in the biological context, it is very common to see a methylated cytosine undergo deamination to become thymine. As this is a transition mutation, it has been shown to increase the expected random ratio from **0.5** to **~2.01**. 
>
>The TiTv Ration for **Novel variants** is around **1.5**.
>
>Furthermore, **CpG islands**, usually found in primer regions, have higher concentrations of methylcytosines. By including these regions, **whole exome sequencing** shows an even stronger lean towards transition mutations, with an expected ratio of **3.0**-**3.3**.



### 5.1 VariantEval 

```
$GATK VariantEval \
   -O $SampleVariants_Evaluation.eval.grp \
   --eval $genotypeGVCFs.filter.biallelic.vcf.gz
```
--------------

### 6 Base Quality Score Recalibration (BQSR)
Since we do not have ```known_sites```, we use the vcf generated above as the ```known_sites``` to adjust the data, and run the variant call again.  

>"Note,	if	you	wish	to	do	BQSR	on	non-human	samples,	you	can	use	the	above	
filtered	file	(but	generated	from	the	whole	genome)	as	the	“known”	variant	
input.			This	set	does	not	need	to	be	precise	since	the	amount	of	error	in	the	
reads	usually	far	exceeds	the	number	of	variants	that	are	called,	and	true	
positives	should	not	generally	exhibit	typical	BQSR	captured	patterns." https://qcb.ucla.edu/wp-content/uploads/sites/14/2017/08/gatkDay3.pdf

### 6.1 Analyze patterns of covariation in the dataset, model the error modes and compute adjustments

```
$GATK BaseRecalibrator \
    -R $genome_fasta \
    -I $sample.sort.dedup.bam \
    -O $recalibration_pre.table \
    --known-sites $genotypeGVCFs.filter.biallelic.vcf.gz 
```

### 6.2. Apply recalibration adjustments to BAM

```
$GATK ApplyBQSR \
    -R $genome_fasta \
    -I $sample.sort.dedup.bam \
    -O $sample.sort.dedup.rec.bam \
    -bqsr $recalibration_pre.table \
    --add-output-sam-program-record 
```

### 6.3 Visualizing the quality of a recalibration run
### 6.3.1 Get recalibration table again 
Using _sample.sort.dedup.rec.bam_ as input, using _sample.sort.dedup.rec.bam_  in the downstream analysis if all looks good with the AnalyzeCovariates plots

```
$GATK BaseRecalibrator \
    -R $genome_fasta \
    -I $sample.sort.dedup.rec.bam \
    -O $recalibration_post.table  \
    --known-sites $genotypeGVCFs.filter.biallelic.vcf.gz \
```

### 6.3.2 Generate before/after plots, and csv

```
$GATK AnalyzeCovariates \
    -before $recalibration_pre.table \
    -after $recalibration_post.table \
    -plots $AnalyzeCovariates.pdf
```
The analysis of figures are in https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
**NOTE**: the recalibration may change the quality score (such as from ~40-50 to ~20-30), so don't set the in -stand-call-conf in haplotypeCaller. 

--------------------

### 7 HaplotypeCaller using recalibrated data

**NOTE**: run step 7 - 9 for each sample

```
#heterozygosity is different in different species
$GATK HaplotypeCaller \
  -R $genome_fasta \
  -I $sample.sort.dedup.rec.bam \
  --dbsnp $genotypeGVCFs.filter.biallelic.vcf.gz \
  --heterozygosity 0.015 \
  -O $sample.vcf.gz
```

### 8 Filter Variants (Hard filters)
### 8.1a SNP hard filters
```
# Extract SNP from call set
$GATK SelectVariants \
    -select-type SNP \
    -V $sample.vcf.gz \
    -O $sample.snp.vcf.gz

# Plot the figure, determine parameters for filtering SNPs

# Apply the filter to the SNP call set
$GATK VariantFiltration \
    -V $sample.snp.vcf.gz \
    --filter-expression "DP > 480.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || ExcessHet > 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O $sample.snp.filter.vcf.gz

# Plot the figure, double check the filter performance

```

### 8.1b Indel hard filters

```
# Extract Indel from call set
$GATK SelectVariants \
    -select-type INDEL \
    -V $sample.vcf.gz \
    -O $sample.indel.vcf.gz

# Plot the figure, determine parameters for filtering SNPs

# Apply the filter to the SNP call set
$GATK VariantFiltration \
    -V $sample.indel.vcf.gz \
    --filter-expression "DP > 480.0 || QD < 2.0 || FS > 200.0 || ExcessHet > 40.0 || SOR > 3.0 || ReadPosRankSum < -20.0" \
    --filter-name "Filter" \
    -O $sample.indel.filter.vcf.gz

# Plot the figure, double check the filter performance

```

### 8.2 Exclude filter sites
Since VariantFilteration only mark the read as "Filter" (the name is defined by yourself) and PASS (the name is defined by GATK) in fields[6], get the "Pass reads"

```
$GATK SelectVariants \
    --exclude-filtered \
    -V $sample.XXX.filter.vcf.gz \
    -O $sample.XXX.removeFilter.vcf.gz
```

### 8.3 Merge SNP and Indel

```
$GATK MergeVcfs \
    -I $sample.snp.removeFilter.vcf.gz \
    -I $sample.indel.removeFilter.vcf.gz \
    -O $sample.filter.vcf.gz
```

### 8.4 Remove multiallelic sites

```
#get the biallelic sites
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    -o $sample.filter.biallelic.vcf \
    $sample.filter.vcf.gz


bgzip $sample.filter.biallelic.vcf

$GATK IndexFeatureFile \
    -F $sample.filter.biallelic.vcf.gz
```
### 8.5 Remove sites in the repeat regions

```
custom python script
```

### 9 Evaluating the quality of  variant callset
```
$GATK VariantEval \
   -O $SampleVariants_Evaluation.eval.grp \
   --eval $sample.filter.biallelic.vcf.gz
```

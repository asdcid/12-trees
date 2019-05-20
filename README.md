# 12-trees
##DATA

Reference:

```
/home/rob/work/Eucalyptus_pauciflora/genome/data/Eucalyptus_pauciflora/masurca_35kb_pilon5.fa
```

Reads (after trim)

```
 /home/raymond/work/Eucalyptus_pauciflora/genome/bin/raw_read_check/8_branches/short_read/2_bbduk/428_result
 /home/raymond/work/Eucalyptus_pauciflora/genome/bin/raw_read_check/8_branches/short_read/2_bbduk/431_result
```

Working Dir

```
/home/rob/work/Eucalyptus_pauciflora/genome/bin/somatic_variation
```

Final result

```
/home/rob/work/Eucalyptus_pauciflora/genome/bin/somatic_variation/result/5_filter/*/*indel/snp.filter.biall.vcf.gz
```

--------------
##METHOD
###1 Map reads to reference (ngm)
*NOTE*: should add ```--rg-id --rg-sm --rg-lb --rg-pb``` for GATK (unique group ID) 

3 bio replicates: such as RL41, RL49, RL57
2 tech replicates: 428, 431
name: 428_RL41: A1, 428_RL49: A2, 428_RL57: A3, 431_RL41: A4, 431_RL49: A5, 431_RL57: A6

8 branches: A,B,C,D,E,F,G,H 

--------------
###2 Mark duplicates
###2.1 Merge different runs according to sample (ABCDEFGH): A, B
***Q: only merge tech replicates (428, 431), or even merge bio replicates (RL41, RL49, RL57)***

*NOTE*: ```picard``` can keep @PG, but ```samtools``` will lose it when merge

```
# sample.aligned.bam is the merge output
$GATK MergeSamFiles\ 
 -I input.aligned.1.bam \
 -I input.aligned.2.bam\
 ...
 -O sample.aligned.bam
```

###2.2 Sort BAM files

```
$GATK SortSam \ 
    -I sample.aligned.sam \ 
    -O sample.sort.bam \ 
    -SO coordinate

#or

samtools sort \
    -@ $threads \
    -O bam \
    -o sample.sort.bam \
    sample.aligned.bam
``` 

###2.3 Make duplicateds

```
$GATK MarkDuplicates \ 
    -I sample.sort.bam \ 
    -O sample.sort.dedup.bam \
    -M metrics_sample.txt
```

###2.4 Build index

```
samtools index sample.sort.dedup.bam

#or

$GATK BuildBamIndex \ 
    -I sample.sort.dedup.bam
```

###2.5 Plot the insert size figure

```
$GATK CollectInsertSizeMetrics \
    -I sample.sort.dedup.bam  \
    -O insert_size_metrics \
    -H insert_size_histogram.pdf \
#M=0.5
```
--------------
###NOTE1
The next step should be "Base Quality Score Recalibration" (BQSR), but it requires ***knownSites***. According to GATK suggests, I run the first round without BSQR, then use the SNP result to call the recalibrated data, and run SNP again. Repeat this step a few times until convergence.
>**Q: I'm working on a genome that doesn't really have a good SNP database yet. I'm wondering if it still makes sense to run base quality score recalibration without known SNPs.**

>A: The base quality score recalibrator treats every reference mismatch as indicative of machine error. True polymorphisms are legitimate mismatches to the reference and shouldn't be counted against the quality of a base. We use a database of known polymorphisms to skip over most polymorphic sites. Unfortunately without this information the data becomes almost completely unusable since the quality of the bases will be inferred to be much much lower than it actually is as a result of the reference-mismatching SNP sites.

>However, all is not lost if you are willing to experiment a bit. You can bootstrap a database of known SNPs. Here's how it works:

>First do an initial round of SNP calling on your original, unrecalibrated data.
Then take the SNPs that you have the highest confidence in and use that set as the database of known SNPs by feeding it as a VCF file to the base quality score recalibrator.
>Finally, do a real round of SNP calling with the recalibrated data. These steps could be repeated several times until convergence.

--------------

###3. Somatic variatance call
###3.1 Create dict for genome reference

```
$GATK CreateSequenceDictionary \
    -R ${ref_fasta} \
    -O ${ref/.fasta/.dict}
```
--------------

###3.2 Call variants per-sample with HaplotypeCaller

```
$GATK HaplotypeCaller \
  -R ${ref_fasta} \
  --emit-ref-confidence GVCF \
  -I sample.sort.dedup.bam \
  -stand-call-conf 50 \
  --heterozygosity $heterozygosity \
  -O sample.g.vcf.gz
```

###3.3 Consolidate GVCFs (joint genotype)
This step may not necessary for 8 branches somatic variation call if we merge the 3 bio replicates together previously. 

###3.3.1 Combine all sample's `gvcf` files
```
$GATK CombineGVCFs \
      -R ${ref_fasta} \
      -V sample1.g.vcf.gz \
      -V sample2.g.vcf.gz \ 
      -O combined.g.vcf.gz

#or 
for sample in $samples
do
    sample_gvcfs=${sample_gvcfs}"-V ${sample_gvcfs_outputDir}/${sample}.g.vcf.gz"
done

$GATK CombineGVCFs \
      -R ${ref_fasta} \
      ${sample_gvcfs} \
      -O combined.g.vcf.gz
```

###3.3.2 Using GenotypeGVCFs to call joint genotype

```
gatk GenotypeGVCFs \
    -R ${ref_fasta} \
    --heterozygosity 0.015 \
    --tmp-dir $TMP \
    -A ExcessHet \
    -A FisherStrand \
    -A QualByDepth \
    -A StrandOddsRatio \
    -A RMSMappingQuality \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest \
    -V combined.g.vcf.gz \
    -O genotypeGVCFs.vcf.gz
```

--------------

###4 Filter Variants by Variant (Quality Score) Recalibration

**NOTE2**

Two filter methods:
>**VQSR**:  recalibrate variant quality scores and produce a callset filtered for the desired levels of sensitivity and specificity

>**Hard-filters**: consists of choosing specific thresholds for one or more annotations and throwing out any variants that have annotation values above or below the set thresholds 

Introduction about Hard-filters


https://software.broadinstitute.org/gatk/documentation/article.php?id=6925


https://www.jianshu.com/p/ff8204ae7ebf

GATK example


https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set

However, VQSR require ```1. known variation sites``` and ```2. enough variants```. We don't have the ```known variation sites```, so we use the hard-filters. 

We need to extract  ```SNP``` and ```indel```  from call set  (`vcf`) to do the hard-filters.

>**Parametersfor hard-filters**

>**1. QualByDepth (QD)** usually <2.0
>This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-reference samples.

>**2. FisherStrand (FS)**
>Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.

>**3. StrandOddsRatio (SOR)**
>The StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. Higher values indicate more strand bias.

>**4. RMSMappingQuality (MQ)**, not for indel, usually <40
>This is the Root Mean Square of the mapping quality of the reads across all samples.
>The score of *best mapping* for BWA is 60.  

>**5. MappingQualityRankSumTest (MQRankSum)**, not for indel
>This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele). Note that the mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.

>**6. ReadPosRankSumTest (ReadPosRankSum)**
>This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.

>**7. InbreedingCoeff** Optional (must have >10 samples), https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.3.0/org_broadinstitute_hellbender_tools_walkers_annotator_InbreedingCoeff.php
>Describes the heterozygosity of the called samples, though without explicitly taking into account the number of samples
>Likelihood-based test for the consanguinuity among samples
This annotation estimates whether there is evidence of consanguinuity in a population. The higher the score, the higher the chance that some samples are related. If samples are known to be related, a pedigree file can be provided so that the calculation is only performed on founders and offspring are excluded.
> only work with more than 10 samples, diploid genome. 

>**8. ExcessHet** Optional, related annotation: inbreedingCoeff, https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_annotator_ExcessHet.php
>Describes the heterozygosity of the called samples, giving a probability of excess heterozygosity being observed
>Phred-scaled p-value for exact test of excess heterozygosity.
This annotation estimates the probability of the called samples exhibiting excess heterozygosity with respect to the null hypothesis that the samples are unrelated. The higher the score, the higher the chance that the variant is a technical artifact or that there is consanguinuity among the samples. In contrast to Inbreeding Coefficient, there is no minimal number of samples for this annotation. If samples are known to be related, a pedigree file can be provided so that the calculation is only performed on founders and offspring are excluded.

Hard-filters standards for unknown company, BGI, GATK and somatic variation pipeline

```
               DP   QD    FS      SOR    MQ      MQRankSum  ReadPosRankSum InbreedingCoeff
SNP   COMPANY <10   <2.0  >60.0   >4.0   <40.0   <-12.5      <-8.0
      BGI           <2.0  >60.0   >3.0   <40.0   <-12.5      <-8.0
      GATK          <2.0  >60.0   >3.0   <40.0   <-12.5      <-8.0

INDEL COMPANY <10   <2.0  >200.0  >10,0                      <-20.0          <-0.8
      BGI           <2.0  >200.0  >10.0          <-12.5      <-8.0
      GATK          <2.0  >200.0  >10.0                      <-20.0

```

Hard-filters for somatic variation pipeline 

https://github.com/adamjorr/somatic-variation

```
bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' repfiltered.vcf 
```

***All filter need to be plotted from the vcf file first before filtering to ensure the suitable parameters.***
***The filter should be more strict for BQSR, and can be normal for the last round (after get the final BQSR result)***

###4.1a SNP hard-filter

```
# Extract SNP from call set
$GATK SelectVariants \
    -select-type SNP \
    -V genotypeGVCFs.vcf.gz \
    -O genotypeGVCFs.snp.vcf.gz

# Plot the figure, determine parameters for filtering SNPs

# Apply the filter to the SNP call set
$GATK VariantFiltration \
    -V genotypeGVCFs.snp.vcf.gz \
    --filter-expression "DP > 480.0 || QD < 2.0 || MQ < 40.0 || FS > 50.0 || ExcessHet > 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O genotypeGVCFs.snp.filter.vcf.gz

# Plot the figure, double check the filter performance

```

###4.1b Indel hard-filter

```
# Extract Indel from call set
$GATK SelectVariants \
    -select-type INDEL \
    -V genotypeGVCFs.vcf.gz \
    -O genotypeGVCFs.indel.vcf.gz

# Plot the figure, determine parameters for filtering SNPs

# Apply the filter to the SNP call set
$GATK VariantFiltration \
    -V genotypeGVCFs.indel.vcf.gz \
    --filter-expression "DP > 480.0 || QD < 2.0 || FS > 50.0 || ExcessHet > 40.0 || SOR > 3.0 || ReadPosRankSum < -20.0" \
    --filter-name "Filter" \
    -O genotypeGVCFs.indel.filter.vcf.gz

# Plot the figure, double check the filter performance

```

###4.2 Exclude filter sites
Since VariantFilteration only mark the read as "Filter" (the name is defined by yourself) and PASS (the name is defined by GATK) in fields[6], get the "Pass reads"

```
$GATK SelectVariants \
    --exclude-filtered \
    -V input.XXX.filter.vcf.gz \
    -O removeFilter.vcf.gz
```


###4.3 Merge SNP and Indel

```
$GATK MergeVcfs \
    -I genotypeGVCFs.snp.removeFilter.vcf.gz \
    -I genotypeGVCFs.indel.removeFilter.vcf.gz \
    -O genotypeGVCFs.filter.vcf.gz

# remove useless doc
rm -f \
    genotypeGVCFs.snp.vcf.gz* \
    genotypeGVCFs.snp.filter.vcf.gz* \
    genotypeGVCFs.indel.vcf.gz* \
    genotypeGVCFs.indel.filter.vcf.gz*

```


--------------

###5 Evaluating the quality of  variant callset

GATK introduction
Summary: https://gatkforums.broadinstitute.org/gatk/discussion/6308/evaluating-the-quality-of-a-variant-callset
VariantEval: https://gatkforums.broadinstitute.org/gatk/discussion/6211/howto-evaluate-a-callset-with-varianteval
CollectVariantCallingMetrics: https://software.broadinstitute.org/gatk/documentation/article?id=6186
https://www.jianshu.com/p/d760f55dfc6c

Tools: ```CollectVariantCallingMetrics```,  ```VariantEval```
Metrics: ```Number of Indels & SNPs```,  ```Indel Ratio```, ```TiTv Ration```

>**1. Number of Indels & SNPs**
>It counts only biallelic sites and filters out multiallelic sites
>Compare the number of Indels and SNPs between different sample.

>**2. Indel Ratio**
>The indel ratio is determined to be the total number of insertions divided by the total number of deletions; this tool does not include filtered variants in this calculation. Usually, the indel ratio is around 1, as insertions occur typically as frequently as deletions. However, in rare variant studies, indel ratio should be around 0.2-0.5. For example, an indel ratio of ~0.95 indicates that these variants are not likely to have a bias affecting their insertion/deletion ratio.

>**3. TiTv Ratio**: This metric is the ratio of **t**rans**i**tion (Ti, A<->G or C<->T) to **t**rans**v**ersion (Tv, A<->C，A<->T，G<->C和G<->T) SNPs. 

>If the distribution of transition and transversion mutations were random (i.e. without any biological influence) we would expect a ratio of **0.5**. This is simply due to the fact that there are twice as many possible transversion mutations than there are transitions. 

>However, in the biological context, it is very common to see a methylated cytosine undergo deamination to become thymine. As this is a transition mutation, it has been shown to increase the expected random ratio from **0.5** to **~2.01**. 

>The TiTv Ration for **Novel variants** is around **1.5**.

>Furthermore, **CpG islands**, usually found in primer regions, have higher concentrations of methylcytosines. By including these regions, **whole exome sequencing** shows an even stronger lean towards transition mutations, with an expected ratio of **3.0**-**3.3**.



###5.1 VariantEval 

```
# Basic output:
$GATK VariantEval \
   -O SampleVariants_Evaluation.eval.grp \
   --eval genotypeGVCFs.filter.vcf.gz


```

###5.1a CollectVariantCallingMetrics (May not work, we don't have the *True* SNP set)

```
CollectVariantCallingMetrics \
    -I genotypeGVCFs.filter.vcf.gz \
    -O output_metrics \
    --DBSNP dbsnp_138.b37.excluding_sites_after_129.vcf 
```

###5.2 Remove multiallelic sites

```
#get the biallelic sites
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    -o genotypeGVCFs.filter.biallelic.vcf \
    genotypeGVCFs.filter.vcf.gz

#can't be gzip
bgzip genotypeGVCFs.filter.biallelic.vcf

$GATK IndexFeatureFile \
    -F genotypeGVCFs.filter.biallelic.vcf.gz

```
--------------

**Biallelic sites**
Only include the biallelic site, specific locus in a genome that contains two observed alleles (-m2 -M2: remove multiallelic site, a specific locus in a genome that contains three or more observed alleles).

In GATK:
>True multiallelic sites are not observed very frequently unless you look at very large cohorts, so they are often taken as a sign of a noisy region where artifacts are likely.

In somatic variation pipeline (https://github.com/adamjorr/somatic-variation)
>We recommend following this with a depth and heterozygosity filter and removing any non-snps. This can be accomplished by running

--------------


###5.3 Check quality again after removed multiallelic sites

```
# Basic output:
$GATK VariantEval \
   -O SampleVariants_Evaluation.eval.grp \
   --eval genotypeGVCFs.filter.biallelic.vcf.gz
```
--------------


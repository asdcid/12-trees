# encoding: utf-8
"""
8_calculate.py
get sites:1. present in all three bio replicates in the same branch
          2. not present in all eight branch

Created by www on 12:57 am, May 18, 2019
"""
import numpy as np
import sys
import os
import gzip
import copy

def loadFile(inputFile, o, snps, sample, mask):
    sites   = 0
    with gzip.open(inputFile, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            line    = line.strip()
            if not line:
                continue
            #the end of file, can't find another method
            if line.startswith('TBI'):
                break
            #print(line)
            info = line.split()
            ref     = info[0]
            alt     = info[4]
            pos     = int(info[1])
            #remove site in repeat or homoploymer regions
            if mask[ref][pos - 1] != 0:
                continue
            if not pos in snps[ref]:
                snps[ref][pos] = {}
            if not alt in snps[ref][pos]:
                snps[ref][pos][alt] = {}
            snps[ref][pos][alt][sample] = ''
            sites += 1
    o.write('%s\t%d\n' % (sample, sites))
    return snps

def output(samples, o, d):
    snps         = {}
    diffSites    = 0

    o.write('Number of SNP sites in all three replicates in the same branch\n')

    for ID in samples:
        commonThreeRep = 0
        for ref in samples[ID]:
            if not ref in snps:
                snps[ref] = {}
            for pos in samples[ID][ref]:
                #if have more than two genotype, remove. (keep site: A/A, A/T, if the site is A/A, A/T, A/C, remove)
                if len(samples[ID][ref][pos]) != 1:
                    continue
                #get SNP exist in all 3 bio replicates
                for alt in samples[ID][ref][pos]:
                    if len(samples[ID][ref][pos][alt]) == 3:
                        commonThreeRep += 1
                        if not pos in snps[ref]:
                            snps[ref][pos] = {}
                        if not alt in snps[ref][pos]:
                            snps[ref][pos][alt] = {}
                        snps[ref][pos][alt][ID] = ''
        o.write('%s\t%d\n' % (ID, commonThreeRep))


    for ref in snps:
        for pos in snps[ref]:
            #if have more than two genotype, remove. (keep site: A/A, A/T, if the site is A/A, A/T, A/C, remove)
            if len(snps[ref][pos]) != 1:
                continue
            for alt in snps[ref][pos]:
                if len(snps[ref][pos][alt]) == 8:
                    #print('All samples have the same alt')
                    #print(ref, pos, alt, snps[ref][pos][alt])
                    continue
                diffSites += 1
                d.write('%s\t%s\t%s\t' % (ref, alt, pos))
                for ID in snps[ref][pos][alt]:
                    d.write('%s,' % ID)
                d.write('\n')
    o.write('Number of SNP sites do not co-exist in all 8 branches\n%d\n' % diffSites)


def loadRef(refFile, removeFile):
    refs    = {}
    mask    = {}
    with open(refFile) as f:
        for line in f:
            info    = line.split()
            ref     = info[0]
            length  = int(info[1])
            refs[ref] = {}
            mask[ref] = np.zeros(length)
    with open(removeFile) as f:
        for line in f:
            info    = line.split()
            ref     = info[0]
            start   = int(info[1])
            end     = int(info[2])
            mask[ref][start - 1 : end] = 1

    return mask, refs


def main():
    inputDir    = 'result/5_filter/'
    outputFile  = 'result/5_filter/summary'
    outputDetail= 'result/5_filter/detail_site'
    removeFile  = 'ref/masked/Epau.repeat.hmopolymer.sort.bed'
    refFile     = 'ref/Epau.fa.fai'

    o       = open(outputFile, 'w+')
    d       = open(outputDetail, 'w+')

    samples     = {}
    o.write('Number of SNP sites in each sample in the same branch\n')

    f   = os.walk(inputDir)

    mask, refs    = loadRef(refFile, removeFile)

    for root, dirs, names in f:
        for name in names:
            if not 'snp.filter.biall.vcf.gz' in name:
                continue
            if name.split('.')[-1] == 'tbi':
                continue
            inputFile   = os.path.join(root, name)

            #sample: A_1, A_2,...H_3
            sample  = name.split('.snp')[0]
            #ID: A B C D E ... H
            ID      = sample.split('_')[0]
            if not ID in samples:
                samples[ID] = copy.deepcopy(refs)

            samples[ID] = loadFile(inputFile, o, samples[ID], sample, mask)

    output(samples, o, d)

    o.close()
    d.close()

if __name__ == '__main__':
    main()

#!/usr/bin/env python3

# See GitHub repositories https://github.com/banfieldlab/mattolm-public-scripts and 
# https://github.com/christophertbrown/iRep for up-to-date versions and dependencies

# Matt Olm
# mattolm@berkeley.edu

# Version 0.2
# - Fixed bug where the reference sequence has a non-ACTG character, it'll now just skip that base instead of crashing

import sys
import os
import argparse
import iRep
import numpy as np
import genome_variation as gv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time
import pandas as pd

from fasta import iterate_fasta as parse_fasta

'''
Known bugs:
* When calculating the coverage, uses the coverage in the mpileup, which
  includes indels. When looking for SNPs, ingores indels and just uses
  regular bases (so requires the min coverage of bases without indels).
  However, when calculating the number of masked bases, it uses the coverage
  with indels. This should be a very small effect, however

*** Data structure used is dictionary of dictionaries ***

genomes[genome]['masked_bases'] = set()
    contains globally masked bases, such as those in tRNAs / end of scaffolds
genomes[genome]['contig_order'] = []
    this determines the location of each base relative to the genome
    see location_to_key
genomes[genome]['length'] = int
genomes[genome]['samples'][sample]['cov'] = []
genomes[genome]['samples'][sample]['SNPs'] = []
genomes[genome]['samples'][sample]['totals']['unmasked_length'] = int
genomes[genome]['samples'][sample]['totals']['average_coverage'] = float
genomes[genome]['samples'][sample]['totals']['median_coverage'] = float
genomes[genome]['samples'][sample]['totals']['consensus_ANI'] = float

*********************************************************

'''

# Return dictionary of base : count
def pile_baseCounts(line):
    data = line.strip().split('\t')
    bp = data[1]
    bases = data[4].upper()
    ref = data[2].upper()
    types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}
    i = 0
    while i < len(bases):
            base = bases[i]
            if base == '^' or base == '$':
                    i += 1
            elif base == '-':
                    i += 1
            elif base == '*':
                    types['-'] += 1
            elif base == '+':
                    #determine length of insertion
                    j = 1
                    while True:
                        if bases[i+j+1].isdigit():
                            j+=1
                        else:
                            break
                    if j > 1:
                        addNum = int(bases[i+1:i+1+j])
                    else:
                        addNum = int(bases[i+1])
                    i+=len(str(addNum))
                    addSeq = ''
                    for a in range(addNum):
                            i += 1
                            addSeq += bases[i]
                    types['+'].append(addSeq)
            elif base == '.' or base == ',':
                    types[ref] += 1
            else:
                    if base in types:
                            types[base] += 1
                    else:
                            types['X'].append(base)
            i += 1
    adds = '.'
    if len(types['+']) > 0:
            adds = ','.join(types['+'])
    amb = '.'
    if len(types['X']) > 0:
            amb = ','.join(types['X'])
    out = {'A':types['A'],'G':types['G'],'C':types['C'],'T':types['T']}
    return out

def location_to_key(genomes,genome,location):
    tally = 0
    for i,contig in enumerate(genomes[genome]['contig_order']):
        length = genomes[genome]['contig_length'][i]
        if location < tally + length:
            return "{0}:{1}".format(contig,location-tally)
        tally += length
    return("NON-VALID_LOCATION")

# Return the genome and starting location of the contig passed to it
def contig_lookup(genomes,saught_contig):
    for genome in genomes:
        if saught_contig in genomes[genome]['contig_order']:
            start = 0
            for i,contig in enumerate(genomes[genome]['contig_order']):
                if contig == saught_contig:
                    return(genome,start)
                start += genomes[genome]['contig_length'][i]
    raise ValueError('{0} not in fasta input'.format(saught_contig))

def iterate_pileup_special(genomes,pileup):
    # Return all lines of the pileup, using "genomes" to determine the scaffold and genome_location of each line
    reader = open(pileup)
    current_contig = "START"
    current_genome = "START"
    base_location = "START"
    for line in reader.readlines():
        linewords = line.split()
        contig = linewords[0]
        if contig != current_contig:
            try:
                current_genome, base_location = contig_lookup(genomes,contig)
            except:
                continue
            current_contig = contig
        loc = linewords[1]
        yield (current_genome,int(loc) + int(base_location),line)

def setup_genomes(fastas,pileups):
    genomes = {}
    for genome in fastas:
        genomes[genome] =  g = {'samples':{}}
        g['length'] = 0
        g['contig_order'] = []
        g['contig_length'] = []
        g['masked_bases'] = set()
        for seq in parse_fasta(genome):
            ID = seq[0].split('>', 1)[1].split()[0]
            length = len(seq[1])
            g['contig_order'].append(ID)
            g['contig_length'].append(length)
            g['length'] += length
        for sample in pileups:
            g['samples'][sample] = {'SNPs':[]}
            g['samples'][sample]['cov'] = [0 for i in range(0, g['length'])]
            g['samples'][sample]['var'] = [0 for i in range(0, g['length'])]
            g['samples'][sample]['unmasked_length'] = 0
    return genomes

def parse_piluep_line(line, SNP_chriteria):
    linewords = line.split()
    cov = int(linewords[3])
    ref = linewords[2].upper()
    if ref == 'N':
        return(cov,False,True)
    if ref not in ['A','C','T','G']:
        return(cov,False,True)
    min_cov, min_perc = SNP_chriteria
    base_counts = pile_baseCounts(line)
    total_bases = sum(base_counts.values())
    isSNP = False
    if total_bases >= min_cov:
        isSNP = (((total_bases - base_counts[ref]) / total_bases) > min_perc)
    return(cov,isSNP,False)

def read_pileups(genomes,pileups,SNP_chriteria):
    for pileup in pileups:
        for touple in iterate_pileup_special(genomes,pileup):
             genome, location, line = touple
             cov, isSNP, isN = parse_piluep_line(line, SNP_chriteria)
             genomes[genome]['samples'][pileup]['cov'][location-1] = cov
             if isSNP: genomes[genome]['samples'][pileup]['SNPs'].append(int(location))
             if isN: genomes[genome]['masked_bases'].add(int(location))
    return genomes

def get_edges(genomes,genome,edge_mask):
    bases = []
    for i,contig in enumerate(genomes[genome]['contig_order']):
        genome, start = contig_lookup(genomes,contig)
        length = genomes[genome]['contig_length'][i]
        bases += range(start,start + edge_mask)
        bases += range(length-edge_mask - 1,length - 1)
    return bases

def mask_SNPs(genomes,genome,masked_bases):
    for sample in genomes[genome]['samples']:
        genomes[genome]['samples'][sample]['SNPs'] = list(set(genomes[genome]['samples'][sample]['SNPs']))
        snps = genomes[genome]['samples'][sample]['SNPs'].copy()
        for snp in snps:
            if snp in masked_bases:
                genomes[genome]['samples'][sample]['SNPs'].remove(snp)
    return genomes[genome]

def get_locatations_from_file(genomes,genome,file_mask):
    reader = open(file_mask)
    bases = []
    for line in reader.readlines():
        line = line.strip()
        scaffold, start, end = line.split()
        try:
            contig_genome, length = contig_lookup(genomes,scaffold)
        except:
            continue
        if contig_genome != genome:
            continue
        if int(end) + length -1 > genomes[genome]['length']:
            raise ValueError("{0} higher than the length of {1}".format(int(end) + length,genome))
        bases += range(int(start) + length, int(end) + length)
    reader.close()
    return bases

def coverage_mask(genomes,genome,sample,SNP_chriteria):
    min_cov, min_perc = SNP_chriteria
    # number of positions masked by coverage
    C = len([i for i in genomes[genome]['samples'][sample]['cov'] if i < min_cov])
    # number of positions masked on the genome (rRNA, edge of scaffold, ect.)
    M = len(genomes[genome]['masked_bases'])
    # number of positions counted twice
    B = len([i for i in genomes[genome]['masked_bases'] if genomes[genome]['samples'][sample]['cov'][i-1] < min_cov])
    return (C + M - B)


def refine_snps(genomes,file_mask,edge_mask,SNP_chriteria):
    for genome in genomes:
        if edge_mask != None:
            genomes[genome]['masked_bases'] = set(get_edges(genomes,genome,edge_mask)) | genomes[genome]['masked_bases']
        if file_mask != None:
            genomes[genome]['masked_bases'] = set(get_locatations_from_file(genomes,genome,file_mask)) | genomes[genome]['masked_bases']
        genomes[genome] = mask_SNPs(genomes,genome,genomes[genome]['masked_bases'])
        for sample in genomes[genome]['samples']:
            total_masked_bases = coverage_mask(genomes,genome,sample,SNP_chriteria)
            genomes[genome]['samples'][sample]['unmasked_length'] = genomes[genome]['length'] - total_masked_bases
    return genomes

def calculate_totals(genomes):
    for genome in genomes:
        for sample in genomes[genome]['samples']:
            current = genomes[genome]['samples'][sample]['totals'] = {}
            current['average_coverage'] = np.mean(genomes[genome]['samples'][sample]['cov'])
            current['median_coverage'] = np.median(genomes[genome]['samples'][sample]['cov'])
            current['std_coverage'] = np.std(genomes[genome]['samples'][sample]['cov'])
            current['breadth'] = 1 - (list(genomes[genome]['samples'][sample]['cov']).count(0) / genomes[genome]['length'])
            current['unmasked_breadth'] = genomes[genome]['samples'][sample]['unmasked_length'] / genomes[genome]['length']
            try :
                current['consensus_ANI'] = 1 - (len(genomes[genome]['samples'][sample]['SNPs']) / genomes[genome]['samples'][sample]['unmasked_length'])
            except :
                current['consensus_ANI'] = "NaN"
    return(genomes)

def totals_to_csv(genomes,outfile):

    genomeD = {}
    sampleD = {}
    avg_covD = {}
    med_covD = {}
    std_covD = {}
    conc_ANID = {}
    breadth_D = {}
    Unbreadth_D = {}
    snpD = {}

    i = 0
    for genome in genomes:
         for sample in genomes[genome]['samples']:
             genomeD[i] = os.path.basename(genome)
             sampleD[i] = os.path.basename(sample)
             avg_covD[i] = genomes[genome]['samples'][sample]['totals']['average_coverage']
             med_covD[i] = genomes[genome]['samples'][sample]['totals']['median_coverage']
             std_covD[i] = genomes[genome]['samples'][sample]['totals']['std_coverage']
             conc_ANID[i] = genomes[genome]['samples'][sample]['totals']['consensus_ANI']
             breadth_D[i] = genomes[genome]['samples'][sample]['totals']['breadth']
             snpD[i] = len(genomes[genome]['samples'][sample]['SNPs'])
             Unbreadth_D[i] = genomes[genome]['samples'][sample]['totals']['unmasked_breadth']
             i += 1

    data = pd.DataFrame({ 'genome' : pd.Series(genomeD),
                      'sample': pd.Series(sampleD),
                      'average_coverage': pd.Series(avg_covD),
                      'median_coverage': pd.Series(med_covD),
                      'std_coverage': pd.Series(std_covD),
                      'breadth': pd.Series(breadth_D),
                      'unmasked_breadth': pd.Series(Unbreadth_D),
                      'SNPs': pd.Series(snpD),
                      'consensus_ANI': pd.Series(conc_ANID)})

    data.to_csv(outfile + "_totals.csv",index = False)

def totals_to_snps(genomes,outfile):
    o = open(outfile + "_snps.txt",'w')
    for genome in genomes:
        for sample in genomes[genome]['samples']:
            for snp in genomes[genome]['samples'][sample]['SNPs']:
                o.write("{0}\t{1}\t{2}\n".format(os.path.basename(genome),os.path.basename(sample),location_to_key(genomes,genome,snp)))
    o.close()

def status_report(genomes):
    for genome in genomes:
        print("*** genome ***")
        print("{0}".format(genome))
        print("* masked bases *")
        print(genomes[genome]['masked_bases'])
        for sample in genomes[genome]['samples']:
            print("*** sample ***")
            print(sample)
            print("* SNPs *")
            print(genomes[genome]['samples'][sample]['SNPs'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# filter sam file based on mismatches')
    parser.add_argument('-p', nargs = '*', action = 'store', required = True, help = 'pileup(s)')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', required = False, \
            help = 'fasta(s)')
    parser.add_argument(\
            '-o', required = True, type = str, \
            help = 'prefix for output files (table and plots)')
    parser.add_argument(\
            '-c', default = '5', \
            help = 'minimum coverage required to support base call [5]')
    parser.add_argument(\
            '-mS', default = '0.8', \
            help = 'minimum percent to call SNP [0.8]')
    parser.add_argument(\
            '-file_mask', \
            help = 'file of positions to mask in the format: contig start end')
    parser.add_argument(\
            '-edge_mask', default = '50', type = int, \
            help = 'mask this distance from the edge of each contig [50]')
    args = vars(parser.parse_args())

    pileups, fastas, min_cov = args['p'], args['f'], int(args['c'])
    SNP_chriteria = [int(args['c']),float(args['mS'])]

    start = time.clock()
    print("Reading pileups")
    genomes = setup_genomes(fastas,pileups)
    genomes = read_pileups(genomes,pileups, SNP_chriteria)
    end = time.clock()
    print("Run time of {0} seconds".format(end-start))

    genomes = refine_snps(genomes,args['file_mask'],args['edge_mask'],SNP_chriteria)

    start = time.clock()
    print("Calculating totals")
    genomes = calculate_totals(genomes)
    end = time.clock()
    print("Run time of {0} seconds".format(end-start))

    totals_to_csv(genomes,args['o'])
    totals_to_snps(genomes,args['o'])

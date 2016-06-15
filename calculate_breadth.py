#!/usr/bin/env python3

# Adopted heavily from Chris Brown's Calculate_coverage.py and iRep (2.15.16)

import sys
import os
import argparse
import numpy as np
import pandas as pd

from mapped import get_reads as mapped
from fasta import iterate_fasta as parse_fasta

def length_and_bases(coverage, sam):
    for line in sam:
        line = line.strip()
        if line.startswith('@SQ'):
            line = line.strip().split()
            scaffold, length = line[1].split(':', 1)[1], float(line[-1].split(':', 1)[1])
            if scaffold not in coverage:
                coverage[scaffold] = [length, {}]
        elif line.startswith('@') is False:
            line = line.split('\t')
            scaffold, bases = line[2], float(len(line[9]))
            if scaffold not in coverage:
                coverage[scaffold] = [length, {}]
            map = int(line[1])
            if map != 4 and map != 8:
                if scaffold != '*':
                    if sam not in coverage[scaffold][1]:
                        coverage[scaffold][1][sam] = 0
                    coverage[scaffold][1][sam] += bases
    return coverage

def combine_by_sample(coverage):
    combined_coverage = {}
    sams = []
    for scaffold in coverage:
        length, combined = coverage[scaffold][0], {}
        for sam in coverage[scaffold][1]:
            comb = '.'.join(sam.name.split('.')[0:2])
            if comb not in sams:
                sams.append(comb)
            if comb not in combined:
                combined[comb] = 0
            combined[comb] += coverage[scaffold][1][sam]
        combined_coverage[scaffold] = [length, combined]
    return combined_coverage, sams

def calculate_coverage(bases):
    coverage = {}
    for scaffold in bases:
        length, counts = bases[scaffold][0], bases[scaffold][1]
        scaffold_coverage = {}
        for count in counts:
            scaffold_coverage[count] = float(counts[count] / length)
        coverage[scaffold] = [length, scaffold_coverage]
    return coverage

def print_coverage(coverage, sams):
    out = ['# scaffold: length']
    for sam in sorted(sams):
        out.append(sam.name)
    yield out
    for scaffold in coverage:
        length, cov = coverage[scaffold][0], coverage[scaffold][1]
        out = ['%s: %s' % (scaffold, length)]
        for sam in sorted(sams):
            if sam in cov:
                out.append(cov[sam])
            else:
                out.append(0)
        yield out

def open_files(files):
    """
    open files in list, use stdin if first 
    item in list is '-'
    """
    if files is None:
        return files
    if files[0] == '-':
        return [sys.stdin]
    return [open(i) for i in files]

def iterate_sams(sams, combine = False):
    coverage = {}
    for sam in sams:
        coverage = length_and_bases(coverage, sam)
    if combine is True:
        coverage, sams = combine_by_sample(coverage)
    coverage = calculate_coverage(coverage)
    return coverage, sams
    
def filter_mapping(sam, mismatches, sort_sam = False, sbuffer = 100):
    """
    create generator for filtering mapping
    """
    for type, read in mapped(sam, False, mismatches, \
            'both', sort_sam, False, False, sbuffer):
        if type == 10 or type == 20:
            yield read
            
def parse_genomes_fa(fastas, mappings):
    """
    genomes[genome name] = {order: [contig order], samples: {}}
        samples[sample name] = {cov: [coverage by position], contigs: {}}
            contigs[contig name] = [coverage by position]
    """
    id2g = {} # contig ID to genome lookup 
    genomes = {} # dictionary for saving genome info
    for genome in fastas:
        name = genome.name
        samples = {s[0]:{'contigs':{}, 'cov':[]} for s in mappings}
        g = genomes[name] = {'order':[], 'samples':samples}
        g['len'] = 0
        for seq in parse_fasta(genome): 
            ID = seq[0].split('>', 1)[1].split()[0]
            g['order'].append(ID)
            id2g[ID] = name
            length = len(seq[1])
            g['len'] += length
            for sample in list(samples.keys()):
                g['samples'][sample]['contigs'][ID] = \
                    [0 for i in range(0, length)]
    return genomes, id2g

def calc_coverage(genomes, mappings, id2g):
    """
    for each sample:
        calcualte coverage at each position in genome
    # genomes[genome]['samples'][sample]['contigs'][ID] = cov
    """
    print("# parsing mapping files", file=sys.stderr)
    sys.stderr.flush()
    for sample, reads in mappings:
        for read in reads:
            c = read[2] # contig
            # left-most position of read on contig
            start, length = int(read[3]), len(read[9])
            end = start + length - 1
            for i in range(start - 1, end):
                try: 
                    genomes[id2g[c]]['samples'][sample]\
                            ['contigs'][c][i] += 1
                except:
                    continue
    # combine coverage data for contigs
    for genome in list(genomes.values()):
        order, samples = genome['order'], genome['samples'] 
        for sample in list(samples.values()):
            for contig in order:
                try:
                    sample['cov'].extend(sample['contigs'][contig])
                    #del sample['contigs'][contig]
                except:
                    continue
            sample['avg_cov'] = np.average(sample['cov'])
    return genomes
            
def iRep(fastas, id2g, mappings, out, threads):
    """
    est. growth from slope of coverage
     1) calculate coverage over windows

    """
    # get genome info from fastas
    genomes, id2g = parse_genomes_fa(fastas, mappings)
    # get coverage from sam files
    genomes = calc_coverage(genomes, mappings, id2g)
    
    return genomes
    
def toPandas(genomes, outt):
    genome = {}
    samples = {}
    average_coverage = {}
    median_coverage = {}
    std_deviation = {}
    length = {}
    zero_count = {}
    br = {}
    
    i = 0
    keys = {}
    gens = {}
    samps = {}
    for g in genomes:
        for sample in genomes[g]['samples']:
            genome[i] = g
            samples[i] = sample
            median_coverage[i] = np.median(genomes[g]['samples'][sample]['cov'])
            average_coverage[i] = genomes[g]['samples'][sample]['avg_cov']
            std_deviation[i] = np.std(genomes[g]['samples'][sample]['cov'])
            length[i] = genomes[g]['len']
            zero_count[i] = list(genomes[g]['samples'][sample]['cov']).count(0)
            br[i] = 1 - (zero_count[i] /  length[i])
            i += 1
    
    data = pd.DataFrame({ 'genome' : pd.Series(genome),
                            'sample' : pd.Series(samples),
                            'average_cov' : pd.Series(average_coverage),
                            'std_dev' : pd.Series(std_deviation),
                            'length' : pd.Series(length),
                            'median_cov' : pd.Series(median_coverage),
                            'zero_count' : pd.Series(zero_count),
                            'breadth' : pd.Series(br)})
                                
    data.to_csv(outt + ".csv")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# calculate coverage and breadth from sam files')
    parser.add_argument(\
            '-s', nargs = '*', action = 'store', required = True, \
            help = 'sam(s)')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', required = False, \
            help = 'fasta(s)')
    parser.add_argument(\
            '-mm', required = False, default = False, type = int, \
            help = 'max. # of read mismatches allowed (default: no limit)')
    parser.add_argument(\
            '-o', required = True, type = str, \
            help = 'prefix for output files (table and plots)')
    parser.add_argument(\
            '--verbose', action = 'store_true', \
            help = 'optional - print output to screen as well as file')
    parser.add_argument(\
            '-t', required = False, default = 6, type = int, \
            help = 'threads (default: 6)')
    args = vars(parser.parse_args())
    
    fastas = open_files(args['f'])
    sams, mm = args['s'], args['mm']
    mappings = [[s, filter_mapping(s, mm)] for s in sams]
    s2bin = None
    
    genomes = iRep(fastas, s2bin, mappings, args['o'], args['t'])
    
    if args['verbose']:
        for genome in genomes:
            print(genome)
            for sample in genomes[genome]['samples']:
                print(sample)
                print("Median:")
                print((np.median(genomes[genome]['samples'][sample]['cov'])))
                print("Coverage:")
                print((np.median(genomes[genome]['samples'][sample]['avg_cov'])))
                print("Std")
                print((np.std(genomes[genome]['samples'][sample]['cov'])))
                print("")
            
    toPandas(genomes, args['o'])

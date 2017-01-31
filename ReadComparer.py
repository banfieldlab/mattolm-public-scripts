#!/usr/bin/python3

#Version 0.2


#   This program is meant to do various tasks with the output of the program VarScan (calling SNPs with reads)
#   As of now it supports compare_vars. This method takes in a series of VarScan files, then looks at all SNPs that are common between them, as well as SNPs unique to each
#       It then outputs summary data, as well as a list of common SNPs and unique SNPs
#   Functions to add:
#       DONE Putting the SNPs into .gff format for genious
#       Statistics
#       DONE Of the SNPs that aren't common, is it due to coverage? Are any of the other type seen?
#       Have a .gb be imported as well, and retreave tags for SNP genes
#       Use the .gb to give % of nucleotides with enough coverage

import seaborn as sns
from sys import argv
import sys
import argparse
import ntpath
import glob
import Bio
import copy
import pandas as pd
import numpy as np
from Bio import SeqIO

def setup(inputs, bases, ignore):
    for Input in inputs:
        print('Analyzing ' + ntpath.basename(Input) + '\n')
        names.append(ntpath.basename(Input))
        results = getSNPs(Input, bases, ignore)
        covSets.append(results[0])
        snpCounts.append(results[1])
        dics.append(results[2])
        if printSNPs: print_snps(results, Input)
        
    return (names, covSets, snpCounts, dics)
        
def print_snps(touple, Input):
    o = open(outputfile + ntpath.basename(Input) + "_snips.txt", 'w')
    for w in sorted(touple[2], key=touple[2].get, reverse=True):
        o.write( str(w) +  "\t" + str(touple[2][w][0]) + "\t" + str(touple[2][w][1]) + '\n')
    o.close()
    

def matrixComparison(inputs, outputfile, reflength):
    CompResults = [[0 for x in range(len(inputs))] for x in range(len(inputs))]
    for index1, var1 in enumerate(inputs):
        for index2, var2 in enumerate(inputs):
            if (index1 < index2):
                CompResults[index1][index2] = compare_vars(index1,index2,names,dics,covSets,reflength)
                print(CompResults[index1][index2])
            elif (index1 == index2):
                CompResults[index1][index2] = selfCompare(index1)
                print(CompResults[index1][index1])
    writeResultsAll(CompResults, outputfile + '_all.txt')
    if (dendrogram):
        writePandas(CompResults)
        writeDendrogram(CompResults, 1, '_dendrogram.txt')
    writeResults(CompResults, 1, '_ANI.txt')
    writeResults(CompResults, 2, '_Cov.txt')
    writeResults(CompResults, 0, '_ComCover.txt')

def writePandas(rMatrix):
    dist_array = np.empty([len(rMatrix[0]) , len(rMatrix[0])])
    #print(names)
    for r in range(len(rMatrix[0])):
        for c in range(len(rMatrix[0])):
            if (c < r):
                dist = rMatrix[c][r][1]
                dist_array[c][r] = dist * 100
            elif (c == r):
                dist = 1
                dist_array[c][r] = dist * 100
            else:
                dist = rMatrix[r][c][1]
                dist_array[c][r] = dist * 100
    pand_array = pd.DataFrame(dist_array,index = names,columns = names)
    
    # Add reference genome information
    ref_dic = {}
    for r in range(len(rMatrix[0])):
        dist = rMatrix[r][r][1]
        ref_dic[names[r]] = dist * 100
    ref_dic['Reference Genome'] = 100
    pand_array['Reference Genome'] = pd.Series(ref_dic)
    pand_array.loc['Reference Genome'] = pd.Series(ref_dic)
    #cols = pand_array.columns.values
    #cols[0] = 'Project1'
    #pand_array.columns = cols
    
    #print(dist_array)
    o = open(outputfile + '_pandasOut.csv','w')
    o.write(pand_array.to_csv())
    o.close()
    
    graphDendrogram(outputfile + '_pandasOut.csv')

def graphDendrogram(csv):
    dend_fn = csv
    dend_data = pd.read_csv(dend_fn, na_values = 'n/a')
    dend_data = dend_data.rename(columns = {'Unnamed: 0':'Project1'})
    dend_data = pd.melt(dend_data, id_vars=['Project1'])
    dend_data = dend_data.rename(columns = {'variable':'project2','value':'ANI'})
    dend_data = dend_data.pivot("Project1", "project2", "ANI")
    g = sns.clustermap(dend_data)
    g.savefig(outputfile + "_dendrogram.pdf")     

def writeDendrogram(rMatrix, numr, ty): 
 o = open(outputfile + ty, 'w')
 index = 0
 for r in range(len(rMatrix[0])):
  o.write('\n' + names[index] + '\t')
  index += 1
  for c in range(len(rMatrix[0])):
   if (c < r):
    toup = rMatrix[c][r]
    asn = toup[int(numr)]
    o.write(str(1-asn))
   elif (c == r):
    o.write('0')
   else:
    toup = rMatrix[r][c]
    asn = toup[int(numr)]
    o.write(str(1-asn))
   o.write('\t')
  toup = rMatrix[r][r]
  asn = toup[int(numr)]
  o.write(str(1-int(asn)))
 o.write("\nReference Genome\t")
 for r in range(len(rMatrix[0])):
  toup = rMatrix[r][r]
  asn = toup[int(numr)]
  o.write(str(1-asn) + '\t')
 o.write('0')

'''
def commonCoverage(out):
    o = open(out + "_covList.txt", 'w')
    if zero:
        for i in range(0,int(reflength)):
            if str(i) not in coverageDic:
                coverageDic[i] = 0
 for w in sorted(coverageDic, key=coverageDic.get, reverse=True):
  o.write( str(w) +  ": " + str(coverageDic[w]) + '\n')
'''
   

def writeResultsAll(rMatrix, out):
 o = open(out, 'w')
 o.write('Format: % common coverage between the two, ANI amoung covered areas, common coverage / refseq \n')
 for name in names:
  o.write(name + '\t')
 for r in range(len(rMatrix[0])):
  o.write('\n')
  for c in range(len(rMatrix[0])):
   var = rMatrix[c][r]
   o.write(str(var))
   o.write('\t')

def writeResults(rMatrix, numr, ty): 
 o = open(outputfile + ty, 'w')
 o.write('Format: % common coverage between the two, ANI amoung covered areas, common coverage / refseq \n')
 for name in names:
  o.write(name + '\t')
 for r in range(len(rMatrix[0])):
  o.write('\n')
  for c in range(len(rMatrix[0])):
   if (c <= r):
    toup = rMatrix[c][r]
    asn = toup[int(numr)]
    o.write(str(asn))
   else:
    o.write('0')
   o.write('\t')

def selfCompare(index):
    print("Comparing " + names[index] + " to itself")
    if (len(covSets[index]) > 0):
        ANI = (len(covSets[index])-snpCounts[index]) / float(len(covSets[index]))
        coverage = len(covSets[index]) / float(reflength)
        return('NA',ANI , coverage)
    else:
        print("ERRRRR")
        return ('NA',0,0)

# Accept a VarScan file, return the number of bases with coverage (that aren't ignored),
#   the number of SNPs, and a dictionary mapping base number to base call.
def getSNPs(infile, bases, ignore):
    
    #linewords[0] is the scaffold
    #linewords[1] is the base number
    #linewords[2] is the reference genome call
    #linewords[3] is the consensus call
 
    snpCount = 0
    dic = {}
    Cov = set([])
    
    reader = open(infile)
    data = reader.readlines()
    for index, line in enumerate(data):
        if (index == 0):
            continue
        linewords = line.split()
        poly = False
        snp = False
        
        if linewords[3] not in bases:
            poly = True
        if linewords[2] != linewords[3]:
            snp = True
        
        if not (ignore and poly):
            Cov.add(linewords[0] + ':' + linewords[1])
            if snp:
                snpCount += 1
                dic[linewords[0] + ':' + linewords[1]] = (linewords[2],linewords[3])
    
    return(Cov, snpCount, dic)
        
#flag = False
#if (iPos != None): 
#if int(linewords[1]) in iPos: flag = True
   
   #if(ignore):
    #if (linewords[3] in bases):
     #Cov.add(linewords[1])
   #else:
    #Cov.add(linewords[1])
   
   #if linewords[1] in coverageDic:
    #coverageDic[linewords[1]] += 1
   #else:
    #coverageDic[linewords[1]] = 1
   
   #if not flag:
    #if (linewords[2] != linewords[3]):
     #if (ignore):
      #if (linewords[3] in bases):
       #snpCount += 1
       #dic[linewords[1]] = (linewords[2], linewords[3])
    
     #else:
      #snpCount += 1
      #dic[linewords[1]] = (linewords[2], linewords[3])

 #return(Cov, snpCount, dic)   


def compare_vars(var1, var2, names, dics, covSets, reflength):
# Takes the indexes of two inputs and returns (% Common Coverage, ANI, Common Coverage / RefSeq)

    print("I'm going to compare ", names[var1], " and ", names[var2], " now") 
 
    Intersection = (covSets[var1] & covSets[var2])
    
    if (listContention):
    # Calculate ANI between regions of common coverage AND write file
        o = open(names[var1] + '_' + names[var2] + '_contestedBases.txt', 'w')
        diff = 0
        for loc in Intersection:
            if loc in dics[var1]:
                if loc in dics[var2]:
                    if (dics[var1][loc][1] != dics[var2][loc][1]):
                        diff += 1
                        o.write(str(loc) + '\n')
	                    # polymorphic SNP
                else:
                    diff += 1
                    o.write(str(loc) + '\t' + '\n')
            elif loc in dics[var2]:
                diff += 1
                o.write(str(loc) + '\t' + '\n')
        if (len(Intersection) == 0):
            ANI = 0
        else:
            ANI = (len(Intersection) - diff) / float(len(Intersection))
        print(ANI)
 
    else:
        # Calculate ANI between regions of common coverage
        diff = 0
        for loc in Intersection:
            if loc in dics[var1]:
                if loc in dics[var2]:
                    if (dics[var1][loc][1] != dics[var2][loc][1]):
                        diff += 1
                        #polymorphic SNP
                else:
                    diff += 1
            elif loc in dics[var2]:
                diff += 1
        if (len(Intersection) == 0):
            ANI = 0
        else:
            ANI = (len(Intersection) - diff) / float(len(Intersection))

 # Calculate % of common coverage
    ComCovPer = (len(Intersection) * 2) / float((len(covSets[var1]) + len(covSets[var2])))
 
    return (ComCovPer, ANI, (len(Intersection)/float(reflength)))

def getRange(fil, len):
 o = open(fil, 'r')
 data = o.readlines()
 ignores = {}
 for line in data:
  #print(line)
  linewords = line.split()
  r1 = linewords[0]
  r2 = linewords[2]
  for x in range(int(r1),int(r2) + 1):
   ignores[x] = True
 return ignores
  
def process_genome(file):
    genome_dic = {}
    l = 0
    genome = SeqIO.parse(file, "fasta")
    for record in genome:
        l += len(record.seq)
    return l
    
def trimInputs(inp, refL, min_b):
    trimmed = []
    for input in inp:
        keep = (file_len(input) / float(refL)) > min_b
        #print(input + "\t" + str(file_len(input) / float(refL)) + "\t" + str(keep))
        if keep: trimmed.append(input)
    return trimmed
    sys.exit()
    
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i    


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required Arguments
RedArgs = parser.add_argument_group('REQUIRED ARGUMENTS')
RedArgs.add_argument('-h', action="help", help="show this help message and exit")
RedArgs.add_argument("-i", "--input", nargs="+",
    help='input VarScan files seperated by commas, or a folder with --folder)')
RedArgs.add_argument('-o', "--output",
    help='output file name')
RedArgs.add_argument('-g', "--genome",
    help='reference genome')
    
#Optional Argumnets
OpArgs = parser.add_argument_group('OPTIONAL ARGUMENTS')
OpArgs.add_argument('-r', "--range_ignore",
    help='BROKEN    File containing list of ranges of bases to ignore, inclusive (typically 16S). Format: 1 - 3\n90 - 95')
OpArgs.add_argument('--folder', default = False, action = 'store_true',
    help='Input is a folder')
OpArgs.add_argument('--min_breadth', default = 0.8,
    help='The mininum percentage of the genome that must be represented in the VarScan files')
OpArgs.add_argument('--matrix', default = False, action = 'store_true',
    help='Generate a comparison matrix')
OpArgs.add_argument('--list', default = False, action = 'store_true',
    help='List bases of contention between each pair')
OpArgs.add_argument('--ignore', default = True, action = 'store_false',
    help='Exlcude from analysis bases that are considered polymorphic by varscan')
OpArgs.add_argument('--smart_ignore', default = False, action = 'store_true',
    help='Include polymorphisms called by varscan in comparison')
OpArgs.add_argument('--dend', default = False, action = 'store_true',
    help='Generate a data file for constructing a dendrogram')
OpArgs.add_argument('--SNPs', default = False, action = 'store_true',
    help='Print all called SNPs')
OpArgs.add_argument('--CC', default = False, action = 'store_true',
    help='Print the number of projects each base has coverage in')

args = parser.parse_args()

# Establish GLOBAL variables (I DON'T LIKE THIS)
names = []
covSets = []
snpCounts = []
dics = []
coverageDic = {}
iPos = None
if args.range_ignore != None:
    iPos = getRange(args.range_ignore, args.length)
outputfile = args.output


# Establish global flags
listContention = args.list
ignore = args.ignore
dendrogram = args.dend
printSNPs = args.SNPs
listContention = args.list
bases = ['A','G','C','T']
if args.smart_ignore:
    bases = ['A','G','C','T','W','S','M','K','R','Y','B','D','H','V']
    ignore = True

# Should you count 'zero' values when doing commonCoverage?
zero = True


if args.folder:
    inputs = glob.glob(args.input[0] + '*')
else:
    inputs = args.input

# Run the program
reflength = process_genome(args.genome)
inputs = trimInputs(inputs, reflength, float(args.min_breadth))
setup(inputs, bases, args.ignore)
if args.matrix:
    matrixComparison(inputs, args.output, reflength)
if args.CC:
    print("Nope, this isn't supported at the moment in this shitty program!")
    #commonCoverage(args.output)

#!/usr/bin/python3

import argparse
import glob
import pandas as pd
import seaborn as sns
import os

class pBase(object):
    location = ""
    ori = 'A'
    var = 'B'
    projects = {}

    def __init__(self, location, ori, var):
        self.location = location
        self.ori = ori
        self.var = var
        self.projects = {}

    def __str__(self):
        return str(self.location) + " " + self.ori + '/' + self.var + str(self.projects)

    def addProject(self, name, oriNumber, varNumber):
        self.projects[name] = (oriNumber, varNumber)

    def getProject(self, name):
        if name in self.projects:
            return self.projects[name]
        else:
            return(0,0)

    def percentRef(self, name):
        if name in self.projects:
            if self.projects[name][0] + self.projects[name][1] == 0:
                return 'NA'
            else:
                return float(self.projects[name][0]) / (int(self.projects[name][1]) + int(self.projects[name][0]))
        else:
            return 'NA'

def runProject(pileup):
    reader = open(pileup)
    data = reader.readlines()
    for line in data:
        linewords = line.split()
        loc = str(linewords[0]) + ":" + str(linewords[1])
        if loc in locList:
            #print("I found " + loc)
            addData(pileup, linewords, loc)

def addData(name, linewords, location):
    for poly in polys:
        if poly.location == location:
            current = poly

    bases = linewords[4].upper()
        #print("bases: " + str(linewords[4]))
        #print("ref: " + str(linewords[2]))
        #print("var: " + str(current.var))

    ref = linewords[2].upper()
    refC = 0
    var = current.var.upper()
    varC = 0

    if not ref == current.ori:
        print("LIST CALLED THE WRONG REFERENCE BASE AT " + str(current.location))

    #print("var is " + var)
    #print("ref is " + ref)

    for base in bases:
        if base == var:
            varC += 1
        if (base == '.' or base == ','):
            refC += 1
    #print(base)

    current.addProject(name, varC, refC)


def parseBases(infile, pileups):
    reader = open(infile)
    data = reader.readlines()
    tempList = []
    for line in data:
        linewords = line.split()
        location = linewords[0].strip()
        tempList.append(location)
    dic = {}
    ref = {}
    for t in tempList:
        # A C G T
        dic[t] = [0,0,0,0]
    for pileup in pileups:
        reader = open(pileup)
        data = reader.readlines()
        for line in data:
            linewords = line.split()
            code = linewords[0].strip() + ":" + linewords[1].strip()
            if code in tempList:
                bases = parsePile(line)
                dic[code] = [x + y for x, y in zip(bases,dic[code])]
                if code not in ref:
                    ref[code] = linewords[2].strip()
                else:
                    if ref[code] != linewords[2].strip():
                        print("ERROR AT POSITION " + str(code))
    rList = []
    for t in dic:
        var = dic[t].index(max(dic[t]))
        if var is 0: b = 'A'
        elif var is 1: b = 'G'
        elif var is 2: b = 'C'
        elif var is 3: b = 'T'
        if b == ref[t]:
            dic[t][var] = 0
            var = dic[t].index(max(dic[t]))
            if var is 0: b = 'A'
            elif var is 1: b = 'G'
            elif var is 2: b = 'C'
            elif var is 3: b = 'T'
        q = pBase(t,ref[t],b)
        rList.append(q)
    return rList

def parsePile(line):
    #print("line- " + str(line))
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
                    i += 1
                    addNum = int(bases[i])
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
    out = [types['A'],types['G'],types['C'],types['T']]
    return out

def parseList(infile):
    reader = open(infile)
    data = reader.readlines()
    tempList = []
    for line in data:
        linewords = line.split()
        base = pBase(linewords[0],linewords[1],linewords[2])
        tempList.append(base)
    return tempList

def printOutput():
 for poly in polys:
  print('\t' + str(poly.location) + '\t\t', end="")
 print('')
 for poly in polys:
  print('\t ' + str(poly.ori) + '\t' + str(poly.var) +  "\t ref%", end="")
 print('')
 for name in pileups:
  print(name + '\t', end = "")
  for poly in polys:
   toup = poly.getProject(name)
   print(str(toup[0]) + '\t' + str(toup[1]) + '\t' + str(poly.percentRef(name)) + '\t', end="")
  print('\n')

def printeOutput():
    print('Name', end="\t")
    for poly in polys:
        print(str(poly.location), end = "\t")
    print('')
    for name in pileups:
        print(name + '\t', end = "")
        for poly in polys:
            toup = poly.getProject(name)
            print(str(poly.percentRef(name)) + '\t' + str(int(toup[0]) + int(toup[1])) + '\t', end="")
        print('\n')
    for poly in polys:
        print('\t' + str(poly.ori) + ' > ' + str(poly.var), end="\t")

def pandasOut():

    firstIteration = True

    scaff = {}
    pos = {}
    varC = {}
    refC = {}
    per_ref = {}
    proj = {}
    ref = {}
    var = {}
    i = 0
    for name in pileups:
        for poly in polys:
            toup = poly.getProject(name)
            scaffold,position = poly.location.split(':')

            scaff[i] = scaffold
            pos[i] = position
            varC[i] = toup[0]
            refC[i] = toup[1]
            per_ref[i] = poly.percentRef(name)
            proj[i] = name

            if firstIteration:
                ref[i] = poly.ori
                var[i] = poly.var

            i+=1

        firstIteration = False

    data = pd.DataFrame({ 'scaffold' : pd.Series(scaff),
                      'location': pd.Series(pos),
                      'variant_base_count' : pd.Series(varC),
                      'reference_base_count' : pd.Series(refC),
                      'percent_reference' : pd.Series(per_ref),
                      'project': pd.Series(proj)})

    base_data = pd.DataFrame({ 'scaffold' : pd.Series(scaff),
                      'location': pd.Series(pos),
                      'reference_base' : pd.Series(ref),
                      'varient_base': pd.Series(var)})

    wd = os.getcwd()
    data.to_csv(wd + "/polymorpher_pandas.csv")
    base_data.to_csv(wd + "/bases_polymorpher_pandas.csv")


def printRout(name):
 o = open(str(name) + '.csv','w')
 o.write("Sample,Location,PercentReference,BaseChange\n")
 for proj in pileups:
  for poly in polys:
   o.write(proj + ',' + str(poly.location) + ',' + str(poly.percentRef(proj)) + ',' + str(poly.ori) + '>' + str(poly.var) + '\n')


parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mpileup", nargs="+",
    help='input mpileup files seperated by commas (or a folder with -f)')
parser.add_argument('-l', "--list",
    help='list of bases to anlyzed in the format (tab seperated): base number, orig nucl, variant nucl., or just a list of bases (in the format scaffold:position) with the appropriate flag')
parser.add_argument('-f', "--folder_input", default = False, action = "store_true",
    help="folder input")
parser.add_argument("-b", "--bases", default = False, action = "store_true",
    help = "input is just a list of bases, varient nucleotides will be determined [default: False].")
parser.add_argument("-R", "--rOut",
    help = "an output for use in R is made with this name.")
args = parser.parse_args()

pileups = args.mpileup
if args.folder_input: pileups = glob.glob(args.mpileup[0] + '*')

if not args.bases: polys = parseList(args.list)
else:
    polys = parseBases(args.list, pileups)
    print("Done processing base-list")

locList = []
for poy in polys:
    locList.append(str(poy.location))

for pileup in pileups:
    runProject(pileup)

#printeOutput()
pandasOut()
#if args.rOut is not None:
 #   printRout(args.rOut)

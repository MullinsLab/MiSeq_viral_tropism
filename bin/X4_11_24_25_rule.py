#!/usr/bin/python



import sys
import re
import os

def translate(codon):
    codon2aa = {
        "ATT" : "I",
        "ATC" : "I",
        "ATA" : "I",
        "CTT" : "L",
        "CTC" : "L",
        "CTA" : "L",
        "CTG" : "L",
        "TTA" : "L",
        "TTG" : "L",
        "GTT" : "V",
        "GTC" : "V",
        "GTA" : "V",
        "GTG" : "V",
        "TTT" : "F",
        "TTC" : "F",
        "ATG" : "M",
        "TGT" : "C",
        "TGC" : "C",
        "GCT" : "A",
        "GCC" : "A",
        "GCA" : "A",
        "GCG" : "A",
        "GGT" : "G",
        "GGC" : "G",
        "GGA" : "G",
        "GGG" : "G",
        "CCT" : "P",
        "CCC" : "P",
        "CCA" : "P",
        "CCG" : "P",
        "ACT" : "T",
        "ACC" : "T",
        "ACA" : "T",
        "ACG" : "T",
        "TCT" : "S",
        "TCC" : "S",
        "TCA" : "S",
        "TCG" : "S",
        "AGT" : "S",
        "AGC" : "S",
        "TAT" : "Y",
        "TAC" : "Y",
        "TGG" : "W",
        "CAA" : "Q",
        "CAG" : "Q",
        "AAT" : "N",
        "AAC" : "N",
        "CAT" : "H",
        "CAC" : "H",
        "GAA" : "E",
        "GAG" : "E",
        "GAT" : "D",
        "GAC" : "D",
        "AAA" : "K",
        "AAG" : "K",
        "CGT" : "R",
        "CGC" : "R",
        "CGA" : "R",
        "CGG" : "R",
        "AGA" : "R",
        "AGG" : "R",
        "TAA" : "*",
        "TAG" : "*",
        "TGA" : "*"
    }
    return codon2aa[codon]

usage = "usage python X4_11-24-25_rule.py infastafile inreffile outcsvfile outexcludefile"
if len(sys.argv) < 5:
    sys.exit(usage)

infile = sys.argv[1]
reffile = sys.argv[2]
outfile = sys.argv[3]
excludefile = sys.argv[4]
totalcount, count, stopcodoncount, validcount, x4count, countat11, countat24, countat25 = 0, 0, 0, 0, 0, 0, 0, 0
id, inreffile, inrefalignfile, refname, refseq, aaseq = '', '', '', '', '', ''
names, refaas = [], []
nameseq, namecount = ({} for i in range(2))
infilematch = re.search('^(.*?)_(.*).fasta', infile)
if infilematch:
    id = infilematch.group(1)
    inreffile = infilematch.group(1) + "_" + infilematch.group(2) + "_withRef.fasta"
    inrefalignfile = infilematch.group(1) + "_" + infilematch.group(2) + "_withRef_aligned.fasta"
else:
    sys.exit("input sequence file format is not correct, must be .fasta")

with open(reffile, "r") as rfp:
    for line in rfp:
        line = line.strip()
        linematch = re.search("^>(\S+)", line)
        if linematch:
            refname = linematch.group(1)
        else:
            refseq += line

with open(inreffile, "w") as irfp:
    with open(excludefile, "a") as efp:
        irfp.write(">" + refname + "\n" + refseq + "\n")
        with open(infile, "r") as ifp:
            for line in ifp:
                line = line.strip()
                linematch = re.search("^>(\S+)", line)
                if linematch:
                    name = linematch.group(1)
                    aaseq = ""
                    namematch = re.search('(.*)_(\d+)$', name)
                    if namematch:
                        count = namematch.group(2)
                        totalcount += int(count)
                    else:
                        sys.exit("sequence name is not formatted: " + name)
                else:
                    line = line.replace('-', '')
                    nts = list(line)
                    seqlen = len(nts)
                    for i in range(0, seqlen, 3):
                        if (i + 2 <= seqlen - 1):
                            codon = nts[i] + nts[i+1] + nts[i+2]
                            aaseq += translate(codon)

                    aaseqmatch = re.search('\*', aaseq)
                    if aaseqmatch:
                        seqname = id + "_" + name
                        efp.write(">"+seqname+"\n"+aaseq+"\n")
                        print(seqname+", "+aaseq)
                        stopcodoncount += int(count)
                    else:
                        irfp.write(">" + name + "\n" + aaseq + "\n")

musclecommand = "muscle -quiet -in " + inreffile + " -out " + inrefalignfile
os.system(musclecommand)

with open(inrefalignfile) as irafp:
    for line in irafp:
        line = line.strip()
        linematch = re.search('^>(\S+)', line)
        if linematch:
            name = linematch.group(1)
            nameseq[name] = ""
            if name != refname:
                names.append(name)
                namematch = re.search('(.*)_(\d+)$', name)
                if namematch:
                    count = namematch.group(2)
                    validcount += int(count)
                    namecount[name] = int(count)
        else:
            nameseq[name] += line

refaas = list(nameseq[refname])
idx, idx1, idx2, idx3 = 0, 0, 0, 0
for i in range(len(refaas)):
    aamatch = re.match('[A-Z]', refaas[i])
    if aamatch:
        idx += 1
    if idx == 11:
        idx1 = i
    elif idx == 24:
        idx2 = i
    elif idx == 25:
        idx3 = i
        break

for name in names:
    aas = list(nameseq[name])
    if (aas[idx1] == "R" or aas[idx1] == "K" or aas[idx1] == "H"):
        x4count += namecount[name]
        countat11 += namecount[name]
    elif (aas[idx3] == "R" or aas[idx3] == "K" or aas[idx3] == "H"):
        countat25 += namecount[name]
        x4count += namecount[name]
    elif (aas[idx2] == "R" or aas[idx2] == "K" or aas[idx2] == "H"):
        countat24 += namecount[name]
        x4count += namecount[name]
x4percent = float(x4count) / validcount * 100

with open(outfile, "a") as ofp:
    ofp.write(infile+","+str(totalcount)+","+str(stopcodoncount)+","+str(validcount)+","+str(x4count)+","+str(x4percent))
    ofp.write(","+str(countat11)+","+str(countat25)+","+str(countat24)+"\n")

print(infile+","+str(totalcount)+","+str(stopcodoncount)+","+str(validcount)+","+str(x4count)+","+str(x4percent)
      +","+str(countat11)+","+str(countat25)+","+str(countat24))



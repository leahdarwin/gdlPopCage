#!/usr/bin/env python
'''
Author: R. Axel W. Wiberg, Leah Darwin (added function for my experiment)

This script takes a popoolation2 .sync file and runs either a G-test,
Binomial GLM, or QuasiBinomial GLM for each SNP.

The script is modeled after cmh.pl from the Popoolation2 package.

.sync file count part format needs to be:

pop1-treat1      pop-1_treat2     pop-2_treat1    pop-2_treat2  ...\n
1:0:6:155:0:0   2:0:9:239:0:0    1:0:3:127:0:0   0:0:5:155:0:0  ...\n

...

Output format:

[SNP RAW DATA] treatment_line_effect

'''

import sys
import argparse
from argparse import RawTextHelpFormatter
import os
import csv
import operator

def GetMajorAlleles(cnts,mincnt,minc):
    # checks if a SNP is truly biallelic across all pops
    # returns the major and minor alleles across all pops
    # also checks if there are too many INDELS
    Acnt = 0
    Tcnt = 0
    Ccnt = 0
    Gcnt = 0
    Dcnt = 0
    allele_i = {"A":0, "T":1, "C":2, "G":3,"INDEL":5}
    all_cnt = 0
    for i in cnts:
        alleles = i.split(":")
        Acnt = Acnt + int(alleles[allele_i['A']])
        Tcnt = Tcnt + int(alleles[allele_i['T']])
        Ccnt = Ccnt + int(alleles[allele_i['C']])
        Gcnt = Gcnt + int(alleles[allele_i['G']])
        Dcnt = Dcnt + int(alleles[allele_i['INDEL']])
    # only count an allele if the number of reads >= mincnt
    if Acnt >= mincnt:
        all_cnt = all_cnt + 1
    if Tcnt >= mincnt:
        all_cnt = all_cnt + 1
    if Ccnt >= mincnt:
        all_cnt = all_cnt + 1
    if Gcnt >= mincnt:
        all_cnt = all_cnt + 1
    # if there are not enough alleles, the position is not a SNP
    #print Acnt, Tcnt, Ccnt, Gcnt, all_cnt
    if all_cnt <= 1:
        return "not SNP" #not enough alleles
    counts = [("A",Acnt), ("T",Tcnt), ("C",Ccnt), ("G",Gcnt)]
    counts.sort(key=operator.itemgetter(1))
    # last two alleles are the major and minor alleles
    major_alleles = counts[-2:]
    #print major_alleles, mincnt
    # if Dcnt has >= mincnt then allele is "tainted" by INDELS
    if Dcnt >= mincnt:
        return "not SNP"# SNP tainted by INDELS
    #print Dcnt
    # if the *third* most common allele has count >= mincnt then the SNP
    # is not truly biallelic. There are more than 2 alleles.
    #print counts
    if counts[1][1] >= mincnt:
        #print counts[1][1]
        return "not SNP"# SNP is not biallelic
    # if the *second* most common allele has count < mincnt then the
    # allele is fixed
    #print major_alleles
    for i in major_alleles:
        if i[1] < mincnt:
            return "not SNP"# major allele considered fixed
    return major_alleles


def checkSNP(line,maxc,minc,mincnt):
    # checks the SNP line for coverage
    # 1) each population has at least minc
    # 2) no population has > maxc
    # 3) no indels
    line = [i.replace("\n","") for i in line.split("\t")]
    cnts = line[3:]
    #print cnts #script tester line
    #print maxc, minc #script tester line
    for i in range(0,len(cnts)):
        pop = [int(j) for j in cnts[i].split(":")]
        #print pop, sum(pop), sum(pop) > maxc, sum(pop) < minc #script tester line
        # If the coverage within a population is > maxc or < minc
        # then the SNP is considered invalid
        # coverage is counted across As,Ts,Cs and Gs,
        # Ns and INDELs *not* counted
 	if sum(pop[:-2]) > maxc or sum(pop[:-2]) < minc:
            return "not SNP"#: coverage too high or too low
        
    return cnts

def printRlines(cnts,major_alleles,npops,nlevels,n,\
                mincnt,line,rescale,scale,zeroes):
    # prints a line that will run a glm in R and a line that will
    # print the results of that glm
    counts = []
    totcounts = []
    if rescale == 'nr':
    # don't rescale the data
    # don't change zeroes to 1s
        for p in range(0,len(cnts),nlevels):
            pop=cnts[p:p+nlevels]
            #print pop #script tester line
            samples = []
            for sam in range(0,nlevels):
                exec "sample_%s=pop[sam]" %(sam)
            	exec "samples.append(sample_%s)" %(sam)
		#exec "print sample_%s" %(i)
                allele_i = {"A":0, "T":1, "C":2, "G":3}
                #print samples #script tester line
                # Get alleles counts from each sample
            for sam in samples:
                #print sam #script tester line
                alleles = sam.split(":")
                counts.append(alleles[allele_i[major_alleles[1][0]]])
            for sam in samples:
                alleles = sam.split(":")
                counts.append(alleles[allele_i[major_alleles[0][0]]])
                        
    if rescale == 'r':
    # rescale the data
        scale=float(scale)
        #print scale
        for p in range(0,len(cnts),nlevels):
            pop=cnts[p:p+nlevels]
            #print pop #script tester line
            samples = []
            for sam in range(0,nlevels):
                exec "sample_%s=pop[sam]" %(sam)
                exec "samples.append(sample_%s)" %(sam)
                #exec "print sample_%s" %(i)
                allele_i = {"A":0, "T":1, "C":2, "G":3}
                #print samples #script tester line
                # Get alleles counts from each sample
            for sam in samples:
                #print sam #script tester line
                alleles = sam.split(":")
                asum=int(alleles[allele_i[major_alleles[1][0]]])+int(alleles[allele_i[major_alleles[0][0]]])
                afreq=float(alleles[allele_i[major_alleles[1][0]]])/float(asum)
                #print "ASUM: ", asum, "SCALE: ", scale, "ALLELE FREQUENCY: ", afreq # script tester line
                ac=afreq*scale
                counts.append(str(int(round(ac))))        
            for sam in samples:
                alleles=sam.split(":")
                asum=int(alleles[allele_i[major_alleles[1][0]]])+int(alleles[allele_i[major_alleles[0][0]]])
                afreq=float(alleles[allele_i[major_alleles[0][0]]])/float(asum)
                #print "ASUM: ", asum, "SCALE: ", scale, "ALLELE FREQUENCY: ", afreq # script tester line
                ac=afreq*scale
                counts.append(str(int(round(ac))))
    if rescale == 'neff':
    # rescale the data
        #print scale
        for p in range(0,len(cnts),nlevels):
            pop=cnts[p:p+nlevels]
            #print pop #script tester line
            samples = []
            for sam in range(0,nlevels):
                exec "sample_%s=pop[sam]" %(sam)
                exec "samples.append(sample_%s)" %(sam)
                #exec "print sample_%s" %(i)
                allele_i = {"A":0, "T":1, "C":2, "G":3}
                #print samples #script tester line
                # Get alleles counts from each sample
            for sam in samples:
                # Deal with Major Allele
                #print sam #script tester line
                alleles = sam.split(":")
                # Calculate the coverage at the SNP (asum)
                asum=int(alleles[allele_i[major_alleles[1][0]]])+int(alleles[allele_i[major_alleles[0][0]]])
                afreq=float(alleles[allele_i[major_alleles[1][0]]])/float(asum)
                #print "ASUM: ", asum, "SCALE: ", scale, "ALLELE FREQUENCY: ", afreq # script tester line
                # Calculate neff and rescale the counts 
                neff=(asum*(n*2)-1)/((n*2)+asum)
                ac=afreq*neff
                counts.append(str(int(round(ac))))        
            for sam in samples:
                # Deal with Minor Allele
                alleles=sam.split(":")
                asum=int(alleles[allele_i[major_alleles[1][0]]])+int(alleles[allele_i[major_alleles[0][0]]])
                afreq=float(alleles[allele_i[major_alleles[0][0]]])/float(asum)
                #print "ASUM: ", asum, "SCALE: ", scale, "ALLELE FREQUENCY: ", afreq # script tester line
                neff=(asum*(n*2)-1)/((n*2)+asum)
                ac=afreq*neff
                counts.append(str(int(round(ac))))        
                #print "NEFF: ", neff, "AC: ", ac, "ASUM: ", asum, "N(n*2): ", n, "(",n*2,")"

    #print counts #script tester line
    counts=','.join(counts)
    print 'matrix<-array(c('+counts+'),'+\
          'dim=c('+str(nlevels)+',2,'+str(npops)+'))'
    #print line #script tester line
    print 'dat<-get_dat(matrix,zeroes='+zeroes+')'
    #print 'print(dat)'#Script tester line
    print 'res<-glm(cbind(A_Cnt,Tot_Cnt-A_Cnt)~tr_l*time+rep+tr_l:rep,family="quasibinomial",data=dat)'
    print 'mod<-anova(res, test="LRT")'
    #print 'print(summary(res))'#Script tester line
    # The treatment and additionaly variable results start after npops rows of the summary
    #print 'n_rows<-nrow(summary(res)$coefficients)'
    print 'cat(c("'+line.replace('\n','')+'"'+\
          ',mod$`Pr(>Chi)`[2]'+\
          ',mod$`Pr(>Chi)`[3]'+\
          ',mod$`Pr(>Chi)`[4]'+\
          ',mod$`Pr(>Chi)`[5]'+\
          ',mod$`Pr(>Chi)`[6])'+\
          ',sep="\\t","\\n")'

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'],formatter_class=RawTextHelpFormatter)

#MAIN CONTROLS
    parser.add_argument('-filename',
                       help = '.sync file name')
  
    parser.add_argument('-npops',
                        help = 'Number of populations.')

    parser.add_argument('-nlevels',
                        help = 'Number of levels in the population.')

    parser.add_argument('-n',
                        help = 'Number of individuals in the pool.')

    parser.add_argument('-mincnt', default = 16,
                        help = 'Minimum count needed to consider a SNP.'+\
                        ' default = 16')

    parser.add_argument('-minc', default = 30,
                        help = 'Minimum coverage needed to consider a SNP.'+\
                        ' default = 30')

    parser.add_argument('-maxc', default = 400,
                        help = 'Maximum coverage needed to consider a SNP.'+\
                        ' default = 400')

    parser.add_argument('-rescale', default = 'nr',
                        help = 'r = rescale to -scale, '+\
                        'nr = no rescaling, neff = rescale to neff')

    parser.add_argument('-scale',default = 'XX',
                        help = 'The new total count for rescaling')

    parser.add_argument('-zeroes', default = 1,
                        help = 'Add (1) or not (0) 1 to each cell if any'+\
                        'cell has count 0')

    args = vars(parser.parse_args())

    try:
        filnam = args['filename']
        npops = int(args['npops'])
        nlevels = int(args['nlevels'])
        n = int(args['n'])
        mincnt = int(args['mincnt'])
        minc = int(args['minc'])
        maxc = int(args['maxc'])
        rescale = args['rescale']
        scale = args['scale']
        zeroes=str(args['zeroes'])
        lines = open(filnam, 'rb')
        scriptsdir=os.path.dirname(os.path.realpath(__file__))
        print 'suppressWarnings(library(methods))'
        print 'currdir<<-"'+scriptsdir+'"'
        print 'source(paste(currdir,"/poolFreqDiff_AF.R",sep=""))'
        print '#Parameters: ',"npops =",npops,"nlevels =",nlevels,"mincnt =",mincnt,"min coverege =",minc,\
              "max coverage =",maxc,"rescale =",rescale,"scale =",scale,"zeroes =",zeroes
        print '#Pvalues order: tr_l, time, rep, tr_l:time, tr_l:rep'
        for line in lines:
            SNP = checkSNP(line,maxc,minc,mincnt)
            #print SNP #script tester line
            if SNP != "not SNP":
                #print SNP #script tester line
                major_alleles = GetMajorAlleles(SNP,mincnt,minc)
                if major_alleles != "not SNP":
                    #print major_alleles #script tester line
                    #print SNP #script tester line
                    printRlines(SNP,major_alleles,npops,nlevels,n,mincnt,\
                                line,rescale,scale,zeroes)
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass

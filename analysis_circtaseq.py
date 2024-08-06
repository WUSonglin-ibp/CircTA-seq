# Copyright 2024 Institute of Biophysics, Chinese Academy of Sciences
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

##########################
# CircTA-seq data analysis
##########################
# Citation: 
# Weidong An, Yunxiao Yan and Keqiong Ye
# High resolution landscape of ribosomal RNA processing and surveillance. 
# Nucleic Acids Research, 2024,  https://doi.org/10.1093/nar/gkae606
##########################
# Basic steps:
# Remove adaptor 
# Align reads with blastn. Normal reads should be mapped to two different sites on 35S
# Extract coordinates of 5' and 3' end of RNA fragment
# Convert 5' and 3' end coordinate to closest processing sites (array: site_35S_5end, site_35S_3end, =>sites)
# Convert sites to a defined pre-rRNA, otherwise as UNK. (prerrna_dict)
# Resolve overlapped junctions
# Statistics and output
##########################
# Remove adaptor and trim low-quality tail, convert to fasta
#  cutadapt -a TGGAATTCTCGGGTGCC -A AGATCGGAAGAGCGTCG  -e 0.1 -m 25 --overlap 5 -q 20 -o ${x}_trimq_1.fq -p ${x}_trimq_2.fq ${x}_R1.fq.gz ${x}_R2.fq.gz > cutadapt_trimq.out
#  seqtk seq -A ${x}_trimq_1.fq > ${x}_trimq_1.fasta
#  seqtk seq -A ${x}_trimq_2.fq > ${x}_trimq_2.fasta
##########################
# making database for alignment in ../genome
# makeblastdb -in 35Sext -dbtype nucl -parse_seqids
##########################
# Alignment by blastn, input fasta files, output table (fmt 6).
#  blastn -db ../genome/35Sext -task "blastn-short" -evalue 1e-3 -strand plus -outfmt 6 -query ${x}_trimq_1.fasta -out ${x}_trimq_r1_blast.tab
#  blastn -db ../genome/35Sext -task "blastn-short" -evalue 1e-3 -strand minus -outfmt 6 -query ${x}_trimq_2.fasta -out ${x}_trimq_r2_blast.tab
##########################
# Major processing steps
#
#  python ../analysis_circtaseq.py -a blasttab -genome35s '../genome/35Sext' -tab1  ${x}_trimq_r1_blast.tab -tab2 ${x}_trimq_r2_blast.tab -r1 ${x}_trimq_1.fq -r2 ${x}_trimq_2.fq -basename b6 
#
# tab.gz fq.gz files can also be used
# Output files:
# b6_hist.stat: all statistics 
# b6_prerrna_3EXT.dat: 3' extension length, read number, fraction
# b6_prerrna_5EXT.dat: 5' extension length, read number, fraction
# b6_prerrna_PA.dat: polyA length, read number, fraction
# b6_prerrna_PAB.dat: range of polyA length, read number, fraction
# b6_all.stat: results for all aligned and normal read pairs, used for subsequent analysis
# content in b6_all.stat
#27SB	0	7	5	RandomSplit	B1-B2	-2	1.00	TCTGATTTGTtttttat||attaaAAACTTTCAA	E00492:356:HTFYMCCXY:8:1114:5518:18573_1	4100	7897
# (1) pre-rRNA species 
# (2) polyA length 
# (3) 3' end extension (+) or deletion (-) of pre-rRNA (deletion5)
# (4) 5' end extension or deletion of pre-rRNA (deletion3) 
# (5) Junction type 
# (6) Sites 
# (7) dist_err: overlapped nt (<0), original pA length (>0) for single-read, or diff_len for crossread  
# (8) fraction of adenine in polyA region
# (9) Junction sequences, 3' end (extension in lower letters)| polyA | 5' end (extension in lower letters)
# (10) read name. _1: read 1, _2RC: read 2 
# (11) 5' end in 35Sext. 
# (12) 3' end in 35Sext.

##########################
# Analyze polyA length at individual position
# All reads are processed, but only major pre-rRNAs with total reads > 2000 (totalread_cutoff) are reported
# Histogram data are reported only for species >= 20 reads (subtotalread_cutoff) 
#
#  python ../analysis_circtaseq.py -a del3 -genome35s '../genome/35Sext' -mapfile b6_all.stat -basename b6
#
# Output files: 
# b6_prerrna_end3_all.stat, b6_prerrna_hisend3_all.stat  b6_prerrna_end3_C20.stat, b6_prerrna_hisend3_C20.stat (for figure)
# end3: polyA profile for all polyA length
# hisend3: histogram of polyA distribution
# all: all sites are reported
# C20: only sites with read >= 20 are reported
# b6_prerrna_hisend3_C20.stat is used for gnuplot drawing. 
# If data are incomplete, we need manually add empty data at the beginning of the file for drawing figure in gnuplot 
#
# If need to analyze polyA profile for species with a specific 5' extension (such as 21S', 25S')
# first extract reads containing all such 5' extension, then do the above analysis. 
# awk 'BEGIN {OFS="\t"} $1=="21S" && $4==0' < b6_all.stat > b6_21S_5ext0.stat
##########################
# Count pre-rRNA species with same 5' and 3' end coordinates. for 2D plot
# We need only 100k reads for this kind of analysis. 
#
#  head -100000 b6_all.stat > b6_100k_all.stat
#  python ../analysis_circtaseq.py -a unk -genome35s '../genome/35Sext' -mapfile b6_100k_all.stat -basename b6
#  sort -k3,3n b6_UNK_coord.dat > b6_UNK_coord.dat.sort
#
##########################

from __future__ import division
from Bio import SeqIO
import re
import collections
import numpy as np
import random
from random import randint
import argparse
import gzip

# Global parameters 
f_adenine_cutoff = 0.9  # cutoff for A level in polyA region
pident_cutoff = 90 # cutoff for alignment identity score

# read number cutoff for outputing invidiual ending data for a pre-rRNA
totalread_cutoff = 2000
#read number cutoff for analyzing polyA distribution at individual position
subtotalread_cutoff = 20  # output file is labled as '_C20'

########################################
# 
########################################
#crtseqs stores retrieved and derived information for blast alignments. 
#crtseqs is dict of dict
#crtseqs = {'seq_id': {'seq_id':'xxx','polyA':2, 'deletion3':0, 'deletion5':1, 'prerRNA':'20S', 'junction':'XXXaaa|NNNNN' ...} ...}
crtseqs = collections.defaultdict(dict)

# stat is dict of dict, stores polyA distrubtion of each prerRNA
# stat = {'20S': {polyA_length: number of sequence, 0:309, 1:20, ...} ...}
stat = collections.defaultdict(dict)

# stat_perfect {'20S': prefect matched seq}
stat_perfect = collections.defaultdict(dict)

# statistics of deletion3, dict of dict
# stat_del3 = {'20S': {deletion_length:number of sequence, 0:309, 1:20, ...} ...}
stat_del3 = collections.defaultdict(dict)
# statistics of deletion5, dict of dict
stat_del5 = collections.defaultdict(dict)

#################################################
# Yeast pre-rRNA specific setting, Start
#################################################

#pre-rRNA species
prerRNA_list= ["35S", "33S", "32S", "23S", "22S", "21S", "20S", "18S", "27SA2", "27SA3", "27SB", "7S","58S", "26S", "25S","UNK"]

# coordinate =  nt 3' (right) of a processing site, 1-base
# template: 35Sext (8101nt = NTS 1243nt + RND37-1 6858nt)
# 5ETS-A0 610, A0-A1 90, 18S 1800, D-A2 212, A2-A3 72, A3-B1 77, 5.8S/B1-E 158, E-C2 139, C2-C1 95, 25S/C1-B2 3394, 3ETS 210
site = {'START':1, '5ETS':1244, 'A0':1854, 'A1':1944, 'D':3744, 'A2':3956, 'A3':4028, 'B1': 4105, 'E': 4263, 'C2': 4402, 'C1': 4497, 'B2':7891, 'END': 8101}
# template 35S, not used
#site = {'START':-1242, '5ETS':1, 'A0':611, 'A1':701, 'D':2501, 'A2':2713, 'A3':2785, 'B1': 2862, 'E': 3020, 'C2':3159, 'C1': 3254, 'B2':6648, 'END': 6858}

## prerrna_dict must correspond to exactly the prerRNA_list
prerrna_dict = {'5ETS-B2':'35S', 'A0-B2':'33S', 'A1-B2':'32S','5ETS-A3':'23S', 'A0-A3':'22S', 'A1-A3':'21S', 'A1-A2':'20S', 'A1-D':'18S', 'A2-B2':'27SA2', 'A3-B2':'27SA3', 'B1-B2':'27SB', 'B1-E':'58S', 'B1-C2':'7S', 'C2-B2':'26S', 'C1-B2':'25S'}
prerrna_3end = {'35S':'B2', '33S':'B2', '32S':'B2','23S':'A3', '22S':'A3', '21S':'A3', '20S':'A2', '18S':'D', '27SA2':'B2', '27SA3':'B2', '27SB':'B2', '58S':'E', '7S':'C2', '26S':'B2', '25S':'B2'}
prerrna_5end = {'35S':'5ETS', '33S':'A0', '32S':'A1','23S':'5ETS', '22S':'A0', '21S':'A1', '20S':'A1', '18S':'A1', '27SA2':'A2', '27SA3':'A3', '27SB':'B1', '58S':'B1', '7S':'B1', '26S':'C2', '25S':'C1'}

def set_site35S_5end ():
  """ mapping coordinates (1-base from blast) to processing site, site_35S_5end is a list. For 5 end mapping """
  global site_35S_5end
  site_35S_5end = ['NO']*8102  # 0 is not used
  for i in range(0, site['A0']-10):        # 35S, 23S, 5ETS fragment
    site_35S_5end[i] = "5ETS"
  for i in range(site['A0']-10, site['A1']-10):  # 33S, 22S
    site_35S_5end[i] = "A0"  
  for i in range(site['A1']-10, site['A1']+50): # 32S, 21S, 20S, 18S
    site_35S_5end[i] = "A1"
  for i in range(site['A1']+50, site['A2']-10): # degradation, span most 18S region
    site_35S_5end[i] = "D"
  for i in range(site['A2']-10, site['A3']-10): # 27SA2 and its 5' degradation products
    site_35S_5end[i] = "A2"
  for i in range(site['A3']-10, site['A3']+10): # 27SA3
    site_35S_5end[i] = "A3"
  for i in range(site['A3']+10, site['B1']+50): # processing products of 27SA3->27SB, 5.8S, 5' degradation
    site_35S_5end[i] = "B1"
  for i in range(site['B1']+50, site['C2']-10): # degradation, span most 5.8S
    site_35S_5end[i] = "E"
  for i in range(site['C2']-10, site['C1']-20): # 26S
    site_35S_5end[i] = "C2"
  for i in range(site['C1']-20, site['C1']+50): # 25S
    site_35S_5end[i] = "C1"
  for i in range(site['C1']+50, site['END']+1): # Nothing
    site_35S_5end[i] = "B2"

def set_site35S_3end ():
  """ mapping coordinates (1-base from blast) to processing site, site_35S_3end is a list. For 3 end mapping """
  global site_35S_3end
  site_35S_3end = ['NO']*8102  # 0 is not used
  for i in range(0, site['5ETS']+10):        # Nothing
    site_35S_3end[i] = "5ETS"
  for i in range(site['5ETS']+10, site['A0']+10):  # 5ETS-A0 fragment
    site_35S_3end[i] = "A0"  
  for i in range(site['A0']+10, site['D']-50): # degradation, span 18S, or A0-A1 fragment
    site_35S_3end[i] = "A1"
  for i in range(site['D']-50, site['D']+10): # 18S
    site_35S_3end[i] = "D"
  for i in range(site['D']+10, site['A2']+10): # 20S
    site_35S_3end[i] = "A2"
  for i in range(site['A2']+10, site['A3']+10): # 23S
    site_35S_3end[i] = "A3"
  for i in range(site['A3']+10, site['E']-50): # degradtion, span 5.8S
    site_35S_3end[i] = "B1"
  for i in range(site['E']-50, site['E']+50): # 5.8S-6S
    site_35S_3end[i] = "E"
  for i in range(site['E']+50, site['C2']+10): # 7S, arbitary division from 6S
    site_35S_3end[i] = "C2"
  for i in range(site['C2']+10, site['B2']-50): # degradation, span 25S
    site_35S_3end[i] = "C1"
  for i in range(site['B2']-50, site['END']+1): # 35S, 33S, 32S, 27S, 25S
    site_35S_3end[i] = "B2"

#################################################
# Yeast pre-rRNA specific setting, End
#################################################

def fraction_adenine (seq):
  count=0
  if len(seq) == 0:
    return 1
  else: 
    for x in seq:
      if x =="A":
        count +=1
    return count/len(seq)
    
def safe_divide (a,b):
  if b == 0:
    return 0
  else: 
    return a/b

def block_histogram(stat_dict,key1,bin_edges=[0,1,4,10,21,201]):
  # stat_dict[key1][length] is nested dictionary holding histogram of polya length
  # The first key could be certain class of RNA, such as prerRNA, 3' deletion/extension length
  # The second key must be length of pA
  hist = [] 
  for i in range(0,len(bin_edges)-1):
    hist.append(0)
    for length in (set(range(bin_edges[i],bin_edges[i+1],1)) & set(stat_dict[key1].keys())): 
      hist[i] += stat_dict[key1][length]
  return hist, bin_edges

def write_block_histogram(output_handle, hist, bin_edges, total_hit):
  for i in range(0,len(bin_edges)-1):
    output_handle.write ('%s-%s\t%d\t%f\n' % (bin_edges[i],bin_edges[i+1]-1,hist[i],hist[i]/float(total_hit)))

def write_block_histogram2(output_handle, hist, bin_edges, total_hit, deletion_length):
  for i in range(0,len(bin_edges)-1):
    output_handle.write ('%d\t%d\t%d\t%f\t%s-%s\n' % (deletion_length, i, hist[i], hist[i]/float(total_hit), bin_edges[i], bin_edges[i+1]-1))

def output_stat ():
  global crtseqs
  seq_mapped=0
  seq_total=0
  crossread_num=0
  singleread_num=0
  seq_perfectmatch =0
  for prerRNA in prerRNA_list:
    stat_perfect[prerRNA] = 0
  for seq_id, v in crtseqs.items():
    #print v
    seq_total +=1
    if v['prerRNA'] != "NO" and v['frac_a'] >=f_adenine_cutoff:
      seq_mapped +=1
      if v['type_junc'][:9] == "CrossJunc":
        crossread_num +=1
      else:
        singleread_num +=1

      prerRNA=v['prerRNA']
      if v['score5'] == 100 and v['score3'] == 100:
        seq_perfectmatch +=1
        stat_perfect[prerRNA] +=1

      polya_length =int (v['polya'])
      if polya_length in stat[prerRNA]: 
        stat[prerRNA][polya_length] +=1
      else: 
        stat[prerRNA][polya_length] =1

      deletion3_length =int (v['deletion3'])
      if deletion3_length in stat_del3[prerRNA]: 
        stat_del3[prerRNA][deletion3_length] +=1
      else: 
        stat_del3[prerRNA][deletion3_length] =1

      deletion5_length =int (v['deletion5'])
      if deletion5_length in stat_del5[prerRNA]: 
        stat_del5[prerRNA][deletion5_length] +=1
      else: 
        stat_del5[prerRNA][deletion5_length] =1
  print ("mapped seq=", seq_mapped, "total seq=", seq_total)
  
  # histogram of polyA length, by prerRNA
  out_handle = open('%s_hist.stat'%(basename), 'w')
  out_handle.write ('## Input files %s, %s, prefix= %s\n' %(args.r1, args.r2, args.basename))
  out_handle.write ('A fraction cutoff >= %.2f, pident cutoff> %.1f\n' %(f_adenine_cutoff, pident_cutoff))
  out_handle.write ('mapped sequences = %d, total seq =%d, mapped ratio=%.3f, prefect match=%d, prefect ratio=%.3f\n'%(seq_mapped, seq_total, seq_mapped/seq_total, seq_perfectmatch, safe_divide(seq_perfectmatch,seq_mapped)))
  out_handle.write ('crossread = %d, singleread=%d, crossread ratio=%.3f, singleread ratio=%.3f  \n'%(crossread_num, singleread_num, crossread_num/seq_mapped, singleread_num/seq_mapped))
  out_handle.write ('###################################################\n' )
  out_handle.write ('#######  Length of Poly A       ##################\n' )
  out_handle.write ('Length\tNumber\tRatio\n')
  for prerRNA in prerRNA_list:
    total_hit = sum (stat[prerRNA].values())
    if total_hit > 0:
      each_handle = open('%s_%s_PA.dat'%(basename,prerRNA), 'w')
      out_handle.write ('####### PolyA %s Total hit=%d\n' % (prerRNA, total_hit))
      for length in sorted(stat[prerRNA].keys()):
        out_handle.write ('%d\t%d\t%f\n'%(length,stat[prerRNA][length],stat[prerRNA][length]/total_hit))
        each_handle.write ('%d\t%d\t%f\n'%(length,stat[prerRNA][length],stat[prerRNA][length]/total_hit))
      each_handle.close()

      ## block statistics
      if max(stat[prerRNA].keys()) > 1:
        each_handle = open('%s_%s_PAB.dat'%(basename,prerRNA), 'w')
        out_handle.write ('\n>>>>>>> Block statistics of  %s \n' % (prerRNA))
        hist,bin_edges = block_histogram(stat_dict=stat, key1=prerRNA)
        write_block_histogram(each_handle, hist, bin_edges, total_hit)
        write_block_histogram(out_handle, hist, bin_edges, total_hit)
        each_handle.close()

  # histogram of deltion3 length, by prerRNA
  out_handle.write ('\n###################################################\n' )
  out_handle.write ('#######  5 end extension of RNA  ##################\n' )
  for prerRNA in prerRNA_list:
    total_hit = sum (stat_del3[prerRNA].values())
    if total_hit > 0:
      each_handle = open('%s_%s_5EXT.dat'%(basename,prerRNA), 'w')
      out_handle.write ('#######5 EXT %s Total hit=%d\n' % (prerRNA, total_hit))
      for length in sorted(stat_del3[prerRNA].keys()):
        out_handle.write ('%d\t%d\t%f\n'%(length,stat_del3[prerRNA][length],stat_del3[prerRNA][length]/total_hit))
        each_handle.write ('%d\t%d\t%f\n'%(length,stat_del3[prerRNA][length],stat_del3[prerRNA][length]/total_hit))
      each_handle.close()

  # histogram of deltion5 length, by prerRNA
  out_handle.write ('\n###################################################\n' )
  out_handle.write ('#######  3 end extension of RNA  ##################\n' )
  for prerRNA in prerRNA_list:
    total_hit = sum (stat_del5[prerRNA].values())
    if total_hit > 0:
      each_handle = open('%s_%s_3EXT.dat'%(basename,prerRNA), 'w')
      out_handle.write ('#######3 EXT %s Total hit=%d\n' % (prerRNA, total_hit))
      for length in sorted(stat_del5[prerRNA].keys()):
        out_handle.write ('%d\t%d\t%f\n'%(length,stat_del5[prerRNA][length],stat_del5[prerRNA][length]/total_hit))
        each_handle.write ('%d\t%d\t%f\n'%(length,stat_del5[prerRNA][length],stat_del5[prerRNA][length]/total_hit))
      each_handle.close()
    else:
      pass
  out_handle.close()

  # Output all aligned reads
  out_handle = open('%s_all.stat'%(basename), 'w')
  print ("writing %s_all.stat" % basename)
  for seq_id, v in crtseqs.items():
    # some seq_id don't have mapped data
    if v['prerRNA'] != "NO" and v['frac_a'] >=f_adenine_cutoff:
      m = re.match ('(\w*)-(\w*)', v['sites'])
      if m:
        x5=site[m.group(1)]-v['deletion3']
        x3=site[m.group(2)]+v['deletion5']-1
      out_handle.write ('%s\t%d\t%d\t%d\t%s\t%s\t%s\t%.2f\t%s\t%s\t%d\t%d\n'%(v['prerRNA'],v['polya'],v['deletion5'],v['deletion3'],v['type_junc'],v['sites'],v['dist_err'],v['frac_a'],v['junction'],v['seq_id'],x5,x3))
  out_handle.close()

def read_map (inputfile):
#22S     3       0       0       NormalJunc      A0-A3   -2      1.00    GAGGTAACAA|aaa|CTTCTAGCAA       E00513:268:HNCJHCCXY:3:1124:9191:60940_1
  #print inputfile
  with open (inputfile) as f:
    for line in f:
      w = line.split('\t')
      seqid = w[9].rstrip()
      #print line, seqid
      if re.match("SeqID", seqid):
        pass
      else:
        crtseqs[seqid]['prerRNA']=w[0]
        crtseqs[seqid]['polya']=int(w[1])
        crtseqs[seqid]['deletion5']=int(w[2])
        crtseqs[seqid]['deletion3']=int(w[3])
        m = re.match ('(\w*)-(\w*)', w[5])
        if m:
          crtseqs[seqid]['end3']=m.group(2)
          crtseqs[seqid]['end5']=m.group(1)
        m1 = re.match ('(\w*)\|(\w*)\|(\w*)', w[8])
        if m1: 
          crtseqs[seqid]['polyaseq']=m1.group(2)
          #print m1.group(1), m1.group(2), m1.group(3)

def output_unk():
# for gnuplot 2D figure
  """ output 5 and 3 end coordinates of species with at least depth_cutoff reads """
  depth_cutoff = 0
  unkpeak = collections.defaultdict(dict)
  highpeak = collections.defaultdict(dict)
  for seq_id, v in crtseqs.items():
    if True:  
    #if v['prerRNA'] == "UNK":  
      x5=site[v['end5']]-v['deletion3']
      x3=site[v['end3']]+v['deletion5']-1
      if x5 in unkpeak.keys():
        if x3 in unkpeak[x5].keys():
          unkpeak[x5][x3] +=1
        else:
          unkpeak[x5][x3] =1
      else:
        unkpeak[x5][x3] =1
  out_handle = open('%s_%s_coord.dat'%(basename, "UNK"), 'w')
  x3_list = []
  depth_list = []
  for x5 in sorted(unkpeak.keys()):
    for x3 in sorted(unkpeak[x5].keys()):
      if unkpeak[x5][x3] > depth_cutoff:        
        depth_list.append (unkpeak[x5][x3])
        if x3 not in x3_list:
          x3_list.append(x3)
        highpeak[x5][x3] = unkpeak[x5][x3]
        out_handle.write ('%d\t%d\t%d\n'%(x5-1243, x3-1243, unkpeak[x5][x3]))
  peak_max = max(depth_list)
  print ("reads of the most abundant species: ", peak_max)
  

def output_del3_new ():
# for gnuplot drawing
  prerRNA_list.remove("UNK")  # don't count UNK
  for prerRNA in prerRNA_list:
    seq_total =0
    # stat_del is dict of dict {{deletion5: {polya_length: number of sequence}}}
    stat_del = collections.defaultdict(dict)  # re-initialize, important
    for seq_id, v in crtseqs.items():
      if v['prerRNA'] == prerRNA:
        seq_total +=1
        #print v['prerRNA'], prerRNA
        polya_length =int (v['polya'])
        deletion_length =int (v['deletion5'])
        if deletion_length in stat_del.keys():
          if polya_length in stat_del[deletion_length].keys():
            stat_del[deletion_length][polya_length] +=1
          else:
            stat_del[deletion_length][polya_length] =1
        else:
          stat_del[deletion_length][polya_length] =1
    if seq_total > totalread_cutoff:
      #out_handle = open('%s_%s_end3.stat'%(basename, prerRNA), 'w')  # single output
      #hist_handle = open('%s_%s_histend3.stat'%(basename, prerRNA), 'w')
      out_handle = open('%s_%s_end3_C20.stat'%(basename, prerRNA), 'w')  # subtotal cutoff
      hist_handle = open('%s_%s_histend3_C20.stat'%(basename, prerRNA), 'w') # subtotal cutoff
      out2_handle = open('%s_%s_end3_all.stat'%(basename, prerRNA), 'w')  # all positions
      hist2_handle = open('%s_%s_histend3_all.stat'%(basename, prerRNA), 'w') # all positions
      out_handle.write ('# %s \n' % (prerRNA))
      hist_handle.write ('# %s \n' % (prerRNA))
      out2_handle.write ('# %s \n' % (prerRNA))
      hist2_handle.write ('# %s \n' % (prerRNA))
      total_sub = 0
      total_rna = 0
      length_list = sorted(stat_del.keys(), reverse = False)
      for deletion_length in range(length_list[0],length_list[-1]+1,1):   # change the order , what if one length is missing
        if deletion_length in length_list:
          total_sub = sum (stat_del[deletion_length].values())
        else:
          total_sub = 0
        total_rna += total_sub
        if prerRNA == "UNK":
          deleted_nt = "N"
        else:
          coord_end = site[prerrna_3end[prerRNA]]+deletion_length -2  # nt at 3' end exact. any base, random split
          if coord_end < len(RNA35S):   # a safty check
            deleted_nt = RNA35S[coord_end] 
          else:
            deleted_nt = "N"
        combine_label = str(deletion_length) + deleted_nt + "/" + str(total_sub)
        combine_label2 = ('%+d' % deletion_length) + "/" + str(total_sub)
        if deletion_length in length_list:
          if (total_sub >= 0):    # output everything
            out2_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
            out2_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
            for polya_length in sorted(stat_del[deletion_length].keys()):
              out2_handle.write ('%d\t%d\t%d\t%f\n'%(deletion_length, polya_length, stat_del[deletion_length][polya_length], float(stat_del[deletion_length][polya_length])/float(total_sub)))
            out2_handle.write ('\n\n'% ())
  
            hist2_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
            hist2_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
            hist,bin_edges = block_histogram(stat_dict=stat_del, key1=deletion_length)
            write_block_histogram2(hist2_handle, hist, bin_edges, total_sub, deletion_length)
            hist2_handle.write ('\n\n'% ())
          else:   # not enough reads, output blank, total_sub > 0
            out2_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
            out2_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
            out2_handle.write ('\n\n'% ())
            hist2_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
            hist2_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
            hist2_handle.write ('\n\n'% ())
 
          if (total_sub >= subtotalread_cutoff):
            out_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
            out_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
            for polya_length in sorted(stat_del[deletion_length].keys()):
              out_handle.write ('%d\t%d\t%d\t%f\n'%(deletion_length, polya_length, stat_del[deletion_length][polya_length], float(stat_del[deletion_length][polya_length])/float(total_sub)))
            out_handle.write ('\n\n'% ())
  
            hist_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
            hist_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
            hist,bin_edges = block_histogram(stat_dict=stat_del, key1=deletion_length)
            write_block_histogram2(hist_handle, hist, bin_edges, total_sub, deletion_length)
            hist_handle.write ('\n\n'% ())
          else:   # not enough reads, output blank, total_sub > 0
            out_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
            out_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
            out_handle.write ('\n\n'% ())
            hist_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
            hist_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
            hist_handle.write ('\n\n'% ())
        else:   # deletion_length is empty, make fake data, total_sub =0
          out_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
          out_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
          out_handle.write ('\n\n'% ())
          hist_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
          hist_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
          hist_handle.write ('\n\n'% ())
          out2_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
          out2_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
          out2_handle.write ('\n\n'% ())
          hist2_handle.write ('#Deletion=%d, Subtotal=%d, Deleted nt=%s\n'% (deletion_length, total_sub,deleted_nt))
          hist2_handle.write ('%s\t%s\t%s\t%s\t%s\n'%(str(deletion_length),total_sub,deleted_nt, combine_label, combine_label2))
          hist2_handle.write ('\n\n'% ())

      out_handle.write ('#%s  Total=%d\n'% (prerRNA, total_rna))
      hist_handle.write ('#%s  Total=%d\n'% (prerRNA, total_rna))
      out_handle.close()
      hist_handle.close()
      out2_handle.write ('#%s  Total=%d\n'% (prerRNA, total_rna))
      hist2_handle.write ('#%s  Total=%d\n'% (prerRNA, total_rna))
      out2_handle.close()
      hist2_handle.close()


def polya_composition():
  stat_polya = collections.defaultdict(dict)
  stat_sum = collections.defaultdict(dict)
  out_handle = open('%s_polyacomposition.stat'%(basename), 'w')
  count_U,count_G,count_C,count_all=0,0,0,0
  for seqid, v in crtseqs.items():
    if v['polya'] >= 10:
       polyaseq=v['polyaseq'].upper()
       prerRNA=v['prerRNA']
       polya_length=int(v['polya'])
       stat_polya[polya_length]['U'] = stat_polya[polya_length].get('U', 0) + polyaseq.count('T')
       stat_polya[polya_length]['C'] = stat_polya[polya_length].get('C', 0) + polyaseq.count('C')
       stat_polya[polya_length]['G'] = stat_polya[polya_length].get('G', 0) + polyaseq.count('G')
       stat_polya[polya_length]['all'] = stat_polya[polya_length].get('all', 0) + len(polyaseq)
       count_U += polyaseq.count('T')
       count_G += polyaseq.count('G')
       count_C += polyaseq.count('C')
       count_all +=len(polyaseq)
  out_handle.write ('#U=%d\tG=%d\tC=%d\tTotal=%d\n' % (count_U, count_G, count_C, count_all))
  #print ('U=%f\tG=%f\tC=%f' % (100*count_U/count_all, 100*count_G/count_all, 100*count_C/count_all))
  binsize = 20
  for k in range (10,160,binsize):
    stat_sum[k]['U'] = 0
    stat_sum[k]['G'] = 0
    stat_sum[k]['C'] = 0
    stat_sum[k]['all'] = 0
    for i in (set(range(k,k+binsize,1)) & set(stat_polya.keys())):
      stat_sum[k]['U'] +=  stat_polya[i].get('U', 0)
      stat_sum[k]['G'] +=  stat_polya[i].get('G', 0)
      stat_sum[k]['C'] +=  stat_polya[i].get('C', 0)
      stat_sum[k]['all'] +=  stat_polya[i].get('all', 0)
  out_handle.write ('#Bin\tU\tG\tC\tTotMut\tTotRes\n' % ())
  out_handle.write ('#%s\t%f\t%f\t%f\t%f\t%d\n' % ('All', 100*count_U/count_all, 100*count_G/count_all, 100*count_C/count_all,100*(count_C+count_U+count_G)/count_all, count_all))
  for i in sorted(stat_sum.keys()):
    if stat_sum[i]['all'] > 0:
      out_handle.write ('%d\t%f\t%f\t%f\t%f\t%d\n' % (i, 100*stat_sum[i]['U']/stat_sum[i]['all'], 100*stat_sum[i]['G']/stat_sum[i]['all'],100*stat_sum[i]['C']/stat_sum[i]['all'],100*(stat_sum[i]['U']+stat_sum[i]['G']+stat_sum[i]['C'])/stat_sum[i]['all'],stat_sum[i]['all']))
  out_handle.close()
      

def read_35S (infile):
  with open(infile) as f:
    return SeqIO.read(f, "fasta").seq

def read_reads (r1file, r2file):
  if re.search('.gz$', r1file):
    f1=gzip.open(r1file)
  else:
    f1=open(r1file) 
  for k in  SeqIO.parse(f1, "fastq"):
    #print k.letter_annotations
    crtseqs[k.id]['seq1'] = k.seq
    # not store phred1 to save memory
    #crtseqs[k.id]['phred1'] = k.letter_annotations["phred_quality"]
    crtseqs[k.id]['r1hitcount'] = 0
    crtseqs[k.id]['prerRNA'] = 'NO'
  f1.close()

  if re.search('.gz$', r2file):
    f2=gzip.open(r2file)
  else:
    f2=open(r2file) 
  for k in  SeqIO.parse(f2, "fastq"):
    crtseqs[k.id]['seq2'] = k.seq.reverse_complement()
    #crtseqs[k.id]['phred2'] = k.letter_annotations["phred_quality"][::-1]
    crtseqs[k.id]['r2hitcount'] = 0
  f2.close()

## alignment  by blast, 1-base coordinate
# Single-read junction
# Normal: q1 end< q2start
# [q1start------------q1end][AAAAAA][q2start---------------q2end]
#                        ext5       ext3
#
# [s2start------------s2end]XXXXXXXXXXXXXX[s1start]------------------s1end][AAAAAA]
# DEL3                   A2] Del3                                    DEL5
#
# Cross-read junction, only 1 alignment for read 1
# read 1
# [q1start------------q1end][AAAAAA]
# read 2RC
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA][q2start---------------q2end]


def read_blasttable (r1file,r2file):
# fmt 6 of blastn
#E00492:356:HTFYMCCXY:8:1114:26311:18309	35SEXT	100.000	40	0	0	1	40	7860	7899	4.58e-19	79.8
#E00492:356:HTFYMCCXY:8:1114:26311:18309	35SEXT	100.000	28	0	0	41	68	4105	4132	6.63e-12	56.0

  if re.search('.gz$', r1file):
    f1=gzip.open(r1file)
  else:
    f1=open(r1file)
  for line in f1:
    line = line.rstrip('\n')
    w = line.split('\t')
    seq_id = w[0]
    gap = int(w[5])   # 4 mut, 5 gap
    pident = float(w[2])
    crtseqs[seq_id]['seq_id'] = seq_id + "_1"
    if pident > pident_cutoff:  # percentage of identical residue
      if crtseqs[seq_id]['r1hitcount'] == 0:
        crtseqs[seq_id]['q1start']= int(w[6])
        crtseqs[seq_id]['q1end']= int(w[7])
        crtseqs[seq_id]['s1start']= int(w[8])
        crtseqs[seq_id]['s1end']= int(w[9])
        crtseqs[seq_id]['score5']= float(w[2])
        crtseqs[seq_id]['r1hitcount'] += 1
        crtseqs[seq_id]['mutdel1']= "m" + w[4] + "/d" + w[5]
      else:
        crtseqs[seq_id]['r1hitcount'] += 1  # hitcount > 2 is abnormal, will keep first (q1) and last (q2) hit 
        if crtseqs[seq_id]['q1start'] < int(w[6]):    # hit1 is on the left, hit2 on the right of read 1 
          crtseqs[seq_id]['q2start']= int(w[6])
          crtseqs[seq_id]['q2end']= int(w[7])
          crtseqs[seq_id]['s2start']= int(w[8])
          crtseqs[seq_id]['s2end']= int(w[9])
          crtseqs[seq_id]['score3']= float(w[2])
          crtseqs[seq_id]['mutdel2']= "m" + w[4] + "/d" + w[5]
        else:  # exchange 
          crtseqs[seq_id]['q2start']=  crtseqs[seq_id]['q1start']
          crtseqs[seq_id]['q2end']=  crtseqs[seq_id]['q1end']
          crtseqs[seq_id]['s2start']= crtseqs[seq_id]['s1start']
          crtseqs[seq_id]['s2end']= crtseqs[seq_id]['s1end']
          crtseqs[seq_id]['score3']= crtseqs[seq_id]['score5']
          crtseqs[seq_id]['mutdel2']= crtseqs[seq_id]['mutdel1']
          crtseqs[seq_id]['q1start']= int(w[6])
          crtseqs[seq_id]['q1end']= int(w[7])
          crtseqs[seq_id]['s1start']= int(w[8])
          crtseqs[seq_id]['s1end']= int(w[9])
          crtseqs[seq_id]['score5']= float(w[2])
          crtseqs[seq_id]['mutdel1']= "m" + w[4] + "/d" + w[5]
  f1.close()
  if re.search('.gz$', r2file):
    f2=gzip.open(r2file)
  else:
    f2=open(r2file)
  for line in f2:
    line = line.rstrip('\n')
    w = line.split('\t')
    seq_id = w[0]
    gap = int(w[5])
    pident = float(w[2])
    if crtseqs[seq_id]['r1hitcount'] !=2:    # process read 2 only if no single-read junction found in Read 1
      crtseqs[seq_id]['seq_id'] = seq_id + "_2RC"
      #print crtseqs[seq_id]['seq_id'], crtseqs[seq_id]['r1hitcount'], crtseqs[seq_id]['r2hitcount']
      if pident > pident_cutoff:  # 
        len2 = len(crtseqs[seq_id]['seq2'])
        qstart = len2-int(w[7])+1            # convert to coordinates on reverse_complement of query
        qend = len2-int(w[6])+1
        sstart = int(w[9])                   # subject coordinates not changed
        send = int(w[8])
        if crtseqs[seq_id]['r2hitcount'] == 0:
          # store as hit2 first for r2, in case both r1 and r2 has single hits and we will use hit1 of r1 and hit2 of r2 as cross-read junction
          crtseqs[seq_id]['q2start']= qstart  
          crtseqs[seq_id]['q2end']= qend
          crtseqs[seq_id]['s2start']= sstart
          crtseqs[seq_id]['s2end']= send
          crtseqs[seq_id]['r2hitcount'] += 1
          crtseqs[seq_id]['score3']= float(w[2])
          crtseqs[seq_id]['mutdel2']= "m" + w[4] + "/d" + w[5]
        else:    # if the 2nd hit found for read 2, it will override hit 1 of read 1
          crtseqs[seq_id]['r2hitcount'] += 1
          if crtseqs[seq_id]['q2start'] > qstart:    # hit1 is on the left, hit2 on the right of read 1 
            crtseqs[seq_id]['q1start']= qstart
            crtseqs[seq_id]['q1end']= qend
            crtseqs[seq_id]['s1start']= sstart
            crtseqs[seq_id]['s1end']= send
            crtseqs[seq_id]['score5']= float(w[2])
            crtseqs[seq_id]['mutdel1']= "m" + w[4] + "/d" + w[5]
          else:    # exchange
            crtseqs[seq_id]['q1start']=  crtseqs[seq_id]['q2start']
            crtseqs[seq_id]['q1end']=  crtseqs[seq_id]['q2end']
            crtseqs[seq_id]['s1start']= crtseqs[seq_id]['s2start']
            crtseqs[seq_id]['s1end']= crtseqs[seq_id]['s2end']
            crtseqs[seq_id]['score5']=  crtseqs[seq_id]['score3']
            crtseqs[seq_id]['mutdel1']=  crtseqs[seq_id]['mutdel2']
            crtseqs[seq_id]['q2start']= qstart
            crtseqs[seq_id]['q2end']= qend
            crtseqs[seq_id]['s2start']= sstart
            crtseqs[seq_id]['s2end']= send
            crtseqs[seq_id]['score3']= float(w[2])
            crtseqs[seq_id]['mutdel2']= "m" + w[4] + "/d" + w[5]
  f1.close()

def output_blast ():
  #out_handle = open('%s_blast.stat'%(basename), 'w') # for debugging
  #out_handle.write ("#prerRNA\tPolyA\t3'end\t5'end\tType\tsites\tdist_err\tFrac_A\tJunction\tSeqID\n")
  random.seed(100101)
  crosshit_count =0
  for k, v in crtseqs.items():
    #print k, v
    if (v['r1hitcount'] == 2) or (v['r2hitcount'] == 2) :    # single-read junction
      v['type_junc'] = 'NormalJunc'
      sites = site_35S_5end[v['s2start']] + "-" + site_35S_3end[v['s1end']]  # for example 'B1-B2'
      if v['r1hitcount'] == 2:
        inputseq = v['seq1']
      else:
        inputseq = v['seq2']
      #Length of extension/deletion relative to a processing site, coordinates on 35Sext is 1-base from blast
      ext5_len = v['s1end']- (site[site_35S_3end[v['s1end']]]-1)  # 5' half of read or 3' end of pre-rRNA
      ext3_len = site[site_35S_5end[v['s2start']]] - v['s2start']  # 3' half of read or 5' end of pre-rRNA
      # half-open coordiations for polya region: 0-base for polya_start, 1-base for polya_end (like in .bed)
      # only ext5_len, ext3_len will be used in splitting read sequence
      polya_start = v['q1end']
      polya_end = v['q2start'] -1
      dist_err = v['q1end']- v['q2start']+1        	# overlapped nt, = minus original polya length
      #print sites, site[site_35S_5end[v['s2start']]], v['s2start'], ext3_len, "dist_err=", dist_err, "ext3_len, adjust==>"
      # resolve overlapping ends => ext5_len, ext3_len  are always reduced by dist_err
      if dist_err > 0:                          # polya = none,   5' and 3' ends are overlapped by dist_err nt
        # Try to make one intact end (ext_len = 0), ext3_len =0
        if ext3_len <= dist_err and  ext3_len >= 0:   # first make 5' end intact 
          dist_err2 = ext3_len
          dist_err1 = dist_err - dist_err2  # split dist_err into two parts
          polya_start = polya_start - dist_err1
          polya_end = polya_start
          ext5_len = ext5_len - dist_err1 #  reduce extension by part 1 of dist_err
          ext3_len = ext3_len - dist_err2 # reduce extension by part 2 of dist_err
          v['type_junc'] = 'Keep5End'
        elif ext5_len <= dist_err and  ext5_len >= 0:   # make 3' end intact
          dist_err1 = ext5_len
          dist_err2 = dist_err - dist_err1
          polya_start = polya_start - dist_err1
          polya_end = polya_start
          ext5_len = ext5_len - dist_err1
          ext3_len = ext3_len - dist_err2
          v['type_junc'] = 'Keep3End'
        else:                                            # Random split
          dist_err1 = randint(0,dist_err)
          dist_err2 = dist_err - dist_err1
          polya_start = polya_start - dist_err1 
          polya_end = polya_start
          ext5_len = ext5_len - dist_err1
          ext3_len = ext3_len - dist_err2
          v['type_junc'] = 'RandomSplit'

      elif dist_err < 0:   # pA > 0,  resolve conflict between terminal base and pA if pA > 0 and terminal base is A*n. Do nothing if pA=0
        term5 = str(inputseq[polya_start-10:polya_start])  # terminal sequence of 5' half of read or 3' end of prerRNA
        term3 = str(inputseq[polya_end:polya_end+10])  # terminal sequence  of 3' half of read or 5' end of prerRNA
        segment_polya = inputseq[polya_start:polya_end]
        #print ">>>", term5, segment_polya, term3
        m1 = re.search ('(A+)$', term5)
        #print "5 end before adj", ext5_len, polya_start
        if m1:
          term5a_length = len(m1.group(1)) # overlapped As
          if 0 < ext5_len <= term5a_length:
            polya_start = polya_start - ext5_len
            ext5_len = 0                   # make 3' end of pre-rRNA intact if possible
            v['type_junc'] = v['type_junc']+'-K3A' # Keep 3' pA
          elif ext5_len != 0:                           # random split
            splitsize = randint(0,term5a_length)
            ext5_len = ext5_len - splitsize
            polya_start = polya_start - splitsize
            v['type_junc'] = v['type_junc']+'-R3A' + str(splitsize)  # Random 3' pA 
          #print "5 end after adj", term5, m1.group(1), ext5_len, polya_start
        #print "3 end before adj", ext3_len, polya_end
        m2 = re.search ('^(A+)', term3)
        if m2:
          term3a_length = len(m2.group(1))
          if 0 < ext3_len <= term3a_length:
            polya_end = polya_end + ext3_len
            ext3_len = 0                   # make 5' end intact if possible
            v['type_junc'] = v['type_junc'] + '-K5A'  # keep 5' pA
          elif ext3_len != 0:                           # random split
            splitsize = randint(0,term3a_length)
            ext3_len = ext3_len - splitsize
            polya_end = polya_end + splitsize
            v['type_junc'] = v['type_junc'] + '-R5A' + str(splitsize)  # random 5' pA 
          #print "3 end after adj", term3, m2.group(1), ext3_len, polya_end


      # for correct assignment of seed and extension regions
      if ext5_len > 0:
        if polya_start - ext5_len-10 >= 0:           # When inputseq is long enough to write out both 5s (10 nt) and 5e (extension)
          segment_5s = inputseq[polya_start - ext5_len - 10:polya_start - ext5_len]  # seed region, mature rRNA
          segment_5e = inputseq[polya_start - ext5_len:polya_start]  # extension
        elif polya_start - ext5_len >= 0 :           # 5s is truncated, 5e is complete
          segment_5s = inputseq[0:polya_start - ext5_len]  
          segment_5e = inputseq[polya_start - ext5_len:polya_start] 
        else:                                      # 5s is none, 5e is 10 nt
          segment_5s = ''
          segment_5e = inputseq[polya_start - 10:polya_start]    
      else:                                        # no extension for 5s (3 end)
        segment_5s = inputseq[polya_start - 10:polya_start]
        segment_5e = ''

      segment_polya = inputseq[polya_start:polya_end]

      inputseq_len = len(inputseq)
      if ext3_len > 0:
        if polya_end + ext3_len + 10 <= inputseq_len:
          segment_3e = inputseq[polya_end:polya_end + ext3_len]
          segment_3s = inputseq[polya_end + ext3_len:polya_end + ext3_len+10]
        elif polya_end + ext3_len <= inputseq_len:
          segment_3e = inputseq[polya_end:polya_end + ext3_len]
          segment_3s = inputseq[polya_end + ext3_len:inputseq_len]
        else:
          segment_3e = inputseq[polya_end:polya_end + 10]
          segment_3s = ''
      else:
        segment_3e = ''
        segment_3s = inputseq[polya_end:polya_end + 10]

## Further adjustment for special cases
#################################################
# Yeast pre-rRNA specific setting, Start
#################################################
      if segment_5e == 'A': # special treatment for  20S, 23S/22S/21S, move A to pA
        segment_5e = ''
        segment_polya = 'A' + segment_polya
        ext5_len -= 1    
      if segment_5e == "AA":  # special treatment for 18S , move AA to pA
        segment_5e = ""
        segment_polya = "AA" + segment_polya
        ext5_len -= 2
#################################################
# Yeast pre-rRNA specific setting, End
#################################################
#  All adjustments done     
      v['junction'] = segment_5s + segment_5e.lower() + "|" + segment_polya.lower() + "|" + segment_3e.lower() + segment_3s
      v['deletion5'] = ext5_len    # 5' part of read,  3' end of RNA, + extension, - deletion, 0 intact
      v['deletion3'] = ext3_len    # 3' part of read,  5' end of RNA
      v['polya'] = len(segment_polya)
      v['polyaseq'] = segment_polya
      v['frac_a'] = fraction_adenine(segment_polya)
      v['prerRNA'] = prerrna_dict.get(sites, "UNK")
      v['sites'] = sites
      v['dist_err'] = -dist_err   ## change the sign, equal to original pA length before correction
      #output allhits for blast output, for debugging
      #out_handle.write ('%s\t%d\t%d\t%d\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%.1f\t%.1f\n'%(v['prerRNA'],v['polya'],v['deletion5'],v['deletion3'],v['type_junc'],v['sites'],v['dist_err'],v['frac_a'],v['junction'],v['seq_id'], v['mutdel1'], v['mutdel2'], v['score5'], v['score3']))

## very long pA across two reads, or long product with pA =0
    if (v['r1hitcount'] == 1) and (v['r2hitcount'] == 1) :
      if v['s1start'] -  v['s2start'] > 100:  # not the normal RNA, not too close
        v['type_junc'] = 'CrossJunc'
        ext5_len = v['s1end'] - (site[site_35S_3end[v['s1end']]] - 1)
        ext3_len = site[site_35S_5end[v['s2start']]] - v['s2start']
        polya_start = v['q1end']
        polya_end = v['q2start'] -1

        segment_polya1= str(v['seq1'][polya_start:])
        segment_polya2= str(v['seq2'][0:polya_end])

        if len(segment_polya1+segment_polya2) > 0:   # pA > 0,  resolve conflict between terminal base and pA if pA > 0 and terminal base is A*n. Do nothing if pA=0
          term5 = str(v['seq1'][polya_start-10:polya_start])
          term3 = str(v['seq2'][polya_end:polya_end+10])
          #print ">>> Cross-read 1", term5, segment_polya1
          #print ">>> Cross-read 2", segment_polya2, term3
          m1 = re.search ('(A+)$', term5)
          #print "5 end before adj", ext5_len, polya_start
          if m1:
            term5a_length = len(m1.group(1))
            if 0 < ext5_len <= term5a_length:
              polya_start = polya_start - ext5_len
              ext5_len = 0                   # make 3' end intact if possible
              v['type_junc'] = v['type_junc'] + '-K3A'
            elif ext5_len != 0:                           # random split
              splitsize = randint(0,term5a_length)
              ext5_len = ext5_len - splitsize
              polya_start = polya_start - splitsize
              v['type_junc'] = v['type_junc'] + '-R3A' + str(splitsize)
            #print "5 end after adj", term5, m1.group(1), ext5_len, polya_start
          #print "3 end before adj", ext3_len, polya_end
          m2 = re.search ('^(A+)', term3)
          if m2:
            term3a_length = len(m2.group(1))
            if 0 < ext3_len <= term3a_length:
              polya_end = polya_end + ext3_len
              ext3_len = 0                   # make 5' end intact if possible
              v['type_junc'] = v['type_junc'] + '-K5A'
            elif ext3_len != 0:                           # random split
              splitsize = randint(0,term3a_length)
              ext3_len = ext3_len - splitsize
              polya_end = polya_end + splitsize
              v['type_junc'] = v['type_junc'] + '-R5A' + str(splitsize)
            #print "3 end after adj", term3, m2.group(1), ext3_len, polya_end

        # select A track from polya
        segment_polya1= extract_polya_bwa(v['seq1'][polya_start:], is_3end=False)
        segment_polya2= extract_polya_bwa(v['seq2'][0:polya_end], is_3end=True)
        #print segment_polya1
        #print segment_polya2
        diff_len = len(segment_polya1) - len(segment_polya2)
        if diff_len >= 0:  # chose the longer one as polya
          segment_polya = segment_polya1
        else:
          segment_polya = segment_polya2
        # for correct assignment of seed and extension regions
        if ext5_len > 0:
          if polya_start - ext5_len-10 >= 0:
            segment_5s = v['seq1'][polya_start - ext5_len-10:polya_start - ext5_len]  
            segment_5e = v['seq1'][polya_start - ext5_len:polya_start]
          elif polya_start-ext5_len >= 0 :  # ext5_len is very long
            segment_5s = v['seq1'][0:polya_start - ext5_len]  
            segment_5e = v['seq1'][polya_start - ext5_len:polya_start] 
          else:
            segment_5s = ''
            segment_5e = v['seq1'][polya_start - 10:polya_start]    
        else:
          segment_5s = v['seq1'][polya_start - 10:polya_start]
          segment_5e = ''
  
        seq2_len = len(v['seq2'])
        if ext3_len > 0:
          if polya_end + ext3_len + 10 <= seq2_len:
            segment_3e = v['seq2'][polya_end:polya_end + ext3_len]
            segment_3s = v['seq2'][polya_end + ext3_len:polya_end + ext3_len + 10]
          elif polya_end + ext3_len <= seq2_len:
            segment_3e = v['seq2'][polya_end:polya_end + ext3_len]
            segment_3s = v['seq2'][polya_end + ext3_len:seq2_len]
          else:
            segment_3e = v['seq2'][polya_end:polya_end + 10]
            segment_3s = ''
        else:
          segment_3e = ''
          segment_3s = v['seq2'][polya_end:polya_end + 10]

#################################################
# Yeast pre-rRNA specific setting, Start
#################################################
        ## Further adjustment for special cases, even for pA = 0, The pA > 0 case shall already be solved.
        if segment_5e == 'A': # special treatment for 20S, 23S/22S/21S, move A to pA
          segment_5e = ''
          if diff_len > 0:  # Don't adjust segment_polya2, a small bug: after adjustment, diff_len and choice of pA may change.
            segment_polya = 'A' + segment_polya
          ext5_len -= 1    
        if segment_5e == "AA":  # special treatment for 18S , move AA to pA
          segment_5e = ""
          if diff_len > 0:
            segment_polya = "AA" + segment_polya
          ext5_len -= 2
#################################################
# Yeast pre-rRNA specific setting, End
#################################################

        crosshit_count +=1
        v['junction'] = segment_5s + segment_5e.lower() + "|" + segment_polya.lower() + "|" + segment_3e.lower() + segment_3s
        v['deletion5'] = ext5_len
        v['deletion3'] = ext3_len
        v['polya'] = len(segment_polya)
        v['polyaseq'] = segment_polya
        v['frac_a'] = fraction_adenine(segment_polya)
        sites = site_35S_5end[v['s2start']] + "-" + site_35S_3end[v['s1end']] 
        v['prerRNA'] = prerrna_dict.get(sites, "UNK")
        v['sites'] = sites
        v['dist_err'] = diff_len   # Different meaning from diff_err
        #output allhits for blast output, for debugging
        #out_handle.write ('%s\t%d\t%d\t%d\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%.1f\t%.1f\n'%(v['prerRNA'],v['polya'],v['deletion5'],v['deletion3'],v['type_junc'],v['sites'],v['dist_err'],v['frac_a'],v['junction'],v['seq_id'], v['mutdel1'], v['mutdel2'], v['score5'], v['score3']))
        #print ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s' % (v['prerRNA'],  v['polya'], v['deletion5'], v['deletion3'], v['type_junc'], v['sites'],  v['dist_err'], v['frac_a'], v['junction'], v['seq_id']))
  print ("Cross-read hits before filtering:", crosshit_count)
  #out_handle.close()

def extract_polya_bwa (seq, is_3end=True):
  """ mutation is present in polyA track """
  """ BWA-like approach to select A-track residues from 3 or 5 end """
  """ 5-AAAAXAAAAAAAANNNNN- """
  """   score = 1 for A, -3 for other base  """
  """   1234123456789630  """
  if seq == None or seq == "" :
    return ""
  if is_3end:    # reverse
    seq = seq.upper()[::-1]  
  else:
    seq = seq.upper()
  seq_len = len(seq)
  num_list = np.zeros(seq_len)
  sum_list = np.zeros(seq_len)
  total = 0
  for i in range(seq_len):
    if seq[i] == 'A':
      num_list[i] = -1
    else:
      num_list[i] = 3
  for j in range(0,seq_len):
    total += num_list[j]
    sum_list[j] = total
    peak_ind = np.where(sum_list == min(sum_list))[0]
  #print num_list
  #print sum_list
  #print peak_ind[-1], seq[0:peak_ind[-1]+1], seq[peak_ind[0]+1:]
  #print seq[0:peak_ind[-1]+1][::-1]
  if is_3end:
    return seq[0:peak_ind[-1]+1][::-1]
  else:
    return seq[0:peak_ind[-1]+1]


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Processing CircTA-seq data")
  parser.add_argument( '-a', '--action', help="", required=True, choices=['blasttab', 'del3', 'unk', 'test'])
  parser.add_argument( '-r1', '--r1', help="read 1, fasta or fastq")
  parser.add_argument( '-r2', '--r2', help="read 2, fasta or fastq")
  parser.add_argument( '-mapfile', '--mapfile', help="basename_all.stat output")
  parser.add_argument( '-basename', '--basename', help="basename of files")
  parser.add_argument( '-tab1', '--table_blast1', help="blast table fmt 6")
  parser.add_argument( '-tab2', '--table_blast2', help="blast table fmt 6")
  parser.add_argument( '-genome35s', '--genome35s', help="fasta file for genome 35Sext used for blastn alignment")
  args = parser.parse_args()
  basename = args.basename
  genome35s = args.genome35s

  if args.action == "blasttab": 
    RNA35S = read_35S(genome35s)
    read_reads(args.r1, args.r2)
    read_blasttable(args.table_blast1, args.table_blast2)
    set_site35S_5end ()
    set_site35S_3end ()
    output_blast()
    output_stat()
  if args.action == "del3": 
    RNA35S = read_35S(genome35s)
    read_map(args.mapfile)
    output_del3_new()
    #polya_composition()
  if args.action == "unk": 
    RNA35S = read_35S(genome35s)
    set_site35S_5end ()
    set_site35S_3end ()
    read_map(args.mapfile)
    output_unk()
  if args.action == "test":
    pass

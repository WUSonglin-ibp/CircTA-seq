#Shell script for processing CircTA-seq data. Comment out the needed steps

date
# x is a directory to be processed, sequencing fastq file should be named as  x_R1.fq.gz, x_R2.fq.gz
#for x in X11 X12 X13 X14
for x in X19

do
  cd $x
  pwd
  # cutadapt 
  echo "remove adaptor.."
  cutadapt -a TGGAATTCTCGGGTGCC -A AGATCGGAAGAGCGTCG  -e 0.1 -m 25 --overlap 5 -q 20 -o ${x}_trimq_1.fq -p ${x}_trimq_2.fq ${x}_R1.fq.gz ${x}_R2.fq.gz > cutadapt_trimq.out

  # convert fastq to fasta, seqtk 1.3-r106
  echo "to fasta.."
  seqtk seq -A ${x}_trimq_1.fq > ${x}_trimq_1.fasta
  seqtk seq -A ${x}_trimq_2.fq > ${x}_trimq_2.fasta

  # Blast alignment, require fasta
  echo "blastn read 1/2..."
  blastn -db ../genome/35Sext -task "blastn-short" -evalue 1e-3 -strand plus -outfmt 6 -query ${x}_trimq_1.fasta -out ${x}_trimq_r1_blast.tab
  blastn -db ../genome/35Sext -task "blastn-short" -evalue 1e-3 -strand minus -outfmt 6 -query ${x}_trimq_2.fasta -out ${x}_trimq_r2_blast.tab

  #################
  # Major analysis step
  #################
  echo "analyze circular junction"
  python ../analysis_circtaseq.py -a blasttab -genome35s '../genome/35Sext' -tab1  ${x}_trimq_r1_blast.tab -tab2 ${x}_trimq_r2_blast.tab -r1 ${x}_trimq_1.fq -r2 ${x}_trimq_2.fq -basename b6 

  #sort -k1,1 -k3,3n  -k4,4n -k2,2n b6_all.stat > b6_all.stat.sort

  #################
  # analyze poly-A distribution at individual positions of 3' tail
  #################
  echo "analyze poly-A at individual positions"
  python ../analysis_circtaseq.py -a del3 -genome35s '../genome/35Sext' -mapfile b6_all.stat -basename b6
 
  #################
  # output coordinate file for drawing 2D map, using 100k reads
  #################
  echo "output coordinate file for drawing 2D map"
  head -100000 b6_all.stat > b6_100k_all.stat
  python ../analysis_circtaseq.py -a unk -genome35s '../genome/35Sext' -mapfile b6_100k_all.stat -basename b6
  sort -k3,3n b6_UNK_coord.dat > b6_UNK_coord.dat.sort

  cd ../
done
echo "Everything done"
date

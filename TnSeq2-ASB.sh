#!/usr/bin/env bash
#This requires fqgrep (https://github.com/indraniel/fqgrep)
usage () {
  echo "usage: $0 [-p <primer seq>] [-i <IR seq>] [-a <assembly>] <pfx> "
  echo "Required parameters:"
# echo "-p     The sequence of your Tn-seq primer specific to your transposon"
# echo "-i     The sequence of the transposon end sequence remaining (for junction authentication)"
  echo "-a     The name of the assembly you're using (PAO1)"
  echo "-m     The number of mismatches/indels you want to tolerate during search"

  echo ""
  echo "The required parameters must precede the file prefix for your sequence file:"
  echo "  (e.g. if your sequence file is named condition1_R1.fastq,"
  echo "   the prefix is \"condition1\")"
  echo ""
  echo "Example:"
# echo "$0 -p CGTCCAGGACGCTACTTGTG -i TATAAGAGTCAG -a PAO1 -m 1 condition1"
  echo "$0 -a PAO1 -m 1 condition1"
}

# Read in the important options
# while getopts ":p:i:a:m:" option; do
while getopts ":a:m:" option; do
  case "$option" in
# p)  PRIMER="$OPTARG" ;;
# i)  IR="$OPTARG" ;;
  	a)  ASSEMBLY="$OPTARG" ;;
  	m)  MISMATCHES="$OPTARG" ;;
    h)  # it's always useful to provide some help 
        usage
        exit 0 
        ;;
    :)  echo "Error: -$option requires an argument" 
        usage
        exit 1
        ;;
    ?)  echo "Error: unknown option -$option" 
        usage
        exit 1
        ;;
  esac
done    
shift $(( OPTIND - 1 ))

# Do some error checking to make sure parameters are defined
# if [ -z "$PRIMER" ]; then
  # echo "Error: you must specify the primer sequence using -p"
  # usage
  # exit 1
# fi

# if [ -z "$IR" ]; then
  # echo "Error: you must specify the Tn end sequence using -i"
  # usage
  # exit 1
# fi

if [ -z "$ASSEMBLY" ]; then
  echo "Error: you must specify an assembly using -a"
  usage
  exit 1
fi

if [ -z "$MISMATCHES" ]; then
  echo "Error: you must specify a number of mismatches using -m"
  usage
  exit 1
fi

# Give the usage if there aren't enough parameters
if [ $# -lt 1 ] ; then
  echo "you must provide a file prefix for analysis"
  usage
  exit 1
fi

PREFIX=$1
R1=${PREFIX}_R1
BOWTIEREF=$REFGENOME/$ASSEMBLY/$ASSEMBLY
# BOWTIEREF=/Users/papinlab/ref_genome/WKB2/WKB2

echo "Performing TnSeq analysis on $PREFIX..."
echo "TnSeq processing stats for $PREFIX" > $PREFIX-TnSeq.txt
echo "Total sequences: " >> $PREFIX-TnSeq.txt
egrep -c '^@HWI|^@M|^@NS|^@SRR|^@' $R1.fastq >> $PREFIX-TnSeq.txt

# Reads with primer
# echo "$PREFIX: Searching for reads with primer..."
# FQ=$(fqgrep -m $MISMATCHES -C -p $PRIMER $R1.fastq)
# [[ "$FQ" =~ \ ([0-9]+)\  ]]
# PRIMERCOUNT=${BASH_REMATCH[1]}
# echo "Reads with primer:" >> $PREFIX-TnSeq.txt
# echo $PRIMERCOUNT >> $PREFIX-TnSeq.txt

# IRs
# echo "$PREFIX: Searching for reads with an IR in right location..."
# let "MIN = ${#PRIMER} + ${#IR} + 2"
# let "MAX = ${#PRIMER} + ${#IR} + 8"
# fqgrep -m $MISMATCHES -r -p $PRIMER$IR $R1.fastq | awk -v min=$MIN -v max=$MAX -F "\t" '(($8 >= min && $8 <= max) || $1=="read name")' | trimmer --5-prime > $PREFIX-IR-clip.fastq

#I use for "flexbar_executable" for personal reasons- you probably want to switch t whatever your flexbar command is... probably "flexbar"
#
# flexbar_executable -f fastq-i1.8 -n 16 -ao 8 -m 18 -z 35 -ae RIGHT -a ~/adapters/3_adapter_seq.fasta -r $PREFIX-IR-clip.fastq -t $PREFIX-IR-clip.trim >> /dev/null 2>&1
# mv $PREFIX-IR-clip.trim.fastq $PREFIX-IR-clip.fastq
# IRSFOUND=$(egrep -c '^@HWI|^@M|^@' $PREFIX-IR-clip.fastq)
# echo "Molecules with IR in right location:" >> $PREFIX-TnSeq.txt
# echo $IRSFOUND >> $PREFIX-TnSeq.txt

# Map and convert - feel free to change bowtie2 parameters yourself
echo "$PREFIX: Mapping with Bowtie2..."
echo "Bowtie2 report:" >> $PREFIX-TnSeq.txt
# bowtie2 --end-to-end -p 16 -a -x $BOWTIEREF -U $PREFIX-IR-clip.fastq -S $PREFIX.sam 2>> $PREFIX-TnSeq.txt
bowtie2 --end-to-end -p 16 -a -x $BOWTIEREF -U $PRIMER $R1.fastq -S $PREFIX.sam 2>> $PREFIX-TnSeq.txt
grep '^@' $PREFIX.sam > $PREFIX-mapped.sam
# ************ changed awk to gawk
cat $PREFIX.sam | grep -v '^@' | gawk -F "\t" '(and($2, 0x4) != 0x4)' | sort -u -k1,1 >> $PREFIX-mapped.sam
echo "Number of reads mapping at high enough score:" >> $PREFIX-TnSeq.txt
cat $PREFIX-mapped.sam | wc -l >> $PREFIX-TnSeq.txt
echo "$PREFIX: Tallying mapping results..."
# ************ changed awk to gawk
grep -v '^@' $PREFIX-mapped.sam | gawk -F "\t" 'and($2, 0x100) != 0x100 {if (and($2, 0x10) != 0x10) print $4; else print $4+length($10)}' | grep '[0-9]' | sort | uniq -c | sort -n -r > $PREFIX-sites.txt
echo "Number of insertion sites identified:" >> $PREFIX-TnSeq.txt
wc -l $PREFIX-sites.txt >> $PREFIX-TnSeq.txt
echo "Most frequent sites:" >> $PREFIX-TnSeq.txt
head $PREFIX-sites.txt >> $PREFIX-TnSeq.txt

# Generate IGV files
echo "$PREFIX: Generating files for IGV..."
samtools view -b -S $PREFIX-mapped.sam > $PREFIX-mapped.bam 2> /dev/null
# ********
# samtools sort $PREFIX-mapped.bam $PREFIX-mapped.sorted > /dev/null
samtools sort $PREFIX-mapped.bam -T temp -o $PREFIX-mapped.sorted.bam > /dev/null
samtools index $PREFIX-mapped.sorted.bam > /dev/null

# Sort output, cleanup
echo "$PREFIX: Cleaning up..."
mkdir $PREFIX 2> /dev/null
mv $PREFIX-TnSeq.txt $PREFIX/
# mv $PREFIX-IR-clip.fastq $PREFIX/
mv $PREFIX.sam $PREFIX/
mv $PREFIX-mapped.sam $PREFIX/
mv $PREFIX-mapped.bam $PREFIX/
mv $PREFIX-sites.txt $PREFIX/
mkdir $PREFIX/IGV 2> /dev/null
mv $PREFIX-mapped.sorted.bam $PREFIX/IGV
mv $PREFIX-mapped.sorted.bam.bai $PREFIX/IGV

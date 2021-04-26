#!/bin/bash
######################################################################
# required bins:
#   bowtie2     2.3.4.1
#   samtools    1.3
#   bamCoverage 3.0.0
######################################################################
if [[ $# -ne 7 ]]; then
  echo "    "
  echo "  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.- "
  echo " / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / /   "
  echo "'-'   '-'-'   '-'-'   '-'-'   '-'-'   '-'-'   '-'-'   '-'-'    "
  echo "               Orthrus by Robin H. van der Weide               "
  echo "                  an updated ChIP-Rx approach                  "
  echo "                    reference-adjusted RPKM                    "
  echo "  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.- "
  echo " / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / /   "
  echo "'-'   '-'-'   '-'-'   '-'-'   '-'-'   '-'-'   '-'-'   '-'-'    "
  echo "    "
  echo "typical usage:"
  echo "    								"
  echo "  orthrus IP.fq.gz /DATA/references/spikeIn/hg19_mm10 MM10 HG19 30 200 'narrow/broad'"
  echo "    								"
  echo "  IP.fastq.gz     A fastq of ChIP-seq (please no .fq.gz)        "
  echo "  hg19_mm10       The bowtie-index path of the concatenated     "
  echo "                    reference genomes: every chromosome should  "
  echo "                    be prefixed with the appropriate basname    "
  echo "                    (e.g. the chr1 of hg19 will be HG19_chr1).  "
  echo "  HG19            The basename of the subject-genome	        "
  echo "  MM10            The basename of the spiked-genome 	        "
  echo "  30              The number of threads used		        "
  echo "  200             Extend reads by []bp                          "
  echo " 					                        "
  exit
fi

echo "[orthrus]    mapping..." >&2
######################################################################
# bowtie2 mapping
b=$(basename $1 .fastq.gz | sed 's/$/.bam/')
bowtie2 -x $2 -U $1 -p $5 --quiet |\
	samtools view -Sb - |\
        samtools sort -T $(mktemp) - |\
        samtools rmdup -s - $b

######################################################################
echo "[orthrus]    indexing..." >&2
samtools index $b

MM_COUNT=`samtools idxstats $b | cut -f1,3 | grep $4 | awk '{SUM += $2} END{print SUM}'`
SF=`bc <<<"scale=10; 1000000/$MM_COUNT"`

echo "[orthrus]    scale factor= " $SF >&2
# remove MM-chromosomes and rename human chromosomes
bs=$(basename $b | sed 's/.bam/_sampleOnly.bam/')
br=$(basename $b | sed 's/.bam/_referenceOnly.bam/')
samtools view -h $b | grep -v $4 | sed "s/$3//g" | sed "s/_chr/chr/g" | samtools view -Shbq 15 > $bs
samtools view -h $b | grep -v $3 | sed "s/$4//g" | sed "s/_chr/chr/g" | samtools view -Shbq 15 > $br
samtools index $bs
samtools index $br

echo "[orthrus]    generating bigwigs..." >&2
# make bigwig with bamcoverage
if [ $3 == "MM10" ]; then SAGEN=2188620071; else SAGEN=2451960000; fi
if [ $4 == "HG19" ]; then REFGEN=2451960000; else REFGEN=2188620071; fi
echo $SAGEN
echo $REFGEN
READABLESF=$(echo $SF | cut -c 1,2,3,4,5)
bsname=$(basename $1 .fastq.gz | sed 's/$/_SF/' | sed "s/$/$READABLESF/" | sed 's/$/_MQ15_sample.bw/')
brsample=$(basename $1 .fastq.gz | sed 's/$/_SF/' | sed "s/$/$READABLESF/" | sed 's/$/_MQ15_reference.bw/')
bamCoverage -b $bs -o $bsname -p $5 --scaleFactor $SF --ignoreForNormalization chrM --binSize 10 --minMappingQuality 15 --normalizeUsing RPGC --effectiveGenomeSize $SAGEN --extendReads $6
bamCoverage -b $br -o $brsample -p $5 --scaleFactor $SF --ignoreForNormalization chrM --binSize 10 --minMappingQuality 15 --normalizeUsing RPGC --effectiveGenomeSize $REFGEN --extendReads $6

echo "[orthrus]    calling peaks..." >&2
# peak input settings
if [ $3 == "MM10" ]; then REFGEN="mm"; else REFGEN="hs"; fi
if [ $4 == "HG19" ]; then REFREF="hs"; else REFREF="mm"; fi
bsp=$(basename $b | sed 's/.bam/_sampleOnly/')
brp=$(basename $b | sed 's/.bam/_referenceOnly/')
if [ $7 == "narrow" ]
then
macs2 callpeak -t $bs -f BAM -g $REFGEN -n $bsp -q 0.01
macs2 callpeak -t $br -f BAM -g $REFREF -n $brp -q 0.01
elif [ $7 == "broad" ]
then
macs2 callpeak -t $bs -f BAM -g $REFGEN -n $bsp -q 0.01 --broad
macs2 callpeak -t $br -f BAM -g $REFREF -n $brp -q 0.01 --broad
else echo "Not for now!"
fi

# remove tmp-bam
rm $b
rm $(basename $1 .fastq.gz | sed 's/$/.bam.bai/')
rm *.r *.xls *_summits.bed

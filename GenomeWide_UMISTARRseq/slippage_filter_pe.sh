#!/bin/bash
set -o errexit
set -o pipefail

# Stark Lab in house pipeline. Adapted by Bernardo Almeida 2020

################################################################################
# Set default values
################################################################################

MM="2"
L="10"
OFFSET="1"
D="2000"
memory=4G #memory restrictions for shell sort

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
  echo >&2 "
$(basename $0) - Filter putative slippage fragments in paired-end data

USAGE: cat BED file | $(basename $0) - [OPTIONS] --
 -m     Maximum number of mismatches for sequence comparison [default: $MM]
 -l     Sequence length for comparison [default: $L]
 -o     Offset for sequence to compare [default: $OFFSET]
 -d     Maximum distance for filtering within a cluster [default: $D]
 -m    Memory restrictions for shell sort  [default: $memory]

"
  exit 1
fi

################################################################################
# Parse input
################################################################################

while getopts "m:l:o:d:" o
do
  case "$o" in
    m) MM="$OPTARG";;
    l) L="$OPTARG";;
    o) OFFSET="$OPTARG";;
    d) D="$OPTARG";;
    m)  memory="$OPTARG";;
   \?) exit 1;;
  esac
done

################################################################################
# Run program
###############################################################################

# Sort by strand, chr, start (lower on top), end (higher on top), and mismatches
# The longest fragment always appears first

sort -k6,6 -k1,1 -k2,2n -k3,3nr -k5,5n -S $memory| \
  awk -vD=$D -vL=$L -vMM=$MM -vO=$OFFSET '
  function dist(a,b,   i,sa,sb,m)
    {
      m=0
      for(i=1;i<=L;i++){
        sa=substr(a,i,1);sb=substr(b,i,1)
        if(sa!=sb && sa!="N" && sb!="N"){m++}
      }
      return m
    }
  {split($4,seq,"_"); ok=1}
  (NR>1 && $6==strand && $1==chr) {
    if($2==start && $3>=end-D && $3<=end && strand=="+"){
      SEQS[seq_r]=1; for(SEQ in SEQS){DIFF=dist(substr(seq[2],O,L),SEQ);if(DIFF<=MM){ok=0; break}}
    }
    else if($2==start && $3>=end-D && $3<=end && strand=="-"){
      SEQS[seq_l]=1; for(SEQ in SEQS){DIFF=dist(substr(seq[1],O,L),SEQ);if(DIFF<=MM){ok=0; break}}
    }
    else if($2<=start+D && $2>start && $3==end && strand=="+"){
      SEQS[seq_l]=1; for(SEQ in SEQS){DIFF=dist(substr(seq[1],O,L),SEQ);if(DIFF<=MM){ok=0; break}}
    }
    else if($2<=start+D && $2>start && $3==end && strand=="-"){
      SEQS[seq_r]=1; for(SEQ in SEQS){DIFF=dist(substr(seq[2],O,L),SEQ);if(DIFF<=MM){ok=0; break}}
    }
    else{
      for (SEQ in SEQS){delete SEQS[SEQ]}
    }
  }
  (NR>1 && ($6!=strand || $1!=chr)){for (SEQ in SEQS){delete SEQS[SEQ]}}
  {chr=$1;start=$2;end=$3;strand=$6;seq_l=substr(seq[1],O,L);seq_r=substr(seq[2],O,L)}
  {if(ok){print $0}}' | \
  sort -k1,1 -k2,2n -k3,3nr -k6,6 -k5,5n -k4,4 -S $memory| \
  awk '!x[$1" "$2" "$3" "$6]++'

exit 0

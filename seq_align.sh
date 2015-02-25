#!/bin/bash


# A script to extract sequence files from macrogen and align them to a reference

FILEDIR=$1

echo "Looking for files in " $FILEDIR
cd $FILEDIR

shift

echo "Removing $# files:"
echo $@

if [[ "$#" > 0 ]]; then
    IGNORE=$(echo $@ | awk -v FILEDIR=$FILEDIR -v q="\"" '{for(f=1;f <= NF; f++){$f=q""$f"*"q" ! -name"}}{sub(/ ! -name$/,"",$0);print "find "FILEDIR" -name \"*.txt\" ! -name "$0" -print "}')
else
    IGNORE="find $FILEDIR -name '*.txt' -print "
fi

#echo $IGNORE
#eval $IGNORE

cat `eval $IGNORE` > ${FILEDIR}/chosen_files.fa

awk '$1~/>/{gsub(/>/,"\n>",$1)}{print $0}' ${FILEDIR}/chosen_files.fa > ${FILEDIR}/chosen_files.neat.fa

sed '/^$/d' ${FILEDIR}/chosen_files.neat.fa > ${FILEDIR}/chosen_files.fa


muscle -quiet       < ${FILEDIR}/chosen_files.fa > ${FILEDIR}/chosen_files.muscle
muscle -quiet -html < ${FILEDIR}/chosen_files.fa > ${FILEDIR}/chosen_files.aligned.html

rm ${FILEDIR}/chosen_files.neat.fa

open ${FILEDIR}/chosen_files.aligned.html



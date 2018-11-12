#!/bin/bash

#$ -S /bin/bash
#$ -N predict
#$ -P longrun
#$ -pe smp 8
#$ -cwd

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date)"

module load samtools

### define run directory
rundir=$PWD

### define reference
reference="$root/references/$rname/sequence/$rname.fasta"

### SVM install dir
libsvm="$root/libsvm-3.22"

### filter & convert patches
$root/bin/features.pl \
    --class 0 \
    --features "$features" \
    --mtotal "$mTotal" \
    --logfile candidates.log \
    "$patches" > candidates

$libsvm/svm-scale -r $model/range candidates > candidates.scale
$libsvm/svm-predict candidates.scale $model/model candidates.predict

### merge predictions
echo "SVM_Prediction" > candidates.tmp
cat candidates.tmp candidates.predict > candidates.tmp.2
paste -d, candidates.log candidates.tmp.2 > candidates.tmp.3
(head --lines 1 candidates.tmp.3 && (cat candidates.tmp.3 | egrep ",1$")) > patches.csv

### cleanup
rm -f candidates*

### generate summary
$root/bin/patches-summary.pl ../coverage.csv patches.csv > patches-summary.csv

echo ""
echo "Finished on $(date)"

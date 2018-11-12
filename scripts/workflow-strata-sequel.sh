#!/bin/bash

#$ -S /bin/bash
#$ -N strata
#$ -P longrun
#$ -cwd

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date)"

module load smrtlink

### define run directory
rundir=$PWD

### define reference
reference="$root/references/$rname/sequence/$rname.fasta"

### protocol directory
protocoldir="0-$protocol"
mkdir -p "$rundir/$protocoldir/data"
cd "$rundir/$protocoldir/data"

### run smrtpipe
echo ""
cat << EOF > aligned_reads.sh
#\$ -S /bin/bash
#\$ -N aligned_reads
#\$ -P longrun
#\$ -pe smp 8
#\$ -o aligned_reads.log
#\$ -j yes
#\$ -cwd

module load smrtlink

pbalign "$collectionPathUri" "$reference" aligned_reads.bam \
    --algorithmOptions "--minMatch 12 --bestn 10 --minPctSimilarity 70.0 --refineConcordantAlignments" \
    --nproc 7 \
    --minAccuracy 70.0 \
    --minLength 50 \
    --hitPolicy randombest \
    --concordant \
    --tmpDir "$PWD"
EOF
qsub -sync yes aligned_reads.sh
echo "Task 1 completed at $(date)"

cd "$rundir"

### stratification
stratadir="1-strata"

echo ""
"$root/bin/strata.pl" --prefix "$rundir/aligned_reads" "$rundir/$protocoldir/data/aligned_reads.bam"
echo "Task 2 completed at $(date)"

### single-molecule modifications
echo ""
cat << EOF > modifications.sh
#\$ -S /bin/bash
#\$ -N mods
#\$ -P longrun
#\$ -pe smp 8
#\$ -p -1023
#\$ -cwd

module load smrtlink

### run directory
jobdir=\`printf "$rundir/$stratadir/%05i" \$SGE_TASK_ID\`
mkdir -p "\$jobdir"
cd "\$jobdir"

echo ""
cat "$rundir/aligned_reads.strata.csv" | egrep ",\$SGE_TASK_ID\$" | cut -d, -f1,2 --output-delimiter='/' | sort -t, -k2 -g > "\$jobdir/whitelist.txt"
echo "Task 1 completed at \$(date)"

echo ""
"$root/bin/h5tool.pl" \\
    "$rundir/$protocoldir/data/aligned_reads.bam" \\
    "\$jobdir/whitelist.txt" > aligned_reads.sam

samtools view -Sb aligned_reads.sam > aligned_reads.tmp.bam
samtools sort --threads 7 -m 4G -o aligned_reads.bam aligned_reads.tmp.bam
samtools index aligned_reads.bam aligned_reads.bam.bai
pbindex aligned_reads.bam
echo "Task 2 completed at \$(date)"

echo ""
ipdSummary \\
    -v \\
    --methylFraction \\
    --identify m6A,m4C \\
    --numWorkers 7 \\
    --gff modifications.gff \\
    --csv modifications.csv \\
    --reference "$reference" \\
    aligned_reads.bam
echo "Task 3 completed at \$(date)"

echo ""
gzip modifications.gff
gzip modifications.csv
rm -f aligned_reads.cmp.h5
echo "Task 4 completed at \$(date)"
EOF
echo "Task 3 completed at $(date)"

echo ""
mkdir -p "$rundir/$stratadir/logs"
last=`cat "$rundir/aligned_reads.strata.csv" | cut -d, -f5 | sort -g | tail --lines 1`
qsub -sync yes -t 1-$last -o "$rundir/$stratadir/logs/modifications.log.\$TASK_ID" -j yes modifications.sh
echo "Task 4 completed at $(date)"

### collect coverage
echo ""
cd "$rundir"
cat << EOF > coverage.sh
#\$ -S /bin/bash
#\$ -N coverage
#\$ -P longrun
#\$ -o coverage.log
#\$ -j yes
#\$ -cwd
"$root/bin/coverage.pl" "$rundir"/$stratadir/*/modifications.csv.gz > coverage.csv
EOF
qsub -sync yes coverage.sh
echo "Task 5 completed at $(date)"

echo ""
echo "Finished on $(date)"

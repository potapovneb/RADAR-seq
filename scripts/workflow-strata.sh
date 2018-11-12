#!/bin/bash

#$ -S /bin/bash
#$ -N strata
#$ -P longrun
#$ -cwd

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date)"

module load samtools

### define run directory
rundir=$PWD

### SMRT Analysis
smrtwrap="/opt/pacbio/smrtanalysis/current/smrtcmds/bin/smrtwrap"
smrtpipe="/opt/pacbio/smrtanalysis/current/smrtcmds/bin/smrtpipe"

### protocol directory
protocoldir="0-$protocol"
mkdir -p "$rundir/$protocoldir"
cd "$rundir/$protocoldir"

### input.xml
echo ""
find "$collectionPathUri/Analysis_Results" -type f -name "*.bax.h5" | sort -u > input.fofn
"$smrtwrap" fofnToSmrtpipeInput.py input.fofn > input.xml
echo "Task 1 completed at $(date)"

### settings.xml
echo ""
cat "$root/xml/$protocol" | sed -e "s|REFERENCE|common/references/$rname|" > settings.xml
echo "Task 2 completed at $(date)"

### run smrtpipe
echo ""
$smrtpipe \
    --output="$rundir/$protocoldir" \
    -D TMP="$rundir/$protocoldir/temp" \
    --distribute \
    --params=settings.xml xml:input.xml
rm -rf "$rundir/$protocoldir/temp"
echo "Task 3 completed at $(date)"

cd "$rundir"

### stratification
stratadir="1-strata"

echo ""
"$root"/bin/strata.pl --prefix "$rundir"/aligned_reads "$rundir/$protocoldir/data/aligned_reads.bam"
echo "Task 4 completed at $(date)"

### single-molecule modifications
echo ""
cat << EOF > modifications.sh
#\$ -S /bin/bash
#\$ -N mods
#\$ -P longrun
#\$ -pe smp 8
#\$ -p -1023
#\$ -cwd

### run directory
jobdir=\`printf "$rundir/$stratadir/%05i" \$SGE_TASK_ID\`
mkdir -p "\$jobdir"
cd "\$jobdir"

echo ""
cat "$rundir/aligned_reads.strata.csv" | egrep ",\$SGE_TASK_ID\$" | cut -d, -f1,2 --output-delimiter='/' | sort -t, -k2 -g > "\$jobdir/whitelist.txt"
echo "Task 1 completed at \$(date)"

echo ""
$smrtwrap python $root/bin/h5tool.py \\
    --outFile aligned_reads.cmp.h5 \\
    "$rundir/$protocoldir/data/aligned_reads.cmp.h5" \\
    "\$jobdir/whitelist.txt"
echo "Task 2 completed at \$(date)"

echo ""
$smrtwrap ipdSummary.py \\
    -v \\
    --methylFraction \\
    --identify m6A,m4C \\
    --paramsPath /opt/pacbio/smrtanalysis/current/analysis/etc/algorithm_parameters/2015-11/kineticsTools \\
    --numWorkers 7 \\
    --summary_h5 temp_kinetics.h5 \\
    --gff modifications.gff \\
    --csv modifications.csv \\
    --reference /opt/pacbio/smrtanalysis/current/common/references/$rname/sequence/$rname.fasta \\
    aligned_reads.cmp.h5
echo "Task 3 completed at \$(date)"

echo ""
gzip modifications.gff
gzip modifications.csv
rm aligned_reads.cmp.h5
echo "Task 4 completed at \$(date)"
EOF
echo "Task 5 completed at $(date)"

echo ""
mkdir -p "$rundir/$stratadir/logs"
last=`cat "$rundir/aligned_reads.strata.csv" | cut -d, -f5 | sort -g | tail --lines 1`
qsub \
    -sync yes \
    -t 1-$last \
    -o "$rundir"/$stratadir/logs/modifications.log.\$TASK_ID \
    -j yes \
    modifications.sh
echo "Task 6 completed at $(date)"

### collect coverage
echo ""
cat << EOF > coverage.sh
#\$ -S /bin/bash
#\$ -N coverage
#\$ -P longrun
#\$ -o coverage.log
#\$ -j yes
#\$ -cwd
"$root"/bin/coverage.pl "$rundir"/$stratadir/*/modifications.csv.gz > coverage.csv
EOF
qsub -sync yes coverage.sh
echo "Task 7 completed at $(date)"

echo ""
echo "Finished on $(date)"

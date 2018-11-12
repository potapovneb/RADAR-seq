#!/bin/bash

#$ -S /bin/bash
#$ -N detect
#$ -P longrun
#$ -pe smp 8
#$ -cwd

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date)"

### define run directory
rundir=$PWD

### define reference
reference="$root/references/$rname/sequence/$rname.fasta"

### define site files
methylation_sites="methylation-sites.csv"

### methylation sites
case "$rname" in
    TKO)
	$root/bin/find-methylation-sites-tko.pl "$reference" > $methylation_sites
	;;
    NC_000913|ecoli_hgap_26661)
	$root/bin/find-methylation-sites-ecoli.pl "$reference" > $methylation_sites
	;;
    *)
	echo "Reference,Position,Strand,Motif,Offset" > $methylation_sites
	;;
esac

### detect patches
$root/bin/detect.pl \
    --edist "$edist" \
    --ndist "$ndist" \
    --bases "$bases" \
    --nsites "$nsites" \
    --msfile "$methylation_sites" \
    "$reference" \
    $rundir/../1-strata/*/modifications.csv.gz > detect.csv

echo ""
echo "Finished on $(date)"

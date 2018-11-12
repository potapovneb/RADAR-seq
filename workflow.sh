#!/bin/bash

module load smrtlink

samples=$1

root=$PWD

while read line
do
    if [[ $line =~ "SampleID" ]] ; then
	continue
    fi

    if [[ $line =~ "#" ]] ; then
	continue
    fi

    sampleId=`echo $line | cut -d, -f1`
    sample=`printf "%05i" $sampleId`
    collectionPathUri=`echo $line | cut -d, -f2`
    reference=`echo $line | cut -d, -f5`
    enzyme=`echo $line | cut -d, -f6`
    bases=`echo $line | cut -d, -f7`
    barcodes=`echo $line | cut -d, -f9`
    instrument=`echo $line | cut -d, -f10`
    rundir=`printf "%s/samples/%05i.1" $root $sampleId`

    nsites=4
    
    echo ""
    echo "sampleId          : $sampleId ($sample)"
    echo "collectionPathUri : $collectionPathUri"
    echo "reference         : $reference"
    echo "enzyme            : $enzyme"
    echo "bases             : $bases"
    echo "barcodes          : $barcodes"
    echo "instrument        : $instrument"
    echo "rundir            : $rundir"

    ### Complete protocol
    mkdir -p "$rundir"
    cd "$rundir"

    if [ "$barcodes" == "none" ]
    then
	echo ""
	echo "Regular workflow"

	qsub \
    	    -cwd \
    	    -v instrument="$instrument",root="$root",sample="$sample",collectionPathUri="$collectionPathUri",rname="$reference",protocol="RS_Modification_Detection.1",enzyme="$enzyme",bases="$bases",edist="30",ndist="30",nsites="$nsites" \
    	    -o workflow.log \
    	    -j yes \
    	    "$root"/scripts/workflow.sh
    else
	echo ""
	echo "Barcoded workflow"

	### barcode sequences
	pacbio_barcodes="$root/references/pacbio_barcodes_96_revseq.fasta"

	### score barcodes
	echo ""
	echo "Scoring barcodes"
	bam2bam "$collectionPathUri" -o m.barcoded --barcodes "$pacbio_barcodes"
	bamfile="$rundir/m.barcoded.subreads.bam"

	### analyze barcoded samples
	IFS=';' read -r -a array <<< "$barcodes"

	for bc in ${array[@]}
	do
	    echo ""
	    echo "Sample $sampleId :: Barcode $bc"

	    ### make barcode-specific run directory
	    bcdir="$rundir/$bc"
	    echo "  + $bcdir"
	    echo "  + $reference"
	    mkdir -p "$bcdir"

	    ### extract barcoded sequences
	    barcodeId=`printf "S,%i,%i" $((bc-1)) $((bc-1))`
	    echo "  + $barcodeId"
	    $root/bin/barco.pl --bc $barcodeId $bamfile > $bcdir/temp.sam 2>$bcdir/barcodes.log

	    ### convert, sort and index
	    samtools view -Sb $bcdir/temp.sam > $bcdir/temp.bam 2>/dev/null
	    samtools sort $bcdir/temp.bam > $bcdir/movie.subreads.bam 2>/dev/null
	    samtools index $bcdir/movie.subreads.bam 2>/dev/null

	    ### PacBio-specific indexing
	    pbindex $bcdir/movie.subreads.bam 2>/dev/null
	    
	    ### cleanup
	    rm -f "$bcdir/temp.sam"
	    rm -f "$bcdir/temp.bam"
	    rm -f "$bcdir/movie.subreads.sam"

	    ### run Sequel scripts
	    cd "$bcdir"
	    collectionPathUri="$bcdir/movie.subreads.bam"

	    qsub \
    		-cwd \
    		-v instrument="$instrument",root="$root",sample="$sample",collectionPathUri="$collectionPathUri",rname="$reference",protocol="RS_Modification_Detection.1",enzyme="$enzyme",bases="$bases",edist="30",ndist="30",nsites="$nsites" \
    		-o workflow.log \
    		-j yes \
    		"$root"/scripts/workflow.sh
	done
    fi
done < $samples

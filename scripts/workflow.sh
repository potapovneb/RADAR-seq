#!/bin/bash

#$ -S /bin/bash
#$ -N ctl
#$ -P longrun
#$ -cwd

rundir=$PWD

### Step 1. RS_Modification_Detection.1 & Stratification
cd "$rundir"

case "$instrument" in
    RSII)
	qsub \
	    -sync yes \
	    -v root="$root",rname="$rname",protocol="RS_Modification_Detection.1",collectionPathUri="$collectionPathUri" \
	    -o workflow-strata.log \
	    -j yes \
	    $root/scripts/workflow-strata.sh
	;;
    SEQUEL)
	qsub \
	    -sync yes \
	    -v root="$root",rname="$rname",protocol="RS_Modification_Detection.1",collectionPathUri="$collectionPathUri" \
	    -o workflow-strata-sequel.log \
	    -j yes \
	    $root/scripts/workflow-strata-sequel.sh
	;;
    *)
	echo "[ERROR] cannot recognize PacBio instrument name"
	exit
	;;
esac


# ### Step 2. Patch Detection
# cd "$rundir"

mkdir 2-detect
cd 2-detect

qsub \
    -sync yes \
    -v root="$root",rname="$rname",enzyme="$enzyme",bases="$bases",edist="30",ndist="30",nsites="4" \
    -o workflow-detect.log \
    -j yes \
    $root/scripts/workflow-detect.sh


### Step 3. Patch Prediction
cd "$rundir"

mkdir 3-predict
cd 3-predict

qsub \
    -sync yes \
    -v root="$root",rname="$rname",patches="$rundir/2-detect/detect.csv",model="$root/svm_model",features="Length;mA;mC;mTotal;A;C;iA;iC;iTtotal;Modifiable;Fraction",mTotal="5" \
    -o workflow-predict.log \
    -j yes \
    "$root"/scripts/workflow-predict.sh

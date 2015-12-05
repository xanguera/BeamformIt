#!/bin/bash

# Copyright 2015, Mitsubishi Electric Research Laboratories, MERL (Author: Shinji Watanabe)

# Set bash to 'debug' mode, it will exit on :
# -e 'error', -u 'undefined variable', -o ... 'error in pipeline', -x 'print commands',
set -e
set -u
set -o pipefail
#set -x

# CHiME3 data directory. Please specify the CHiME3 data root
chime3_data=/db/laputa1/data/processed/public/CHiME3 # MERL environment
sdir=$chime3_data/data/audio/16kHz/isolated
if [ ! -d $sdir ]; then
  echo "Please specify CHiME3 data root"
fi

# Enhanced data are stored
odir=chime3/enhanced_data
# Logfiles, data list etc.
wdir=chime3/list
ldir=chime3/log

# BeamformIt binary
BEAMFORMIT=../BeamformIt
if [ ! -e $BEAMFORMIT ] ; then
  echo "Please specify BeamformIt binary"
fi

# BeamformIt config
bconf=../cfg-files/chime3.cfg
if [ ! -e $bconf ] ; then
  echo "No config file found at $bconf"
fi

# we use the following channel signals, and remove 2nd channel signal, which located on the back of
# tablet, and behaves very different from the other front channel signals.
bmf="1 3 4 5 6"

echo "Enhanced data would be stored at $odir"
mkdir -p $odir
mkdir -p $wdir
mkdir -p $ldir

echo "Will use the following channels: $bmf"
# number of channels
numch=`echo $bmf | tr ' ' '\n' | wc -l`
echo "the number of channels: $numch"

# wavfiles.list can be used as the name of the output files
output_wavfiles=$wdir/wavfiles.list
find $sdir/*{simu,real} | grep CH1.wav \
  | awk -F '/' '{print $(NF-1) "/" $NF}' | sed -e "s/\.CH1\.wav//" | sort > $output_wavfiles

# this is an input file list of the microphones
# format: 1st_wav 2nd_wav ... nth_wav
input_arrays=$wdir/channels_$numch
for x in `cat $output_wavfiles`; do
  echo -n "$x"
  for ch in $bmf; do
    echo -n " $x.CH$ch.wav"
  done
  echo ""
done > $input_arrays

# making a subdirectory for the output wav files
for x in `awk -F '/' '{print $1}' $output_wavfiles | sort | uniq`; do
  mkdir -p $odir/$x
done

echo "Beamforming"
set -x ## set -x 'print commands',
# making a shell script for each job
while read line; do
  $BEAMFORMIT -s $line -c $input_arrays \
    --config_file $bconf \
    --source_dir $sdir \
    --result_dir $odir
done < $output_wavfiles > $ldir/beamform.log

echo "`basename $0` Done."

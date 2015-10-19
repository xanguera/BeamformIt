#!/bin/bash
#helper script that performs beamforming on all wav files found in a given directory
#usage: $0 audio_source_dir output_name

inputDir=$1
outName=$2

#check the number of parameters
if [ $# -ne 2 ]
then
    echo "Usage is: $0 audio_source_dir output_name"
    exit
fi


rm -Rf output/${outName}
mkdir -p output/${outName}

#create the channels file from the list of audio files in the source directory
echo -n "${outName}" > output/${outName}/channels_file
for file in ${inputDir}/*.wav ${inputDir}/*.sph #add any more extensions as needed
do
    if [ -a ${file} ]; then
        echo -n " ${file}" >> output/${outName}/channels_file  
    fi
done

./BeamformIt \
    --scroll_size 250 \
    --window_size 500 \
    --nbest_amount 4 \
    --do_noise_threshold 1 \
    --noise_percent 10 \
    --trans_weight_multi 25 \
    --trans_weight_nbest 25 \
    --print_features 1 \
    --do_avoid_bad_frames 1 \
    --do_compute_reference 1 \
    --do_use_uem_file 0 \
    --do_adapt_weights 1 \
    --do_write_sph_files 1 \
    --channels_file ./output/${outName}/channels_file \
    --show_id ${outName} \
    --result_dir ./output/${outName}




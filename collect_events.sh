#!/usr/bin/env bash

usage="./collect_events.sh results_folder"
results_folder=$1

if [ -z "$results_folder" ]
then
    echo $usage
    exit 1
fi

collectedFolder=$results_folder\_collected

if [ ! -d "$collectedFolder" ]; then
    mkdir $collectedFolder
fi

for ijob in `ls --color=none $results_folder`
do 
    jobid=`echo $ijob | sed 's/job-//'`
    jobPath=$results_folder/$ijob
    if [ -d "$jobPath/results" ]; then
        if [ -e "$jobPath/results/Charged_eta_integrated_vndata.dat" ]; then
            mv $jobPath/results $collectedFolder/event-$jobid
            IPGlasmaFile=`cat $jobPath/music_input_2 | grep "Initial_Distribution_Filename" | sed 's/Initial_Distribution_Filename  //'`
            IPGlasmaPath=`cat $jobPath/music_input_2 | grep "Initial_Distribution_Filename" | sed 's/Initial_Distribution_Filename  //' | sed 's$/u_field_[0-9]*.dat$$' | cut -f 1 -d " "`
            cp $IPGlasmaFile $collectedFolder/event-$jobid
            cp $IPGlasmaPath/config_file.cfg $collectedFolder/event-$jobid
        else
            echo $ijob "is not done ..."
        fi
    fi
done

echo "finished."

#!/usr/bin/env bash

folder_name=$1
echo Moving all the results into $folder_name ... 
mkdir $folder_name
rm -fr spvn*/EOS
rm -fr spvn*/mpihydro
rm -fr surface.dat
rm -fr *.log
rm -fr *.err
mv *.dat $folder_name
mv spvn* $folder_name
cp music_input_2 $folder_name/music_input
mv $folder_name/known_nuclei.dat ./
mv $folder_name/eps_freeze_list_*.dat ./
rm -fr music_input_*

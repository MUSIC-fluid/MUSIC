#!/usr/bin/env bash

folder_name=$1
echo Moving all the results into $folder_name ... 
mkdir $folder_name
mv *.dat $folder_name
mv *.gp $folder_name
cp input $folder_name
mv $folder_name/known_nuclei.dat ./

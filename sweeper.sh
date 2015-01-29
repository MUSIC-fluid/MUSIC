#!/usr/bin/env bash

folder_name=$1
echo Moving all the results into $folder_name ... 
mkdir $folder_name
mv *.dat $folder_name
mv *.err $folder_name
mv *.log $folder_name
cp input $folder_name
mv $folder_name/known_nuclei.dat ./
rm -fr $folder_name/tmp*.dat

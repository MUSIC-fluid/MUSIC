#!/usr/bin/env bash
# Collect the MUSIC output in the current directory into a new folder.
# usage: bash sweeper.sh <folder_name>

if [ $# -ne 1 ] || [ -z "$1" ]; then
    echo "usage: $0 <folder_name>" >&2
    exit 1
fi
folder_name=$1
if [ -e "${folder_name}" ]; then
    echo "${folder_name} already exists; not overwriting it." >&2
    exit 1
fi
echo "Moving all the results into ${folder_name} ..."

shopt -s nullglob

# combine multiple surface files
for ii in surface_eps_*_*.dat; do
    jj=$(echo "$ii" | cut -f 1-3 -d _)
    cat "$ii" >> "$jj.dat"
done
rm -f surface_eps_*_*.dat

mkdir "${folder_name}" || exit 1
results=(*.dat *.err *.log)
if [ ${#results[@]} -gt 0 ]; then
    mv "${results[@]}" "${folder_name}/"
fi
if [ -f music_input_mode_2 ]; then
    cp music_input_mode_2 "${folder_name}/music_input"
fi

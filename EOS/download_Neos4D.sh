#!/usr/bin/env bash

USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="neos"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0/"

# NEOS version and considered quantities
EOS_TYPE="urqmd"
FILE_NAME_LIST=("cs" "mub" "muq" "mus" "p" "t")

for name in "${FILE_NAME_LIST[@]}"; do
    # Define file name
    FILE_PATH="EoS_UrQMD/neos4d_${EOS_TYPE}_${name}_b.dat" 
    # Set the name of URL
    API_FILE_URL="${API_BASE_URL}repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

    LOCAL_DESTINATION="neos4d_${EOS_TYPE}_${name}_b.dat"
    curl -L -o "$LOCAL_DESTINATION" "$API_FILE_URL"
done

mkdir -p neos4D
mv neos4d_*_b.dat neos4D/

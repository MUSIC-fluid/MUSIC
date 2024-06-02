#!/usr/bin/env bash

USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="neos"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0/"

# NEOS version and considered quantities
# UrQMD or pdg

EOS_TYPE=$1
if [ -z ${EOS_TYPE} ]
then
    EOS_TYPE="UrQMD"
fi

# Check for cs file before downloading
if [[ "$EOS_TYPE" == "UrQMD" ]]; then
    FILE_NAME_LIST=("cs" "mub" "muq" "mus" "p" "t")
else
    FILE_NAME_LIST=("mub" "muq" "mus" "p" "t")
fi

# Download 4D EoS
for name in "${FILE_NAME_LIST[@]}"; do
    # Define file name
    FILE_PATH="EoS4D_${EOS_TYPE}/neos4d_${name}_b.dat" 
    # Set the name of URL
    API_FILE_URL="${API_BASE_URL}repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

    LOCAL_DESTINATION="neos4d_${name}_b.dat"
    curl -L -o "$LOCAL_DESTINATION" "$API_FILE_URL"
done

# Move to correct folder for use in MUSIC
mkdir -p neos4D
mv neos4d_*_b.dat neos4D/

# Keep a label for the EoS tables
touch neos4D/Neos4D_${EOS_TYPE}
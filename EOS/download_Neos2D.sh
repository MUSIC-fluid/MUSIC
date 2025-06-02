#!/usr/bin/env bash

FILETYPE=$1
if [ -z ${FILETYPE} ]
then
    # default download NEoS-BQS
    FILETYPE="bqs"
fi

USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="neos"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0/"

# NEOS version and considered quantities
EOS_TYPE="pdg"

# Define file name
FILE_PATH="EoS_2D_${EOS_TYPE}/neos_${FILETYPE}.tar.gz"
# Set the name of URL
API_FILE_URL="${API_BASE_URL}repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

LOCAL_DESTINATION="neos_${FILETYPE}.tar.gz"
curl -L -o "$LOCAL_DESTINATION" "$API_FILE_URL"

tar -xf neos_${FILETYPE}.tar.gz
rm -fr neos_${FILETYPE}.tar.gz


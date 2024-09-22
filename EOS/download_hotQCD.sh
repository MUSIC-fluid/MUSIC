#!/usr/bin/env bash

FILETYPE=$1
if [ -z ${FILETYPE} ]
then
    # options: binary (match with UrQCD, chemical equilibrium) [default]
    #          SMASH_binary (match with SMASH, chemical equilibrium)
    FILETYPE="binary"
fi

USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="hotqcd"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0"

# Define file name
FILE_PATH="hrg_hotqcd_eos_${FILETYPE}.dat"
# Set the name of URL
API_FILE_URL="${API_BASE_URL}/repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

curl -L -o "${FILE_PATH}" "${API_FILE_URL}"

mkdir -p hotQCD
mv ${FILE_PATH} hotQCD/${FILE_PATH}

#!/usr/bin/env bash

FILETYPE=$1
if [ -z ${FILETYPE} ]
then
    # default download v1.2 (match with UrQMD, chemical equilibrium)
    # options: v1 (match with pdg05, chemical equilibrium)
    #          PCE-v1 (match with pdg05, PCE at T = 150 MeV)
    #          PCE155 (match with pdg05, PCE at T = 155 MeV)
    #          PCE160 (match with pdg05, PCE at T = 160 MeV)
    #          PCE165-v0 (match with pdg05, PCE at T = 165 MeV)
    FILETYPE="v1.2"
fi

USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="s95p"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0"

# Define file name
FILE_PATH="s95p-${FILETYPE}.tar.gz"
# Set the name of URL
API_FILE_URL="${API_BASE_URL}/repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

LOCAL_DESTINATION="s95p-${FILETYPE}.tar.gz"
curl -L -o "$LOCAL_DESTINATION" "$API_FILE_URL"

tar -xf s95p-${FILETYPE}.tar.gz
rm -fr s95p-${FILETYPE}.tar.gz


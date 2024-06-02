#!/usr/bin/env bash

USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="hotqcd"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0"

# Define file name
FILE_PATH="EOS-Q.tar.gz"
# Set the name of URL
API_FILE_URL="${API_BASE_URL}/repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

curl -L -o "${FILE_PATH}" "${API_FILE_URL}"

tar -xf ${FILE_PATH}
rm -fr ${FILE_PATH}

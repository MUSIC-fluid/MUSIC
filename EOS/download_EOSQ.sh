#!/usr/bin/env bash
set -euo pipefail
# download into this script's directory (EOS/), where MUSIC looks
cd "$(dirname "${BASH_SOURCE[0]}")"

USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="eos_q"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0"

# Define file name
FILE_PATH="EOS-Q.tar.gz"
# Set the name of URL
API_FILE_URL="${API_BASE_URL}/repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

curl -fL -o "${FILE_PATH}" "${API_FILE_URL}"

tar -xf ${FILE_PATH}
rm -fr ${FILE_PATH}

#!/usr/bin/env bash

# NOTE: The environment variables are present during the CI job.

URL_RELEASES=https://api.github.com/repos/$GITHUB_USER/$GITHUB_REPONAME/releases

# create the release
API_JSON=$(printf '{"tag_name":"%s",
                    "target_commitish":"main",
                    "name":"%s",
                    "body":"",
                    "draft":false,
                    "prerelease":false}' $CI_COMMIT_TAG $CI_COMMIT_TAG)
AUTH_HEADER="Authorization: token $GITHUB_RELEASE_TOKEN"
curl $URL_RELEASES --data "$API_JSON" -H "$AUTH_HEADER"

# sleep 5 seconds, otherwise the check below may fail
sleep 5

# get latest release tag name
LATEST_RELEASE=$(curl -s "$URL_RELEASES"/latest | grep -oP '"tag_name": "\K(.*)(?=")')

if [ "$LATEST_RELEASE" != "$CI_COMMIT_TAG" ]; then
 echo "Error: Could not properly create the release!" 1>&2
 exit 1
fi

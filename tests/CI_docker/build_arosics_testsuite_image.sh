#!/usr/bin/env bash

pkgname="arosics"
repourl="https://git.gfz-potsdam.de/danschef/arosics"

context_dir="./context"
dockerfile="${pkgname}_ci.docker"
python_script='
version = {}
with open("../../'${pkgname}'/version.py") as version_file:
    exec(version_file.read(), version)
print(version["__version__"])
'
version=`python -c "$python_script"`
tag="ds__${pkgname}_ci:$version"
gitlab_runner="${pkgname}_gitlab_CI_runner"
runnername_remote="${pkgname}_ci_runner__v${version}__${HOSTNAME}"
taglist="${pkgname}_ci_client"

echo "#### Build runner docker image"
docker rmi ${tag}
docker build ${context_dir} \
    --no-cache \
    -f ${context_dir}/${dockerfile} \
    -m 20G \
    -t ${tag}

# create the gitlab-runner docker container for the current project
# NOTE: The 'gitlab-runner' and 'gitlab-ci-multi-runner' services will run within this container.
#       The runner uses a 'config.toml' configuration file at /etc/gitlab-runner within the container which can be
#       modified through additional parameters of the 'gitlab-runner register' command.
echo "#### Create gitlab-runner (daemon) container with tag; ${tag}"
docker stop ${gitlab_runner}
docker rm ${gitlab_runner}
docker run \
    -d \
    --name ${gitlab_runner} \
    --restart always \
    -v /var/run/docker.sock:/var/run/docker.sock \
    gitlab/gitlab-runner:latest

# register the runner at the corresponding GitLab repository via a registration-token
echo "
--------------------------------------------------------------------------
To register the runner at GitLab, go to ${repourl}/-/runners,
click on 'New project runner' and use the following settings:

Tags:                       ${taglist}
Run untagged jobs:          Yes
Runner description:         ${runnername_remote}
Paused:                     No
Protected:                  No
Lock to current projects:   Yes
Maximum job timeout:        7200

Then click 'Create runner'!
--------------------------------------------------------------------------"
read -p "Please enter the GitLab runner authentification token (should start with 'glrt-'): " token
# NOTE: In case of locally stored images (like here), the docker pull policy 'never' must be used
#       (see https://docs.gitlab.com/runner/executors/docker.html#how-pull-policies-work).
docker exec -it ${gitlab_runner} /bin/bash -c "\
gitlab-ci-multi-runner register \
  --non-interactive \
  --executor 'docker' \
  --docker-image '${tag}' \
  --url 'https://git.gfz-potsdam.de' \
  --token '${token}' \
  --description '${runnername_remote}' \
  --docker-pull-policy='never'
  "
ls

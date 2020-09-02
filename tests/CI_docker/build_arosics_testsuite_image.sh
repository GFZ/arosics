#!/usr/bin/env bash

context_dir="./context"
dockerfile="arosics_ci.docker"
tag="arosics_ci:0.9.20"
gitlab_runner="arosics_gitlab_CI_runner"

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

# register the runner at the corresponding GitLab repository via a registration-tok
echo "#### Register container at gitlab, get token here https://gitext.gfz-potsdam.de/danschef/arosics/settings/ci_cd"
read -p "Please enter gitlab token: " token
echo ""
read -p "Please enter gitlab runner name: " runner_name
echo "New gitlab runner image will named  ${gitlab_runner}"
# NOTE: In case of locally stored images (like here), the docker pull policy 'never' must be used
#       (see https://docs.gitlab.com/runner/executors/docker.html#how-pull-policies-work).
docker exec -it ${gitlab_runner} /bin/bash -c "\
export RUNNER_EXECUTOR=docker && \
gitlab-ci-multi-runner register \
  --non-interactive \
  --executor 'docker' \
  --docker-image '${tag}' \
  --url 'https://gitext.gfz-potsdam.de/ci' \
  --registration-token '${token}' \
  --description '${runner_name}' \
  --tag-list arosics_ci_client \
  --run-untagged='true' \
  --locked='true' \
  --access-level='not_protected' \
  --docker-pull-policy='never'
ls

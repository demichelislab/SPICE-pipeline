#!/usr/bin/env bash

usage () {
  echo "USAGE: $0 step_name"
}

required_args=1
. $(dirname $0)/parse_args_bash.sh

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
step_name="${1%%+(/)}"

for pref in logs data;
do
  cur="${project_dir}/${pref}/$step_name"
  [[ -d "$cur" ]] || mkdir -p "$cur"
done


#!/usr/bin/env bash

if [[ ! -f "WDL/cromwell-57.jar" ]]; then
  wget https://github.com/broadinstitute/cromwell/releases/download/57/cromwell-57.jar -o WDL/cromwell-57.jar
fi

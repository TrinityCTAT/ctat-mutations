#!/usr/bin/env bash

if [[ ! -f "WDL/cromwell-58.jar" ]]; then
  wget https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O WDL/cromwell-58.jar
fi

#!/bin/bash

set -e

VERSION=`cat VERSION.txt`

docker build -t trinityctat/ctat_mutations:$VERSION .
docker build -t trinityctat/ctat_mutations:latest .



#!/bin/bash

set -e

VERSION=`cat VERSION.txt`

docker push trinityctat/ctat_mutations:$VERSION 
docker push trinityctat/ctat_mutations:latest 



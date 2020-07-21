#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`

singularity build --disable-cache  ctat_mutations.v${VERSION}.simg docker://trinityctat/ctat_mutations:$VERSION

singularity exec -e ctat_mutations.v${VERSION}.simg env

ln -sf  ctat_mutations.v${VERSION}.simg  ctat_mutations.vLATEST.simg  #for local testing


#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`

singularity build ctat_mutations.v${VERSION}.simg docker://trinityctat/ctat_mutations:$VERSION

singularity exec -e ctat_mutations.v${VERSION}.simg env

cp ctat_mutations.v${VERSION}.simg ctat_mutations.vLATEST.simg


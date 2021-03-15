#!/usr/bin/env bash

pytest

docker run -e CTAT_GENOME_LIB='/genome' -v $(pwd):/data -v "$CTAT_GENOME_LIB":/genome --rm trinityctat/ctat_mutations:2.6.0-alpha.1 \
  /bin/bash -c "cd /usr/local/src/ctat-mutations && pip install pytest && pytest"

singularity build -F ctat_mutations.simg docker://trinityctat/ctat_mutations:2.6.0-alpha.1

singularity exec -B "$CTAT_GENOME_LIB":"$CTAT_GENOME_LIB" ctat_mutations.simg \
  /bin/bash -c "cd /usr/local/src/ctat-mutations && python -m pytest"

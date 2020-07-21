#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals
import os, sys, re
import logging
import argparse
import subprocess

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="add rvboost score annotations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--input_vcf", type=str, required=True, help="input vcf file")

    parser.add_argument(
        "--scores", type=str, required=True, help="RVB-like-R variant scores table"
    )

    parser.add_argument(
        "--output_vcf",
        type=str,
        required=True,
        help="output vcf file including score annotations",
    )

    args = parser.parse_args()

    ## get orig score annotations.

    scores_file = args.scores

    if not os.path.exists(scores_file):
        logger.critical("Error, cannot find file: {}".format(scores_file))
        sys.exit(1)

    chrpos_to_score_annots = dict()
    with open(scores_file, "rt") as fh:
        line = next(fh)
        for line in fh:
            line = line.rstrip()
            (chrpos, RVBLR, RVBLRQ) = line.split("\t")
            annot_text = "RVBLR={:0.3f};RVBLRQ={:0.3f}".format(
                float(RVBLR), float(RVBLRQ)
            )
            chrpos_to_score_annots[chrpos] = annot_text

    ## write annotated vcf
    counter = -1

    ofh = open(args.output_vcf, "wt", encoding='utf-8')

    open_file = open
    if re.search("\.gz$", args.input_vcf):
        import gzip

        open_file = gzip.open

    with open_file(args.input_vcf, "rt", encoding='utf-8') as fh:

        for line in fh:
            line = line.rstrip()
            if line[0] == "#":
                # header line

                if re.match("#CHROM\tPOS", line):
                    print(
                        '##INFO=<ID=RVBLR,Number=1,Type=Float,Description="RVBLR Boosting Original score">\n'
                        + '##INFO=<ID=RVBLRQ,Number=1,Type=Float,Description="RVBLR Boosting algorithm Q-Score according to ecdf of RVBLR on common variants">',
                        file=ofh,
                    )

                print(line, file=ofh)
                continue

            counter += 1
            vals = line.split("\t")
            chrom = vals[0]
            pos = vals[1]
            refbase = vals[3]
            altbase = vals[4]

            chrpos = ":".join([chrom, pos, refbase, altbase])
            annot = chrpos_to_score_annots[chrpos]

            vals[7] += ";{}".format(annot)

            print("\t".join(vals), file=ofh)

    print("Done.", file=sys.stderr)

    sys.exit(0)


####################

if __name__ == "__main__":
    main()

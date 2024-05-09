#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging
sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"]))
from Pipeliner import Pipeliner, Command

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)




def main():

    parser = argparse.ArgumentParser(description="normalize bam in a strand-specific manner", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--input_bam", type=str, required=True, help="input bam filename")
    parser.add_argument("--output_bam", type=str, required=True, help="output for normalied bam file")
    parser.add_argument("--normalize_max_cov_level", type=int, default=1000, help="normalize to max read coverage level before assembly (default: 1000)")


    args = parser.parse_args()

    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam
    normalize_max_cov_level = args.normalize_max_cov_level

    pipeliner = Pipeliner("__chckpts")
    
    # first separate input bam into strand-specific files
    SS_output_prefix = os.path.basename(input_bam_filename) + ".SS"

    scriptdir = os.path.abspath(os.path.dirname(__file__))
    cmd = " ".join([os.path.join(scriptdir, "separate_bam_by_strand.py"),
                    "--bam {}".format(input_bam_filename),
                    "--output_prefix {}".format(SS_output_prefix)])
    
    
    pipeliner.add_commands([Command(cmd, "sep_by_strand.ok")])

    ## run normalizations
    bamsifter_prog = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../plugins/bamsifter/bamsifter")
    SS_bam_files = [SS_output_prefix + x for x in (".+.bam", ".-.bam") ]

    SS_norm_bam_files = list()
    
    for SS_bam_file in SS_bam_files:

        norm_bam_filename = f"{SS_bam_file}.{normalize_max_cov_level}.bam"
        cmd = " ".join([bamsifter_prog,
                        f" -c {normalize_max_cov_level} ",
                        f" -o {norm_bam_filename} ",
                        SS_bam_file])
        
        norm_bam_checkpoint = norm_bam_filename + ".ok"

        pipeliner.add_commands([Command(cmd, norm_bam_checkpoint)])
        
        SS_norm_bam_files.append(norm_bam_filename)
        

    # merge the norm SS bam filenames into the final output file

    cmd = f"samtools merge {output_bam_filename} " + " ".join(SS_norm_bam_files)
    pipeliner.add_commands([Command(cmd, "SS_merge.ok")])


    cmd = f"samtools index {output_bam_filename}"
    pipeliner.add_commands([Command(cmd, "index_merged.ok")])

    pipeliner.run()

    logger.info("Done.  See SS-normalized bam: {}".format(output_bam_filename))

    sys.exit(0)

    
if __name__=='__main__':
    main()

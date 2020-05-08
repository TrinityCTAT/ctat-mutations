#!/usr/bin/env python
import argparse
import os, sys
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate with Cravat")
    parser.add_argument("input", help="Input vcf file.")
    parser.add_argument("output", help="Output base name")

    parser.add_argument("--genome", help="Genome", default='hg38')
    parser.add_argument("--classifier", help="Classifier (ignored)")

    args = parser.parse_args()
    input_vcf = args.input
    output = args.output

    

    genome = args.genome
    output_dir = os.path.dirname(input_vcf)

    # run open-cravat

    cravat_cmd = ['oc', 'run', '--cs', '{"run":{"vcfreporter":{"type":"separate"}}}', input_vcf, '-t', 'vcf', '-l',
                  args.genome, '-d', output_dir, '-n', output]
    
    sys.stderr.write("annotate_with_cravat: {}\n".format(" ".join(cravat_cmd)))
    subprocess.check_call(cravat_cmd)

    sys.exit(0)

    

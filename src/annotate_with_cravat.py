#!/usr/bin/env python
import argparse
import os, sys, re
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate with Cravat")
    parser.add_argument("--input_vcf", help="Input vcf file.", required=True)
    parser.add_argument(
        "--output_vcf", help="Output cravat-annotated vcf file", required=True
    )

    parser.add_argument("--genome", help="genome (hg19 or hg38)", default="hg38")
    parser.add_argument(
        "--cravat_lib_dir", help="cravat resources lib dir", required=True
    )

    args = parser.parse_args()

    input_vcf = args.input_vcf
    output_vcf = args.output_vcf
    cravat_lib_dir = args.cravat_lib_dir

    genome = args.genome
    output_dir = os.path.dirname(output_vcf)
    output_vcf = os.path.basename(output_vcf)

    if re.search(".vcf$", output_vcf):
        output_vcf, count = re.subn(".vcf$", "", output_vcf)
        if count != 1:
            raise RuntimeError(
                "Error, couldnt replace .vcf extension of {}".format(output_vcf)
            )

    # set cravat resource lib:
    # subprocess.check_call(["oc", "config", "md", cravat_lib_dir])

    # run open-cravat
    cravat_cmd = [
        "oc",
        "run",
        "--cs",
        '{"run":{"vcfreporter":{"type":"separate"}}}',
        input_vcf,
        "--system-option",
        "modules_dir={}".format(cravat_lib_dir),
        "-t",
        "vcf",
        "-l",
        args.genome,
        "-d",
        output_dir,
        "-n",
        output_vcf,
    ]

    sys.stderr.write("annotate_with_cravat: {}\n".format(" ".join(cravat_cmd)))
    subprocess.check_call(cravat_cmd)

    sys.exit(0)

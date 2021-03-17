import glob
import os
import shutil
from subprocess import check_call


def check_output(path):
    output_path = os.path.join(path, '**', 'call-VariantFiltration', '*', 'test.*.vcf.gz')
    filtered_vcfs = list(glob.iglob(output_path, recursive=True))
    assert len(filtered_vcfs) == 1, "{} VCFs found".format(len(filtered_vcfs))
    shutil.rmtree(path)


# test all pairwise combinations of boosting method and alg type to ensure they run successfully
def test_boosting(boosting_method, boosting_alg_type, tmp_path):
    tmp_path = str(tmp_path)
    check_call(
        [os.path.abspath('./ctat_mutations'),
         '--left', os.path.abspath('testing/reads_1.fastq.gz'), '--right', os.path.abspath('testing/reads_2.fastq.gz'),
         '--sample_id', 'test', '--boosting_method', boosting_method, '--boosting_alg_type', boosting_alg_type],
        cwd=tmp_path)
    check_output(tmp_path)


# test the full pipeline with fastq input
def test_fastqs(tmp_path):
    tmp_path = str(tmp_path)
    check_call(
        [os.path.abspath('./ctat_mutations'), '--left',
         os.path.abspath('testing/reads_1.fastq.gz'), '--right', os.path.abspath('testing/reads_2.fastq.gz'),
         '--sample_id', 'test', '--boosting_method', 'none'], cwd=tmp_path)
    check_output(tmp_path)


# test the full pipeline with single fastq input
def test_single_fastq(tmp_path):
    tmp_path = str(tmp_path)
    check_call(
        [os.path.abspath('./ctat_mutations'), '--left',
         os.path.abspath('testing/reads_1.fastq.gz'),
         '--sample_id', 'test', '--boosting_method', 'none'], cwd=tmp_path)
    check_output(tmp_path)

#
# # test the pipeline with pre-aligned bam input and diff with known good output
# def test_bam(tmp_path):
#     tmp_path = str(tmp_path)
#     check_call(
#         [os.path.abspath('./ctat_mutations'), '--bam',
#          os.path.abspath('testing/reads_output/test.star.Aligned.sortedByCoord.out.bam'),
#          '--sample_id', 'test', '--boosting_method',
#          'none'], cwd=tmp_path)
#     check_output(tmp_path)
#
#
# # test the pipeline with VCF input, annotate and filter VCF, and diff with known good output
# def test_vcf(tmp_path):
#     tmp_path = str(tmp_path)
#     check_call(
#         [os.path.abspath('./ctat_mutations'), '--bam', os.path.abspath('testing/reads_output/test.bqsr.bam'), '--vcf',
#          os.path.abspath('testing/reads_output/test.vcf.gz'), '--sample_id', 'test', '--boosting_method',
#          'none'], cwd=tmp_path)
#     check_output(tmp_path)

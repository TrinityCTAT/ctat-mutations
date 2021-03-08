import glob
import os
import shutil
from subprocess import check_call


def check_output(path):
    filtered_vcfs = list(glob.iglob(path + '/**/call-VariantFiltration/*/test.filtered.vcf.gz', recursive=True))
    assert len(filtered_vcfs) == 1, "VCF"
    shutil.rmtree(path)


def test_benchmarking(tmp_path):
    tmp_path = str(tmp_path)
    if not os.path.exists('../CTAT-benchmarking'):
        raise ValueError(
            'Please clone https://github.com/broadinstitute/CTAT-benchmarking.git to ' + os.path.abspath('..'))
    check_call([os.path.abspath('../CTAT-benchmarking/CTAT-mutation-benchmarking/BENCHMARK_variant_calling.py'),
                '--pred_vcf', os.path.abspath('testing/reads_output/test.vcf.gz'), '--truth_vcf',
                os.path.abspath('testing/reads_output/test.filtered.vcf.gz'), '--pred_bam',
                os.path.abspath('testing/reads_output/test.bqsr.bam'), '--output_dir', 'out'], cwd=tmp_path)
    check_output(tmp_path)


# test all pairwise combinations of boosting method and alg type to ensure they run successfully
def test_boosting(boosting_method, boosting_alg_type, tmp_path):
    tmp_path = str(tmp_path)
    check_call(
        [os.path.abspath('./ctat_mutations'), '--bam', os.path.abspath('testing/reads_output/test.bqsr.bam'), '--vcf',
         os.path.abspath('testing/reads_output/test.vcf.gz'), '--sample_id', 'test', '--boosting_method',
         boosting_method, '--boosting_alg_type', boosting_alg_type], cwd=tmp_path)
    check_output(tmp_path)


# test the full pipeline with fastq input and diff with known good output
def test_fastqs(tmp_path):
    tmp_path = str(tmp_path)
    check_call(
        [os.path.abspath('./ctat_mutations'), '--left',
         os.path.abspath('testing/reads_1.fastq.gz'), '--right', os.path.abspath('testing/reads_2.fastq.gz'),
         '--sample_id', 'test', '--boosting_method', 'none', '--read_var_pos_annotation_cpu', '1'], cwd=tmp_path)
    check_output(tmp_path)


# test the full pipeline with single fastq input and diff with known good output
def test_single_fastq(tmp_path):
    tmp_path = str(tmp_path)
    check_call(
        [os.path.abspath('./ctat_mutations'), '--left',
         os.path.abspath('testing/reads_1.fastq.gz'),
         '--sample_id', 'test', '--boosting_method', 'none', '--read_var_pos_annotation_cpu', '1'], cwd=tmp_path)
    check_output(tmp_path)


# test the pipeline with pre-aligned bam input and diff with known good output
def test_bam(tmp_path):
    tmp_path = str(tmp_path)
    check_call(
        [os.path.abspath('./ctat_mutations'), '--bam',
         os.path.abspath('testing/reads_output/test.star.Aligned.sortedByCoord.out.bam'),
         '--sample_id', 'test', '--boosting_method',
         'none', '--read_var_pos_annotation_cpu', '1'], cwd=tmp_path)
    check_output(tmp_path)


# test the pipeline with VCF input, annotate and filter VCF, and diff with known good output
def test_vcf(tmp_path):
    tmp_path = str(tmp_path)
    check_call(
        [os.path.abspath('./ctat_mutations'), '--bam', os.path.abspath('testing/reads_output/test.bqsr.bam'), '--vcf',
         os.path.abspath('testing/reads_output/test.vcf.gz'), '--sample_id', 'test', '--boosting_method',
         'none', '--read_var_pos_annotation_cpu', '1'], cwd=tmp_path)
    check_output(tmp_path)

from unittest import TestCase
from genotype_state_and_count_filter import apply_filter_criteria
from tempfile import NamedTemporaryFile


test_vcf_content = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##reference=/reference.fa
##contig=<ID=Chr_00>
##contig=<ID=Chr_01>
##contig=<ID=Chr_02>
##contig=<ID=Chr_03>
##contig=<ID=Chr_04>
##contig=<ID=Chr_05>
##FILTER=<ID=FBQualDepth,Description="Set if true: expression">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_01\tSAMPLE_02\tSAMPLE_03\tSAMPLE_04\tSAMPLE_05
Chr_00\t280\t.\tA\tC\t117.22\tPASS\tAC=3;AN=10\tGT\t0|1\t0|1\t0|0\t0|0\t0/1
Chr_00\t284\t.\tA\tG\t117.22\tPASS\tAC=6;AN=10\tGT\t0|0\t0|0\t1|1\t1|1\t1/1
Chr_00\t424\t.\tC\tT\t219.92\tPASS\tAC=7;AN=10\tGT\t0/1\t0/1\t1/1\t1/1\t0/1
Chr_00\t464\t.\tA\tC\t625.87\tPASS\tAC=6;AN=10\tGT\t0/1\t0/1\t0/1\t0/1\t1/1
Chr_00\t464\t.\tA\tC\t625.87\tPASS\tAC=4;AN=10\tGT\t1/1\t1/1\t0/0\t0/0\t0/0"""




class TestApply_filter_criteria(TestCase):

    def test_nref_2_required_het(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="nref",
                              samples_required_het_arg="SAMPLE_01,SAMPLE_02",
                              samples_required_hom_arg=None,
                              min_het_genotype_count_arg=None,
                              max_het_genotype_count_arg=None,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 3, "Number of pass variants is not as expected")

    def test_nref_2_required_het_max_het_count_3(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="nref",
                              samples_required_het_arg="SAMPLE_01,SAMPLE_02",
                              samples_required_hom_arg=None,
                              min_het_genotype_count_arg=None,
                              max_het_genotype_count_arg=3,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 2, "Number of pass variants is not as expected")

    def test_nref_2_required_het_min_het_count_4(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="nref",
                              samples_required_het_arg="SAMPLE_01,SAMPLE_02",
                              samples_required_hom_arg=None,
                              min_het_genotype_count_arg=4,
                              max_het_genotype_count_arg=None,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 1, "Number of pass variants is not as expected")

    def test_nref_2_required_hom(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="nref",
                              samples_required_het_arg=None,
                              samples_required_hom_arg="SAMPLE_03,SAMPLE_04",
                              min_het_genotype_count_arg=None,
                              max_het_genotype_count_arg=None,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 2, "Number of pass variants is not as expected")

    def test_nref_2_required_het_2_required_hom(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="nref",
                              samples_required_het_arg="SAMPLE_01,SAMPLE_02",
                              samples_required_hom_arg="SAMPLE_03,SAMPLE_04",
                              min_het_genotype_count_arg=None,
                              max_het_genotype_count_arg=None,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 1, "Number of pass variants is not as expected")

    def test_major_2_required_hom(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              "major",
                              None,
                              "SAMPLE_03,SAMPLE_04",
                              None,
                              None,
                              None,
                              None)

        self.assertTrue(count_variant_pass == 4, "Number of pass variants is not as expected")

    def test_major_2_required_hom_max_count_het_0(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="major",
                              samples_required_het_arg=None,
                              samples_required_hom_arg="SAMPLE_03,SAMPLE_04",
                              min_het_genotype_count_arg=None,
                              max_het_genotype_count_arg=0,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 2, "Number of pass variants is not as expected")

    def test_minor_2_required_hom(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="minor",
                              samples_required_het_arg=None,
                              samples_required_hom_arg="SAMPLE_01,SAMPLE_02",
                              min_het_genotype_count_arg=None,
                              max_het_genotype_count_arg=None,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 2, "Number of pass variants is not as expected")

    def test_nonmajor_2_required_het(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="nonmajor",
                              samples_required_het_arg="SAMPLE_01,SAMPLE_02",
                              samples_required_hom_arg=None,
                              min_het_genotype_count_arg=None,
                              max_het_genotype_count_arg=None,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 3, "Number of pass variants is not as expected")

    def test_nonmajor_2_required_het_max_count_het_3(self):

        vcf_file = NamedTemporaryFile(mode="w")
        vcf_file.write(test_vcf_content)
        vcf_file.flush()

        count_variant_pass = apply_filter_criteria(vcf_file,
                              allele_type="nonmajor",
                              samples_required_het_arg="SAMPLE_01,SAMPLE_02",
                              samples_required_hom_arg=None,
                              min_het_genotype_count_arg=None,
                              max_het_genotype_count_arg=3,
                              min_hom_genotype_count_arg=None,
                              max_hom_genotype_count_arg=None)

        self.assertTrue(count_variant_pass == 2, "Number of pass variants is not as expected")




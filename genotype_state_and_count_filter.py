from cyvcf2 import VCF
import argparse
import sys
import numpy as np
import subprocess


def create_allele_count_dictionary(variant):
    """Create a dictionary with the allele count per allele"""

    # Dictionary that stores the count per allele
    allele_count_dictionary = {}

    # Determine the reference allele count
    total_number_of_alleles = variant.INFO["AN"]
    total_number_of_alternative_alleles = sum(wrap_in_list(variant.INFO["AC"]))
    reference_allele_count = total_number_of_alleles - total_number_of_alternative_alleles
    allele_count_dictionary[variant.REF] = reference_allele_count

    # add counts of the alternative alleles to the allele count dictionary
    # loop over the alternative alleles, lookup the allele count from the AC attribute and add dictionary entry
    alternative_allele_index = 0
    for alternative_allele in variant.ALT:
        allele_count_dictionary[alternative_allele] = wrap_in_list(variant.INFO["AC"])[alternative_allele_index]
        alternative_allele_index += 1

    return allele_count_dictionary


def determine_minor_allele(allele_count_dictionary):
    """Determine the minor (least frequent) allele.
    If there are multiple alleles with the same lowest count
    then one of these alleles is chosen at random based on position in the dictionary loop.
    There is no way to chose 1 minor allele over the other"""

    minor_allele = None

    # Loop over all alleles and their counts
    for allele, allele_count in allele_count_dictionary.items():

        # Set the minor allele if:
        # A) minor allele is none
        # B) allele count of current allele is smaller than the allele count for previous determined minor allele
        if not minor_allele or allele_count < allele_count_dictionary[minor_allele]:
            minor_allele = allele

    return minor_allele


def determine_major_allele(allele_count_dictionary):
    """Determine the major (most frequent) allele.
    If there are multiple alleles with the same highest count
    then one of these alleles is chosen at random based position in the dictionary loop.
    There is no way to chose 1 major allele over the other"""

    major_allele = None

    # Loop over all alleles and their counts
    for allele, allele_count in allele_count_dictionary.items():

        # Set the major allele if:
        # A) major allele is none
        # B) allele count of current allele is higher than the allele count for previous determined major allele
        if not major_allele or allele_count > allele_count_dictionary[major_allele]:
            major_allele = allele

    return major_allele


def determine_non_major_alleles(allele_count_dictionary):
    """Determine the non-major (least frequent) allele(s).
    All alleles that are not the major allele are returned"""

    major_allele = determine_major_allele(allele_count_dictionary)

    return set(allele_count_dictionary.keys()) - set(major_allele)


def determine_alleles_of_interest(variant, allele_type):
    """Determine the allele(s) of interest based on the allele type provided in the CLI"""

    alleles_of_interest = None

    # Create a dictionary with the allele count for all alleles
    allele_count_dictionary = create_allele_count_dictionary(variant)

    # If the allele type is non-reference, the alleles of interest are all alternative (non-reference) alleles
    if allele_type == "nref":
        alleles_of_interest = variant.ALT
    # If the allele type alt1, the alleles of interest is the 1st alternative allele
    if allele_type == "alt1":
        alleles_of_interest = variant.ALT[0]
    # If the allele type minor, the alleles of interest is least frequent allele
    if allele_type == "minor":
        alleles_of_interest = list(determine_minor_allele(allele_count_dictionary))
    # If the allele type major, the alleles of interest is most frequent allele
    if allele_type == "major":
        alleles_of_interest = list(determine_major_allele(allele_count_dictionary))
    # If the allele type non-major, the alleles of interest are all alleles except the major allele
    if allele_type == "nonmajor":
        alleles_of_interest = list(determine_non_major_alleles(allele_count_dictionary))

    return alleles_of_interest


def samples_required_pass(variant, samples_required_boolean_index, alleles_of_interest, het_or_hom):
    """Return True if the samples that are required to be heterozygous or homozygous are indeed heterozygous or homozygous.
    If any of the required samples is missing return False.

    The het_or_hom arguments specifies if this function check that the required samples are heterozygous or homozygous.
    The alleles_of_interest specifies the allele for which the required samples should be heterozygous or homozygous.

    If there are multiple alleles of interest (for multi-allelic variants and allele type non-reference or non-major)
    the required samples can be heterozygous or homozygous for any of those alleles for this function to return True."""

    # Number of unique alleles expected
    if het_or_hom == "het":
        expected_nr_unique_alleles = 2
    else:
        expected_nr_unique_alleles = 1

    # Get the genotypes of the required samples by applying boolean index
    genotypes_samples_of_interest = variant.gt_bases[samples_required_boolean_index]

    # Return False if a genotype is missing
    if any(is_missing_genotype(gt) for gt in genotypes_samples_of_interest):
        return False

    # Count how many of the required samples have the allele of interest in the desired het or hom state
    number_of_genotypes_match_desired_state = 0

    # Loop over the genotypes that should be het or hom with an allele of interest
    for genotype in genotypes_samples_of_interest:

        # Check that the genotype
        # A) matches the desired heterozygous or homozygous genotype state
        # (by checking the length of the set of unique alleles of the genotype (het=2, hom=1))
        # B) the genotype allele set contains at least 1 of the alleles of interest.
        # Alleles of interest can be nref, or nonmajor, which can be multiple alleles. (e.g. 1/2 genotype)

        # Convert the genotype strings ("A/T") to genotype allele sets {A,T}
        genotype_set = genotype_to_allele_set(genotype)
        if len(genotype_set) == expected_nr_unique_alleles and any(allele_of_interest in genotype_set for allele_of_interest in alleles_of_interest):
            number_of_genotypes_match_desired_state += 1

    # If all the required samples had the desired het or hom genotype return True
    if number_of_genotypes_match_desired_state == genotypes_samples_of_interest.size:
        return True
    else:
        return False


def count_genotypes_with_allele_of_interest(allele_sets, alleles_of_interest):
    """ Count the number of heterozygous and homozygous genotypes that contain at least 1 allele of interest

    allele_sets argument is a list with the set of alleles for alle samples
    alleles_of_interest is the alleles of interest for the allele type specified in the CLI"""

    count_het_genotype_with_allele_of_interest = 0
    count_hom_genotype_with_allele_of_interest = 0

    # Loop over the allele sets of all samples
    for allele_set in allele_sets:

        # If there is only 1 unique allele, the genotype is homozygous
        if len(allele_set) == 1:
            # if any allele of interest is in the allele set increase the hom genotype counter
            if any(allele_of_interest in allele_set for allele_of_interest in alleles_of_interest):
                count_hom_genotype_with_allele_of_interest += 1
        # Otherwise there are multiple distinct alleles and the genotype is heterozygous
        else:
            # if any allele of interest is in the allele set increase the het genotype counter
            if any(allele_of_interest in allele_set for allele_of_interest in alleles_of_interest):
                count_het_genotype_with_allele_of_interest += 1

    # Return a tuple with the number of genotypes that is heterozygous and homozygous for the alleles of interest
    return count_het_genotype_with_allele_of_interest, count_hom_genotype_with_allele_of_interest


def is_missing_genotype(genotype_string):
    """Return True is a genotype is missing"""
    if genotype_string == "./." or genotype_string == ".|." or genotype_string == ".":
        return True
    else:
        return False


def genotype_to_allele_set(genotype_string):
    """Convert a genotype string to a genotype set
    I.e. 'A/T' to {'A','T'}"""

    if '/' in genotype_string:
        separate_char = "/"
    else:
        separate_char = "|"

    return set(genotype_string.split(separate_char))


def parse_arguments(args):

    parser = argparse.ArgumentParser(
        description="Filter a BCF/VCF input stream on genotype count and state.")
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), required=False, dest="input_file",
                        help="Input BCF/VCF file")
    parser.add_argument('--allele_type', action='store', required=False, dest="allele_type", choices=["nref", "alt1", "minor", "major", "nonmajor" ], default="nref",
                        help="The type of allele on which the MIN and MAX genotype count and required samples list should be applied.")
    parser.add_argument('--samples_required_het', action='store', type=str, required=False, dest='samples_required_het',
                        help='Comma-separated list of samples that are required to be heterozygous for the allele type of interest')
    parser.add_argument('--samples_required_hom', action='store', type=str, required=False, dest='samples_required_hom',
                        help='Comma-separated list of samples that are required to be homozygous for the allele type of interest')
    parser.add_argument('--min_het_genotype_count', action='store', type=int, required=False, dest='min_het_genotype_count',
                        help='Minimum heterozygous genotype count for the allele type of interest')
    parser.add_argument('--max_het_genotype_count', action='store', type=int, required=False, dest='max_het_genotype_count',
                        help='Maximum heterozygous genotype count for the allele type of interest')
    parser.add_argument('--min_hom_genotype_count', action='store', type=int, required=False, dest='min_hom_genotype_count',
                        help='Minimum homozygous genotype count for the allele type of interest')
    parser.add_argument('--max_hom_genotype_count', action='store', type=int, required=False, dest='max_hom_genotype_count',
                        help='Maximum homozygous genotype count for the allele type of interest')

    return parser.parse_args(args)


def print_variant(variant):
    """Output a variant to standard out in the VCF format"""
    print(str(variant), end="")


def apply_filter_criteria(input_vcf,
                          allele_type,
                          samples_required_het_arg,
                          samples_required_hom_arg,
                          min_het_genotype_count_arg,
                          max_het_genotype_count_arg,
                          min_hom_genotype_count_arg,
                          max_hom_genotype_count_arg):
    """Apply the filter criteria proved in the CLI to VCF input file or VCF input stream.
    Only output variants that pass all the filter criteria to standard out formatted as VCF"""

    # Create a reader for reading from file or stdin
    if input_vcf:
        reader = VCF(input_vcf.name)
    else:
        reader = VCF('-')

    # update and print header
    command_line_invocation = subprocess.list2cmdline(sys.argv)
    command_line_invocation = command_line_invocation.replace('"', '')
    new_header_entry = {"ID": "GenotypeStateAndCount", "Description": command_line_invocation}
    reader.add_filter_to_header(new_header_entry)
    print(reader.raw_header, end="")

    # Create boolean index for the required het and hom samples
    # This boolean index is used later to efficiently
    # retrieve the genotypes of required samples from all the variant genotypes
    samples_required_het_boolean_index = None
    samples_required_hom_boolean_index = None

    # If the samples that are required to be heterozygous have been specified in the CLI
    # Create the boolean index for the heterozygous samples
    if samples_required_het_arg:
        samples_required_het = samples_required_het_arg.split(",")
        samples_required_het_boolean_index = np.array([sample in samples_required_het for sample in reader.samples])
    # If the samples that are required to be homozygous have been specified in the CLI
    # Create the boolean index for the homozygous samples
    if samples_required_hom_arg:
        samples_required_hom = samples_required_hom_arg.split(",")
        samples_required_hom_boolean_index = np.array([sample in samples_required_hom for sample in reader.samples])

    # Numpy vectorize function to convert genotype strings to allele sets
    # Is applied later on to "map" all genotype strings to genotype allele sets
    vfunc=np.vectorize(genotype_to_allele_set)

    # Count variants that pass the filter criteria and are outputted, mainly for unit testing
    count_variants_pass = 0

    # Loop over all variants
    for variant in reader:

        # Determine the alleles of intererest for the allele type specified in the CLI
        alleles_of_interest = determine_alleles_of_interest(variant, allele_type)

        # If samples are required to be het and they aren't skip this variant
        if samples_required_het_arg and not samples_required_pass(variant, samples_required_het_boolean_index, alleles_of_interest, "het"):
            continue
        # If samples are required to be hom and they aren't skip this variant
        if samples_required_hom_arg and not samples_required_pass(variant, samples_required_hom_boolean_index, alleles_of_interest, "hom"):
            continue

        # If any MIN or MAX het/hom count is not none
        # Thus a there at least 1 filter criteria genotype count
        if min_het_genotype_count_arg is not None \
                or max_het_genotype_count_arg is not None \
                or min_hom_genotype_count_arg is not None\
                or max_hom_genotype_count_arg is not None:

            # Convert genotype string array to genotype allele sets array
            allele_sets = vfunc(variant.gt_bases)

            # Count the number of genotypes that are heterozygous or homozygous for the alleles of interest
            count_het_genotype_with_allele_of_interest, count_hom_genotype_with_allele_of_interest = count_genotypes_with_allele_of_interest(allele_sets,
                                                                                                                                             alleles_of_interest)

            # Skip variant if count of het genotypes is lower than the minimum
            if min_het_genotype_count_arg is not None and count_het_genotype_with_allele_of_interest < min_het_genotype_count_arg:
                continue
            # Skip variant if count of het genotypes is higher than the maximum
            if max_het_genotype_count_arg is not None and count_het_genotype_with_allele_of_interest > max_het_genotype_count_arg:
                continue
            # Skip variant if count of hom genotypes is lower than the minimum
            if min_hom_genotype_count_arg is not None and count_hom_genotype_with_allele_of_interest < min_hom_genotype_count_arg:
                continue
            # Skip variant if count of hom genotypes is higher than the maximum
            if max_hom_genotype_count_arg is not None and count_hom_genotype_with_allele_of_interest > max_hom_genotype_count_arg:
                continue

        print_variant(variant)

        # Increment the number of variants that passed the filter criteria and are outputted
        count_variants_pass += 1

    return count_variants_pass


def main():

    args = parse_arguments(sys.argv[1:])

    apply_filter_criteria(args.input_file,
                          args.allele_type,
                          args.samples_required_het,
                          args.samples_required_hom,
                          args.min_het_genotype_count,
                          args.max_het_genotype_count,
                          args.min_hom_genotype_count,
                          args.max_hom_genotype_count)


def wrap_in_list(object):
    """If the object is not a list or tuple wrap the object (most often a string or integer) in a list"""

    if isinstance(object, list) or isinstance(object, tuple):
        return object
    else:
        return [object]

if __name__ == "__main__":
    main()



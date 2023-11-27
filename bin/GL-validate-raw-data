#!/usr/bin/env python

"""
This is a program for validating GeneLab raw data
"""

import os
import sys
import argparse
import textwrap
import pandas as pd
import glob
import re
import subprocess
import shutil

parser = argparse.ArgumentParser(description = "This program validates GeneLab raw datasets, renames files if needed \
                                                to follow GL convention, checks and/or creates md5s, and \
                                                runs fastqc and combines its outputs with multiqc. It is meant to be \
                                                executed in the directory holding the read files. \
                                                For version info, run `GL-version`.",
                                 epilog = "Ex. usage: GL-validate-raw-data -g GLDS-480 -n 3")

required = parser.add_argument_group('required arguments')

required.add_argument("-g", "--GLDS-ID", help = 'GLDS ID (e.g. "GLDS-480")', action = "store", required = True)
required.add_argument("-n", "--number-of-samples", help = "The expected number of samples here", action = "store", type = int, required = True)

parser.add_argument("-m", "--md5-file", help = "Include a file holding all md5sum results if avaiable (one will be created if not, or if the filenames are changed)", \
                    action = "store", default = "")
parser.add_argument("-s", "--single-ended", help = "Add this flag if data are single-end sequencing", action = "store_true")
parser.add_argument("--ATACseq", help = "Add this flag if the data are ATACseq (in which case R3 is expected to be the reverse read, and R2 holds barcodes)", action = "store_true")
parser.add_argument("-k", "--keep-fastqc-files", help = "Add this flag if wanting to keep individual-sample fastqc files", action = "store_true")
parser.add_argument("--slurm", help = "Add this flag to submit the job through slurm", action = "store_true")
parser.add_argument("--nodename", help = "Add this flag to specify slurm node to use", action = "store", default = "")
parser.add_argument("--mem", help = "Add this flag to specify the amount of memory to use", action = "store", default = "")
parser.add_argument("--exclude_nodename", help = "Add this flag to specify slurm nodes to exclude. This overrides any --nodename supplied.", action = "store", default = "")
parser.add_argument("--queue", help = "Add this flag to specify slurm queue to use", action = "store", default = "normal")
parser.add_argument('--HRremoved-suffix', action='store_true', 
                    help='Set this flag to indicate that "HRremoved" suffix is expected.')

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

# currently hard-coded things we may want to adjust
R1_designations = ["_R1_", "_R1.", "-R1.", "-R1-", ".R1.", "_1."]
R2_designations = ["_R2_", "_R2.", "-R2.", "-R2-", ".R2.", "_2."]
extensions = [".fq", ".fastq"]
standard_GL_R1_suffix = "_R1_raw.fastq.gz" if not args.HRremoved_suffix else "_R1_HRremoved_raw.fastq.gz"
standard_GL_R2_suffix = "_R2_raw.fastq.gz" if not args.HRremoved_suffix else "_R2_HRremoved_raw.fastq.gz"
standard_GL_SE_suffix = "_raw.fastq.gz" if not args.HRremoved_suffix else "_HRremoved_raw.fastq.gz"

# additional ones for ATACseq
R3_designations = ["_R3_", "_R3.", "-R3.", "-R3-", ".R3.", "_3."]
standard_GL_R3_suffix = "_R3_raw.fastq.gz" # ADD THIS IF ATACseq also needs to allow `HRremoved` in suffix at user request: if not args.HRremoved_suffix else "_R3_HRremoved_raw.fastq.gz"
print(f"Standard GL Suffixes: R1: '{standard_GL_R1_suffix}'; R2: '{standard_GL_R2_suffix}'; R_SE: '{standard_GL_SE_suffix}'; R3 (ATACseq) '{standard_GL_R3_suffix}'")

md5_output_file = "raw_md5sum.txt"
md5_check_output_file = str(args.GLDS_ID) + "-raw-md5-check-results.txt"
fastqc_threads = 8 # allocates 250 MB memory per thread, so 8 is 2 GB
multiqc_output_prefix = "raw_multiqc"
multiqc_data_dir = multiqc_output_prefix + "_data"
multiqc_html = multiqc_output_prefix + ".html"
# multiqc_html = multiqc_output_prefix + "_report.html"
multiqc_dir = multiqc_output_prefix + "_report"

versions_used_file = str(args.GLDS_ID) + "-raw-validation-tool-versions.txt"

log_out = str(args.GLDS_ID) + "-raw-fastqc-multiqc-log.txt"

output_table_filename = str(args.GLDS_ID) + "-raw-validation-summary.tsv"

gzip_or_fastq_format_check_failed_message = "Failed gzip test and/or fastq-format check."

if args.slurm:

    sbatch_file = str(args.GLDS_ID) + "-raw-validation.slurm"
    slurm_out_file = str(args.GLDS_ID) + "-raw-validation-slurm.out"
    slurm_job_name = str(args.GLDS_ID) + "-raw-validation"
    slurm_queue = str(args.queue)

    # This section makes memory allocated via the --mem option override the default of 20GB
    if args.mem:
        slurm_mem = str(args.mem)
    else:
        slurm_mem = 20000

    # This sections makes `exclusions` override any supplied nodelist!
    if args.exclude_nodename:
        slurm_node_exclusions = str(args.exclude_nodename)
        slurm_node = False
    else:
        slurm_node = str(args.nodename)
        slurm_node_exclusions = False


################################################################################

def main():

    # checking that single-end and ATACseq aren't specified together, as they are mutually exclusive
    if args.single_ended and args.ATACseq:

        print_notification("Arguments '--single-ended' and '--ATACseq' can't both be specified. Exiting for now.")
        sys.exit(1)

    if args.slurm:

        # pre-flight check so things aren't passed to slurm before finding a problem
        # so folks will be notified immediately by this script rather than only in the slurm output log
        map_dict = check_files_and_gen_unique_ID_to_file_map(args.number_of_samples)

        # making command and submitting as slurm job if --slurm was specified
        submit_slurm()

    else:

        # creating mapping dictionary of files to unique sample IDs
        map_dict = check_files_and_gen_unique_ID_to_file_map(args.number_of_samples)

        # starting validation output table
        out_tab = start_validation_output_table(map_dict)

        # if an md5 file was provided, checking md5s before changing their names
        if args.md5_file != "":

            md5_df = check_md5s(map_dict, out_tab)

        # renaming files to be in our convention (e.g., _R1_raw.fastq.gz for forward, or "_raw.fastq.gz" 
        # if single-end), also adding to output table
        map_dict, out_tab, files_renamed = rename_files(map_dict, out_tab)

        # checking gzip integrity
        map_dict, out_tab, problem_files_list = check_gzip_integrity(map_dict, out_tab)

        # check fastq format with fastq_utils
        map_dict, out_tab, problem_files_list = check_fastq_format(map_dict, out_tab, problem_files_list)

        # stop to ensure we don't move on if no files remain after gzip and fastq-format checks
        check_any_read_files_remain(map_dict, out_tab, problem_files_list)

        # generating md5s if none were provided
        if args.md5_file == "":

            map_dict, out_tab = gen_md5s(map_dict, out_tab, problem_files_list)

        else:
            # adding them to the output summary table now if they were provided
            out_tab = pd.concat([out_tab, md5_df], axis = 1)

            # writing out new md5s file with renamed filenames
            write_new_md5s_file(map_dict, md5_df, problem_files_list, files_renamed)

        # run fastqc and multiqc
        run_fastqc_and_multiqc(map_dict, problem_files_list)

        # add read counts to summary output table
        out_tab = parse_and_add_read_counts_and_lengths(map_dict, out_tab, problem_files_list)

        # writing output summary table
        out_tab.reset_index(inplace = True)
        out_tab = out_tab.rename(columns = {'index': 'unique_ID'})
        out_tab.to_csv(output_table_filename, sep = "\t", index = False)

        # # zipping up multiqc data dir
        # zip_multiqc_data_dir()

        # making directory for multiqc outputs and zipping
        package_and_zip_multiqc_outputs()

        gen_and_add_multiqc_zip_md5()
        # # making md5 for zipped qc dir and html and adding to md5s file
        # gen_and_add_multiqc_md5s()

        # writing out versions used and protocol text for Curation
        write_versions_used_file()

        # finish message
        finish_message(problem_files_list)


################################################################################

# setting some colors
tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}


### functions ###

def color_text(text, color='green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text):
    """ print wrapper """

    print(textwrap.fill(text, width=80, initial_indent="  ",
          subsequent_indent="  ", break_on_hyphens=False))


def report_failure(message, color = "red"):
    print("")
    wprint(color_text(message, color))
    print("\nRaw data V+V is exiting without completing.\n")
    sys.exit(1)


def print_notification(message, color = "yellow"):
    print("")
    wprint(color_text(message, color))
    print("")


def check_files_and_gen_unique_ID_to_file_map(expected_count):
    """ 
    checks there are the appropriate number of files based on expected number of samples, and
    returns a dictionary with unique IDs as keys and list of read files as values 
    """

    # getting all files in current directory that end with .gz
    gz_files = glob.glob("*.gz")
    # keeping only those that have .fq or .fastq in them
    gz_files = [file for file in gz_files if any(extension in file for extension in extensions)]

    num_files = len(gz_files)

    # checking the expected number of files were detected
    # for single-end
    if args.single_ended:

        if num_files != args.number_of_samples:
            report_failure("We expected " + str(args.number_of_samples) + " read files, but " + str(num_files) + " were detected.")

    # for ATACseq
    elif args.ATACseq:

        if num_files != args.number_of_samples * 3 and num_files != args.number_of_samples * 4:
            report_failure("We expected " + str(args.number_of_samples * 3) + " or " + str(args.number_of_samples * 4) + " read files, but " + str(num_files) + " were detected.")

    # for typical paired-end
    else:

        if num_files != args.number_of_samples * 2:
            report_failure("We expected " + str(args.number_of_samples * 2) + " read files, but " + str(num_files) + " were detected.")


    # making sure they have the expected R1/R2 designations if typical paired-end data
    if not args.single_ended and not args.ATACseq:

        forward_read_files_list = [file for file in gz_files if any(designation in file for designation in R1_designations)]
        reverse_read_files_list = [file for file in gz_files if any(designation in file for designation in R2_designations)]

        num_forward_read_files = len(forward_read_files_list)
        num_reverse_read_files = len(reverse_read_files_list)

        if args.number_of_samples != num_forward_read_files:

            report_failure("We expected " + str(args.number_of_samples) + " forward read files, but " + str(num_forward_read_files) + " were found based on currently allowed read-pair designations (e.g. these formats: " + ', '.join(R1_designations) + "). If it's not clear what the problem is here, share this use-case with Mike so we can expand functionality if needed.")

        if args.number_of_samples != num_reverse_read_files:

            report_failure("We expected " + str(args.number_of_samples) + " reverse read files, but " + str(num_forward_read_files) + " were found based on currently allowed read-pair designations (e.g. these formats: " + ', '.join(R2_designations) + "). If it's not clear what the problem is here, share this use-case with Mike so we can expand functionality if needed")


    # making a dictionary with unique IDs as keys and a list of read files as values if typical paired-end data

        # finding which read-pair designation we are using here (these are set in a list at top of script if wanting to add more)
        # putting a check that only one is found (one mechanism to hopefully prevent if these patterns are just part of a 
        # sample name, and not a read-pair designation)

        num_forward_designations_found = 0
        num_reverse_designations_found = 0

        for designation in R1_designations:

            if designation in forward_read_files_list[0]:
                R1_designation = designation
                num_forward_designations_found += 1

        for designation in R2_designations:

            if designation in reverse_read_files_list[0]:
                R2_designation = designation
                num_reverse_designations_found += 1

        if num_forward_designations_found != 1 or num_reverse_designations_found != 1:
            report_failure("We're having trouble delineating read-pairs based on the currently allowed formats (e.g. these formats: " + ', '.join(R1_designations) + "). Share this use-case with Mike to expand functionality.")


        # making dictionary with unique IDs as keys and list of read-files as values
        map_dict = {}
        for file in forward_read_files_list:

            # getting unique portion of names prior to _R?_ designation
            target_string_pattern = R1_designation + ".*$"
            unique_ID = re.sub(target_string_pattern, "", file)

            map_dict[unique_ID] = [file]

        for file in reverse_read_files_list:

            # getting unique portion of names prior to _R?_ designation
            target_string_pattern = R2_designation + ".*$"
            unique_ID = re.sub(target_string_pattern, "", file)

            map_dict[unique_ID].append(file)

        return(map_dict)


    # now the case if single-end data
    elif args.single_ended:

        # find out which fq/fastq extension is being used here so we can get unique names prior to that (or to be able to rename in our convention if needed)
        # making sure not more than one is detected
        num_extensions_found = 0
        for possible_extension in extensions:

            if possible_extension in gz_files[0]:
                used_extension = possible_extension
                num_extensions_found += 1

        if num_extensions_found != 1:
            report_failure("We're not seeing the expected fastq extensions (e.g. any of these: " + ', '.join(extensions) + "). Share this use-case with Mike to expand functionality.")

        # making dictionary with unique IDs as keys and read file as value
        map_dict = {}
        for file in gz_files:

            # getting unique portion of names prior to fasta extension
            target_string_pattern = used_extension + ".*$"
            unique_ID = re.sub(target_string_pattern, "", file)
            # here we are removing "_raw" from the end if it is present, so we don't inadvertently double it at the end of the unique portion of the sample ID
            unique_ID = re.sub("_raw$", "", unique_ID)

            map_dict[unique_ID] = [file]

        return(map_dict)

    # now the case if ATACseq
    elif args.ATACseq:

        # making sure they have the expected R1/R2/R3
        R1_read_files_list = [file for file in gz_files if any(designation in file for designation in R1_designations)]
        R2_read_files_list = [file for file in gz_files if any(designation in file for designation in R2_designations)]
        R3_read_files_list = [file for file in gz_files if any(designation in file for designation in R3_designations)]

        num_R1_read_files = len(R1_read_files_list)
        num_R2_read_files = len(R2_read_files_list)
        num_R3_read_files = len(R3_read_files_list)

        if args.number_of_samples != num_R1_read_files:

            report_failure("We expected " + str(args.number_of_samples) + " R1 read files, but " + str(num_R1_read_files) + " were found based on currently allowed read-pair designations (e.g. these formats: " + ', '.join(R1_designations) + "). If it's not clear what the problem is here, share this use-case with Mike so we can expand functionality if needed.")

        if args.number_of_samples != num_R2_read_files:

            report_failure("We expected " + str(args.number_of_samples) + " R2 read files, but " + str(num_R2_read_files) + " were found based on currently allowed read-pair designations (e.g. these formats: " + ', '.join(R2_designations) + "). If it's not clear what the problem is here, share this use-case with Mike so we can expand functionality if needed.")

        if args.number_of_samples != num_R3_read_files:

            report_failure("We expected " + str(args.number_of_samples) + " R3 read files, but " + str(num_R3_read_files) + " were found based on currently allowed read-pair designations (e.g. these formats: " + ', '.join(R3_designations) + "). If it's not clear what the problem is here, share this use-case with Mike so we can expand functionality if needed.")

        # making a dictionary with unique IDs as keys and a list of read files as values if typical paired-end data

        # finding which read-pair designation we are using here (these are set in a list at top of script if wanting to add more)
        # putting a check that only one is found (one mechanism to hopefully prevent if these patterns are just part of a 
        # sample name, and not a read-pair designation)

        num_R1_designations_found = 0
        num_R2_designations_found = 0
        num_R3_designations_found = 0

        for designation in R1_designations:

            if designation in R1_read_files_list[0]:
                R1_designation = designation
                num_R1_designations_found += 1

        for designation in R2_designations:

            if designation in R2_read_files_list[0]:
                R2_designation = designation
                num_R2_designations_found += 1

        for designation in R3_designations:

            if designation in R3_read_files_list[0]:
                R3_designation = designation
                num_R3_designations_found += 1

        if num_R1_designations_found != 1 or num_R2_designations_found != 1 or num_R3_designations_found != 1:
            report_failure("We're having trouble delineating R1, R2, and R3 read files based on the currently allowed formats (e.g., formatted like this: " + ', '.join(R1_designations) + "). Share this use-case with Mike to expand functionality if needed.")


        # making dictionary with unique IDs as keys and list of read-files as values
        map_dict = {}
        for file in R1_read_files_list:

            # getting unique portion of names prior to _R?_ designation
            target_string_pattern = R1_designation + ".*$"
            unique_ID = re.sub(target_string_pattern, "", file)

            map_dict[unique_ID] = [file]

        for file in R2_read_files_list:

            # getting unique portion of names prior to _R?_ designation
            target_string_pattern = R2_designation + ".*$"
            unique_ID = re.sub(target_string_pattern, "", file)

            map_dict[unique_ID].append(file)

        for file in R3_read_files_list:

            # getting unique portion of names prior to _R?_ designation
            target_string_pattern = R3_designation + ".*$"
            unique_ID = re.sub(target_string_pattern, "", file)

            map_dict[unique_ID].append(file)

        return(map_dict)


def start_validation_output_table(map_dict):
    """ starts output table and adds unique IDs and file names """

    # getting count of files for each entry and putting into dictionary to add to summary output table
    files_count_dict = {}

    for key in map_dict.keys():

        files_count_dict[key] = len(map_dict[key])

    # starting output summary table
    out_tab = pd.DataFrame.from_dict(files_count_dict, orient = "index", columns = ['num_read_files'])

    # adding info on original file names
    # for typical paired-end case
    if not args.single_ended and not args.ATACseq:

        new_tab = pd.DataFrame.from_dict(map_dict, orient = "index", columns = ['orig_R1_filename', 'orig_R2_filename'])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    # for single-end case
    elif args.single_ended:

        new_tab = pd.DataFrame.from_dict(map_dict, orient = "index", columns = ['orig_read_filename'])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    # for ATACseq case
    elif args.ATACseq:

        new_tab = pd.DataFrame.from_dict(map_dict, orient = "index", columns = ['orig_R1_filename', 'orig_R2_filename', 'orig_R3_filename'])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    return(out_tab)


def rename_files(map_dict, out_tab):
    """ 
    this renames the files to have the conventional GeneLab suffixes if needed
    e.g., _R1_raw.fastq.gz
    """

    # value storing if we do rename anything or not
    files_renamed = False

    # dictionary of current samples-to-files map
    current_map_dict = {}
    # dictionary of renamed files (or not), for output table summary
    renamed_map_dict = {}

    # case for typical paired-end
    if not args.single_ended and not args.ATACseq:

        for key in map_dict.keys():

            R1_orig_name = map_dict[key][0]
            R2_orig_name = map_dict[key][1]

            # only renaming if they don't fit our convention
            if not R1_orig_name.endswith(standard_GL_R1_suffix):

                files_renamed = True

                R1_new_name = str(key) + str(standard_GL_R1_suffix)
                os.rename(R1_orig_name, R1_new_name)
                current_map_dict[key] = [R1_new_name]
                renamed_map_dict[key] = [R1_new_name]

            else:
                current_map_dict[key] = [R1_orig_name]
                renamed_map_dict[key] = ["not-renamed"]

            if not R2_orig_name.endswith(standard_GL_R2_suffix):

                files_renamed = True

                R2_new_name = str(key) + str(standard_GL_R2_suffix)
                os.rename(R2_orig_name, R2_new_name)
                current_map_dict[key].append(R2_new_name)
                renamed_map_dict[key].append(R2_new_name)

            else:
                current_map_dict[key].append(R2_orig_name)
                renamed_map_dict[key].append("not-renamed")

        # adding new names to output table
        new_tab = pd.DataFrame.from_dict(renamed_map_dict, orient = "index", columns = ["new_R1_filename", "new_R2_filename"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)


    # case for ATACseq
    elif args.ATACseq:

        for key in map_dict.keys():

            R1_orig_name = map_dict[key][0]
            R2_orig_name = map_dict[key][1]
            R3_orig_name = map_dict[key][2]

            # only renaming if they don't fit our convention
            if not R1_orig_name.endswith(standard_GL_R1_suffix):

                files_renamed = True

                R1_new_name = str(key) + str(standard_GL_R1_suffix)
                os.rename(R1_orig_name, R1_new_name)
                current_map_dict[key] = [R1_new_name]
                renamed_map_dict[key] = [R1_new_name]

            else:
                current_map_dict[key] = [R1_orig_name]
                renamed_map_dict[key] = ["not-renamed"]

            if not R2_orig_name.endswith(standard_GL_R2_suffix):

                files_renamed = True

                R2_new_name = str(key) + str(standard_GL_R2_suffix)
                os.rename(R2_orig_name, R2_new_name)
                current_map_dict[key].append(R2_new_name)
                renamed_map_dict[key].append(R2_new_name)

            else:
                current_map_dict[key].append(R2_orig_name)
                renamed_map_dict[key].append("not-renamed")

            if not R3_orig_name.endswith(standard_GL_R3_suffix):

                files_renamed = True

                R3_new_name = str(key) + str(standard_GL_R3_suffix)
                os.rename(R3_orig_name, R3_new_name)
                current_map_dict[key].append(R3_new_name)
                renamed_map_dict[key].append(R3_new_name)


            else:
                current_map_dict[key].append(R3_orig_name)
                renamed_map_dict[key].append("not-renamed")

        # adding new names to output table
        new_tab = pd.DataFrame.from_dict(renamed_map_dict, orient = "index", columns = ["new_R1_filename", "new_R2_filename", "new_R3_filename"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    # case for single-end
    elif args.single_ended:

        for key in map_dict.keys():

            orig_name = map_dict[key][0]

            # only renaming if doesn't fit our convention
            if not orig_name.endswith(standard_GL_SE_suffix):
    
                new_name = str(key) + str(standard_GL_SE_suffix)
                os.rename(orig_name, new_name)
                current_map_dict[key] = [new_name]
                renamed_map_dict[key] = [new_name]

            else:
                current_map_dict[key] = [orig_name]
                renamed_map_dict[key] = ["not-renamed"]


        # adding new names to output table
        new_tab = pd.DataFrame.from_dict(renamed_map_dict, orient = "index", columns = ["new_filename"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    return(current_map_dict, out_tab, files_renamed)


def check_gzip_integrity(map_dict, out_tab):
    """ runs gzip -t <file> on each to check and report integrity """

    print_notification("Checking gzip integrity...")
    # making dictionary holding gzip test results (PASS/FAIL)
    gzip_test_dict = {}

    # making list of problematic files
    problem_files_list = []

    # case for typical paired-end
    if not args.single_ended and not args.ATACseq:

        for key in map_dict.keys():

            forward_read_file = map_dict[key][0]
            reverse_read_file = map_dict[key][1]

            forward_gzip_test_out = subprocess.run(['gzip', '-t', forward_read_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            reverse_gzip_test_out = subprocess.run(['gzip', '-t', reverse_read_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

            if forward_gzip_test_out.returncode != 0:
                forward_gzip_test_result = "FAIL"
                problem_files_list.append(forward_read_file)
            else:
                forward_gzip_test_result = "PASS"

            if reverse_gzip_test_out.returncode != 0:
                reverse_gzip_test_result = "FAIL"
                problem_files_list.append(reverse_read_file)
            else:
                reverse_gzip_test_result = "PASS"

            gzip_test_dict[key] = [forward_gzip_test_result]
            gzip_test_dict[key].append(reverse_gzip_test_result)

        ## adding results to output summary table
        new_tab = pd.DataFrame.from_dict(gzip_test_dict, orient = "index", columns = ["R1_gzip_test", "R2_gzip_test"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    # case if ATACseq
    elif args.ATACseq:

        for key in map_dict.keys():

            R1_read_file = map_dict[key][0]
            R2_read_file = map_dict[key][1]
            R3_read_file = map_dict[key][2]

            R1_gzip_test_out = subprocess.run(['gzip', '-t', R1_read_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            R2_gzip_test_out = subprocess.run(['gzip', '-t', R2_read_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            R3_gzip_test_out = subprocess.run(['gzip', '-t', R3_read_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

            if R1_gzip_test_out.returncode != 0:
                R1_gzip_test_result = "FAIL"
                problem_files_list.append(R1_read_file)
            else:
                R1_gzip_test_result = "PASS"

            if R2_gzip_test_out.returncode != 0:
                R2_gzip_test_result = "FAIL"
                problem_files_list.append(R2_read_file)
            else:
                R2_gzip_test_result = "PASS"

            if R3_gzip_test_out.returncode != 0:
                R3_gzip_test_result = "FAIL"
                problem_files_list.append(R3_read_file)
            else:
                R3_gzip_test_result = "PASS"

            gzip_test_dict[key] = [R1_gzip_test_result]
            gzip_test_dict[key].append(R2_gzip_test_result)
            gzip_test_dict[key].append(R3_gzip_test_result)

        ## adding results to output summary table
        new_tab = pd.DataFrame.from_dict(gzip_test_dict, orient = "index", columns = ["R1_gzip_test", "R2_gzip_test", "R3_gzip_test"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    # case if single-end
    elif args.single_ended:

        for key in map_dict.keys():

            read_file = map_dict[key][0]

            gzip_test_out = subprocess.run(['gzip', '-t', read_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

            if gzip_test_out.returncode != 0:
                gzip_test_result = "FAIL"
                problem_files_list.append(read_file)
            else:
                gzip_test_result = "PASS"

            gzip_test_dict[key] = [gzip_test_result]

        # adding results to output summary table
        new_tab = pd.DataFrame.from_dict(gzip_test_dict, orient = "index", columns = ["read_file_gzip_test"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)


    return(map_dict, out_tab, problem_files_list)

def run_fastq_check(first_file, second_file = "", problem_files_list = []):

    """
    this runs the fastq_info subprocess and returns needed outputs, used in check_fastq_format function
    """

    if second_file == "":
        
        check_out = subprocess.run(['fastq_info', first_file], capture_output = True, text = True)

    else:

        check_out = subprocess.run(['fastq_info', first_file, second_file], capture_output = True, text = True)

    if check_out.returncode != 0:

        # adding to problem files list if not already present and if there was only one file provided here
        if second_file == "" and first_file not in problem_files_list:
            problem_files_list.append(first_file)

        # getting error line to put in output table (it stops after a single error, so will always only report one)
        for line in iter(check_out.stderr.splitlines()):

            if line.startswith("ERROR"):

                error = re.sub(r".*line", "line", line)

        # catch in case one fails, but does not have an "ERROR" line like searched above
        if 'error' not in locals():

            error = "Fastq check had non-zero exit, but no error reported. Needs a closer look."

        # making report column hold the error message too
        check_result = f"FAIL - {error}"

    else:

        check_result = "PASS"

    return(check_result, problem_files_list)


def check_fastq_format(map_dict, out_tab, problem_files_list):

    """ 
    checks fastq format with fastq_utils 
        - if paired, it checks them individually, and then together (to give a little more info in some cases then if just done on one)
        - if ATACseq, does all 3, then the paired are R1 and R3    
    """

    print_notification("Checking fastq format...")
    # making dictionary holding fastq format validation results
    fastq_validation_dict = {}

    # case for typical paired-end
    if not args.single_ended and not args.ATACseq:

        for key in map_dict.keys():

            forward_read_file = map_dict[key][0]
            reverse_read_file = map_dict[key][1]

            # forward check
            fwd_fastq_check_result, problem_files_list = run_fastq_check(forward_read_file, problem_files_list = problem_files_list)

            # reverse check
            rev_fastq_check_result, problem_files_list = run_fastq_check(reverse_read_file, problem_files_list = problem_files_list)

            # paired check
            paired_fastq_check_result, problem_files_list = run_fastq_check(forward_read_file, reverse_read_file, problem_files_list = problem_files_list)
            
            # adding result to dictionary
            fastq_validation_dict[key] = [fwd_fastq_check_result, rev_fastq_check_result, paired_fastq_check_result]

        ## adding results to output summary table
        new_tab = pd.DataFrame.from_dict(fastq_validation_dict, orient = "index", columns = ["R1_fastq_format_check", "R2_fastq_format_check", "paired_fastq_format_check"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    # case if ATACseq
    elif args.ATACseq:

        for key in map_dict.keys():

            R1_read_file = map_dict[key][0]
            R2_read_file = map_dict[key][1]
            R3_read_file = map_dict[key][2]

            # R1 check
            R1_fastq_check_result, problem_files_list = run_fastq_check(R1_read_file, problem_files_list = problem_files_list)

            # R2 check
            R2_fastq_check_result, problem_files_list = run_fastq_check(R2_read_file, problem_files_list = problem_files_list)

            # R3 check
            R3_fastq_check_result, problem_files_list = run_fastq_check(R3_read_file, problem_files_list = problem_files_list)

            # paired R1-R3 check (because ATACseq)
            paired_R1_R3_fastq_check_result, problem_files_list = run_fastq_check(R1_read_file, R3_read_file, problem_files_list = problem_files_list)

            # adding result to dictionary
            fastq_validation_dict[key] = [R1_fastq_check_result, R2_fastq_check_result, R3_fastq_check_result, paired_R1_R3_fastq_check_result]

        ## adding results to output summary table
        new_tab = pd.DataFrame.from_dict(fastq_validation_dict, orient = "index", columns = ["R1_fastq_format_check", "R2_fastq_format_check", "R3_fastq_format_check", "paired_R1_R3_fastq_format_check"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)


    # case if single-end
    elif args.single_ended:

        for key in map_dict.keys():

            read_file = map_dict[key][0]

            # running check
            fastq_check_result, problem_files_list = run_fastq_check(read_file, problem_files_list = problem_files_list)
            
            # adding result to dictionary
            fastq_validation_dict[key] = [fastq_check_result]

        ## adding results to output summary table
        new_tab = pd.DataFrame.from_dict(fastq_validation_dict, orient = "index", columns = ["fastq_format_check"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    return(map_dict, out_tab, problem_files_list)


def check_md5s(map_dict, out_tab):
    """ checks (if provided) or generates md5sums for all read files """

    print_notification("Checking provided md5s...")

    md5_in_df = pd.read_csv(args.md5_file, sep = "\s+", names = ["md5", "filename"])
    md5_in_df.set_index('filename', inplace = True)

    md5_check_command = "md5sum -c " + str(args.md5_file) + " > " + str(md5_check_output_file)

    md5_check_out = subprocess.run([md5_check_command], shell = True, stderr = subprocess.PIPE)

    # reading in md5sum check results
    md5_check_results_df = pd.read_csv(md5_check_output_file, sep = ": ", names = ["filename", "status"], engine = "python")
    md5_check_results_df.set_index('filename', inplace = True)

    # getting md5 values into dictionary with our map_dict keys
    md5_dict = {}

    # case for typical paired-end
    if not args.single_ended and not args.ATACseq:

        for key in map_dict.keys():

            forward_read_file = map_dict[key][0]
            reverse_read_file = map_dict[key][1]

            forward_file_result = md5_check_results_df.loc[forward_read_file].status
            reverse_file_result = md5_check_results_df.loc[reverse_read_file].status

            forward_file_md5 = md5_in_df.loc[forward_read_file].md5
            reverse_file_md5 = md5_in_df.loc[reverse_read_file].md5

            if forward_file_result == "FAILED":
                md5_dict[key] = ["FAIL"]
            else:
                md5_dict[key] = ["PASS"]

            if reverse_file_result == "FAILED":
                md5_dict[key].append("FAIL")
            else:
                md5_dict[key].append("PASS")

            if forward_file_result == "FAILED":
                md5_dict[key].append("failed-md5-check")
            else:
                md5_dict[key].append(forward_file_md5)

            if reverse_file_result == "FAILED":
                md5_dict[key].append("failed-md5-check")
            else:
                md5_dict[key].append(reverse_file_md5)

        # putting into table to be combined with output summary table later
        md5_df = pd.DataFrame.from_dict(md5_dict, orient = "index", columns = ["R1_md5_check", "R2_md5_check", "R1_md5", "R2_md5"])

    # case for ATACseq
    elif args.ATACseq:

        for key in map_dict.keys():

            R1_read_file = map_dict[key][0]
            R2_read_file = map_dict[key][1]
            R3_read_file = map_dict[key][1]

            R1_file_result = md5_check_results_df.loc[R1_read_file].status
            R2_file_result = md5_check_results_df.loc[R2_read_file].status
            R3_file_result = md5_check_results_df.loc[R3_read_file].status

            R1_file_md5 = md5_in_df.loc[R1_read_file].md5
            R2_file_md5 = md5_in_df.loc[R2_read_file].md5
            R3_file_md5 = md5_in_df.loc[R3_read_file].md5

            if R1_file_result == "FAILED":
                md5_dict[key] = ["FAIL"]
            else:
                md5_dict[key] = ["PASS"]

            if R2_file_result == "FAILED":
                md5_dict[key].append("FAIL")
            else:
                md5_dict[key].append("PASS")

            if R3_file_result == "FAILED":
                md5_dict[key].append("FAIL")
            else:
                md5_dict[key].append("PASS")

            if R1_file_result == "FAILED":
                md5_dict[key].append("failed-md5-check")
            else:
                md5_dict[key].append(R1_file_md5)

            if R2_file_result == "FAILED":
                md5_dict[key].append("failed-md5-check")
            else:
                md5_dict[key].append(R2_file_md5)

            if R3_file_result == "FAILED":
                md5_dict[key].append("failed-md5-check")
            else:
                md5_dict[key].append(R3_file_md5)

        # putting into table to be combined with output summary table later
        md5_df = pd.DataFrame.from_dict(md5_dict, orient = "index", columns = ["R1_md5_check", "R2_md5_check", "R3_md5_check", "R1_md5", "R2_md5", "R3_md5"])

    # case for single-ended
    elif args.single_ended:

        for key in map_dict.keys():

            read_file = map_dict[key][0]
            file_result = md5_check_results_df.loc[read_file].status
            file_md5 = md5_in_df.loc[read_file].md5

            if file_result == "FAILED":
                md5_dict[key] = ["FAIL", "failed-md5-check"]

            else:
                md5_dict[key] = ["PASS", file_md5]

        # putting into table to be combined with output summary table later
        md5_df = pd.DataFrame.from_dict(md5_dict, orient = "index", columns = ["md5_check", "md5"])

    # removing md5_check_output_file
    os.remove(md5_check_output_file)

    return(md5_df)


def gen_md5s(map_dict, out_tab, problem_files_list):
    """ generates md5 checksum file for all read files """

    print_notification("Generating md5s...")
    # making dictionary holding md5 values
    md5_dict = {}

    # case for typical paired-end
    if not args.single_ended and not args.ATACseq:

        for key in map_dict.keys():

            forward_read_file = map_dict[key][0]
            reverse_read_file = map_dict[key][1]

            # only running if file passed the gzip test
            if forward_read_file not in problem_files_list:

                forward_md5_out = subprocess.run(['md5sum', forward_read_file], stdout = subprocess.PIPE).stdout.decode('utf-8')
    
                # getting rid of filename (these come like this '3s9ckfjie  file.txt' from command-line md5sum program)
                forward_md5_out = re.sub(" .*$", "", forward_md5_out).strip()
                md5_dict[key] = [forward_md5_out]

            else:
                md5_dict[key] = [gzip_or_fastq_format_check_failed_message]


            if reverse_read_file not in problem_files_list:

                reverse_md5_out = subprocess.run(['md5sum', reverse_read_file], stdout = subprocess.PIPE).stdout.decode('utf-8')
                reverse_md5_out = re.sub(" .*$", "", reverse_md5_out).strip()
                md5_dict[key].append(reverse_md5_out)

            else:

                md5_dict[key].append(gzip_or_fastq_format_check_failed_message)


        # adding md5 values to output table
        new_tab = pd.DataFrame.from_dict(md5_dict, orient = "index", columns = ["R1_md5", "R2_md5"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

        # writing out md5 file (only including those that are not problematic)
        with open(md5_output_file, "w") as md5_out_file:
            
            for key in map_dict.keys():
                
                forward_read_file = map_dict[key][0]
                reverse_read_file = map_dict[key][1]

                if forward_read_file not in problem_files_list:

                    forward_file_md5 = md5_dict[key][0]
                    md5_out_file.write(str(forward_file_md5) + "  " + str(forward_read_file) + "\n")

                if reverse_read_file not in problem_files_list:

                    reverse_file_md5 = md5_dict[key][1]
                    md5_out_file.write(str(reverse_file_md5) + "  " + str(reverse_read_file) + "\n")


    # case for ATACseq
    elif args.ATACseq:

        for key in map_dict.keys():

            R1_read_file = map_dict[key][0]
            R2_read_file = map_dict[key][1]
            R3_read_file = map_dict[key][2]

            # only running if file passed the gzip test
            if R1_read_file not in problem_files_list:

                R1_md5_out = subprocess.run(['md5sum', R1_read_file], stdout = subprocess.PIPE).stdout.decode('utf-8')
    
                # getting rid of filename (these come like this '3s9ckfjie  file.txt' from command-line md5sum program)
                R1_md5_out = re.sub(" .*$", "", R1_md5_out).strip()
                md5_dict[key] = [R1_md5_out]

            else:
                md5_dict[key] = [gzip_or_fastq_format_check_failed_message]


            if R2_read_file not in problem_files_list:

                R2_md5_out = subprocess.run(['md5sum', R2_read_file], stdout = subprocess.PIPE).stdout.decode('utf-8')
                R2_md5_out = re.sub(" .*$", "", R2_md5_out).strip()
                md5_dict[key].append(R2_md5_out)

            else:

                md5_dict[key].append(gzip_or_fastq_format_check_failed_message)


            if R3_read_file not in problem_files_list:

                R3_md5_out = subprocess.run(['md5sum', R3_read_file], stdout = subprocess.PIPE).stdout.decode('utf-8')
                R3_md5_out = re.sub(" .*$", "", R3_md5_out).strip()
                md5_dict[key].append(R3_md5_out)

            else:

                md5_dict[key].append(gzip_or_fastq_format_check_failed_message)


        # adding md5 values to output table
        new_tab = pd.DataFrame.from_dict(md5_dict, orient = "index", columns = ["R1_md5", "R2_md5", "R3_md5"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

        # writing out md5 file (only including those that are not problematic)
        with open(md5_output_file, "w") as md5_out_file:
            
            for key in map_dict.keys():
                
                R1_read_file = map_dict[key][0]
                R2_read_file = map_dict[key][1]
                R3_read_file = map_dict[key][2]

                if R1_read_file not in problem_files_list:

                    R1_file_md5 = md5_dict[key][0]
                    md5_out_file.write(str(R1_file_md5) + "  " + str(R1_read_file) + "\n")

                if R2_read_file not in problem_files_list:

                    R2_file_md5 = md5_dict[key][1]
                    md5_out_file.write(str(R2_file_md5) + "  " + str(R2_read_file) + "\n")

                if R3_read_file not in problem_files_list:

                    R3_file_md5 = md5_dict[key][2]
                    md5_out_file.write(str(R3_file_md5) + "  " + str(R3_read_file) + "\n")


    # case for single-end
    elif args.single_ended:

        for key in map_dict.keys():

            read_file = map_dict[key][0]

            # only running if file passed the gzip test
            if read_file not in problem_files_list:

                md5_out = subprocess.run(['md5sum', read_file], stdout = subprocess.PIPE).stdout.decode('utf-8')

                # getting rid of filename (these come like this '3s9ckfjie  file.txt' from command-line md5sum program)
                md5_out = re.sub(" .*$", "", md5_out).strip()
                md5_dict[key] = [md5_out]

            else:

                md5_dict[key] = [gzip_or_fastq_format_check_failed_message]                

        # adding md5 values to output table
        new_tab = pd.DataFrame.from_dict(md5_dict, orient = "index", columns = ["read_file_md5"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

        # writing out md5 file
        with open(md5_output_file, "w") as md5_out_file:
            
            for key in map_dict.keys():
                
                read_file = map_dict[key][0]

                if read_file not in problem_files_list:

                    file_md5 = md5_dict[key][0]
                    md5_out_file.write(str(file_md5) + "  " + str(read_file) + "\n")


    return(map_dict, out_tab)


def write_new_md5s_file(map_dict, md5_df, problem_files_list, files_renamed):
    """ writes out a new md5 file with renamed files """

    # just making a copy of the original if we didn't rename things
    if not files_renamed:
        shutil.copy(args.md5_file, md5_output_file)

    md5_dict = {}

    # case if paired-end
    if not args.single_ended and not args.ATACseq:

        with open(md5_output_file, "w") as md5_out_file:

            for key in map_dict.keys():

                forward_read_file = map_dict[key][0]
                reverse_read_file = map_dict[key][1]

                # only writing out if they passed the gzip test
                if forward_read_file not in problem_files_list:

                    forward_file_md5 = md5_df.loc[key].R1_md5
                    md5_out_file.write(str(forward_file_md5) + "  " + str(forward_read_file) + "\n")

                if reverse_read_file not in problem_files_list:

                    reverse_file_md5 = md5_df.loc[key].R2_md5
                    md5_out_file.write(str(reverse_file_md5) + "  " + str(reverse_read_file) + "\n")

    # case if ATACseq
    elif args.ATACseq:

        with open(md5_output_file, "w") as md5_out_file:

            for key in map_dict.keys():

                R1_read_file = map_dict[key][0]
                R2_read_file = map_dict[key][1]
                R3_read_file = map_dict[key][2]

                # only writing out if they passed the gzip test
                if R1_read_file not in problem_files_list:

                    R1_file_md5 = md5_df.loc[key].R1_md5
                    md5_out_file.write(str(R1_file_md5) + "  " + str(R1_read_file) + "\n")

                if R2_read_file not in problem_files_list:

                    R2_file_md5 = md5_df.loc[key].R2_md5
                    md5_out_file.write(str(R2_file_md5) + "  " + str(R2_read_file) + "\n")

                if R3_read_file not in problem_files_list:

                    R3_file_md5 = md5_df.loc[key].R3_md5
                    md5_out_file.write(str(R3_file_md5) + "  " + str(R3_read_file) + "\n")

    # case if single-end
    elif args.single_ended:

        with open(md5_output_file, "w") as md5_out_file:

            for key in map_dict.keys():

                read_file = map_dict[key][0]

                # only writing out if they passed the gzip test
                if read_file not in problem_files_list:

                    file_md5 = md5_df.loc[key].md5
                    md5_out_file.write(str(file_md5) + "  " + str(read_file) + "\n")


def check_any_read_files_remain(map_dict, out_tab, problem_files_list):

    """
    checks if any files remain after gzip and fastq-format tests,
    if not, writes output summary table with info we have so far and reports message to user
    """

    list_of_remaining_files = []

    for key in map_dict.keys():

        # only adding if passed the gzip test and fastq-format check (and therefore filename isn't in the problem_files_list)
        for file in map_dict[key]:
            if file not in problem_files_list:
                list_of_remaining_files += [file]

    if len(list_of_remaining_files) == 0:

        print_notification(f"NOTICE: All input files failed the gzip test and/or fastq-format check.", "red")

        print("    You can find what was the cause for each in the generated summary table.\n\n")

        # writing output summary table
        out_tab.reset_index(inplace = True)
        out_tab = out_tab.rename(columns = {'index': 'unique_ID'})
        out_tab.to_csv(output_table_filename, sep = "\t", index = False)

        exit()


def write_versions_used_file():
    
    """ writes the versions used and protocol text for Curation """

    # getting fastqc version
    fastqc_version = subprocess.run(['fastqc', '--version'], stdout = subprocess.PIPE).stdout.decode('utf-8')
    # getting multiqc version
    multiqc_version = subprocess.run(['multiqc', '--version'], stdout = subprocess.PIPE).stdout.decode('utf-8').replace(", version ", " v")
    # not as easy to get the fastq_utils version...
    fastq_info_out = subprocess.run(['fastq_info', '-h'], capture_output = True, text = True)
    for line in iter(fastq_info_out.stderr.splitlines()):

        if line.startswith("fastq_utils"):

            fastq_utils_version = line

    # writing out versions used to file
    with open(versions_used_file, "w") as out_file:

        out_file.write(str(fastq_utils_version) + str(fastqc_version) + str(multiqc_version))

        out_file.write("\nProtocol text:\n\n")
        out_file.write("Fastq format was checked with fastq_utils v" + str(fastq_utils_version).replace("fastq_utils ", "").replace("\n", "") + ". " + \
                       "Quality assessment of reads was performed with FastQC " + str(fastqc_version).replace("FastQC ", "").replace("\n", "") + \
                       " and reports were combined with MultiQC " + str(multiqc_version).replace("multiqc ", "").replace("\n", "") + ".\n\n")


def run_fastqc_and_multiqc(map_dict, problem_files_list):
    """ runs fastqc on all read files and writes out versions to file """

    print_notification("Running fastqc...")

    # set JAVA options to allow fastqc to use up to 8GB of heap memory
    os.environ["_JAVA_OPTIONS"] = "-Xmx8g"

    # making list of all files
    list_of_all_files = []

    for key in map_dict.keys():

        # only adding if passed the gzip test and fastq-format check (and therefore filename isn't in the problem_files_list)
        for file in map_dict[key]:
            if file not in problem_files_list:
                list_of_all_files += [file]

    # converting list to space-delimited string for passing to subprocess
    list_for_subprocess = " ".join(list_of_all_files)

    # getting fastqc version
    fastqc_version = subprocess.run(['fastqc', '--version'], stdout = subprocess.PIPE).stdout.decode('utf-8')

    # building command since i can't figure out how to pass a list or space-delimited string of positional
    # arguments to subprocess.run
    with open(log_out, "w") as log:
        log.write("\n\n    " + str(fastqc_version).strip() + " log:\n\n")

    fastqc_command = "fastqc -t " + str(fastqc_threads) + " " + str(list_for_subprocess) + " >> " + str(log_out) + " 2>&1"

    subprocess.run([fastqc_command], shell = True)

    # getting multiqc version
    multiqc_version = subprocess.run(['multiqc', '--version'], stdout = subprocess.PIPE).stdout.decode('utf-8').replace(", version ", " v")

    # running multiqc
    with open(log_out, "a") as log:
        log.write("\n\n    " + str(multiqc_version).strip() + " log:\n\n")

    multiqc_command = "multiqc ./ --force --cl-config 'max_table_rows: 99999999' --interactive --filename " + str(multiqc_output_prefix) + " >> " + str(log_out) + " 2>&1"

    subprocess.run([multiqc_command], shell = True)

    # removing fastqc files unless specified not to
    if not args.keep_fastqc_files:
        
        for file in glob.glob("*fastqc*"):

            if file != log_out:
                os.remove(file)


def parse_and_add_read_counts_and_lengths(map_dict, out_tab, problem_files_list):
    """ this gets the read counts from the multiqc output and adds them to the summary table """

    multiqc_general_stats_df = pd.read_csv(str(multiqc_output_prefix) + "_data/multiqc_general_stats.txt", sep = "\t", usecols = [0,5])
    multiqc_general_stats_df.columns = ["sample", "counts"]
    multiqc_general_stats_df.set_index("sample", inplace = True)

    multiqc_fastqc_df = pd.read_csv(str(multiqc_output_prefix) + "_data/multiqc_fastqc.txt", sep = "\t", usecols = [0,6,9])
    multiqc_fastqc_df.columns = ["sample", "length_range", "avg_length"]
    multiqc_fastqc_df.set_index("sample", inplace = True)

    # making dictionary of counts, length ranges, average lengths, and number of fastqc files present in the multiqc report
    counts_dict, length_range_dict, avg_length_dict, fastqc_files_count = {}, {}, {}, {}

    # case for typical paired-end
    if not args.single_ended and not args.ATACseq:

        for key in map_dict.keys():

            # only getting info if the file passed the gzip test
            if map_dict[key][0] not in problem_files_list:

                forward_read_multiqc_ID = re.sub(".fastq.gz", "", map_dict[key][0])

                num_forward_reads = round(multiqc_general_stats_df.loc[forward_read_multiqc_ID, 'counts'])
                forward_read_length_range = multiqc_fastqc_df.loc[forward_read_multiqc_ID, 'length_range']
                forward_avg_read_length = multiqc_fastqc_df.loc[forward_read_multiqc_ID, 'avg_length']

                counts_dict[key] = [num_forward_reads]
                length_range_dict[key] = [forward_read_length_range]
                avg_length_dict[key] = [forward_avg_read_length]

                fastqc_files_count[key] = 1

            else: 

                num_forward_reads = gzip_or_fastq_format_check_failed_message
                counts_dict[key] = [gzip_or_fastq_format_check_failed_message]
                length_range_dict[key] = [gzip_or_fastq_format_check_failed_message]
                avg_length_dict[key] = [gzip_or_fastq_format_check_failed_message]

                forward_avg_read_length = gzip_or_fastq_format_check_failed_message
                forward_read_length_range = gzip_or_fastq_format_check_failed_message

                fastqc_files_count[key] = 0


            # only getting info if the file passed the gzip test
            if map_dict[key][1] not in problem_files_list:

                reverse_read_multiqc_ID = re.sub(".fastq.gz", "", map_dict[key][1])

                num_reverse_reads = round(multiqc_general_stats_df.loc[reverse_read_multiqc_ID, 'counts'])
                reverse_read_length_range = multiqc_fastqc_df.loc[reverse_read_multiqc_ID, 'length_range']
                reverse_avg_read_length = multiqc_fastqc_df.loc[reverse_read_multiqc_ID, 'avg_length']

                counts_dict[key].append(num_reverse_reads)
                length_range_dict[key].append(reverse_read_length_range)
                avg_length_dict[key].append(reverse_avg_read_length)

                fastqc_files_count[key] += 1

            else:

                num_reverse_reads = gzip_or_fastq_format_check_failed_message
                counts_dict[key].append(gzip_or_fastq_format_check_failed_message)
                length_range_dict[key].append(gzip_or_fastq_format_check_failed_message)
                avg_length_dict[key].append(gzip_or_fastq_format_check_failed_message)

                reverse_avg_read_length = gzip_or_fastq_format_check_failed_message
                reverse_read_length_range = gzip_or_fastq_format_check_failed_message

                fastqc_files_count[key] += 0

            # adding value for if counts are equal
            if num_forward_reads == num_reverse_reads:

                counts_dict[key].append("YES")

            else:

                counts_dict[key].append("NO")

            # adding value for if read lengths are equal
            if forward_read_length_range == reverse_read_length_range and forward_avg_read_length == reverse_avg_read_length:

                avg_length_dict[key].append("YES")

            else:

                avg_length_dict[key].append("NO")


        # adding read counts to output table
        new_tab = pd.DataFrame.from_dict(counts_dict, orient = "index", columns = ["R1_num_reads", "R2_num_reads", "R1_and_R2_num_reads_equal"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

        # adding read length ranges to output table
        new_tab = pd.DataFrame.from_dict(length_range_dict, orient = "index", columns = ["R1_read_length_range", "R2_read_length_range"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)        

        # adding read avg lengths to output table
        new_tab = pd.DataFrame.from_dict(avg_length_dict, orient = "index", columns = ["R1_avg_read_length", "R2_avg_read_length", "R1_and_R2_read_lengths_equal"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

        # adding num fastqc files in multiqc report
        new_tab = pd.DataFrame.from_dict(fastqc_files_count, orient = "index", columns = ["num_fastqc_reports_in_multiqc_report"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    # case for typical paired-end
    elif args.ATACseq:

        for key in map_dict.keys():

            # doing R1, only getting info if the file passed the gzip test
            if map_dict[key][0] not in problem_files_list:

                R1_read_multiqc_ID = re.sub(".fastq.gz", "", map_dict[key][0])

                num_R1_reads = round(multiqc_general_stats_df.loc[R1_read_multiqc_ID, 'counts'])
                R1_read_length_range = multiqc_fastqc_df.loc[R1_read_multiqc_ID, 'length_range']
                R1_avg_read_length = multiqc_fastqc_df.loc[R1_read_multiqc_ID, 'avg_length']

                counts_dict[key] = [num_R1_reads]
                length_range_dict[key] = [R1_read_length_range]
                avg_length_dict[key] = [R1_avg_read_length]

                fastqc_files_count[key] = 1

            else: 

                num_R1_reads = gzip_or_fastq_format_check_failed_message
                counts_dict[key] = [gzip_or_fastq_format_check_failed_message]
                length_range_dict[key] = [gzip_or_fastq_format_check_failed_message]
                avg_length_dict[key] = [gzip_or_fastq_format_check_failed_message]

                R1_avg_read_length = gzip_or_fastq_format_check_failed_message
                R1_read_length_range = gzip_or_fastq_format_check_failed_message

                fastqc_files_count[key] = 0


            # doing R2, only getting info if the file passed the gzip test
            if map_dict[key][1] not in problem_files_list:

                R2_read_multiqc_ID = re.sub(".fastq.gz", "", map_dict[key][1])

                num_R2_reads = round(multiqc_general_stats_df.loc[R2_read_multiqc_ID, 'counts'])
                R2_read_length_range = multiqc_fastqc_df.loc[R2_read_multiqc_ID, 'length_range']
                R2_avg_read_length = multiqc_fastqc_df.loc[R2_read_multiqc_ID, 'avg_length']

                counts_dict[key].append(num_R2_reads)
                length_range_dict[key].append(R2_read_length_range)
                avg_length_dict[key].append(R2_avg_read_length)

                fastqc_files_count[key] += 1

            else:

                num_R2_reads = gzip_or_fastq_format_check_failed_message
                counts_dict[key].append(gzip_or_fastq_format_check_failed_message)
                length_range_dict[key].append(gzip_or_fastq_format_check_failed_message)
                avg_length_dict[key].append(gzip_or_fastq_format_check_failed_message)

                R2_avg_read_length = gzip_or_fastq_format_check_failed_message
                R2_read_length_range = gzip_or_fastq_format_check_failed_message

                fastqc_files_count[key] += 0

            # doing R3, only getting info if the file passed the gzip test
            if map_dict[key][2] not in problem_files_list:

                R3_read_multiqc_ID = re.sub(".fastq.gz", "", map_dict[key][2])

                num_R3_reads = round(multiqc_general_stats_df.loc[R3_read_multiqc_ID, 'counts'])
                R3_read_length_range = multiqc_fastqc_df.loc[R3_read_multiqc_ID, 'length_range']
                R3_avg_read_length = multiqc_fastqc_df.loc[R3_read_multiqc_ID, 'avg_length']

                counts_dict[key].append(num_R3_reads)
                length_range_dict[key].append(R3_read_length_range)
                avg_length_dict[key].append(R3_avg_read_length)

                fastqc_files_count[key] += 1

            else:

                num_R3_reads = gzip_or_fastq_format_check_failed_message
                counts_dict[key].append(gzip_or_fastq_format_check_failed_message)
                length_range_dict[key].append(gzip_or_fastq_format_check_failed_message)
                avg_length_dict[key].append(gzip_or_fastq_format_check_failed_message)

                R3_avg_read_length = gzip_or_fastq_format_check_failed_message
                R3_read_length_range = gzip_or_fastq_format_check_failed_message

                fastqc_files_count[key] += 0


            # adding value for if counts are equal in R1, R2, and R3 for ATACseq
            if num_R1_reads == num_R2_reads and num_R1_reads == num_R3_reads:

                counts_dict[key].append("YES")

            else:

                counts_dict[key].append("NO")

            # adding value for if read lengths are equal in R1 and R3 for ATACseq
            if R1_read_length_range == R3_read_length_range and R1_avg_read_length == R3_avg_read_length:

                avg_length_dict[key].append("YES")

            else:

                avg_length_dict[key].append("NO")


        # adding read counts to output table
        new_tab = pd.DataFrame.from_dict(counts_dict, orient = "index", columns = ["R1_num_reads", "R2_num_reads", "R3_num_reads", "R1_R2_and_R3_num_reads_equal"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

        # adding read length ranges to output table
        new_tab = pd.DataFrame.from_dict(length_range_dict, orient = "index", columns = ["R1_read_length_range", "R2_read_length_range", "R3_read_length_range"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

        # adding read avg lengths to output table
        new_tab = pd.DataFrame.from_dict(avg_length_dict, orient = "index", columns = ["R1_avg_read_length", "R2_avg_read_length", "R3_avg_read_length", "R1_and_R3_read_lengths_equal"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

        # adding num fastqc files in multiqc report
        new_tab = pd.DataFrame.from_dict(fastqc_files_count, orient = "index", columns = ["num_fastqc_reports_in_multiqc_report"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)


    # case if single-end
    elif args.single_ended:

        for key in map_dict.keys():

            # only working in if the file passed the gzip test
            if map_dict[key][0] not in problem_files_list:

                read_multiqc_ID = re.sub(".fastq.gz", "", map_dict[key][0])

                num_reads = round(multiqc_general_stats_df.loc[read_multiqc_ID, 'counts'])
                read_length_range = multiqc_fastqc_df.loc[read_multiqc_ID, 'length_range']
                avg_read_length = multiqc_fastqc_df.loc[read_multiqc_ID, 'avg_length']

                counts_dict[key] = [num_reads]
                length_range_dict[key] = [read_length_range]
                avg_length_dict[key] = [avg_read_length]

                fastqc_files_count[key] = 1

            else:

                counts_dict[key] = [gzip_or_fastq_format_check_failed_message]
                length_range_dict[key] = [gzip_or_fastq_format_check_failed_message]
                avg_length_dict[key] = [gzip_or_fastq_format_check_failed_message]

                fastqc_files_count[key] = 0


        # adding read counts to output table
        new_tab = pd.DataFrame.from_dict(counts_dict, orient = "index", columns = ["read_counts"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)


        # adding read length ranges to output table
        new_tab = pd.DataFrame.from_dict(length_range_dict, orient = "index", columns = ["read_length_range"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)        

        # adding read avg lengths to output table
        new_tab = pd.DataFrame.from_dict(avg_length_dict, orient = "index", columns = ["avg_read_length"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)        

        # adding num fastqc files in multiqc report
        new_tab = pd.DataFrame.from_dict(fastqc_files_count, orient = "index", columns = ["num_fastqc_reports_in_multiqc_report"])
        out_tab = pd.concat([out_tab, new_tab], axis = 1)

    return(out_tab)


# def zip_multiqc_data_dir():
#     """ zips up multiqc data dir """

#     # zipping multiqc dir
#     zip_command = "zip -q -m -r " + str(multiqc_data_dir) + ".zip " + str(multiqc_data_dir)

#     subprocess.run([zip_command], shell = True)

#     # renaming output html
#     shutil.move(str(multiqc_output_prefix) + ".html", multiqc_html)


def package_and_zip_multiqc_outputs():
    """ moves multiqc outputs into their own dir and zips it """

    try:
        os.mkdir(multiqc_dir)
    except:
        pass

    shutil.move(multiqc_data_dir, multiqc_dir)
    shutil.move(multiqc_html, multiqc_dir)

    # zipping multiqc dir
    zip_command = "zip -q -m -r " + str(multiqc_dir) + ".zip " + str(multiqc_dir)

    subprocess.run([zip_command], shell = True)


# def gen_and_add_multiqc_md5s():
#     """ generates and adds multiqc md5s to md5 file """

#     zip_file = str(multiqc_data_dir) + ".zip"
#     zip_md5_out = subprocess.run(['md5sum', zip_file], stdout = subprocess.PIPE).stdout.decode('utf-8')

#     html_md5_out = subprocess.run(['md5sum', multiqc_html], stdout = subprocess.PIPE).stdout.decode('utf-8')

#     with open(md5_output_file, "a") as out_file:

#         out_file.write(zip_md5_out)
#         out_file.write(html_md5_out)

def gen_and_add_multiqc_zip_md5():
    """ generates and adds multiqc zip md5 to md5 file """

    zip_file = str(multiqc_dir) + ".zip"
    zip_md5_out = subprocess.run(['md5sum', zip_file], stdout = subprocess.PIPE).stdout.decode('utf-8')

    with open(md5_output_file, "a") as out_file:

        out_file.write(zip_md5_out)


def finish_message(problem_files_list):
    """ reports when done, and notes if there were any files that did not pass the gzip test """

    print_notification("Program complete :)", "green")

    if len(problem_files_list) > 0:

        print_notification("NOTICE: " + str(len(problem_files_list)) + " file(s) did not pass the gzip integrity test and/or the fastq-format check.", "red")

        print("    These are therefore currently missing from the multiqc output as well as the")
        print("    'raw_md5sum.txt'. You can find which these are in the generated summary table.\n\n")


def submit_slurm():
    """ making sbatch file and submitting as slurm job if `--slurm` was specified """

    command_list = sys.argv
    
    # removing --slurm from list for command that will be placed in slurm batch script
    command_list.remove("--slurm")

    command_for_slurm_script = " ".join(command_list)

    # building and writing out slurm sbatch file
    with open(sbatch_file, "w") as out:

        out.write('#!/bin/bash\n\n')

        out.write('#SBATCH --job-name=' + str(slurm_job_name) + "\n")
        out.write('#SBATCH --output=' + str(slurm_out_file) + "\n")
        if slurm_node_exclusions:
            out.write('#SBATCH --exclude=' + str(slurm_node_exclusions) + "\n")
        if slurm_node:
            out.write('#SBATCH --nodelist=' + str(slurm_node) + "\n")
        out.write('#SBATCH --partition=' + str(slurm_queue) + "\n")
        out.write('#SBATCH --mem=' + str(slurm_mem) + "\n\n")
	
        out.write('# sourcing user profile\n')
        out.write('. ~/.profile\n\n')

        out.write('# enables regular conda activate robustly (had trouble with source activate)\n')
        out.write('eval "$(conda shell.bash hook)"\n\n')

        out.write('# tracking time\n')
        out.write('start=$(date +%s)\n\n')

        out.write('# printing out node that is running the job\n')
        out.write('echo -e "Job running on node: ${HOSTNAME}\\n"\n\n')
        
        out.write('# activating sra-tools conda environment\n')
        out.write('conda activate genelab-utils\n\n')

        out.write('call="' + str(command_for_slurm_script) + '"\n\n')

        out.write('# printing out command to log file\n')
        out.write('echo -e "${call}\\n"\n\n')

        out.write('# running the command\n')
        out.write('eval ${call}\n')
        out.write('echo ""\n\n')

        out.write('# storing end time and reporting time taken\n')
        out.write('end=$(date +%s)\n\n')

        out.write('runtime_s=$(echo $(( end - start )))\n')
        out.write('echo "total run time(s): $runtime_s"\n\n')

        out.write('runtime_m=$(echo "scale=2; $runtime_s / 60;" | bc)\n')
        out.write('runtime_h=$(echo "scale=2; $runtime_s / 3600;" | bc)\n\n')

        out.write('echo "total run time(m): $runtime_m"\n')
        out.write('echo -e "total run time(h): $runtime_h\\n"\n\n')

        out.write('echo "slurm job ID: ${SLURM_JOB_ID}"\n')

    # submitting the slurm job
    os.system("sbatch " + str(sbatch_file))

    sys.exit(0)

if __name__ == "__main__":
    main()

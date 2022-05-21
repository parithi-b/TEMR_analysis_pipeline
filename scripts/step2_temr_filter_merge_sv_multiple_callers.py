"""
Script Information

Title :
STEP 1 : Merging SV calls based on Ranking approach
STEP 2 : Removing calls present completely within Gaps and Centromere regions

Input :
BED FILES from the callers containing (DELs,DUPs and INVs), also SV_LENGTH >= 50bp

BED Files should contain 9 columns:
CHR START END SVTYPE CALLER SAMPLEID PR SR DHFFC DHBFC

Also make sure to mention the full paths for the input files and UCSC tracks (check the main function)

Output :
Single bed file containing the following information:
All input bedfiles should have 5 columns if not change the script accordingly
chr start   end svtype  SAMPLEID callers(seperated by ;)

Description :

Identifying SV Calls called by atleast 2 out of the 3 SV callers used.
Once the calls are identified they are merged into a single SV call using a Ranking approach to retain the coordinates
predicted by the Highest Ranked Caller.

Remove calls present within Gaps and Centromeres

Merging is done in two parts 80% RO for SVs under 50Kbp

Search RANKING to change the ranking of the callers

"""


########################################################################################################################
########################################################################################################################
########################################################################################################################

# General use functions


def check_folder_exist(root_folder, trio_id):
    """
    Checks if the folder exists if not create one to save all the files

    :param root_folder: path where the folder will exist
    :param trio_id: Individual ID

    :return: path to the 'Ensemble' folder
    """

    ensemble_folder = root_folder + "/Ensemble"

    if os.path.isdir(ensemble_folder) is False:
        os.makedirs(ensemble_folder)

    if os.path.isdir(ensemble_folder + "/" + trio_id) is False:
        os.makedirs(ensemble_folder + "/" + trio_id)

    return ensemble_folder


def check_if_bed_exists(root_folder, c1_bed, c2_bed, c3_bed):
    """
    Checks if all the bed files exists

    :param root_folder: Main Folder (everything will be saved within this folder)
    :param c1_bed:  path to the manta bed file
    :param c2_bed: path to the lumpy bed file
    :param c3: path to the delly bed file

    :return: True if all files exists else False
    """

    status_flag = True

    if os.path.isdir(root_folder) is False:
        print("Root folder location is incorrect !!!")
        status_flag = False
    if os.path.isfile(c1_bed) is False:
        print("Caller1 bed location is incorrect  or missing !!!")
        status_flag = False
    if os.path.isfile(c2_bed) is False:
        print("Caller2 bed location is incorrect  or missing !!!")
        status_flag = False
    if os.path.isfile(c3_bed) is False:
        print("Caller3 bed location is incorrect  or missing !!!")
        status_flag = False

    return status_flag


def total_number_of_lines_in_a_file(file_path):
    """
    prints the number of lines in the given file

    :param file_path:
    :return:
    """

    count = 0

    for i in open(file_path, 'r'):
        count += 1

    print(file_path.split("/")[-1])
    print(count, " SV Calls")


def splitting_output_files_based_on_svtype(output_file):
    """
    Given a BED file, will be split it based on the SV types

    :param output_file: full path to the BED file that needs to be split
    :return:
    """

    del_bed_write = open(output_file.split(".")[0] + "_DEL.bed", 'w')
    dup_bed_write = open(output_file.split(".")[0] + "_DUP.bed", 'w')
    inv_bed_write = open(output_file.split(".")[0] + "_INV.bed", 'w')

    cat_stage = subprocess.Popen(['cat', output_file], stdout=subprocess.PIPE)
    subprocess.call(['grep', 'DEL'], stdin=cat_stage.stdout, stdout=del_bed_write)
    cat_stage = subprocess.Popen(['cat', output_file], stdout=subprocess.PIPE)
    subprocess.call(['grep', 'DUP'], stdin=cat_stage.stdout, stdout=dup_bed_write)
    cat_stage = subprocess.Popen(['cat', output_file], stdout=subprocess.PIPE)
    subprocess.call(['grep', 'INV'], stdin=cat_stage.stdout, stdout=inv_bed_write)


########################################################################################################################
########################################################################################################################
########################################################################################################################

# Part 1 of the pipeline - Merging


def union_of_caller_files(root_folder, c1, c2, c3, trio_id, technology):
    """
    Create a single file with all the SV calls from each caller and then sort them

    :param root_folder: Folder where all files will be saved
    :param c1:  path to the manta bed file
    :param c2: path to the lumpy bed file
    :param c3: path to the delly bed file
    :param trio_id: Individual ID
    :param technology: short-read or long-read

    :return: path to the union file
    """

    print("Unionzing and filtering step....")
    print("Sample ID ", trio_id)

    support_info = trio_id.split("_")

    print(support_info)

    if technology == "long-read":
        re_support = support_info[1]
        rd_flag = support_info[-1]

        if support_info[1] == "default":
            re_support = '5'

        awk_command = '(' + '$7>=' + re_support + ')'

        if rd_flag == "RD":
            awk_command += '&& (($4=="DEL" && $8<0.7) || ($4=="DUP" && $9>1.3) || ($4=="INV"))'

    elif technology == "short-read":

        if len(support_info) == 3:
            pr_support = '0'
            sr_support = '0'
            rd_info = support_info[-1]
            awk_command = '(' + '$7>=' + pr_support + '|| $8>=' + sr_support + ')'
            if rd_info == "RD":
                awk_command += '&& (($4=="DEL" && $9<0.7) || ($4=="DUP" && $10>1.3) || ($4=="INV"))'
        else:
            pr_support = support_info[1]
            sr_support = support_info[2]
            rd_info = support_info[-1]

            awk_command = '(' + '$7>=' + pr_support + '|| $8>=' + sr_support + ')'
            if rd_info == "RD":
                awk_command += '&& (($4=="DEL" && $9<0.7) || ($4=="DUP" && $10>1.3) || ($4=="INV" && ($7>=20 || $8>=10)))'

    union_file = root_folder + "/" + trio_id + "_union_of_caller.tsv"

    print("FILTER query : \n \t", awk_command)

    write_file = open(union_file, "w")

    cat_process = subprocess.Popen(['cat', c1, c2, c3], stdout=subprocess.PIPE)

    awk_process = subprocess.Popen(['awk', '-F', '\t', awk_command], stdin=cat_process.stdout, stdout=subprocess.PIPE)

    subprocess.call(['sort', '-k1,1', '-k2,2n', '-'], stdin=awk_process.stdout, stdout=write_file)

    return union_file


def create_a_dictionary(input_bed_file):
    """
    convert the bedtools output into a dict()

    :param input_bed_file: bedtools output
    :return: two copies of the dict()
    """

    sv_dict = dict()

    count_cols = 0  # to find the number of columns
    for i in open(input_bed_file):
        if "chr" in i:
            count_cols = len(i.split('\t'))
        break

    print("No of columns in dict() :  ", count_cols)

    key_col_count = count_cols // 2

    for sv in open(input_bed_file, "r"):

        if len(sv) == 0:  # ignoring empty lines
            continue

        y = sv.split()

        ### Key

        key_chromosome = y[0]
        key_start_pos = y[1]
        key_end_pos = y[2]
        key_sv_type = y[3]
        key_caller_used = y[5]

        key_temp = "__".join(str(k) for k in y[:key_col_count])

        if key_temp not in sv_dict:
            sv_dict[key_temp] = []

        ### Value

        value_chromosome = y[0 + key_col_count]
        value_start_pos = y[1 + key_col_count]
        value_end_pos = y[2 + key_col_count]
        value_sv_type = y[3 + key_col_count]
        value_caller_used = y[5 + key_col_count]

        value_temp = "__".join(str(k) for k in y[key_col_count:])

        if value_sv_type == key_sv_type:  # same SV type
            if value_temp not in sv_dict[key_temp] and key_temp != value_temp:
                sv_dict[key_temp].append(value_temp)

    # remove calls that were RO with calls with different SV type

    for i in sv_dict.copy():
        if len(sv_dict[i]) == 0:
            del sv_dict[i]

    # print("NO of KEYs in Dictionary  :  ", len(sv_dict))

    return sv_dict


def remove_unique_sv(root_folder, bed_file, trio_id, RO_value):
    """
    use bedtools intersect and run a 80% RO between the same file, also use -c (count) parameter
    Remove calls that have c = 1 , these are the calls that were called by exactly 1 caller

    :param root_folder: Folder where all files will be saved
    :param bed_file: Union file with all SVs from 3 callers
    :param trio_id: Individual ID

    :return: path to the file which contains non unique calls
    """

    non_unique_file = root_folder + "/" + trio_id + "_non_unique_temporary.bed"

    write_file = open(non_unique_file, "w")

    bedtools_unique = subprocess.Popen(
        ['bedtools', 'intersect', '-f', RO_value, '-r', '-a', bed_file, '-b', bed_file, '-c'], stdout=subprocess.PIPE)

    # when using bedtools intersect keep file a and file b , same.
    #  Remove calls that have count of 1 , as they will be the unique calls

    a = 0  # to find the number of columns
    out, err = bedtools_unique.communicate()
    for i in out.decode().split("\n"):
        if "chr" in i:
            a = len(i.split('\t'))

    print("No of columns :  ", a)

    bedtools_unique = subprocess.Popen(
        ['bedtools', 'intersect', '-f', RO_value, '-r', '-a', bed_file, '-b', bed_file, '-c'], stdout=subprocess.PIPE)

    filter_query = "$" + str(a) + "!=1"

    awk_process = subprocess.Popen(['awk', '-F', '\t', filter_query], stdin=bedtools_unique.stdout,
                                   stdout=subprocess.PIPE)

    # if same caller different sample, this is the place to wrk

    # subprocess.call(['awk', '-v','OFS=\t', '{print $1,$2,$3,$4,$5,$6}'], stdin=awk_process.stdout, stdout=write_file)

    subprocess.call(['cut', '-f', '1,2,3,4,5,6'], stdin=awk_process.stdout, stdout=write_file)

    return non_unique_file


def RO_with_self(root_folder, bed_file, trio_id, RO_value):
    """
    Using bedtools to perform given % RO with the calls present within the same file and sort them

    :param root_folder: Folder where all files will be saved
    :param bed_file: Union file with all SVs from 3 callers
    :param trio_id: Individual ID

    :return: path to the file containing the output
    """

    self_RO_file = root_folder + "/" + trio_id + "_3_caller_intersect_temporary.bed"

    write_file = open(self_RO_file, "w")

    subprocess.call(['bedtools', 'intersect', '-f', RO_value, '-r', '-a', bed_file, '-b', bed_file, '-wa', '-wb'],
                    stdout=write_file)

    return self_RO_file


def selecting_calls_present_in_2_out_of_3_caller(bed_file, technology):
    """
    Select calls that have support from at least 2 callers (2 out of 4)

    Create a dictionary

    :param bed_file: path to file with bedtools intersect with self and '-wa -wb' parameters

    :return: final list of SVs that are called by at least 2 calls

    """

    ########################################################################
    ### RANKING
    ########################################################################

    if technology == "short-read":
        caller_rank = {"lumpy": 3, "delly": 2, "manta": 1}
    elif technology == "long-read":
        caller_rank = {"svim": 3, "sniffles": 2, "pbsv": 1}

    ########################################################################
    ### Creating the dictionary
    ########################################################################

    sv_database = create_a_dictionary(bed_file)

    ########################################################################
    ### Filtering
    ########################################################################

    for key_sv in sv_database.copy():

        key_string = key_sv

        y = key_string.split("__")

        key_caller_used = y[-1]

        caller_flag = 0  # this flag is used to remove index, were both index and key are called by same caller

        if key_sv in sv_database:

            for value_sv in sv_database.copy()[key_string]:

                value_string = value_sv

                y = value_string.split("__")

                value_caller_used = y[-1]

                # if the rank of the value calls lower than index calls
                # then remove the corresponding value calls from the database index

                if key_caller_used != value_caller_used:
                    if caller_rank[key_caller_used] < caller_rank[value_caller_used]:
                        if value_string in sv_database:
                            del sv_database[value_string]
                    caller_flag = 1
                else:
                    # if it is same caller remove the key calls from database index
                    if value_string in sv_database:
                        del sv_database[value_string]

            if caller_flag == 0:
                # if both key and value calls are called by same caller remove them from database index
                if key_string in sv_database:
                    del sv_database[key_string]

    ########################################################################
    ### FINAL REARRANGEMENTS
    ########################################################################
    # exit()

    sv_skipped = 0

    consolidated_call_set = []

    for i in sv_database:

        key_caller_used = i.split("__")[-1]

        sv_considered = i

        called_by_other_caller = []  # add name of the other callers that called this SVs

        for j in sv_database[i]:

            value_caller_used = j.split("__")[-1]

            if value_caller_used not in called_by_other_caller and caller_rank[key_caller_used] < caller_rank[
                value_caller_used]:
                # add lower ranked caller in this group
                # because higher rank would have been deleted because they would be RO with a different call in the database
                called_by_other_caller.append(value_caller_used)

        # sort the callers name based on the rank assume the index call is the highest ranked in the list
        called_by_other_caller = sorted(called_by_other_caller, key=caller_rank.get)
        # merge the callers name to the sv call's information
        sv_considered += ";" + ";".join(called_by_other_caller)

        if len(called_by_other_caller) == 0:
            sv_skipped += 1
            continue
        # write the ones that are called by at least 2 or more callers
        consolidated_call_set.append(sv_considered)

    # print("SV Skipping", sv_skipped)
    return consolidated_call_set


########################################################################################################################
########################################################################################################################
########################################################################################################################

## Function that calls part 1 and part 2

def break_and_run(root_folder, trio_id, bed_file, RO_percentage, technology):
    """

    Function that calls all the other function in an orderly manner based on the input values and obtains the
    merged output

    :param root_folder: main folder
    :param trio_id: ID of the individual
    :param bed_file: sorted union of all calls in a bed file
    :param RO_percentage: Reciprocal Overlap percentage

    :return: list of Merged SV

    """

    # number of SV calls after removing unique SVs
    non_unique_file_path = remove_unique_sv(root_folder, bed_file, trio_id, RO_percentage)

    # 2 out of 3 callers
    self_RO_file = RO_with_self(root_folder, non_unique_file_path, trio_id, RO_percentage)
    atleast_2_out_of_3 = selecting_calls_present_in_2_out_of_3_caller(self_RO_file, technology)

    print("Number of SVs called by atleast 2 out of 3 callers  : ", len(atleast_2_out_of_3))

    return atleast_2_out_of_3


def main():
    #####################################################################
    ## INPUT file paths
    #####################################################################

    print("#### Input commands ##"+"#" * 50)

    input_arg = sys.argv

    technology = input_arg[8].lower()

    if technology not in ["short-read" , "long-read"]:
        print("Check parameters, make sure to mention short-read or long-read")
        exit()
    else:
        if technology == "short-read":
            parameter_length = 13
        elif technology == "long-read":
            parameter_length = 12

    if len(input_arg) != parameter_length:
        print("Incorrect number of input parameter...")
        print("Please check the query, there should be 13 entries including the script and vcf file...")
        print("#" * 72, "\n")
        exit()

    input_files = input_arg[1:4]
    input_caller_order = input_arg[4:7]

    for i in range(len(input_caller_order)):
        input_caller_order[i] = input_caller_order[i].lower()

    sample_id = input_arg[7]

    if technology == "short-read":
        output_folder = os.path.abspath(input_arg[12])

        caller1_index = input_caller_order.index("manta")
        caller2_index = input_caller_order.index("delly")
        caller3_index = input_caller_order.index("lumpy")

        caller1_bed = input_files[caller1_index]
        caller2_bed = input_files[caller2_index]
        caller3_bed = input_files[caller3_index]

        paired_read_filter = int(input_arg[9])  # 10
        split_read_filter = int(input_arg[10])  # 5
        read_depth = input_arg[11]  # True

        if paired_read_filter == 0 and split_read_filter == 0:
            if read_depth is False:
                trio_id = sample_id + "_default_noRD"
            else:
                trio_id = sample_id + "_default_RD"
        else:
            if read_depth is False:
                trio_id = sample_id + "_" + str(paired_read_filter) + "_" + str(split_read_filter) + "_noRD"
            else:
                trio_id = sample_id + "_" + str(paired_read_filter) + "_" + str(split_read_filter) + "_RD"

    elif technology == "long-read":

        output_folder = os.path.abspath(input_arg[11])
        caller1_index = input_caller_order.index("pbsv")
        caller2_index = input_caller_order.index("sniffles")
        caller3_index = input_caller_order.index("svim")

        caller1_bed = input_files[caller1_index]
        caller2_bed = input_files[caller2_index]
        caller3_bed = input_files[caller3_index]

        read_filter = int(input_arg[9])  # 5
        read_depth = input_arg[10]  # True

        if read_filter <= 5:
            if read_depth is False:
                trio_id = sample_id + "_default_noRD"
            else:
                trio_id = sample_id + "_default_RD"
        else:
            if read_depth is False:
                trio_id = sample_id + "_" + str(read_filter) + "_noRD"
            else:
                trio_id = sample_id + "_" + str(read_filter) + "_RD"

    else:
        print("Please enter the correct parameter name....")
        print("acceptable parameter : short-read or long-read (not case sensitive).....")

    print("caller1 file : ",caller1_bed)
    print("caller2 file : ",caller2_bed)
    print("caller3 file : ",caller3_bed)
    print("output_folder : ", output_folder)
    if technology == "short-read":
        print(sample_id, technology, paired_read_filter, split_read_filter, read_depth)
    else:
        print(sample_id, technology, read_filter, read_depth)
    """
    se=5 means se>=5, so calls with se = [0,1,2,3,4] will be removed

    if no filter needs to be applied 
    mention the value as 0
    """

    data_root_folder = check_folder_exist(output_folder, sample_id)
    data_root_folder += "/" + sample_id

    print("Main folder \n", data_root_folder)

    if check_if_bed_exists(data_root_folder, caller1_bed, caller2_bed, caller3_bed) is False:
        print(" Please !! check the bed file locations")
        exit()

    ######################################################################
    ### Unionize all SVs from multiple caller
    ######################################################################

    union_file_path = union_of_caller_files(data_root_folder, caller1_bed, caller2_bed, caller3_bed, trio_id,
                                            technology)

    ######################################################################
    ### Mergeing SVs from multiple caller
    ######################################################################

    print("*" * 50)
    print("*" * 50)

    RO_value = str(0.798)
    calls_under = break_and_run(data_root_folder, trio_id, union_file_path, RO_value, technology)

    print("*" * 50)
    print("*" * 50)

    ######################################################################
    ### Merge both the above callset
    ######################################################################

    final_output_file = data_root_folder + "/" + trio_id + "_merged_temporary.bed"
    merged_file = open(final_output_file, "w")
    for i in calls_under:
        merged_file.write("\t".join(i.split("__")) + "\n")
    merged_file.close()

    sorted_output = final_output_file.split("temporary")[0] + "sorted.tsv"
    write_file = open(sorted_output, 'w')
    subprocess.call(['sort', '-k1,1', '-k2,2n', final_output_file], stdout=write_file)

    total_number_of_lines_in_a_file(sorted_output)

    # #####################################################################
    # ## Remove TEMP files used in this pipeline
    # #####################################################################
    print("\ntemporary files removed")
    for files in os.listdir(data_root_folder + "/"):
        if "temporary" in files:
            # print(data_root_folder + "/" + files)
            print(files)
            os.remove(data_root_folder + "/" + files)

    print("\nProcess completed .......")

if __name__ == "__main__":

    import subprocess, sys, os
    main()


def union_the_bedfile(file1, file2, file3, union_merge):
    """
    This function creates a files containing all the SV from all 3 ensemble files
    :param file1: ensemble file from sample 1
    :param file2: ensemble file from sample 2
    :param file3: ensemble file from sample 3
    :param union_merge: union of all the ensemble files
    :return:

    """

    write_file = open(union_merge, 'w')

    cat_command = subprocess.Popen(['cat', file1, file2, file3], stdout=subprocess.PIPE)

    grep_command = subprocess.Popen(['grep', '-v', 'chrY'], stdin=cat_command.stdout, stdout=subprocess.PIPE)

    subprocess.call(['sort', '-k1,1', '-k2,2n'], stdin=grep_command.stdout, stdout=write_file)

    print("Unionising the calls .. completed ....")


def create_a_dictionary(subprocess_command):
    """
        create a dictionary to analyze the data
    :param subprocess_command: output from bedtools
    :return: the dict
    """

    sv_dict = dict()

    for sv in subprocess_command.decode().split("\n"):
        if "chr" in sv:
            total_columns = len(sv.split('\t'))
            break

    print("Total columns in files  : ", total_columns)

    key_col_count = total_columns // 2

    for sv in subprocess_command.decode().split("\n"):

        y = sv.split()

        if len(y) < 2:  # ignoring empty lines
            continue

        ### Key
        key_sv_type = y[3]
        key_temp = "__".join(k for k in y[:key_col_count])

        if key_temp not in sv_dict:
            sv_dict[key_temp] = []

        ### Value
        value_sv_type = y[3 + key_col_count]
        value_temp = "__".join(k for k in y[key_col_count:])

        if value_sv_type == key_sv_type:  # same SV type
            if value_temp not in sv_dict[key_temp] and key_temp != value_temp:
                sv_dict[key_temp].append(value_temp)

    return sv_dict


def overlapping_individuals(input_bed_file):
    """
        function to overlap calls from multiple samples
    :param input_bed_file:
    :return:
    """

    bedtools_command = subprocess.Popen(
        ['bedtools', 'intersect', '-f', '0.798', '-r', '-a', input_bed_file, '-b', input_bed_file, '-wa', '-wb'],
        stdout=subprocess.PIPE)

    out, err = bedtools_command.communicate()

    combine_dict = create_a_dictionary(out)

    output_file = input_bed_file.split("union")[0] + "merged.tsv"
    write_output = open(output_file, 'w')

    for i in combine_dict.copy():
        if i in combine_dict:
            for j in combine_dict.copy()[i]:
                if j in combine_dict:
                    del combine_dict[j]

    for i in combine_dict:
        temp = []
        for j in combine_dict[i]:
            temp.append(j.split("__")[4])
            # print("\t", j)

        temp += [i.split("__")[4]]
        temp = ";".join(k for k in sorted(temp))
        if len(temp) != 0:
            y = i.split("__")

            y[1] = str(int(y[1]))
            y[2] = str(int(y[2]))

            y[4] = y[4] + ":" + temp

        write_output.write("\t".join(str(k) for k in y) + "\n")
        # print()

    write_output.close()
    print("NO of KEYs in Dictionary  :  ", len(combine_dict))


def make_dict(output_bedtools, key_col_count):

    sv_dict = dict()

    for sv in output_bedtools.decode().split("\n"):

        y = sv.split()

        if len(y) < 2:  # ignoring empty lines
            continue

        ### Key
        key_sv_type = y[3]
        key_temp = "__".join(k for k in y[:key_col_count])

        if key_temp not in sv_dict:
            sv_dict[key_temp] = []

        ### Value
        value_sv_type = y[3 + key_col_count]
        value_temp = "__".join(k for k in y[key_col_count:])

        if value_sv_type == key_sv_type:  # same SV type
            if value_temp not in sv_dict[key_temp] and key_temp != value_temp:
                sv_dict[key_temp].append(value_temp)

    return sv_dict


def main():

    #####################################################################
    ## INPUT file paths
    #####################################################################

    print("#### Input commands ##"+"#" * 50)
    input_arg = sys.argv

    if len(input_arg) != 5:
        print("Incorrect number of input parameter...")
        print("Please check the query, there should be 5 entries including the script and vcf file...")
        print("#" * 72, "\n")
        exit()

    input_file1 = input_arg[1]
    input_file2 = input_arg[2]
    input_file3 = input_arg[3]

    filter_flag = input_file3.split("/")[-1].split("_")
    print(filter_flag)
    print(len(filter_flag))

    if len(filter_flag) == 6:
        filter_flag = filter_flag[-5:-2]
    else:
        filter_flag = filter_flag[-4:-2]

    filter_flag = "_".join(filter_flag)

    print(filter_flag)

    output_folder = input_arg[4]

    output_path = os.path.abspath(output_folder) + "/All_samples_" + filter_flag + "_union.tsv"
    union_the_bedfile(input_file1, input_file2, input_file3, output_path)
    overlapping_individuals(output_path)
    os.remove(output_path)


if __name__ == "__main__":

    import subprocess, os, sys

    main()
"""
This script will merge long-read and short-read calls used in this study
"""


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


def illumina_and_pacbio_merge(illumina_bed, pacbio_bed, output_file):

    final_list = list()
    col_count = 0
    
    for i in open(illumina_bed,"r"):
        if "chr" in i:
            col_count = len(i.split("\t"))
        break

    bedtools_unique = subprocess.Popen(
        ['bedtools', 'intersect', '-f', '0.798', '-r', '-a', illumina_bed, '-b', pacbio_bed, '-loj', ],
        stdout=subprocess.PIPE)

    out, err = bedtools_unique.communicate()

    illumina_dict = make_dict(out, col_count)

    print(len(illumina_dict))
    count = 0

    for i in illumina_dict:
        temp_flag = "shared"
        if len(illumina_dict[i]) == 0:
            temp_flag = "illumina"
            count += 1
        final_list.append(i + "__" + temp_flag)

    print("shared : {}, unique : {} .... total : {}".format(len(illumina_dict) - count, count, len(illumina_dict)))

    # pacbio unique

    bedtools_unique = subprocess.Popen(
        ['bedtools', 'intersect', '-f', '0.798', '-r', '-a', pacbio_bed, '-b', illumina_bed, '-loj', ],
        stdout=subprocess.PIPE)

    out, err = bedtools_unique.communicate()

    pacbio_dict = make_dict(out, col_count)

    print(len(pacbio_dict))
    count = 0
    for i in pacbio_dict:
        if len(pacbio_dict[i]) == 0:
            final_list.append(i + "__pacbio")
            count += 1

    print("shared : {}, unique : {} .... total : {}".format(len(pacbio_dict) - count, count, len(pacbio_dict)))
    print("\n Final")
    print(len(final_list))

    output_write = open(output_file, "w")

    for i in final_list:
        output_write.write("\t".join(str(j) for j in i.split("__")) + "\n")

    output_write.close()


def main():

    #####################################################################
    ## INPUT file paths
    #####################################################################

    print("#### Input commands ##"+"#" * 50)
    input_arg = sys.argv

    if len(input_arg) != 4:
        print("Incorrect number of input parameter...")
        print("Please check the query, there should be 4 entries including the script and vcf file...")
        print("#" * 72, "\n")
        exit()

    input_file1 = input_arg[1]
    input_file2 = input_arg[2]

    output_folder = input_arg[3]
    output_file_path = os.path.abspath(output_folder) + "/All_samples_shortRead_longRead_merged.tsv"

    illumina_and_pacbio_merge(input_file1, input_file2, output_file_path)


if __name__ == "__main__":

    import subprocess, os, sys

    main()

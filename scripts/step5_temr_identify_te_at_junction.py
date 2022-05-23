

def make_dict(output_bedtools, col_count):

    output_file = dict()

    out, err = output_bedtools.communicate()

    for i in out.decode().split("\n"):

        sv_data = i.split()

        if len(sv_data) < 4:
            continue

        a = "__".join(str(j) for j in sv_data[:col_count])
        b = "__".join(str(j) for j in sv_data[col_count:])

        if a not in output_file:
            output_file[a] = []

        if b.split("__")[1] == "-1":
            output_file[a].append("NA")
        else:
            output_file[a].append(b)

    # print(len(output_file))

    return output_file


def repeat_overlapping(sv_data, repeats_within_SV):

    """

    For a given SV region , certain # of repeats would be captured around the breakpoint using our capture window
    (25 to 50bp), with those captured repeats calculate what % of the SV do these repeats contribute
    :param sv_data: information about a SV - chr start end
    :param repeats_within_SV: all repeats captured by our window criteria

    :return:
    return the % of overlapping regions between each repeat and SV
    """

    sv = sv_data

    start_sv = int(sv[1]) + mini_window
    stop_sv = int(sv[2]) - mini_window
    sv_len = stop_sv - start_sv

    percentage_dict = dict()

    for i in repeats_within_SV:

        te_data = i.split("__")

        start_repeat = int(te_data[1])
        stop_repeat = int(te_data[2])

        te_family = te_data[-2]
        bases_overlap = int(te_data[-1])

        # if start_repeat <= start_sv <= stop_repeat:
        #     bases_overlap -= 2
        #
        # if start_repeat <= stop_sv <= stop_repeat:
        #     bases_overlap -= 2

        if te_family not in percentage_dict:
            percentage_dict[te_family] = 0

        percentage_dict[te_family] += bases_overlap

    for i in percentage_dict:
        # print(i, percentage_dict[i], round(percentage_dict[i] / sv_len * 100, 2))
        percentage_dict[i] = round(percentage_dict[i] / sv_len * 100, 2)

    return percentage_dict


def check_TEM(sv_data, te_list):

    # print("Inside check TEM\n")
    # print("#"*50)

    sv = sv_data.split("__")
    # print("SV : ", "\t".join(sv))

    sv_start = int(sv[1]) + mini_window
    sv_end = int(sv[2]) - mini_window

    # print("\tTE in this region" , len(te_list))

    count_te = len(te_list)

    # print("SV START, end : ", sv_start, sv_end)

    flag = "Others"

    if count_te == 1:
        if te_list[0] == "NA":
            flag = "No_Repeat"
        else:
            te_data = te_list[0].split("__")
            start_te = int(te_data[1])
            end_te = int(te_data[2])

            flag_5 = ""
            flag_3 = ""

            if start_te <= sv_start <= end_te:
                flag_5 = "TE_SE_5"

            if start_te <= sv_end <= end_te:
                flag_3 = "TE_SE_3"

            if flag_5 == flag_3:
                flag = "No_Repeat"
            else:
                if flag_5 != "" and flag_3 != "":
                    flag = "No_Repeat"
                else:
                    if flag_3 == "":
                        flag = flag_5
                    else:
                        flag = flag_3
            # print("flag 5 : {}, flag 3 : {}, flag : {}".format(flag_5,flag_3, flag))
            # print(flag)

    else:
        overlap_percent = repeat_overlapping(sv, te_list)
        flag_5 = False
        flag_3 = False
        start_list = []
        end_list = []

        for j in te_list:
            # print("\t",j)
            te_data = j.split("__")

            start_te = int(te_data[1])
            end_te = int(te_data[2])

            if start_te <= sv_start <= end_te:
                flag_5 = True
                start_list.append(j)

            if start_te <= sv_end <= end_te:
                flag_3 = True
                end_list.append(j)

        if flag_3 is True and flag_5 is True:
            if len(start_list) != 1 or len(end_list) != 1:
                flag = "crowded_repeat"
            else:

                te_5 = start_list[0].split("__")
                te_3 = end_list[0].split("__")

                sign_a = te_5[-3]
                sign_b = te_3[-3]

                name_a = te_5[-2]
                name_b = te_3[-2]

                family_a = te_5[-4]
                family_b = te_3[-4]

                orientation = "OPP"

                if sign_a == sign_b:
                    orientation = "SAME"

                if name_a == name_b:
                    if family_a not in ["SINE", "LINE", "LTR", "DNA", "Retroposon"]:
                        flag = "Others"
                    else:
                        broken_flag = False
                        for i in overlap_percent:
                            if overlap_percent[i] > 80:
                                broken_flag = i + ";" + str(overlap_percent[i])
                        # print(overlap_percent, broken_flag, len(te_list))
                        # print(te_list)
                        if broken_flag != False:
                            if len(te_list) == 3:
                                flag = "Broken;" + broken_flag
                            elif len(te_list) < 3:
                                if float(broken_flag.split(";")[1]) < 97:
                                    # flag = "TEM_broken_by_2_" + orientation + ";" + start_list[0] + ";" + end_list[0]
                                    flag = "TEMR_" + orientation + ";" + start_list[0] + ";" + end_list[0]
                                else:
                                    flag = "Broken;" + broken_flag
                            else:
                                flag = "TEMR_" + orientation + ";" + start_list[0] + ";" + end_list[0]
                        else:
                            flag = "TEMR_" + orientation + ";" + start_list[0] + ";" + end_list[0]
                else:
                    flag = "nonTEM"

    return flag


def main():

    input_arg = sys.argv

    if len(input_arg) != 3:
        print("Incorrect number of input parameter...")
        print("Please check the query, there should be 3 entries including the script and vcf file...")
        print("Please unzip and merge thr provided repeatMasker track or download from UCSC table browser")
        print("#" * 72, "\n")
        exit()

    input_file1 = input_arg[1]
    TE_file = input_arg[2]

    tem_file = input_file1.split(".tsv")[0] + "_TEMR.tsv"

    print("#### Input commands ##"+"#" * 50)
    print("\ninput file : {}\nTE repeatMasker file :{}\noutput file : {}\n".format(input_file1, TE_file, tem_file))
    print("#" * 72)
    print("#" * 72)
    print()
    write_output = open(tem_file, "w")

    col_in_bed = 0

    for i in open(input_file1, "r"):
        if "chr" in i:
            col_in_bed = len(i.split("\t"))
        break

    filter_flag = []

    window = 2

    for i in range(1, col_in_bed+1):
        x = "$" + str(i)
        if x == "$2":
            x = "$2-" + str(window)
        elif x == "$3":
            x = "$3+" + str(window)
        filter_flag.append(x)

    awk_filter = ",".join(filter_flag)
    awk_filter = "{print " + awk_filter + "}"
    # print(awk_filter)

    awk_command = subprocess.Popen(['awk', '-v', 'OFS=\t', awk_filter, input_file1], stdout=subprocess.PIPE)
    bedtools_command = subprocess.Popen(['bedtools', 'intersect', '-a', '-', '-b', TE_file, '-wao'],
                                        stdin=awk_command.stdout, stdout=subprocess.PIPE)

    result = make_dict(bedtools_command, col_in_bed)

    for i in result:
        x = check_TEM(i, sorted(result[i]))
        # print(i)
        # if "TE_SE" in x[:5]:
        #     print("$$$","#"*80)
        #     temp = i.split("__")
        #
        #     temp[1] = (int(temp[1]) + 2)
        #     temp[2] = (int(temp[2]) - 2)
        #
        #     write_output.write("\t".join(str(j) for j in temp + [x]) + "\n")
        if "TEMR" in x[:4]:
            y1 = x.split(";")[0].split("_")[-1]
            y2 = x.split(";")[1].split("__")[-2]
            temp = i.split("__")
            temp[1] = str(int(temp[1]) + 2)
            temp[2] = str(int(temp[2]) - 2)
            y = x.split(";")
            te_left = ";".join(y[1].split("__"))
            te_right = ";".join(y[2].split("__"))
            te_flag = y[0] + ";" + y2

            temp += [te_left, te_right, te_flag]
            # print("\t final", temp[-1])
            if "Alu" in temp[-1] or "L1" in temp[-1]:
                write_output.write("\t".join(str(j) for j in temp) + "\n")

    write_output.close()
    print("TEMR (Alu & L1 driven) identification process completed.....\n")


if __name__ == "__main__":

    import subprocess
    import sys
    mini_window = 0
    main()

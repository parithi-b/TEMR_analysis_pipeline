

def make_dict(output_bedtools, col_count):

    output_file = dict()

    out, err = output_bedtools.communicate()

    for i in out.decode().split("\n"):

        sv_data = i.split()
        # print(sv_data)

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

    print(len(output_file))

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
    # print(repeats_within_SV)
    # print(sv)

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

    # exit()
    return percentage_dict


def check_TEM(sv_data, te_list):

    # print("Inside check TEM\n")

    print("#"*50)

    sv  = sv_data.split("__")
    # print("SV : ", "\t".join(sv))

    sv_start = int(sv[1]) + mini_window
    sv_end = int(sv[2]) - mini_window

    print("\tTE in this region" , len(te_list))

    # for j in te_list:
    #     print("\t",j)

    count_te = len(te_list)

    print("SV START, end : ", sv_start, sv_end)

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
                if flag_5!="" and flag_3!="":
                    print("WARNING>........")
                else:
                    if flag_3 =="":
                        flag = flag_5
                    else:
                        flag = flag_3
            # print("flag 5 : {}, flag 3 : {}, flag : {}".format(flag_5,flag_3, flag))
            print(flag)

    else:
        # print("###########TEST TEM")
        overlap_percent = repeat_overlapping(sv,te_list)
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

        if flag_3 == True and flag_5 == True:
            # for j in te_list:
            #     print("\t", j)
            if len(start_list) != 1 or len(end_list)!=1:
                # print("WARNING TEM")
                # print(start_list)
                # print(end_list)
                # print()
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

                # if "ERV" in name_a:
                #     name_a = LTR

                if name_a == name_b:
                    if family_a not in ["SINE","LINE","LTR","DNA","Retroposon"]:
                        flag = "Others"
                    else:
                        broken_flag = False


                        for i in overlap_percent:
                            if overlap_percent[i] > 80:
                                broken_flag = i + ";" + str(overlap_percent[i])
                        # print(overlap_percent, broken_flag, len(te_list))
                        # print(te_list)

                        if broken_flag != False:
                            if len(te_list)==3:
                                flag = "Broken;" + broken_flag
                            # elif len(te_list)==1:
                            #     flag = "Br0ken_mix"
                            elif len(te_list)<3:
                                if float(broken_flag.split(";")[1])<97:
                                    flag = "TEM_broken_by_2_" + orientation + ";" + start_list[0] + ";" + end_list[0]
                                else:
                                    flag = "Broken;" + broken_flag
                            else:
                                flag = "TEM_" + orientation + ";" + start_list[0] + ";" + end_list[0]
                        else:
                            # if len(te_list) == 3 and sum(overlap_percent[i4] for i4 in overlap_percent)>95:
                            #     print("\t#####", sv)
                            #     print("\t#####", overlap_percent)
                            #     flag = "Broken;" + "TE_mix"
                            #     print(flag,sum(overlap_percent[i4] for i4 in overlap_percent))
                            #     print()
                            # else:
                                flag = "TEM_" + orientation + ";" + start_list[0] + ";" + end_list[0]
                else:
                    flag = "nonTEM"

        # print(flag_5, flag_3)
        else:
            if flag_3 ==True or flag_5==True:
                # print(start_list)
                # print(end_list)
                # print("", len(start_list), len(end_list),overlap_percent)
                # print()
                return "TE_SE"
    #
    # print("Flag value : ",flag)
    # print()
    return flag



def main():

    TE_file = "/Users/balacp/Desktop/1000GP_WGS/HGSV/UCSC_files/RepeatMasker/RepeatMasker_canonical_clean_sorted_TE.bed"

    TE_file = "/Users/balacp/Desktop/1000GP_WGS/HGSV/UCSC_files/RepeatMasker/RepeatMasker_canonical_clean_sorted.bed"

    root_folder = "/Users/balacp/Desktop/1000GP_WGS/HGSV/"

    ensemble_file = root_folder + "WGS/Ensemble/Final_phase/All_SV_under50k_illumina_pacbio_combined_minus_SR_MEI_filtered_DEL_homSEQ_with_PAV_corrected.bed"
    # ensemble_file = root_folder+ "WGS/Ensemble/Final_phase/All_SV_under50k_illumina_pacbio_combined_minus_SR_MEI_filtered_DUP_homSEQ_with_PAV_corrected.bed"
    ensemble_file = root_folder + "WGS/Ensemble/Final_phase/All_SV_under50k_illumina_pacbio_combined_minus_SR_MEI_filtered_INV_homSEQ_with_PAV.bed"

    ensemble_file = root_folder + "WGS/Ensemble/Final_phase/INV_post_correction.bed"

    ensemble_file = root_folder + "WGS/Ensemble/All_SV_over50k_default_noRD_3children_pav_coordinates.bed"

    # ensemble_file = root_folder + "WGS/Ensemble/batch2_bad_SVs.tsv"
    #
    # ensemble_file = "/Users/balacp/Desktop/1000GP_WGS/HGSV/WGS/MSI/New_data_April_2021/MSH2KO/Ensemble/Passage5_minus_simpleR_RM.bed"

    # ensemble_file = root_folder + "WGS/Ensemble/All_SV_under50k_3children_merged_pacbio_80RO_RD_minus_SR_filtered_TEM_corrected_breakpoint.bed"
    #
    # ensemble_file = root_folder + "WGS/Ensemble/validation/junction_reconstructed.bed"

    # pacbio_ensemble = root_folder + "PacBio/Ensemble/All_SV_pacbio_under50k_3children_with_RD_merged_minus_SR_filtered.bed"
    # ensemble_file = pacbio_ensemble

    tem_file = ensemble_file.split(".bed")[0] + "_TEM.bed"

    write_output = open(tem_file,"w")

    col_in_bed = 0

    for i in open(ensemble_file,"r"):
        if "chr" in i:
            col_in_bed = len(i.split("\t"))
        break

    filter_flag = []

    window = 2

    for i in range(1,col_in_bed+1):
        x="$" + str(i)
        if x =="$2":
            x = "$2-"+str(window)
        elif x =="$3":
            x = "$3+"+str(window)

        filter_flag.append(x)

    awk_filter = ",".join(filter_flag)

    awk_filter = "{print " + awk_filter + "}"

    # print(awk_filter)

    awk_command = subprocess.Popen(['awk','-v','OFS=\t', awk_filter, ensemble_file],stdout=subprocess.PIPE)

    bedtools_command = subprocess.Popen(['bedtools', 'intersect', '-a', '-', '-b', TE_file, '-wao'], stdin=awk_command.stdout,
                                        stdout=subprocess.PIPE)
    count = 0

    result = make_dict(bedtools_command,col_in_bed)

    tem_count = dict()

    for i in result:
        print("$$$", i)
        # count+=1
        x = check_TEM(i,sorted(result[i]))
        print("####", x)
        # print(i)
        # if "TE_SE" in x[:5]:
        #     print("$$$","#"*80)
        #     temp = i.split("__")
        #
        #     temp[1] = (int(temp[1]) + 2)
        #     temp[2] = (int(temp[2]) - 2)
        #
        #     write_output.write("\t".join(str(j) for j in temp + [x]) + "\n")

        if "TEM" in x[:3]:
            # if "broken_by_2" in x:
            #     continue
            count +=1
            print(i)
            print("\t", x)

            y1 = x.split(";")[0].split("_")[-1]
            y2 = x.split(";")[1].split("__")[-2]
            # print(y1,y2)

            if y1 not in tem_count:
                tem_count[y1] = dict()

            if y2 not in tem_count[y1]:
                tem_count[y1][y2] = 0

            tem_count[y1][y2] += 1

            temp = i.split("__")
            temp[1] = str(int(temp[1]) + 2)
            temp[2] = str(int(temp[2]) - 2)

            y = x.split(";")
            te_left = ";".join(y[1].split("__"))
            te_right = ";".join(y[2].split("__"))
            te_flag = y[0]+";"+y2

            temp += [te_left,te_right,te_flag]
            print("\t final",temp)

            write_output.write("\t".join(str(j) for j in temp) + "\n")



    print("\nRESULTS")
    print(count)
    maybe = dict()
    for i in tem_count:
        print(i)
        for j in tem_count[i]:
            if j not in maybe:
                maybe[j] = 0
            maybe[j] += tem_count[i][j]
            print("\t", j , tem_count[i][j])

    print()
    for i in maybe:
        print(i, maybe[i])

    write_output.close()


if __name__ == "__main__":

    import os, subprocess

    mini_window = 10

    os.environ['PATH'] += ":/Users/balacp/anaconda3/envs/pyth38_plot/bin"

    main()
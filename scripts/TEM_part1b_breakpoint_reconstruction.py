

def make_dict(bedtools_output):
    """

    create a dictonary with the parents and child data using bedtools intersect and '-loj' parameter
    parents should have the same sv type as the child

    :param bedtools_output: output from bedtools

    :return: dict, with key as SVs from child and value being Parent name and corresponding callers

    """

    sv_count = 0

    sv_database = dict()

    key_col_count = 6   # (number of column in main bed + 1)

    ### Creating the dictionary

    for sv in bedtools_output.decode().split("\n"):


        if len(sv) == 0:
            continue

        sv_count += 1

        y = sv.split('\t')
        ### Key

        key_start_pos = int(y[1])
        key_end_pos = int(y[2])

        key_temp = "$$".join(y[:key_col_count])

        if key_temp not in sv_database:
            sv_database[key_temp] = []

        ### Value
        value_start_pos = int(y[1 + key_col_count])

        value_temp = "$$".join(y[key_col_count:])

        # if "chr5" in key_temp and "40848" in key_temp:
        #     print(key_temp, value_temp)

        if value_temp not in sv_database[key_temp]:
            if (key_start_pos - 500 < value_start_pos < key_start_pos + 500) or (
                        key_end_pos - 500 < value_start_pos < key_end_pos + 500):

                sv_database[key_temp].append(value_temp)

    print("SV in child : ", len(sv_database))

    return sv_database


def junction_reconstruction(bed_file, svaba, write_output):

    if write_output is True:
            final_output_file = bed_file.split(".")[0]+ '_svaba_temp.bed'
            merged_file = open(final_output_file, "w")

    bedtools_unique = subprocess.Popen(
        ['bedtools', 'window', '-w', '500', '-a', bed_file, '-b', svaba], stdout=subprocess.PIPE)

    out, err = bedtools_unique.communicate()

    bedtools_dict = make_dict(out)

    counter = 0

    new_counter = 0

    final_list = []

    for i in bedtools_dict:
        # if "chr5" in i and "40848" in i:
        #     print(i , bedtools_dict[i])

        # if "TEM_HOM" not in i:
        #     continue
        #

        caller_data = i.split("$$")
        print(caller_data)

        caller_start = int(caller_data[1])
        caller_end = int(caller_data[2])

        if caller_end - caller_start > 50000:
            continue

        counter += 1

        testing = dict()

        sv_type = caller_data[3]

        for j in bedtools_dict[i]:

            svaba_data = j.split("$$")

            chr_start = svaba_data[0]
            start_pos = int(svaba_data[1])-1
            info_vcf = svaba_data[7]

            qual = svaba_data[5]

            for sub_info in info_vcf.split(";"):
                if "EVDNC=" in sub_info[:6]:
                    evidence = (sub_info.split("=")[-1])
                elif "SPAN=" in sub_info[:5]:
                    sv_len = (sub_info.split("=")[-1])

            if sv_len not in testing:
                testing[sv_len] = []

            testing[sv_len].append((chr_start, start_pos, qual, evidence))

        """
        ##INFO=<ID=EVDNC,Number=1,Type=String,Description="Evidence for variant. 
        ASSMB assembly only, 
        ASDIS assembly+discordant. 
        DSCRD discordant only, 
        TSI_L templated-sequence insertion (local, e.g. AB or BC of an ABC), 
        TSI_G global (e.g. AC of ABC)">
        """
        # else: # when only one BND is present  --> MISSING
        #     print(testing)

        svaba_start = 0
        svaba_end = 0

        if len(testing) == 1:
            for k in testing:
                if len(testing[k]) == 2:
                    new_counter += 1
                    svaba_start = testing[k][0][1]
                    svaba_end = testing[k][1][1]
                    if int(svaba_start) > int(svaba_end):
                        svaba_start = testing[k][1][1]
                        svaba_end = testing[k][0][1]
                else: ### needs better conditions

                    if len(testing[k]) % 2 == 0:
                        temp = dict()
                        new_counter += 1

                        for k1 in testing[k]:
                            temp_q = int(k1[2])
                            if temp_q not in temp:
                                temp[temp_q] = []
                            temp[temp_q].append(k1[1])

                        if len(temp) == 1:
                            size_diff = 5000
                            key_temp = list(temp.keys())[0]
                            half_size = sum(map(len, temp.values()))//2
                            for c in range(0,half_size):
                                a1 = int(temp[key_temp][c])
                                b1 = int(temp[key_temp][c+half_size])

                                if abs(caller_start-a1) + abs(caller_end-b1) < size_diff:
                                    a = a1
                                    b = b1
                                    size_diff = abs(caller_start-a1) + abs(caller_end-b1)
                        else:
                            a = temp[sorted(temp, reverse=True)[0]][0]
                            b = temp[sorted(temp, reverse=True)[0]][1]

                        svaba_start = a
                        svaba_end = b
        else:

            flag = 0
            for k in testing:
                if len(testing[k]) == 2:
                    flag += 1
                    temp = testing[k]

            if flag == 1:
                    new_counter += 1
                    svaba_start = temp[0][1]
                    svaba_end = temp[1][1]
                    if int(svaba_start) > int(svaba_end):
                        svaba_start = temp[1][1]
                        svaba_end = temp[0][1]

            elif flag > 1:

                evd_flag = 0 # evidence flag
                sv_size = int(caller_end) - int(caller_start)
                size_diff = 5000

                for k in testing:

                    if len(testing[k]) == 2:
                        svaba_start = testing[k][0][1]
                        svaba_end = testing[k][1][1]

                        if int(svaba_start) > int(svaba_end):
                            svaba_start = testing[k][1][1]
                            svaba_end = testing[k][0][1]

                        svaba_size = (int(svaba_end) - int(svaba_start))

                        if  svaba_size < (sv_size * .85) or svaba_size > (sv_size * 1.15):
                            # print("asdfasdfsad",sv_size * .85, sv_size * 1.15)
                            # print(svaba_size)
                            continue

                        # print("PASD", svaba_size, testing[k][0][-1],testing[k][1][-1])

                        if ("ASDIS" in testing[k][0][-1] or "ASDIS" in testing[k][1][-1]):
                            evd_flag = 1
                            if abs(sv_size - svaba_size) < size_diff:
                                size_diff = abs(sv_size - svaba_size)
                                a = svaba_start
                                b = svaba_end

                        elif ("ASSMB" in testing[k][0][-1] or "ASSMB" in testing[k][1][-1]) and evd_flag !=1 :
                            if evd_flag == 3:
                                a = svaba_start
                                b = svaba_end
                            else:
                                if abs(sv_size - svaba_size) < size_diff:
                                    size_diff = abs(sv_size - svaba_size)
                                    a = svaba_start
                                    b = svaba_end
                            evd_flag = 2
                        else:
                            if abs(sv_size - svaba_size) < size_diff and (evd_flag==0 or evd_flag==3):
                                size_diff = abs(sv_size - svaba_size)
                                a = svaba_start
                                b = svaba_end
                            evd_flag = 3
                if evd_flag > 0:
                    new_counter += 1
                    svaba_start = a
                    svaba_end = b

        final_sv = caller_data[:6]

        if svaba_start !=0:
            final_sv[1] = str(svaba_start)
            final_sv[2] = str(svaba_end)
            final_sv[-1] += ";svaba"
        #
        # if "chr5" in final_sv and "40848" in final_sv[1]:
        #     print("Adf", final_sv)
        final_list.append(final_sv)

    bedtools_unique = subprocess.Popen(
        ['bedtools', 'window', '-w', '500', '-a', bed_file, '-b', svaba,'-v'], stdout=subprocess.PIPE)

    out, err = bedtools_unique.communicate()

    for i in out.decode().split("\n"):

        y = i.split("\t")
        if len(y) <3:
            continue
        final_sv = y[:6]
        final_list.append(final_sv)

    if write_output is True:
        for i in final_list:
            merged_file.write("\t".join(i) + "\n")

    merged_file.close()

    if write_output is True:
            final_output = bed_file.split(".")[0]+ '_svaba.bed'
            merged_file = open(final_output, "w")

    subprocess.call(['sort', '-k1,1', '-k2,2n', final_output_file], stdout=merged_file)

    subprocess.call(['rm',final_output_file])

    print(counter, new_counter)
    print(len(final_list))


def main():

    write_output_flag = True

    trio_id = "HG00733"

    root_folder = "/Users/balacp/Desktop/1000GP_WGS/HGSV/WGS/"

    input_file = root_folder + "Ensemble/" + trio_id + "/" + trio_id + "_10_5_RD_merged_sorted_minus_UCSC_under50k.bed"

    svaba_bed = root_folder + "svaba/" + trio_id + "/" + trio_id + "_svaba_with_support.bed"

    junction_reconstruction(input_file, svaba_bed, write_output_flag)


if __name__ == "__main__":

    import subprocess
    import pandas as pd
    import os

    os.environ['PATH'] += ":/Users/balacp/anaconda3/envs/pyth38_plot/bin"


    main()
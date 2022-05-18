
"""
This script is used to convert vcf to bed for the pacbio callers used in Beck Lab
input : .vcf file ,  sample name and caller used
output : 3 bed files
    1. containing all DEL, DUP and INV
    2. containing INS
    3. Other types of SVs if any mentioned in the VCF (BND and complex)

Fitlering : Variants under 50bp are filtered out

Columns in bed file so far

CHR
START
END
SVTYPE
CALLER NAME
SAMPLE ID
SV ID from vcf (col 3 in vcf -> unique ID per event in vcf)
READ SUPPORT (based on caller) # the high quality variant in the ALT
FLAG (based on caller --> 7th column ,
            if a caller only reports PASS, this can be substituted with other values,
                like for sniffles)
GT

search for CHECK_READ_SUPPORT will take you to the location to edit the read support


For Peter's project the minimum support for variant caller set was 3
Beck lab project support = 5

sniffles default is 10
svim recommends 10
pbsv --> suggest all caller use a support of 3 in their standard practice for SV calling

"""


def sv_type_categorizing(sv_bed_db):

    """

    :param sv_bed_db: list of SVs are provided
    :return: None :  prints # of SVs in each type
    """
    print("SV type categorizing")

    sv_type_db = dict()

    for i in sv_bed_db:
        sv_data = i.split("_")

        pos = sv_data[1]
        end = sv_data[2]
        given_sv = sv_data[3]

        if given_sv not in sv_type_db:
            sv_type_db[given_sv] = 0

        if given_sv in ["DEL", "DUP", "INV"]:
            if abs(int(end) - int(pos)) >= 50:
                sv_type_db[given_sv] += 1
        else:
            sv_type_db[given_sv] += 1

    for i in sv_type_db:
        print(i, sv_type_db[i])


def sniffles_vcf_to_bed(vcf_path, writing_output_flag, sample):

    print("Converting VCF to BED  -  SNIFFLES ")

    write_output_file = writing_output_flag  # set true if u want to write your output in a bed file

    if write_output_file is True:

        bed_write = open(vcf_path.split(".")[0] + "_DEL_DUP_INV.bed", 'w')
        ins_bed_write = open(vcf_path.split(".")[0] + "_INS.bed", 'w')
        others_bed_write = open(vcf_path.split(".")[0] + "_OTHERS.bed", 'w')

    total_sv = 0

    sv_db = []

    for sv in open(vcf_path, 'r'):

        if "#" in sv:   # skip header
            continue

        if "SVTYPE=BND;" in sv:     # skip BND
            continue

        # if "SVTYPE=DEL;" not in sv and "SVTYPE=DUP;" not in sv and "SVTYPE=INV;" not in sv:
        #     continue

        total_sv += 1

        y = sv.split()

        chr_start = y[0]
        start_pos = int(y[1])
        if start_pos <= 1:
            start_pos = '1'
        else:
            start_pos = str(start_pos)

        alt_seq = y[4]
        info_vcf = y[7]

        sv_len = -1

        for sub_info in info_vcf.split(";"):
            if "SVLEN=" in sub_info[:6]:
                sv_len = (sub_info.split("=")[-1])
                if sv_len == ".":
                    # print("ASDF", info_vcf)
                    # print(sv)
                    continue
                else:
                    sv_len = abs(int(sv_len))
            elif "SVTYPE=" in sub_info[:7]:
                sv_type = (sub_info.split("=")[-1])
            elif "END=" in sub_info[:4]:
                end_pos = str(int(sub_info.split("=")[-1]))
            elif "RE=" in sub_info[:3]:
                pe = sub_info.split("=")[-1]

        if sv_len == ".":
            continue

        len_temp = abs(int(start_pos) - int(end_pos))

        if len_temp != sv_len and sv_type in ["DEL", "DUP", "INV"]:
            sv_len = len_temp

        gt_field_label = y[8].split(':')
        gt_field = y[9].split(':')

        dhffc_location = -1
        dhbfc_location = -1

        for i in range(len(gt_field_label)):
            if gt_field_label[i] == "DHFFC":
                dhffc_location = i
            if gt_field_label[i] == "DHBFC":
                dhbfc_location = i

        dhffc_value = gt_field[dhffc_location]
        dhbfc_value = gt_field[dhbfc_location]

        key_string = chr_start + "_" + start_pos + "_" + end_pos + "_" + sv_type

        temp = list()
        temp.append(chr_start)
        temp.append(start_pos)
        temp.append(end_pos)
        temp.append(sv_type)
        temp.append(sample)
        temp.append("sniffles")
        temp.append(pe)
        temp.append(dhffc_value)
        temp.append(dhbfc_value)

        if ("SVTYPE=DEL" not in sv) and ("SVTYPE=INV" not in sv) and ("SVTYPE=DUP" not in sv):
            temp.append(sv_len)
            temp.append(alt_seq)

        if key_string not in sv_db:
            sv_db.append(key_string)

            if write_output_file is True:
                if sv_type in ["DEL", "DUP", "INV"]:
                    if abs(sv_len) >= 50:
                            bed_write.write("\t".join(str(i) for i in temp) + "\n")
                elif sv_type == "INS":
                    if abs(sv_len) >= 50:
                        ins_bed_write.write("\t".join(str(i) for i in temp) + "\n")
                else:
                    others_bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants = {} ".format(total_sv))
    sv_type_categorizing(sv_db)

    if write_output_file is True:
        bed_write.close()
        ins_bed_write.close()
        others_bed_write.close()

    print("*" * 50)


def pbsv_vcf_to_bed(vcf_path, writing_output_flag, sample):

    print("\nConverting VCF to BED  -  PBSV ")

    write_output_file = writing_output_flag  # set true if u want to write your output in a bed file

    if write_output_file is True:

        bed_write = open(vcf_path.split(".")[0] + "_DEL_DUP_INV.bed", 'w')
        ins_bed_write = open(vcf_path.split(".")[0] + "_INS.bed", 'w')
        others_bed_write = open(vcf_path.split(".")[0] + "_OTHERS.bed", 'w')

    total_sv = 0
    sv_db = []

    for sv in open(vcf_path, 'r'):

        if "#" in sv:
            continue

        if "SVTYPE=BND;" in sv:
            continue

        # if "SVTYPE=DEL;" not in sv and "SVTYPE=DUP" not in sv and "SVTYPE=INV" not in sv:
        #     continue

        total_sv += 1

        y = sv.split()

        chr_start = (y[0])
        start_pos = str(int(y[1]))
        alt_seq = y[4]

        info_vcf = y[7]

        sv_len = -1

        for sub_info in info_vcf.split(";"):

            if "SVLEN=" in sub_info[:6]:
                sv_len = (sub_info.split("=")[-1])
                sv_len = abs(int(sv_len))
            if "SVTYPE=" in sub_info[:7]:
                sv_type = (sub_info.split("=")[-1])
            elif "END=" in sub_info[:4]:
                end_pos = str(int(sub_info.split("=")[-1]))

        re = y[9]  # read evidence

        if sv_type == "cnv":
            # print(re, sv_type)
            # print(chr_start, start_pos, sv_len)
            # print()
            total_reads = int(re.split(":")[0])
        else:
            # ref_support = int(re.split(":")[1].split(',')[0])
            alt_support = int(re.split(":")[1].split(',')[1])

            total_reads = alt_support

        key_string = chr_start + "_" + start_pos + "_" + end_pos + "_" + sv_type

        gt_field_label = y[8].split(':')
        gt_field = y[9].split(':')

        dhffc_location = -1
        dhbfc_location = -1

        for i in range(len(gt_field_label)):
            if gt_field_label[i] == "DHFFC":
                dhffc_location = i
            if gt_field_label[i] == "DHBFC":
                dhbfc_location = i

        dhffc_value = gt_field[dhffc_location]
        dhbfc_value = gt_field[dhbfc_location]

        len_temp = abs(int(start_pos) - int(end_pos))

        if len_temp != sv_len and sv_type in ["DEL", "DUP", "INV"]:
            sv_len = len_temp

        temp = list()

        temp.append(chr_start)
        temp.append(start_pos)
        temp.append(end_pos)
        temp.append(sv_type)
        temp.append(sample)
        temp.append("pbsv")
        temp.append(total_reads)
        temp.append(dhffc_value)
        temp.append(dhbfc_value)

        if sv_type not in ["DEL", "DUP", "INV"]:
            temp.append(sv_len)
            temp.append(alt_seq)

        if key_string not in sv_db:
            sv_db.append(key_string)

            if write_output_file is True:

                if sv_type in ["DEL", "DUP", "INV"]:
                    if abs(sv_len) >= 50:
                        bed_write.write("\t".join(str(i) for i in temp) + "\n")
                elif sv_type == "INS":
                    if abs(sv_len) >= 50:
                        ins_bed_write.write("\t".join(str(i) for i in temp) + "\n")
                else:
                    # print("TEMP::, ", temp)
                    if abs(sv_len) >= 50:

                        others_bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants = {} ".format(total_sv))
    sv_type_categorizing(sv_db)

    if write_output_file is True:
        ins_bed_write.close()
        bed_write.close()
        others_bed_write.close()

    print("*" * 50)


def svim_vcf_to_bed(vcf_path, writing_output_flag, sample):

    print("\nConverting VCF to BED  -  svim ")

    write_output_file = writing_output_flag  # set true if u want to write your output in a bed file

    if write_output_file is True:

        bed_write = open(vcf_path.split(".")[0] + "_DEL_DUP_INV.bed", 'w')
        ins_bed_write = open(vcf_path.split(".")[0] + "_INS.bed", 'w')
        others_bed_write = open(vcf_path.split(".")[0] + "_OTHERS.bed", 'w')

    total_sv = 0
    sv_db = []

    for sv in open(vcf_path, 'r'):

            if "#" in sv:
                continue

            if "SVTYPE=BND;" in sv:
                continue

            total_sv += 1

            y = sv.split()

            chr_start = y[0]
            start_pos = str(int(y[1]))
            alt_seq = y[4]

            info_vcf = y[7]

            sv_len = -1

            for sub_info in info_vcf.split(";"):

                if "SVLEN=" in sub_info[:6]:
                    sv_len = (sub_info.split("=")[-1])
                    sv_len = abs(int(sv_len))
                if "SVTYPE=" in sub_info[:7]:
                    sv_type = (sub_info.split("=")[-1])
                elif "END=" in sub_info[:4]:
                    end_pos = str(int(sub_info.split("=")[-1]))
                elif "SUPPORT=" in sub_info[:8]:
                    re = (sub_info.split("=")[-1])

            if "DUP:TA" in sv_type:
                sv_type = "DUP"

            len_temp = abs(int(start_pos) - int(end_pos))

            if len_temp != sv_len and sv_type in ["DEL", "DUP", "INV"]:
                sv_len = len_temp

            gt_field = y[9].split(':')

            gt_field_label = y[8].split(':')

            dhffc_location = -1
            dhbfc_location = -1

            for i in range(len(gt_field_label)):
                if gt_field_label[i] == "DHFFC":
                    dhffc_location = i
                if gt_field_label[i] == "DHBFC":
                    dhbfc_location = i

            dhffc_value = gt_field[dhffc_location]
            dhbfc_value = gt_field[dhbfc_location]

            key_string = chr_start + "_" + start_pos + "_" + end_pos + "_" + sv_type

            temp = list()

            temp.append(chr_start)
            temp.append(start_pos)
            temp.append(end_pos)
            temp.append(sv_type)
            temp.append(sample)
            temp.append("svim")
            temp.append(re)
            temp.append(dhffc_value)
            temp.append(dhbfc_value)

            if sv_type not in ["DEL", "DUP", "INV"]:
                temp.append(sv_len)
                temp.append(alt_seq)

            if key_string not in sv_db:
                sv_db.append(key_string)
                if write_output_file is True:

                    if sv_type in ["DEL", "DUP", "INV"]:
                        if abs(sv_len) >= 50:
                            bed_write.write("\t".join(str(i) for i in temp) + "\n")
                    elif sv_type == "INS":
                        ins_bed_write.write("\t".join(str(i) for i in temp) + "\n")
                    else:
                        others_bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants = {} ".format(total_sv))
    sv_type_categorizing(sv_db)

    if write_output_file is True:
        bed_write.close()
        ins_bed_write.close()
        others_bed_write.close()

    print("*" * 50)


def main():

    # root_folder = "/Users/balacp/Desktop/Christine/Robinson/svCalling/"
    #
    # # bed_file = root_folder + "PID-1048/PID_1048_CLR_3runs_s3_svim.vcf"
    # # bed_file = root_folder + "PID-1048/PID_1048_CLR_3runs_pbsv.vcf"
    # bed_file = root_folder + "PID-1048/PID_1048_CLR_3runs_sniffles.vcf"
    #
    # caller_used = "sniffles"  # sniffles or pbsv or svim
    # sample_name = "PID1048"  # please avoid -/;/:/_
    # print(sample_name)
    #
    # if caller_used == "sniffles":
    #     sniffles_vcf_to_bed(bed_file, write_output, sample_name)
    # elif caller_used == "pbsv":
    #     pbsv_vcf_to_bed(bed_file, write_output, sample_name)
    # elif caller_used == "svim":
    #     svim_vcf_to_bed(bed_file, write_output, sample_name)
    # else:
    #     print("Please mention the proper name of the tool used")

    data_root_folder = "/Users/balacp/Desktop/1000GP_WGS/HGSV/PacBio/"

    for sample_name in ["HG00514", "HG00733", "NA19240"]:
    # for sample_name in ["HG00514"]:

        for caller_used in ["sniffles", "pbsv","SVIM"]:

            write_output = True  # do not put True or False under quotes
            print(sample_name, caller_used)

            vcf_file = data_root_folder + caller_used + "/" + sample_name
            vcf_file += "/" + sample_name + "_" + caller_used.lower() + "_duphold.vcf"

            if caller_used == "sniffles":
                sniffles_vcf_to_bed(vcf_file, write_output, sample_name)
            elif caller_used == "pbsv":
                pbsv_vcf_to_bed(vcf_file, write_output, sample_name)
            elif caller_used == "SVIM":
                svim_vcf_to_bed(vcf_file, write_output, sample_name)
            else:
                print("Please mention the proper name of the tool used")


if __name__ == "__main__":

    main()

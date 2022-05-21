"""
Make sure duphold is run before using this script
only DEL & DUP & INV are considered

This program is used to convert vcf to bed/tsv file for the following tools

short-read : manta/delly/lumpy
long-read : pbsv/sniffles/svim

BED fILE contains

CHR START END SVTYPE CALLER SAMPLEID PAIRED-READ(PR) SPLIT-READ(SR) DHFFC DHBFC
CHR START END SVTYPE CALLER SAMPLEID READ-SUPPORT(RS) DHFFC DHBFC

"""


def manta_vcf_to_bed(input_bed_file, output_bed_file, sample_id, lower_sv, upper_sv):

    print("Converting VCF to TSV  -  MANTA ")

    bed_write = open(output_bed_file, 'w')

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "#" in sv:
            continue

        if "chrY" in sv: # all female samples
            continue

        if "MinQUAL" in sv:
            continue

        total_sv += 1

        if "SVTYPE=INV" in sv or "SVTYPE=DEL" in sv or "SVTYPE=DUP" in sv:

            y = sv.split()

            chr_start = y[0]
            start_pos = str(int(y[1]))
            info_vcf = y[7]

            for sub_info in info_vcf.split(";"):
                if "SVTYPE=" in sub_info[:7]:
                    sv_type = (sub_info.split("=")[-1])
                elif "END=" in sub_info[:4]:
                    end_pos = str(int(sub_info.split("=")[-1]))

            gt_field = y[9].split(':')

            gt_field_label = y[8].split(':')

            dhffc_location = -1
            dhbfc_location = -1
            pe_location = -1
            sr_location = -1

            sr_support = 0

            for i in range(len(gt_field_label)):
                if gt_field_label[i] == "DHFFC":
                    dhffc_location = i
                elif gt_field_label[i] == "DHBFC":
                    dhbfc_location = i
                elif gt_field_label[i] == "PR":
                    pe_location = i
                elif gt_field_label[i] == "SR":
                    sr_location = i

            dhffc_value = gt_field[dhffc_location]
            dhbfc_value = gt_field[dhbfc_location]
            pe_support = gt_field[pe_location].split(",")[1]

            if sr_location != -1:
                sr_support = gt_field[sr_location].split(",")[1]

            precision_support = gt_field[1]

            if precision_support != "PASS":
                continue

            key_temp = chr_start + "_" + start_pos + "_" + end_pos + "_" + sv_type

            temp = list()
            temp.append(chr_start)
            temp.append(start_pos)
            temp.append(end_pos)
            temp.append(sv_type)
            temp.append(sample_id)
            temp.append("manta")
            temp.append(pe_support)
            temp.append(sr_support)
            temp.append(dhffc_value)
            temp.append(dhbfc_value)

            if key_temp not in sv_db:
                sv_db.append(key_temp)
                sv_length = abs(int(end_pos) - int(start_pos))
                if sv_length >= lower_sv and sv_length < upper_sv:
                    sv_count += 1
                    bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("DEL, DUP, & INV only and between {}bp and {}bp".format(lower_sv, upper_sv))
    print("Total SVs considered  = {} ".format(sv_count))

    bed_write.close()
    print("*" * 50)


def delly_vcf_to_bed(input_bed_file, output_bed_file, sample_id, lower_sv, upper_sv):

    print("Converting VCF to TSV  -  DELLY ")

    bed_write = open(output_bed_file, 'w')

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "#" in sv:
            continue

        if "chrY" in sv: # all female samples
            continue

        if "SVTYPE=INV;" in sv or "SVTYPE=DEL;" in sv or "SVTYPE=DUP;" in sv:

            y = sv.split()

            if y[6] == "LowQual":  # remove lowQual reads
                continue

            total_sv += 1

            chr_start = y[0]
            start_pos = str(int(y[1]))

            info_vcf = y[7]

            sr_support = 0
            for sub_info in info_vcf.split(";"):
                if "SVTYPE=" in sub_info[:7]:
                    sv_type = (sub_info.split("=")[-1])
                elif "END=" in sub_info[:4]:
                    end_pos = str(int(sub_info.split("=")[-1]))
                elif "PE=" in sub_info[:3]:
                    pe_support = int(sub_info.split("=")[-1])
                elif "SR=" in sub_info[:3]:
                    sr_support = int(sub_info.split("=")[-1])

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

            key_temp = chr_start + "_" + start_pos + "_" + end_pos + "_" + sv_type

            temp = list()
            temp.append(chr_start)
            temp.append(start_pos)
            temp.append(end_pos)
            temp.append(sv_type)
            temp.append(sample_id)
            temp.append("delly")
            temp.append(pe_support)
            temp.append(sr_support)
            temp.append(dhffc_value)
            temp.append(dhbfc_value)

            if key_temp not in sv_db:
                sv_db.append(key_temp)
                sv_length = abs(int(end_pos) - int(start_pos))
                if sv_length >= lower_sv and sv_length < upper_sv:
                    sv_count += 1
                    bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("DEL, DUP, & INV only and between {}bp and {}bp".format(lower_sv, upper_sv))
    print("Total SVs considered  = {} ".format(sv_count))

    print("*" * 50)


def lumpy_vcf_to_bed(input_bed_file, output_bed_file, sample_id, lower_sv, upper_sv):

    print("Converting VCF to TSV  -  LUMPY ")

    bed_write = open(output_bed_file, 'w')

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "#" in sv:
            continue

        if "chrY" in sv: # all female samples
            continue

        total_sv += 1

        if "SVTYPE=INV" in sv or "SVTYPE=DEL" in sv or "SVTYPE=DUP" in sv:

            y = sv.split()

            chr_start = y[0]
            start_pos = str(int(y[1]))

            info_vcf = y[7]

            for sub_info in info_vcf.split(";"):

                if "SVTYPE=" in sub_info[:7]:
                    sv_type = (sub_info.split("=")[-1])
                elif "END=" in sub_info[:4]:
                    end_pos = str(int(sub_info.split("=")[-1]))
                elif "PE=" in sub_info[:3]:
                    pe_support = (sub_info.split("=")[-1])
                elif "SR=" in sub_info[:3]:
                    sr_support = (sub_info.split("=")[-1])

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

            key_temp = chr_start + "_" + start_pos + "_" + end_pos + "_" + sv_type

            temp = list()
            temp.append(chr_start)
            temp.append(start_pos)
            temp.append(end_pos)
            temp.append(sv_type)
            temp.append(sample_id)
            temp.append("lumpy")
            temp.append(pe_support)
            temp.append(sr_support)
            temp.append(dhffc_value)
            temp.append(dhbfc_value)

            if key_temp not in sv_db:
                sv_db.append(key_temp)
                sv_length = abs(int(end_pos) - int(start_pos))
                if sv_length >= lower_sv and sv_length < upper_sv:
                    sv_count += 1
                    bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("DEL, DUP, & INV only and between {}bp and {}bp".format(lower_sv, upper_sv))
    print("Total SVs considered  = {} ".format(sv_count))

    print("*" * 50)


def sniffles_vcf_to_bed(input_bed_file, output_bed_file, sample_id, lower_sv, upper_sv):

    print("Converting VCF to TSV  -  SNIFFLES ")

    bed_write = open(output_bed_file, 'w')

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "#" in sv:   # skip header
            continue

        if "chrY" in sv: # all female samples
            continue

        if "SVTYPE=BND;" in sv:     # skip BND
            continue

        if "SVTYPE=DEL;" not in sv and "SVTYPE=DUP;" not in sv and "SVTYPE=INV;" not in sv:
            continue

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
        temp.append(sample_id)
        temp.append("sniffles")
        temp.append(pe)
        temp.append(dhffc_value)
        temp.append(dhbfc_value)

        if key_string not in sv_db:
            sv_db.append(key_string)
            sv_length = abs(sv_len)
            if sv_length >= lower_sv and sv_length < upper_sv:
                sv_count += 1
                bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("DEL, DUP, & INV only and between {}bp and {}bp".format(lower_sv, upper_sv))
    print("Total SVs considered  = {} ".format(sv_count))

    bed_write.close()
    print("*" * 50)


def pbsv_vcf_to_bed(input_bed_file, output_bed_file, sample_id, lower_sv, upper_sv):

    print("\nConverting VCF to TSV  -  pbsv ")

    bed_write = open(output_bed_file, 'w')

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "#" in sv:
            continue

        if "chrY" in sv: # all female samples
            continue

        if "SVTYPE=BND;" in sv:
            continue

        if "SVTYPE=DEL;" not in sv and "SVTYPE=DUP" not in sv and "SVTYPE=INV" not in sv:
            continue

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
        temp.append(sample_id)
        temp.append("pbsv")
        temp.append(total_reads)
        temp.append(dhffc_value)
        temp.append(dhbfc_value)

        if key_string not in sv_db:
            sv_db.append(key_string)
            sv_length = abs(sv_len)
            if sv_length >= lower_sv and sv_length < upper_sv:
                sv_count += 1
                bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("DEL, DUP, & INV only and between {}bp and {}bp".format(lower_sv, upper_sv))
    print("Total SVs considered  = {} ".format(sv_count))

    bed_write.close()
    print("*" * 50)

def svim_vcf_to_bed(input_bed_file, output_bed_file, sample_id, lower_sv, upper_sv):

    print("\nConverting VCF to BED  -  svim ")

    bed_write = open(output_bed_file, 'w')

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "chrY" in sv: # all female samples
            continue

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

        if sv_type not in ["DEL", "DUP", "INV"]:
            continue

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
        temp.append(sample_id)
        temp.append("svim")
        temp.append(re)
        temp.append(dhffc_value)
        temp.append(dhbfc_value)

        if key_string not in sv_db:
            sv_db.append(key_string)
            sv_length = abs(sv_len)
            if sv_length >= lower_sv and sv_length < upper_sv:
                sv_count += 1
                bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("DEL, DUP, & INV only and between {}bp and {}bp".format(lower_sv, upper_sv))
    print("Total SVs considered  = {} ".format(sv_count))

    bed_write.close()
    print("*" * 50)


def filter_gap_centromere_simple_repeats(input_file, gap_cent, simple_rep):
    """

    :param input_file: unfiltered tsv file
    :param gap_cent: gap and cent location file
    :param simple_rep: simple repeat file
    :return:
    """

    print("input: ", input_file)
    output_file_temp = input_file.split("unfiltered.tsv")[0] + 'filtered_temptemr.tsv'
    write_file_temp = open(output_file_temp, "w")

    ### remove gap and cent
    subprocess.call(['bedtools', 'window', '-w', '500', '-a', input_file, '-b', gap_cent, '-v'],
                    stdout=write_file_temp)

    ### remove simpleR
    output_file = input_file.split("unfiltered.tsv")[0] + 'filtered.tsv'
    write_file = open(output_file, "w")

    subprocess.call(['bedtools', 'intersect', '-f', '.50', '-a', output_file_temp, '-b', simple_rep, '-v'],
                   stdout=write_file)

    os.remove(output_file_temp)

    count = 0
    for i in open(output_file, "r"):
        if len(i.split()) > 5:
            count += 1

    print("SVs after filtering : {}".format(count))

    return output_file


def main():

    print("#### Input commands ##"+"#" * 50)
    input_arg = sys.argv

    if len(input_arg) != 8:
        print("Incorrect number of input parameter...")
        print("Please check the query, there should be 8 entries including the script and vcf file...")
        print("Make sure to provide the gap, centromere, simple repeat files provided with the script...")
        print("#" * 72, "\n")
        exit()

    input_file = input_arg[1]
    gap_and_centromere = input_arg[2]
    simple_repeat = input_arg[3]
    sample_id = input_arg[4]
    caller = input_arg[5].lower()
    lower_limit = int(input_arg[6])
    upper_limit = int(input_arg[7])

    if input_file.split(".")[-1] != "vcf":
        print("Please provide a vcf file and make sure the file is unzipped...")
        print("#" * 72, "\n")
        exit()

    if gap_and_centromere.split(".")[-1] == "gz" or simple_repeat.split(".")[-1] == "gz":
        print("Please make sure the files are unzipped...")
        print("#" * 72, "\n")
        exit()

    print(
        "File Name : {} \nSAMPLE ID : {}  \nCALLER : {} \nSV length range : {}bp to {}bp \n"
        "gap&centromere file : {} \nsimple repeat file : {}".format(
            input_file, sample_id, caller, lower_limit, upper_limit, gap_and_centromere, simple_repeat))

    print("#" * 72, "\n")

    output_file = os.path.abspath(input_file).rsplit("/", 1)[
                      0] + "/" + sample_id + "_" + caller + "_duphold_sv_unfiltered.tsv"

    if caller == "manta":
        manta_vcf_to_bed(input_file, output_file, sample_id, lower_limit, upper_limit)
        filter_gap_centromere_simple_repeats(output_file, gap_and_centromere, simple_repeat)
    elif caller == "delly":
        delly_vcf_to_bed(input_file, output_file, sample_id, lower_limit, upper_limit)
        filter_gap_centromere_simple_repeats(output_file, gap_and_centromere, simple_repeat)
    elif caller == "lumpy":
        lumpy_vcf_to_bed(input_file, output_file, sample_id, lower_limit, upper_limit)
        filter_gap_centromere_simple_repeats(output_file, gap_and_centromere, simple_repeat)
    elif caller == "pbsv":
        pbsv_vcf_to_bed(input_file, output_file, sample_id, lower_limit, upper_limit)
        filter_gap_centromere_simple_repeats(output_file, gap_and_centromere, simple_repeat)
    elif caller == "sniffles":
        sniffles_vcf_to_bed(input_file, output_file, sample_id, lower_limit, upper_limit)
        filter_gap_centromere_simple_repeats(output_file, gap_and_centromere, simple_repeat)
    elif caller == "svim":
        svim_vcf_to_bed(input_file, output_file, sample_id, lower_limit, upper_limit)
        filter_gap_centromere_simple_repeats(output_file, gap_and_centromere, simple_repeat)

    print("Conversion completed....")


if __name__ == "__main__":

    import os, sys, subprocess
    import sys

    main()

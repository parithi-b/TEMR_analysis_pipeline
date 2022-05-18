"""
This progrom is used to convert vcf to bed file for the following tools
Manta/Delly/Lumpy

Make sure Duphold is run before using this function

make sure the vcf file names are in proper structure
the bed file just adds "_with_support.bed" to the vcf file name

vcf_file name = sampleID_toolName_duphold.vcf

conversion for each tool is called inside a seperate function

BED fILE contains

CHR START END SVTYPE CALLER SAMPLEID PAIRED-READ(PR) SPLIT-READ(SR) DHFFC DHBFC

"""


def sv_type_categorizing(sv_bed_db):
    # print()
    print("SV type categorizing")

    sv_type_db = dict()

    for i in ["DEL", "DUP", "INV"]:
        sv_type_db[i] = 0

    for i in sv_bed_db:
        sv_data = i.split("_")
        pos = sv_data[1]
        end = sv_data[2]
        given_sv = sv_data[3]

        if abs(int(end) - int(pos)) >= 50:
            sv_type_db[given_sv] += 1

    for i in sv_type_db:
        print(i, sv_type_db[i])


def splitting_output_files_based_on_svtype(output_file):
    print("FILE NAME : ", output_file)

    del_bed_write = open(output_file.split(".")[0] + "_DEL.bed", 'w')
    dup_bed_write = open(output_file.split(".")[0] + "_DUP.bed", 'w')
    inv_bed_write = open(output_file.split(".")[0] + "_INV.bed", 'w')

    cat_stage = subprocess.Popen(['cat', output_file], stdout=subprocess.PIPE)
    subprocess.call(['grep', 'DEL'], stdin=cat_stage.stdout, stdout=del_bed_write)
    cat_stage = subprocess.Popen(['cat', output_file], stdout=subprocess.PIPE)
    subprocess.call(['grep', 'DUP'], stdin=cat_stage.stdout, stdout=dup_bed_write)
    cat_stage = subprocess.Popen(['cat', output_file], stdout=subprocess.PIPE)
    subprocess.call(['grep', 'INV'], stdin=cat_stage.stdout, stdout=inv_bed_write)


def delly_vcf_to_bed(input_bed_file, sample_id, write_output_file):

    print("Converting VCF to BED  -  DELLY ")

    if write_output_file is True:
        bed_write = open(input_bed_file.split(".")[0] + "_with_support.bed", 'w')

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "#" in sv:
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

                if abs(int(end_pos) - int(start_pos)-1) >= 50:
                    sv_count += 1
                    if write_output_file is True:
                        bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("Total SV Considered   = {} ".format(sv_count))
    print("Delly")
    sv_type_categorizing(sv_db)

    if write_output_file is True:
        bed_write.close()
        # splitting_output_files_based_on_svtype(bed_write.name)

    print("*" * 50)


def lumpy_vcf_to_bed(input_bed_file, sample_id, write_output_file):

    print("Converting VCF to BED  -  LUMPY ")

    if write_output_file is True:
        bed_write = open(input_bed_file.split(".")[0] + "_with_support.bed", 'w')  # for duphold

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "#" in sv:
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
                if abs(int(end_pos) - int(start_pos) -1) >= 50:
                    sv_count += 1
                    if write_output_file is True:
                        bed_write.write("\t".join(i for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("Total SV Considered   = {} ".format(sv_count))

    sv_type_categorizing(sv_db)

    if write_output_file is True:
        bed_write.close()
        # splitting_output_files_based_on_svtype(bed_write.name)

    print("*" * 50)


def manta_vcf_to_bed(input_bed_file, sample_id, write_output_file):

    print("Converting VCF to BED  -  MANTA ")

    if write_output_file is True:
        bed_write = open(input_bed_file.split(".")[0] + "_with_support.bed", 'w')  # for duphold

    sv_count = 0
    total_sv = 0
    sv_db = []

    for sv in open(input_bed_file, 'r'):

        if "#" in sv:
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

                if abs(int(end_pos) - int(start_pos) -1) >= 50:
                    sv_count += 1

                    if write_output_file is True:
                        bed_write.write("\t".join(str(i) for i in temp) + "\n")

    print("Total Variants called = {} ".format(total_sv))
    print("Total SV Considered   = {} ".format(sv_count))

    sv_type_categorizing(sv_db)

    if write_output_file is True:
        bed_write.close()
        # splitting_output_files_based_on_svtype(bed_write.name)

    print("*" * 50)


def chaisson_dataset(trio_id, writing_output_flag):

    print("Converting VCF to BED  -  Chaisson ")

    ### only vcf is in the folder else will throw error # fix: if file already exists delete it

    write_output_file = writing_output_flag  # set true if u want to write your output in a bed file

    # duphold

    vcf_path = "/Users/balacp/Desktop/1000GP_WGS/HGSV/PacBio/Duphold/" + trio_id + "_1000GP_duphold.vcf"

    if write_output_file is True:
        bed_write = open(vcf_path.split(".")[0] + "_with_support.bed", 'w')

    total_sv = 0
    sv_db = []
    sv_50 = 0

    for sv in open(vcf_path, 'r'):

        if "#" in sv:
            continue

            # if "SVtpe" in sv or "<DEL>" in sv or "<DUP>" in sv:

        if "<DEL>" in sv:
            total_sv += 1

            y = sv.split()

            chr_start = y[0]
            start_pos = y[1]

            info_vcf = y[7]

            for sub_info in info_vcf.split(";"):
                if "SVLEN=" in sub_info[:6]:
                    sv_len = (sub_info.split("=")[-1])
                    if sv_len.isalpha():
                        sv_len = 0
                    else:
                        sv_len = int(sv_len)
                elif "SVTYPE=" in sub_info[:7]:
                    sv_type = (sub_info.split("=")[-1])
                elif "END=" in sub_info[:4]:
                    end_pos = (sub_info.split("=")[-1])
                elif "SVCLASS=" in sub_info[:8]:
                    class_support = sub_info.split("=")[-1]
                elif "UNION=" in sub_info[:6]:
                    union_support = sub_info.split("=")[-1]

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

            ayz = chr_start + "_" + start_pos + "_" + end_pos + "_" + sv_type

            # print(y)
            # print(info_vcf)

            temp = list()
            temp.append(chr_start)
            temp.append(start_pos)
            temp.append(end_pos)
            temp.append(sv_type)
            temp.append("hgsv")
            temp.append(trio_id)
            temp.append(class_support)
            temp.append(union_support)
            temp.append(dhffc_value)
            temp.append(dhbfc_value)

            print(temp)

            # if int(pe_support) < 5:
            #     continue

            if ayz not in sv_db:
                sv_db.append(ayz)
                if abs(sv_len) >= 50:
                    sv_50 += 1
                    if write_output_file is True:
                        bed_write.write("\t".join(i for i in temp) + "\n")
    print("Total Variants                    = {} ".format(total_sv))
    print("Total SV Considered               = {} ".format(sv_50))

    if write_output_file is True:
        bed_write.close()
        # h = remove_calls_from_gaps_centromeres(bed_after_filter_write.name)
        # splitting_output_files_based_on_svtype(bed_write.name)

    print("*" * 50)


def main():

    write_output_flag = True  # do not put True or False under quotes

    # trio_id_db = ["HG00731", "HG00732", "HG00733","HG00512", "HG00513", "HG00514", "NA19239", "NA19238", "NA19240"]
    #
    # trio_id_db = ["NA19240"]
    trio_id_db = ["HG00733", "HG00514", "NA19240"]

    # for trio_id in trio_id_db:
    #     for caller_used in ["delly", "lumpy", "manta"]:
    #
    #         print(trio_id)
    #
    #         if caller_used == "delly":
    #             folder_path = "/Users/balacp/Desktop/1000GP_WGS/HGSV/WGS/Delly/" + trio_id + "/" + trio_id + "_delly_duphold.vcf"
    #             delly_vcf_to_bed(folder_path, trio_id, write_output_flag)
    #         elif caller_used == "lumpy":
    #             folder_path = "/Users/balacp/Desktop/1000GP_WGS/HGSV/WGS/Lumpy/" + trio_id + "/" + trio_id + "_lumpy_duphold.vcf"
    #             lumpy_vcf_to_bed(folder_path, trio_id, write_output_flag)
    #         elif caller_used == "manta":
    #             folder_path = "/Users/balacp/Desktop/1000GP_WGS/HGSV/WGS/Manta/" + trio_id + "/" + trio_id + "_manta_duphold.vcf"
    #             manta_vcf_to_bed(folder_path, trio_id, write_output_flag)

    msi_path = "/Users/balacp/Desktop/1000GP_WGS/HGSV/WGS/MSI/New_data_April_2021/MSH2KO/"
    for caller in os.listdir(msi_path):
        print(caller)
        if caller == "Manta" :
            for i in os.listdir(msi_path + "/Manta/"):
                trio_id = i.split("_")[0] + "-" + i.split("_")[1]
                if ".vcf" in i:
                    folder_path = msi_path + "/Manta/" + i
                    print(trio_id, folder_path)
                    manta_vcf_to_bed(folder_path, trio_id, write_output_flag)
        elif caller == "Lumpy":
            for i in os.listdir(msi_path + "/Lumpy/"):
                trio_id = i.split("_")[0] + "-" + i.split("_")[1]
                if ".vcf" in i:
                    folder_path = msi_path + "/Lumpy/" + i
                    print(trio_id, folder_path)
                    lumpy_vcf_to_bed(folder_path, trio_id, write_output_flag)
        elif caller == "Delly":
            for i in os.listdir(msi_path + "/Delly/"):
                trio_id = i.split("_")[0] + "-" + i.split("_")[1]
                if ".vcf" in i:
                    folder_path = msi_path + "/Delly/" + i
                    print(trio_id, folder_path)
                    delly_vcf_to_bed(folder_path, trio_id, write_output_flag)

    # chaisson_dataset(Trio_ID,write_output_flag)


if __name__ == "__main__":

    import subprocess,os

    main()

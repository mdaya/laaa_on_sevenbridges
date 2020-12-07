import sys

allele_file_name=sys.argv[1]
outfile_prefix="plink_in"

allele_file=open(allele_file_name)
tfam_file=open(outfile_prefix + ".tfam", "w")
tped_file=open(outfile_prefix + ".tped", "w")
head_line=allele_file.readline()
head_el=head_line.strip().split()
for i in range(3, len(head_el), 1):
    tfam_file.write(head_el[i] + " " + head_el[i] + " 0 0 0 0\n")
for line in allele_file:
    el=line.strip().split()
    tped_file.write("1 " + el[0] + ":" + el[1] + ":" + el[2] + " 0 " + el[0])
    for j in range(3, len(el), 1):
        dose = el[j]
        if dose == "0":
            tped_file.write(" 1 1")
        if dose == "1":
            tped_file.write(" 1 2")
        if dose == "2":
            tped_file.write(" 2 2")
    tped_file.write("\n")
allele_file.close()
tfam_file.close()
tped_file.close()



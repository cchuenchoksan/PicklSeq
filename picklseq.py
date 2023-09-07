import subprocess
import re
import pickle
import sys
import os

current_file = os.path.abspath(__file__)
current_directory = os.path.dirname(current_file)

# Take arguments
file_name = sys.argv[1]
try:
    seq_type = sys.argv[2]
except IndexError:
    print(f"WARNING: A sequence type was not given. Defaulting to CRT.")
    seq_type = "CRT"

if seq_type == "DHPS":
    fasta_file = f"{current_directory}/utils/DHPS.fasta"
    seq_len = 642
elif seq_type == "DHFR":
    fasta_file = f"{current_directory}/utils/DHFR.fasta"
    seq_len = 491
elif seq_type == "CRT":
    fasta_file = f"{current_directory}/utils/CRT.fasta"
    seq_len = 178
else:
    print(
        f"WARNING: The sequence type {seq_type} is not recognised by the script. Defaulting to CRT.")
    fasta_file = f"{current_directory}/utils/CRT.fasta"
    seq_len = 178

print("Running subprocess...")

process = subprocess.Popen(['bash', f"{current_directory}/utils/commands.sh", file_name,
                           fasta_file, current_directory], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = process.communicate()

print("#####Process Completed#####")
print("Subprocess Output:")
print(stdout.decode())
print(stderr.decode())

print("Start picklseq trimming...")

cigar_num_pattern = r"\d+"
cigar_alpha_pattern = r"[A-Z]"


def get_match_count(cigar_nums, cigar_alphas):
    output = 0
    for num, alpha in zip(cigar_nums, cigar_alphas):
        if alpha in {"M"}:
            output += int(num)
    return output


data = []
with open(f"{current_directory}/sorted_alignment.txt", "r") as f:
    for line in f:
        line_array = line.split('\t')
        flag = int(line_array[1])
        dna_seq = line_array[9]
        cigar_str = line_array[5]
        pos = int(line_array[3])

        cigar_nums = re.findall(cigar_num_pattern, line_array[5])
        cigar_alphas = re.findall(cigar_alpha_pattern, line_array[5])
        if len(cigar_nums) > 0:
            if cigar_alphas[0] == 'S':
                dna_seq = dna_seq[int(cigar_nums[0]):]
            if len(dna_seq) < seq_len:
                continue

            interested = [line_array[2], cigar_str, dna_seq[:seq_len]]
            if flag in [0, 1, 16, 2048, 2064]:
                if get_match_count(cigar_nums, cigar_alphas) > 0:
                    if pos == 1:  # will need to pad
                        data.append(interested)

os.remove(f"{current_directory}/sorted_alignment.txt")

print("Trimming done")
print("Length of data: ", len(data))

with open("output.pkl", "wb") as file:
    pickle.dump(data, file)

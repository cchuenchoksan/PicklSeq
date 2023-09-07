import subprocess
import re
import pickle
import sys

# Take arguments
file_name = sys.argv[1]
seq_type = sys.argv[2]

if seq_type == "DHPS":
    fasta_file = "./utils/DHPS.fasta"
elif seq_type == "DHFR":
    fasta_file = "./utils/DHFR.fasta"
elif seq_type == "CRT":
    fasta_file = "./utils/CRT.fasta"

print(file_name)
print(fasta_file)

process = subprocess.Popen(['bash', "./utils/commands.sh", file_name, fasta_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = process.communicate()

print("Standard Output:")
print(stdout.decode())
print(stderr.decode())

process.wait()
print("Process Completed ################################################################")

cigar_num_pattern = r"\d+"
cigar_alpha_pattern = r"[A-Z]"
seq_len = 178

def get_match_count(cigar_nums, cigar_alphas):
    output = 0
    for num, alpha in zip(cigar_nums, cigar_alphas):
        if alpha in {"M"}:
            output += int(num)
    return output

data = []
with open('./sorted_example_alignment.txt', 'r') as f:
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

print(len(data))

with open("output.pkl", "wb") as file:
    pickle.dump(data, file)

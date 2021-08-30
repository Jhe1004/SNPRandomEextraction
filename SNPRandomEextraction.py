import os
from Bio import SeqIO
import random

outgroup = "Clematis_ochotensis_2019051801" #外类群的样品名称，仅可有1个外类群
gap = "-" #文件中的gap用什么符号表示
snp_num_per_gene = 1 #每个基因中随机选择几个SNP




#函数get_file_list：获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:
        if ".fasta" in each:
            file_list.append(each)
    return file_list

def get_seq_names(fasta_file_list):
    '''
    遍历所有文件，获取所有序列的名称
    '''
    seq_name = []
    for each_file in fasta_file_list:
        with open(each_file, "r") as read_file:
            for each_line in read_file:
                if ">" in each_line:
                    if each_line.replace("\n", "") not in seq_name:
                        seq_name.append(each_line.replace("\n", ""))
    return seq_name


def make_snp_matrix(each_file, seq_names, snp_list):
    file_dict = SeqIO.to_dict(SeqIO.parse(each_file, "fasta"))
    alignmet_len = len(file_dict[outgroup].seq)
    temp_list = []
    for each_index in range(alignmet_len):
        judge_list = []
        for each_seq in seq_names:
            if each_seq[1:] != outgroup:
                judge_list.append(str(file_dict[each_seq[1:]]
                                 .seq[each_index:each_index+1]))
            else:
                judge_list.insert(0, str(file_dict[each_seq[1:]]
                                 .seq[each_index:each_index+1]))
        judge_list_raw = judge_list.copy()
        judge_list2 = judge_list[1:]
        judge_list2.sort()
        if gap not in judge_list:
            if judge_list2[0] != judge_list2[-1]:
                temp_list.append(judge_list_raw)
    if len(temp_list) >= snp_num_per_gene:
        for each_list in random.sample(temp_list,snp_num_per_gene):
            snp_list.append(each_list)
    


def main():
    fasta_file_list = get_file_list() 
    seq_names = get_seq_names(fasta_file_list)
    snp_list = []
    for each_file in fasta_file_list:
        with open(each_file, "r") as read_file:
            n = 0
            for each_line in read_file:
                if ">" in each_line:
                    n = n + 1
            if n == len(seq_names):
                make_snp_matrix(each_file, seq_names, snp_list)
    with open("result.fasta", "a") as write_file:
        n = 0
        for each_name in seq_names:
            write_file.write(">" + each_name + "\n")
            for each in snp_list:
                write_file.write(each[n])
            write_file.write("\n")
            n = n + 1
            




if __name__ == "__main__":
    main()

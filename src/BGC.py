import os
from Bio import SeqIO

class BGC:
    '''
    It takes a BGC fasta file and translate it into proteins in six open reading frames.
    '''
    def __init__(self, mode, out_files_path, BGC_file_path, gene_len_thresh=300):
        self.BGC_file_path = BGC_file_path
        self.BGC_name = os.path.basename(self.BGC_file_path).split('.')[0]
        self.BGC_DNA_fasta = self.read_sequence()
        self.six_orf = self.find_six_open_reading_frame()
        self.BGC_gene = self.orf_gene_screening(gene_len_thresh)
        self.BGC_protein = self.translation()  # generate a protein fasta file_

        if mode == "debug":
            gene_file_name_total = os.path.join(out_files_path, self.BGC_name + "_genes.fasta")
            self.write_sequence_into_file("gene", gene_file_name_total)

        protein_file_name_total = os.path.join(out_files_path,  self.BGC_name + "_translated_protein.fasta")
        self.write_sequence_into_file("protein", protein_file_name_total)
    

    def read_sequence(self):
        # BGC fasta should contain only one sequence
        if len([1 for line in open(self.BGC_file_path) if line.startswith(">")]) != 1:
            raise ValueError("Please select the genome file contained 1 sequence.")

        inputfile = open(self.BGC_file_path)
        record = SeqIO.read(inputfile, "fasta")

        return record

    def find_six_open_reading_frame(self):
        DNA_sequence = self.BGC_DNA_fasta.seq  #this line only works for the fna file contained only 1 squence.
        complement_DNA_sequence = DNA_sequence.reverse_complement()
        reading_frame = {}
        reading_frame[1] = DNA_sequence
        reading_frame[2] = DNA_sequence[1:]
        reading_frame[3] = DNA_sequence[2:]
        reading_frame[4] = complement_DNA_sequence
        reading_frame[5] = complement_DNA_sequence[1:]
        reading_frame[6] = complement_DNA_sequence[2:]

        return reading_frame

    def orf_gene_screening(self, gene_len_thresh):
        start_codon = "ATG"
        stop_codon = ["TAA", "TAG", "TGA"]
        genes = {}
        count = 1
        for orf in self.six_orf.keys():
            current_sequence = self.six_orf[orf]
            start_codon_index = []
            stop_codon_index = []
            for i in range(0, len(current_sequence) - 3, 3):
                if current_sequence[i:i + 3] == start_codon:
                    start_codon_index.append(i)
                if current_sequence[i:i + 3] in stop_codon:
                    stop_codon_index.append(i)

            for i in start_codon_index:
                for j in stop_codon_index:
                    if j > i: # to make sure we eliminate shorter sequences.
                        if len(current_sequence[i:j+3]) > gene_len_thresh:
                            genes[count] = {"sequence": current_sequence[i:j + 3], "orf": orf, "start position": str(i),
                                            "stop position": str(j + 2)}
                            if orf <= 3:
                                genes[count]["strand"] = "+"
                                count += 1
                            else:
                                genes[count]["strand"] = "-"
                                count += 1
                            break
                        else:
                            break
        return genes

    def write_sequence_into_file(self, sequence_type, fileName):
        with open(fileName, "w") as fw:
            for i in self.BGC_gene.keys():
                if sequence_type == "gene":
                    fw.write(
                        ">" + str(i) + " " + self.BGC_gene[i]["strand"] + " " + str(self.BGC_gene[i]["orf"]) + " " +
                        self.BGC_gene[i]["start position"] + " " + self.BGC_gene[i][
                            "stop position"] + "\n" + str(self.BGC_gene[i]["sequence"]) + "\n")
                elif sequence_type == "protein":
                    fw.write(
                        ">" + str(i) + " " + self.BGC_gene[i]["strand"] + " " + str(self.BGC_gene[i]["orf"]) + " " +
                        self.BGC_gene[i]["start position"] + " " + self.BGC_gene[i][
                            "stop position"] + "\n" + str(self.BGC_protein[i]) + "\n")

    def translation(self):
        protein_frame = {}
        for i in self.BGC_gene.keys():
            protein_frame[i] = str(self.BGC_gene[i]["sequence"].translate())
        return protein_frame
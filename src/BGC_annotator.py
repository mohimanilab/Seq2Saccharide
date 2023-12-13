import shutil
import pandas as pd
import os
from Bio import SearchIO


class BGC_annotator:
    '''
    This class uses HMMer to annotate monomer synthesis genes and modification genes in the BGC
    using the pre-built HMM files.
    '''

    def __init__(self, primary_metabolite, BGC_name, mode, primary_mod_file, output_folder="", 
                 monomer_hmm_cutoff="../data/monomer_hmm_cutoff.csv",
                 modification_hmm_cutoff="../data/modification_hmm_cutoff.csv",
                 path_monomer_hmm_folder="../data/HmmModel_version4_with_predicted_threshold",
                 path_modification_hmm_folder="../data/Hmm_modification_version5_with_predicted_threshold"
                 ):
        self.BGC = BGC_name
        self.out_put_folder = output_folder
        self.protseq_file_name = os.path.join(output_folder, BGC_name + "_translated_protein.fasta")
        self.primary_metabolite = primary_metabolite
        self.monomer_hmm = self.read_monomer_hmm(monomer_hmm_cutoff)
        self.monomer_hmm_list = self.get_monomer_hmm_list()
        self.modification_hmm, self.mod_hmm = self.read_modification_hmm(modification_hmm_cutoff)
        self.primary_mod_file = primary_mod_file
        self.primary_mod = self.read_primary_mod()
        self.threshold_dict_monomer = self.create_threshold_dict(monomer_hmm_cutoff)
        self.threshold_dict_mod = self.create_threshold_dict(modification_hmm_cutoff)
        self.matched_monomer_gene, self.result_monomer_path = self.hmm_search("monomer", path_monomer_hmm_folder)
        self.matched_modificaion_gene, self.result_modification_path = self.hmm_search("modification", path_modification_hmm_folder)
        self.generated_monomers = self.find_monomers()  # return a dictionary, where the key is the monomer name and the value is hmm used for pathways.
        print("matched monomer hmm:", self.matched_monomer_gene)
        print("modification hmm: ", self.matched_modificaion_gene)
        print("generated monomers: ", list(self.generated_monomers.keys()))

        if mode == "run":
            shutil.rmtree(self.result_monomer_path)
            shutil.rmtree(self.result_modification_path)

    def read_primary_mod(self):
      primary_mod = []
      with open(self.primary_mod_file) as f:
        lines =f.readlines()
        for line in lines:
          primary_mod.append(line.strip("\n"))
        f.close()

      return primary_mod

    def read_monomer_hmm(self, monomer_hmm_cutoff):
        hmm_dict = {}
        monomer_hmm_sheet = pd.read_csv(monomer_hmm_cutoff)
        for i in range(0, len(monomer_hmm_sheet["monomer"])):
            if monomer_hmm_sheet["monomer"][i] not in self.primary_metabolite:
                if monomer_hmm_sheet["monomer"][i] not in hmm_dict.keys():
                    hmm_dict[monomer_hmm_sheet["monomer"][i]] = {1: [monomer_hmm_sheet["hmm_name"][i]]}
                else:
                    if monomer_hmm_sheet["pathway"][i] in hmm_dict[monomer_hmm_sheet["monomer"][i]].keys():
                        hmm_dict[monomer_hmm_sheet["monomer"][i]][monomer_hmm_sheet["pathway"][i]].append(
                            monomer_hmm_sheet["hmm_name"][i])
                    else:
                        hmm_dict[monomer_hmm_sheet["monomer"][i]][monomer_hmm_sheet["pathway"][i]] = [
                            monomer_hmm_sheet["hmm_name"][i]]

        return hmm_dict

    def get_monomer_hmm_list(self):
        hmmlist = []
        for monomer in self.monomer_hmm.keys():
            sub_dict = self.monomer_hmm[monomer]
            for sub_list in sub_dict.values():
                for hmmname in sub_list:
                    if hmmname not in hmmlist:
                        hmmlist.append(hmmname)
        return hmmlist

    def read_modification_hmm(self, modification_hmm_cutoff):
        modification_hmm = pd.read_csv(modification_hmm_cutoff)
        modification_hmm_list = []
        mod_hmm = {}
        for i in range(len(modification_hmm["hmm_name"])):
            modification_hmm_list.append(modification_hmm["hmm_name"][i])
            mod_hmm[modification_hmm["hmm_name"][i]] = modification_hmm["modification_name"][i]

        return modification_hmm_list, mod_hmm

    def hmm_search(self, hmm_type, folder_path):
        matched_monomer_hmm = {}
        hmm_list = []
        if hmm_type == "monomer":
            hmm_list = self.monomer_hmm_list
        elif hmm_type == "modification":
            hmm_list = self.modification_hmm
        for hmm in hmm_list:
            search_result_accession = []
            search_result_score = []
            if_threshold = False
            if hmm_type == "monomer":
                if_threshold = self.check_if_threshold(str(hmm), "monomer")
            elif hmm_type == "modification":
                if_threshold = self.check_if_threshold(str(hmm), "modification")

            path = os.path.join(self.out_put_folder, self.BGC + "hmm_search_" + hmm_type + "result")
            isExists=os.path.exists(path)
            if not isExists:
                os.makedirs(path)

            output_file_name = os.path.join(path, self.BGC + "_" + str(hmm) + ".txt")

            if if_threshold:
                os.system("/home/yanjing2/aminoglycoside/hmmer/src/./hmmsearch --cut_ga " + os.path.join(folder_path, str(
                    hmm)) + " " + self.protseq_file_name + " > " + output_file_name)
            else:
                os.system(
                    "/home/yanjing2/aminoglycoside/hmmer/src/./hmmsearch " + os.path.join(folder_path, str(hmm)) + " " + self.protseq_file_name + " > " + output_file_name)

            handle = open(output_file_name)
            qresults = SearchIO.parse(handle, "hmmer3-text")

            for qresult in qresults:
                for Hit in qresult:
                    if Hit.bitscore > 0:
                        search_result_accession.append(qresult.id)
                        search_result_score.append(Hit.bitscore)
                        break

            if len(search_result_accession) != 0:
                matched_monomer_hmm[hmm] = search_result_accession[search_result_score.index(max(search_result_score))]

            handle.close()

        return matched_monomer_hmm, path

    def find_monomers(self):
        monomer_dict = {}
        for hmm in self.matched_monomer_gene:
            for monomer in self.monomer_hmm.keys():
                for hmmvalues in self.monomer_hmm[monomer].values():
                    if hmm in hmmvalues:
                        if monomer not in monomer_dict.keys():
                            for pathway_name in self.monomer_hmm[monomer].keys():
                                if hmmvalues == self.monomer_hmm[monomer][pathway_name]:
                                    monomer_dict[monomer] = {pathway_name: hmmvalues}
                        else:
                            if hmmvalues not in monomer_dict[monomer].values():
                                for pathway_name in self.monomer_hmm[monomer].keys():
                                    if hmmvalues == self.monomer_hmm[monomer][pathway_name]:
                                        monomer_dict[monomer][pathway_name] = hmmvalues

        return monomer_dict

    def create_threshold_dict(self, hmm_file):
        hmm_sheet = pd.read_csv(hmm_file)
        threshold_dict = {}
        for i in range(0, len(hmm_sheet["hmm_name"])):
            if hmm_sheet["hmm_name"][i] not in threshold_dict.keys():
                threshold_dict[hmm_sheet["hmm_name"][i]] = hmm_sheet["cutoff"][i]
        return threshold_dict

    def check_if_threshold(self, hmm, hmm_type):

        threshold_dict = None

        if hmm_type == "monomer":
            threshold_dict = self.threshold_dict_monomer
        elif hmm_type == "modification":
            threshold_dict = self.threshold_dict_mod

        if threshold_dict[hmm] != 0:
            return True
        else:
            return False
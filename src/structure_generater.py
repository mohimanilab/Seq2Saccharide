import os
import json
import itertools
import copy
from scipy.stats import fisher_exact
from rdkit import Chem
from rdkit.Chem import Descriptors
import csv
import pandas as pd

from BGC import BGC
from BGC_annotator import BGC_annotator

class structure_generater:
  '''
  This is the main class of the tool.
  It contains three main parts:
  (1) Select most confident monomer sets by Fisher or MLE
  (2) Translate the monomer sets into chemical-valid saccharide bakcbones in the dictionary format
  (3) Translate the dictionary into product SMILES
  '''
  def __init__(self, mode, BGC_file, model, monomer_hmm_cutoff, modification_hmm_cutoff, path_monomer_hmm_folder, path_modification_hmm_folder, out_files, connectivity, train_list, train_BGC_dir, BGC_backbone, bond_list, mod_sif, primary_mod_file, length, mod_depth, gene_length, bond, primary_opt, set_num):
    self.BGCfile = BGC_file
    self.primary_metabolite = ["D-ribose", "D-mannose", "D-glucose", "D-glucosamine", "N-acetylglucosamine", "D-xylose"]
    self.out_files = out_files
    self.model = model
    self.connectivity = connectivity
    self.train_list = train_list
    self.train_BGC_dir = train_BGC_dir
    self.BGC_backbone = BGC_backbone
    self.bond_list = bond_list
    self.mod_sif = mod_sif
    self.length = length
    self.gene_length = gene_length
    self.mod_depth = mod_depth
    self.bond = bond
    self.primary_opt = primary_opt
    self.set_num = set_num
    self.primary_mod_file = primary_mod_file
    self.BGC = BGC(mode, out_files, BGC_file, gene_length)
    self.annotator = BGC_annotator(self.primary_metabolite, self.BGC.BGC_name, mode, primary_mod_file, out_files, monomer_hmm_cutoff, modification_hmm_cutoff, path_monomer_hmm_folder, path_modification_hmm_folder)

    self.all_monomer_set = self.monomer_sets()  #a list containing all monomer sets
    self.monomerbase = self.read_monomer()
    self.residuebase = self.read_residue()
    
    if model == "Fisher":
      self.selected_monomer_sets = self.selected_monomer_sets_fisher() #top 100 monomer sets dictionary
      with open(os.path.join(self.out_files, os.path.basename(self.BGCfile) + "Fisher_selected_monomer_sets.txt"), "w") as f:
        for i in self.selected_monomer_sets:
          f.write(str(i) + "\n")
        f.close()


    elif model == "MLE":
      self.backbone_dict = self.read_known_backbone()
      self.parameter = self.parameter_estimation(monomer_hmm_cutoff, modification_hmm_cutoff, path_monomer_hmm_folder, path_modification_hmm_folder)
      self.BGC_pattern = self.BGCpattern(self.annotator.matched_monomer_gene, self.annotator.monomer_hmm_list, self.primary_metabolite)
      self.selected_monomer_sets = self.selected_monomer_sets_mle()
      with open(os.path.join(self.out_files, os.path.basename(self.BGCfile) + "MLE_selected_monomer_sets.txt"), "w") as f:
        for i in self.selected_monomer_sets:
          f.write(str(i) + "\n")
        f.close()


    self.backbone = self.generate_backbone() #backbone dictionary
    self.smiles, count = self.backbone_to_smiles() #generate backbone smiles
    self.products = self.add_modification() #generate final products smiles
    self.generate_chemical_database()



  def read_monomer(self):
    sugarDBsheet = pd.read_csv(self.connectivity)

    monomerbase = {}
    for i in range(0,len(sugarDBsheet["monomer"])):
        monomer = sugarDBsheet["monomer"][i]
        residuename = sugarDBsheet["residue"][i]
        connectionpoint = json.loads(str(sugarDBsheet["connection_point"][i]))
        residuedict = {"residuesmiles": sugarDBsheet["smiles"][i], "category": sugarDBsheet["category"][i], "connectionpoint": connectionpoint}

        if monomer in monomerbase.keys():
          monomerbase[monomer][residuename] = residuedict
        else:
          monomerbase[monomer] = {}
          monomerbase[monomer][residuename] = residuedict
          
    return monomerbase

  
  def read_residue(self):
    sugarDBsheet = pd.read_csv(self.connectivity)

    residuebase = {}
    for i in range(0, len(sugarDBsheet["residue"])):
      residue = str(sugarDBsheet["residue"][i])
      residuesmiles = str(sugarDBsheet["smiles"][i])
      connectionpint = json.loads(str(sugarDBsheet["connection_point"][i]))
      residuebase[residue] = {"smiles": residuesmiles, "connectionpoint": connectionpint}

    return residuebase


  # Generate all possible monomer sets
  def monomer_sets(self):
    monomer_dict = self.annotator.generated_monomers
    monomer = list(monomer_dict.keys())

    monomer_sets = []
    min_length = min(len(monomer), self.length)
    for i in range(1, min_length+1):
      comb = list(itertools.combinations(monomer, i))
      monomer_sets += comb   #[(),(),()] tuple in a list

    return monomer_sets

  
  def selected_monomer_sets_fisher(self):
    all_gene = self.BGC.BGC_gene  #gene dictionary
    all_matched_gene = self.annotator.matched_monomer_gene  #list of matched hmm
    needed_genes = self.annotator.generated_monomers #return a dictionary, where the key is the monomer name and the value is possibly used hmm used in different pathways.
    
    #For each monomer set, find all possible hmm combinations (for different pathways)
    p_value_summary = {}
    for monomer_set in self.all_monomer_set:
      needed_hmm = []
      for monomer in monomer_set:
        needed_hmm.append(list(needed_genes[monomer].values()))
      possible_hmm_sets = list(itertools.product(*needed_hmm))

      concatenated_hmm_sets = []
      for hmm_set in possible_hmm_sets:
        concatenation = []
        for sub_set in hmm_set:
            concatenation += sub_set
        concatenated_hmm_sets.append(set(concatenation))
        
      #fisher exact test
      fisher_p_value = []
      for hmm_set in concatenated_hmm_sets:
        all_gene_number = len(self.annotator.monomer_hmm_list)
        needed_gene_number = len(hmm_set)
        needed_and_annotated_number = len(set(hmm_set) & set(all_matched_gene))
        annotated_number = len(all_matched_gene)
        table = [[needed_and_annotated_number, needed_gene_number - needed_and_annotated_number], [annotated_number - needed_and_annotated_number, all_gene_number - needed_gene_number - annotated_number + needed_and_annotated_number]]
        oddsr, p = fisher_exact(table, alternative='greater')
        fisher_p_value.append(p)
      
      p_value_summary[tuple(monomer_set)] = max(fisher_p_value)

    sorted_set = sorted(p_value_summary.items(), key = lambda kv:(kv[1], kv[0]))
    
    ranked_set = [sorted_set[i][0] for i in range(len(sorted_set))]

    with open(os.path.join(self.out_files, str(self.mod_depth) + os.path.basename(self.BGCfile) + self.model + "_monomer_set_rank.csv"), "w") as f:
      for i in range(len(sorted_set)):
        f.write("{},{}\n".format(sorted_set[i][0],sorted_set[i][1]))
      f.close()

    selected_set = {}
    if len(sorted_set) <= self.set_num:
      for i in range(len(sorted_set)):
        selected_set[sorted_set[i][0]] = sorted_set[i][1]
    else:
      for i in range(len(sorted_set)):
        if i <= self.set_num-1:
          selected_set[sorted_set[i][0]] = sorted_set[i][1]
        elif sorted_set[i][1] == sorted_set[self.set_num-1][1]:
          selected_set[sorted_set[i][0]] = sorted_set[i][1]

    primary_metabolite_combinations = []
    final_monomer_set = {}
    for i in self.primary_metabolite:
      primary_metabolite_combinations.append([i])
    for i in range(2, len(self.primary_metabolite)+1):
      temp_combination = list(itertools.combinations(self.primary_metabolite, i))
      for j in temp_combination:
        primary_metabolite_combinations.append(list(j))
    primary_metabolite_combinations.append([])

    
    for i in selected_set.keys():
      for j in primary_metabolite_combinations:
        if len(i) + len(j) <= self.length:
          final_monomer_set[tuple(list(i) + list(j))] = selected_set[i]

    return final_monomer_set


  # Estimate the MLE parameters
  def parameter_estimation(self, monomer_hmm_cutoff, modification_hmm_cutoff, path_monomer_hmm_folder, path_modification_hmm_folder):
    filelist = pd.read_csv(self.train_list)["BGC_list"]
    parameter_list = []

    for filename in filelist:
        train_BGC = BGC("run", self.out_files, os.path.join(self.train_BGC_dir, filename), self.gene_length)
        train_BGC_annotator = BGC_annotator(self.primary_metabolite, train_BGC.BGC_name, "run", self.primary_mod_file, self.out_files, monomer_hmm_cutoff, modification_hmm_cutoff, path_monomer_hmm_folder, path_modification_hmm_folder)
        partitioned_filename = filename.partition(".")
        BGC_id = partitioned_filename[0]
        matched_hmm = list(train_BGC_annotator.matched_monomer_gene.keys())
        BGC_pattern = self.BGCpattern(matched_hmm, self.annotator.monomer_hmm_list, self.primary_metabolite)
        backbone_node = self.backbone_dict[BGC_id]

        exist_monomer_dict = {}
        for element in backbone_node:
          for node in element:
            if node in self.annotator.monomer_hmm.keys() and node not in exist_monomer_dict.keys():  #check here
              exist_monomer_dict[node] = self.annotator.monomer_hmm[node]

        product = []
        for element in backbone_node:
          sub_product = []
          for monomer in element:
            if monomer in self.annotator.monomer_hmm.keys() or monomer in self.primary_metabolite:  #ensure that the monomer is either in the database or is primary metabolite
                sub_product.append(monomer)
          
          if sub_product not in product:
              product.append(sub_product)
        
        product_pattern = self.ProductPattern(product, train_BGC_annotator.monomer_hmm_list, exist_monomer_dict)

        parameter_list = self.number_count(BGC_pattern, product_pattern, parameter_list)

    accu_parameter_list = []

    for i in parameter_list:
      accu_parameter_list.append(sorted(i)[0])

    alpha_number = 0
    non_alpha_number = 0
    beta_number = 0
    non_beta_number = 0
    primary_metabolite_exist_number = 0
    primary_metabolite_non_exist_number = 0
    for i in accu_parameter_list:
      alpha_number += i[0]
      non_alpha_number += i[3]
      beta_number += i[1]
      non_beta_number += i[4]
      primary_metabolite_exist_number += i[2]
      primary_metabolite_non_exist_number += i[5]

    alpha = float(alpha_number / (alpha_number + non_alpha_number))
    beta = float(beta_number / (beta_number + non_beta_number))
    gamma = float(primary_metabolite_exist_number / (primary_metabolite_exist_number + primary_metabolite_non_exist_number))
    
    return [alpha, beta, gamma]
  
  # read training data
  def read_known_backbone(self):
    backbone_dict = {}
    backbone_sheet = pd.read_csv(self.BGC_backbone)
    for i in range(0, len(backbone_sheet["BGC"])):
        monomerlist = []
        converted_dict = json.loads(backbone_sheet["backbone_graph"][i])
        if converted_dict != {}:
          nodelist = list(converted_dict["node"].values())
          for node in nodelist:
            partitioned_node = node.partition("_")
            monomerlist.append(partitioned_node[0])
          monomerlist = list(set(monomerlist))
          if backbone_sheet["BGC"][i] not in backbone_dict.keys():
            backbone_dict[backbone_sheet["BGC"][i]] = [monomerlist]
          else:
            if monomerlist not in backbone_dict[backbone_sheet["BGC"][i]]:
              backbone_dict[backbone_sheet["BGC"][i]].append(monomerlist)

    return backbone_dict

  # Check if required gene is annotated in the BGC
  def BGCpattern(self, matched_hmm, hmmlist, primary_metabolite):
    BGC_pattern = []
    for gene in hmmlist:
        if gene in matched_hmm:
            BGC_pattern.append(1)
        else:
            BGC_pattern.append(0)

    for i in range(0, len(primary_metabolite)):
        BGC_pattern.append(2)

    return BGC_pattern

  # Check if genes required by the BGC product are annotated in the BGC
  def ProductPattern(self, product, hmmlist, monomer_dict):
    required_backbone_gene = self.required_by_backbone(product, monomer_dict)
    product_pattern = []

    for i in range(0, len(required_backbone_gene)):
      backbone = required_backbone_gene[i]
      for sub_required_gene in backbone:
        sub_pattern = []
        for gene in hmmlist:
          if gene in sub_required_gene:
              sub_pattern.append(1)
          else:
              sub_pattern.append(0)
        for monomer in self.primary_metabolite:
          if monomer in product[i]:
            sub_pattern.append(2)
          else:
            sub_pattern.append(3)

        product_pattern.append(sub_pattern)    

    return product_pattern

  # Get the genes required by the BGC product
  def required_by_backbone(self, comb, monomer_dict):
    genelist = []
    for sub_comb in comb:
      current_comb_gene = []
      gene_required_by_backbone = []
      backbone_gene = []
      for i in range(0, len(sub_comb)):
        temp = []
        if sub_comb[i] in monomer_dict.keys():
          for j in monomer_dict[sub_comb[i]].values():
            temp.append(j)
          backbone_gene.append(temp)

      if len(sub_comb) == 1:
        gene_required_by_backbone = [backbone_gene[0]] 
      else:
        gene_required_by_backbone = list(itertools.product(*backbone_gene))
      for i in gene_required_by_backbone:
        concate = []
        for j in i:
          concate += j
        concate = list(set(concate))
        current_comb_gene.append(concate)
   
      genelist.append(current_comb_gene)

    return genelist


  def number_count(self, BGC_pattern, product_pattern, parameter_list):
    combine_sub_parameter_list = []
    for sub_pattern in product_pattern:
      sub_parameter_list = [0,0,0,0,0,0]
      for i in range(0, len(BGC_pattern)):
          if BGC_pattern[i] == 1 and sub_pattern[i] == 1:
              sub_parameter_list[0] += 1
          elif BGC_pattern[i] == 0 and sub_pattern[i] == 0:
              sub_parameter_list[1] += 1
          elif BGC_pattern[i] == 2 and sub_pattern[i] == 2:
              sub_parameter_list[2] += 1
          elif BGC_pattern[i] == 1 and sub_pattern[i] == 0:
              sub_parameter_list[3] += 1
          elif BGC_pattern[i] == 0 and sub_pattern[i] == 1:
              sub_parameter_list[4] += 1
          elif BGC_pattern[i] == 2 and sub_pattern[i] == 3:
              sub_parameter_list[5] += 1 

      combine_sub_parameter_list.append(sub_parameter_list)
    parameter_list.append(combine_sub_parameter_list)

    return parameter_list


  def selected_monomer_sets_mle(self):
    final_combinations = []
    primary_metabolite_combinations = []
    for i in self.primary_metabolite:
      primary_metabolite_combinations.append([i])
    for i in range(2, len(self.primary_metabolite)+1):
      temp_combination = list(itertools.combinations(self.primary_metabolite, i))
      for j in temp_combination:
        primary_metabolite_combinations.append(list(j))
    primary_metabolite_combinations.append([])
    
    for i in self.all_monomer_set:
      for j in primary_metabolite_combinations:
        if len(j) <= self.primary_opt:
          final_combinations.append(tuple(list(i) + j))

    monomer_set_p = {}
    for monomer_set in final_combinations:
      parameter_list = []
      product_pattern = self.ProductPattern([monomer_set], self.annotator.monomer_hmm_list, self.annotator.generated_monomers)
      parameter_list = self.number_count(self.BGC_pattern, product_pattern, parameter_list)
      probability = self.probability_calculation(parameter_list)
      monomer_set_p[tuple(monomer_set)] = probability

    sorted_set = sorted(monomer_set_p.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
    
    ranked_set = [sorted_set[i][0] for i in range(len(sorted_set))]

    with open(os.path.join(self.out_files, str(self.mod_depth) + os.path.basename(self.BGCfile) + self.model + "_monomer_set_rank.csv"), "w") as f:
      for i in range(len(sorted_set)):
        f.write("{},{}\n".format(sorted_set[i][0],sorted_set[i][1]))
      f.close()


    selected_set = {}
    if len(sorted_set) <= self.set_num:
      for i in range(len(sorted_set)):
        selected_set[sorted_set[i][0]] = sorted_set[i][1]
    else:
      for i in range(len(sorted_set)):
        if i <= self.set_num-1:
          selected_set[sorted_set[i][0]] = sorted_set[i][1]
        elif sorted_set[i][1] == sorted_set[self.set_num-1][1]:
          selected_set[sorted_set[i][0]] = sorted_set[i][1]
        
    return selected_set
    

  def probability_calculation(self, parameter_list):
    alpha = self.parameter[0]
    beta = self.parameter[1]
    gamma = self.parameter[2]
    probability = 0.0
    for sub_parameter_list in parameter_list[0]:
      probability += ((alpha ** sub_parameter_list[0]) * (beta ** sub_parameter_list[1]) * (gamma ** sub_parameter_list[2]) * ((1 - alpha) ** sub_parameter_list[3]) * ((1 - beta) ** sub_parameter_list[4]) * ((1 - gamma) ** sub_parameter_list[5]))
    
    return probability  




  def generate_backbone(self):
    bond_data = []
    with open(self.bond_list) as f:
      lines = f.readlines()
      for line in lines:
        bond_data.append({line.strip("\n").split(" ")[0], line.strip("\n").split(" ")[1]})

    #generate backbone dictionary
    backbone_dict = {}
    for i in self.selected_monomer_sets.keys():
      all_comb = []
      if len(i) < self.length: 
        all_comb.append(i)
        add_part = []
        possibleproduct = [[]]

        for l in i:
            possibleproduct.append([l])
        possibleproduct = possibleproduct[1:]

        currlength = len(possibleproduct)

        for j in range(1, self.length - len(i)):
            for k in i:   
                for m in range(0, currlength):
                    newsubset = possibleproduct[m].copy()
                    newsubset.append(k)
                    if newsubset not in possibleproduct:
                        possibleproduct.append(newsubset)
            currlength = len(possibleproduct)

        count_dict_list = []
        for l in possibleproduct:
            count_dict = {}
            for n in l:
                if n not in count_dict.keys():
                    count_dict[n] = 1
                else:
                    count_dict[n] += 1
            if count_dict not in count_dict_list:
                add_part.append(l)
                count_dict_list.append(count_dict)

        for l in add_part:
          all_comb.append(list(i) + list(l))
        
        for l in all_comb:
          backbone_perm = list(set(list(itertools.permutations(l))))
          for k in backbone_perm:
            bond_not_in_data = 0
            for index in range(len(k)-1):
              if {k[index], k[index+1]} not in bond_data:
                bond_not_in_data += 1
            if bond_not_in_data <= self.bond:
              backbone_dict[tuple(k)] = self.selected_monomer_sets[i]

    return backbone_dict
  

  def backbone_to_smiles(self):
    backbone_smiles = {}
    count = 0
    if self.model == "Fisher":
        outname = os.path.join(self.out_files, os.path.basename(self.BGCfile) + "Fisher_backbone_dict.txt")
    else:
        outname = os.path.join(self.out_files, os.path.basename(self.BGCfile) + "MLE_backbone_dict.txt")
    with open(outname, "w") as f:     
        for i in self.backbone.keys():
            j = 1
            dictlist = []
            productdict = {}
            self.structure_dictionary_generater(i, j, productdict, dictlist, self.monomerbase)
            f.write(str(i) + "\n")
            count += len(dictlist)

            for k in dictlist:
              backbone_smiles[str(k)] = {"monomers": i, "p": self.backbone[i], "smiles": self.DictToSmiles(k)}
        f.close()
    
    return backbone_smiles, count


  def structure_dictionary_generater(self, product, j, productdict, dictlist, monomerbase):
    if len(product) == 1:
      productresidue = list(monomerbase[product[0]].keys())[0]
      productdict = {"node": {"1": productresidue}, "bond":{}}
      dictlist.append(productdict)
      return dictlist
    else:
      if j == 1: 
        for residue1 in monomerbase[product[0]].keys():
          if product[0] + "_OH" in residue1:
            productdict = {"node": {"1": residue1}, "bond":{"1":["1","OH"]}}
            self.structure_dictionary_generater(product, j+1, productdict, dictlist, monomerbase)
      elif j != len(product):
        monomerpattern = product[j-1] + "_inner"
        for residue2 in monomerbase[product[j-1]].keys():
          if monomerpattern in residue2:
            productdict = copy.deepcopy(productdict)
            if productdict["bond"][str(j-1)][1] == "OH":
              if "OH" not in monomerbase[product[j-1]][residue2]["connectionpoint"].keys():
                productdict["node"][str(j)] = residue2
                if len(productdict["bond"][str(j-1)]) == 2:
                  productdict["bond"][str(j-1)].extend((str(j),"H1"))
                productdict["bond"][str(j)] = [str(j), "H2"]
                self.structure_dictionary_generater(product, j+1, productdict, dictlist, monomerbase)
              else:
                productdict["node"][str(j)] = residue2
                if len(productdict["bond"][str(j-1)]) == 2:
                  productdict["bond"][str(j-1)].extend((str(j),"H"))
                productdict["bond"][str(j)] = [str(j), "OH"]
                self.structure_dictionary_generater(product, j+1, productdict, dictlist, monomerbase)
            else:
              if "OH" in monomerbase[product[j-1]][residue2]["connectionpoint"].keys():
                productdict["node"][str(j)] = residue2
                productdict["bond"][str(j-1)].extend((str(j),"OH"))
                productdict["bond"][str(j)] = [str(j), "H"]
                self.structure_dictionary_generater(product, j+1, productdict, dictlist, monomerbase)
      if j == len(product):
        if productdict["bond"][str(j-1)][1] == "OH":
          monomerpattern = product[j-1] + "_H"
          for residue3 in monomerbase[product[j-1]].keys():
            if monomerpattern in residue3:
              productdict = copy.deepcopy(productdict)
              productdict["node"][str(j)] = residue3
              if len(productdict["bond"][str(j-1)]) == 2:
                productdict["bond"][str(j-1)].extend((str(j), "H"))
              dictlist.append(productdict)
        else:
          monomerpattern = product[j-1] + "_OH"
          for residue3 in monomerbase[product[j-1]].keys():
            if monomerpattern in residue3:
              productdict = copy.deepcopy(productdict)
              productdict["node"][str(j)] = residue3
              if len(productdict["bond"][str(j-1)]) == 2:
                productdict["bond"][str(j-1)].extend((str(j), "OH"))
              dictlist.append(productdict)


  def DictToSmiles(self,backbonedict):
    productsmiles = ""
    if len(backbonedict["node"]) == 1:
        productsmiles = self.residuebase[backbonedict["node"]["1"]]["smiles"]
        for connect in self.residuebase[backbonedict["node"]["1"]]["connectionpoint"].keys():
          if "OH" in connect:
            productsmiles = productsmiles.replace("[*:"+str(self.residuebase[backbonedict["node"]["1"]]["connectionpoint"][connect])+"]", "O")
          if "H" in connect:
            productsmiles = productsmiles.replace("([*:"+str(self.residuebase[backbonedict["node"]["1"]]["connectionpoint"][connect])+"])", "")
    else:
        for i in backbonedict["node"].keys():
            if int(i) == 1:
                    group = backbonedict["bond"][i][1]
                    monomername = backbonedict["node"][i]
                    residuesmiles = self.residuebase[monomername]["smiles"]
                    connectat = self.residuebase[monomername]["connectionpoint"][group]
                    productsmiles = residuesmiles.replace("([*:"+str(connectat)+"])", "%99")
                  
            elif int(i) != len(backbonedict["node"]):
                    monomername = backbonedict["node"][i]
                    firstgroup = backbonedict["bond"][str(int(i)-1)][3]
                    secondgroup = backbonedict["bond"][i][1]
                    connectat1 = self.residuebase[monomername]["connectionpoint"][firstgroup]
                    connectat2 = self.residuebase[monomername]["connectionpoint"][secondgroup]
                    monomersmiles = self.residuebase[monomername]["smiles"]
                    productsmiles = productsmiles + "." + monomersmiles.replace("([*:"+str(connectat1)+"])","%" + str(99-int(i)+2))
                    productsmiles = productsmiles.replace("([*:"+str(connectat2)+"])","%" + str(99-int(i)+1))

            elif int(i) == len(backbonedict["node"]):
                    lastmonomer = backbonedict["node"][i]
                    group = backbonedict["bond"][str(int(i)-1)][3]
                    connectat = self.residuebase[lastmonomer]["connectionpoint"][group]
                    residuesmiles = self.residuebase[lastmonomer]["smiles"]
                    residuesmiles = residuesmiles.replace("([*:"+str(connectat)+"])","%" + str(99-int(i)+2))
                    productsmiles = productsmiles + "." + residuesmiles

    mol = Chem.MolFromSmiles(productsmiles)
    smi = Chem.MolToSmiles(mol)

    return smi


  def add_modification(self):
    primary_mod = self.annotator.primary_mod
    modification_hmm = list(self.annotator.matched_modificaion_gene.keys())
    backbone_smiles = self.smiles
    mod = []
    search_string = ""

    mod_to_add = set([self.annotator.mod_hmm[i] for i in modification_hmm])

    for i in mod_to_add:
      if (i == "OH-dehydration" and "OH-dehydration_amination" not in mod_to_add and "OH-dehydration_amination_methylation" not in mod_to_add) or (i == "OH-dehydration_amination" and "OH-dehydration_amination_methylation" not in mod_to_add) or (i == "D-glucosamine_N-acetylglucosamine_glucose_dehydration" and "D-glucosamine_N-acetylglucosamine_glucose_dehydration_amination" not in mod_to_add and "D-glucosamine_N-acetylglucosamine_glucose_dehydration_amination_N-methylation" not in mod_to_add) or (i == "D-glucosamine_N-acetylglucosamine_glucose_dehydration_amination" and "D-glucosamine_N-acetylglucosamine_glucose_dehydration_amination_N-methylation" not in mod_to_add):
        search_string += i
        search_string += " "

    for i in primary_mod:
        search_string += i
        search_string += " "

    with open(os.path.join(self.out_files, str(self.mod_depth) + os.path.basename(self.BGCfile) + self.model + "_modification_results.txt"), "w") as f:
      for i in backbone_smiles.keys():
        f.write(backbone_smiles[i]["smiles"] + "\n")

      f.close()

    for i in backbone_smiles.keys():
      #for cluster running
      mol = Chem.MolFromSmiles(backbone_smiles[i]["smiles"])
      Chem.Kekulize(mol)
      non_aromatic_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
      filename = os.path.join(self.out_files, str(self.mod_depth) + os.path.basename(self.BGCfile) + self.model + "_modification_results.txt")
      command = "singularity run {} core2pks  -c '{}' -d {} -m {} -o 20 >> {}".format(self.mod_sif, non_aromatic_smiles, str(self.mod_depth), search_string, filename)
      print(command)
      os.system(command)

    with open(os.path.join(self.out_files, str(self.mod_depth) + os.path.basename(self.BGCfile) + self.model + "_modification_results.txt")) as f:
      for line in f: 
        line.strip("\n")
        mod.append(line) 
      for i in backbone_smiles.keys():
        mod.append(backbone_smiles[i]["smiles"])
      f.close()
    
    return mod


  def generate_chemical_database(self):
    number = 1
    rows = []

    filename = os.path.join(self.out_files, "{}_{}_{}_compound_database.csv".format(self.bond, self.mod_depth, os.path.basename(self.BGCfile)))
    print("save file: ", filename)
    print("number of final products: ", len(self.products))
    with open(filename, "w") as f:
        writer = csv.writer(f)
        header = ["name", "smiles", "inchikey14", "exact_mass"]
        writer.writerow(header)
        for smi in self.products:
          rows.append("compound_{},{},{},{}".format(number, smi.strip("\n"), Chem.MolToInchiKey(Chem.MolFromSmiles(smi))[0:14], Descriptors.ExactMolWt(Chem.MolFromSmiles(smi))).split(","))
          number += 1
        writer.writerows(rows)
        f.close()
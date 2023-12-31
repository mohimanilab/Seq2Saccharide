import argparse
import os
from structure_generater import structure_generater
from config import args

structure_generater(args.mode, args.BGC, args.model, args.mon_cutoff, args.mod_cutoff, args.mon_hmm, args.mod_hmm, args.out_files, args.connectivity, args.BGClist, args.train_BGC_dir, args.BGCbackbone, args.bond_list, args.mod_sif, args.primary_mod_file,
                   args.length, args.mod_depth, args.gene_length, args.bond, args.primary_opt, args.set_num)

if args.if_spec:
    if args.mass_spec is None or args.prob is None:
        raise Exception("Please give both the mass spectrum data and probability json file for analysis.")
    else:
        args.mod_sif = os.path.abspath(args.mod_sif)
        command = "singularity run -B {},{} {} dereplicator_plus score --fragments {} --database {} --spectra-files {} -P {} > {}_outputs.txt".format(args.mod_sif, args.mass_spec, args.mod_sif, args.prob, os.path.join(args.out_files, "{}_{}_{}_compound_database.csv".format(args.bond, args.mod_depth, os.path.basename(args.BGC))), args.mass_spec, args.mass_error, os.path.basename(args.BGC))
        print(command)
        os.system(command)










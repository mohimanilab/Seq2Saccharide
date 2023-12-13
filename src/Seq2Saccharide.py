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
        if args.out_files != "":
            command = "dereplicator dereplicator_plus explain -f {} -d {}/{}_compound_database.csv -s {} -P 0.01 > {}_outputs.txt".format(args.prob, args.out_files, os.path.basename(args.BGC), args.mass_spec, os.path.basename(args.BGC))
        else:
            command = "dereplicator dereplicator_plus explain -f {} -d {}_compound_database.csv -s {} -P 0.01 > {}_outputs.txt".format(args.prob, os.path.basename(args.BGC), args.mass_spec, os.path.basename(args.BGC))

        # all_command = (
        #         'module use /projects/mohimanilab/mustafa/tools/modulefiles && ' 
        #         'module load dereplicator && ' 
        #         f'{command}'
        # )
        os.system(command)










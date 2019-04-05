
import os, sys
proj_dir = os.getcwd()
msprime_dir = proj_dir +"/output/msprime/"
archie_src_dir = proj_dir + "/src/sprime/"
archie_out_dir = proj_dir + "/output/sprime/"


sys.path.insert(0, proj_dir + '/src/stdpopsim')

#from stdpopsim import homo_sapiens, models
import homo_sapiens, models
import msprime
import itertools
import random
import numpy as np
import pandas as pd
import math

def add_mnms(ts, model_label, output_dir="./", rep_label=0, mnm_dist=100, mnm_frac=0.015):
    
    prefix = model_label + "_mnm" + str(mnm_dist) + "-" + str(mnm_frac) + "_" + str(rep_label)
    
    mnm_dict = {}

    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    l_IDs = []
    for popID in ["AFR", "EUR", "EA"]:
        for indID in range(50):
            l_IDs.append("%s_ind%s" % (popID, indID))

    header.extend(l_IDs)

    if rep_label == str(0):
        with open(output_dir + prefix + ".indID", "w") as text_file:
            for indID in l_IDs:
                text_file.write(indID + "\n")

    with open(output_dir + prefix + ".vcf", "w") as text_file:
        text_file.write("\t".join(header) + "\n") 

        for variant in ts.variants():
            pos = round(variant.site.position)
            # output current variant to VCF
            snv = ["1", str(pos), "1_%s" % pos, "A", "G", "50", "PASS", "VT=SNP", "GT"] 
            snv.extend([("|").join([str(g) for g in variant.genotypes[i: i+2]]) for i in range(0, len(variant.genotypes), 2)])
            text_file.write("\t".join(snv) + "\n") 

            random.seed(variant.site.position)
            if random.random() < mnm_frac:
                dist = random.randint(1,mnm_dist)
                mnm_cand = variant.site.position + dist 
                mnm_cand_r = str(round(mnm_cand))
                snv = ["1", str(mnm_cand_r), "1_%s" % mnm_cand_r, "A", "G", "50", "PASS", "VT=MNM", "GT"] 
                snv.extend([("|").join([str(g) for g in variant.genotypes[i: i+2]]) for i in range(0, len(variant.genotypes), 2)])
                text_file.write("\t".join(snv) + "\n") 



# ## Simulate models
# 
# Simulate 150 samples (50 for each of African, European, and EA ancestry) under each of the specified models (with 1 sample)



# coalescent simulation parameters
sample_size = 100 #haploid, each pop
length = 50000
mu = 1.15e-8
rr = 1e-8
seed = 30

# Gutenkunst 3-population model
GutenkunstThreePop_model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
GutenkunstThreePop_ts = msprime.simulate(
    # first 100 haploid samples from AFR, next 100 from EUR, followed by another 100 from EA
    samples=[msprime.Sample(0, 0)]*sample_size + [msprime.Sample(1, 0)]*sample_size + [msprime.Sample(2, 0)]*sample_size,
    length=length, 
    mutation_rate=mu, 
    recombination_rate=rr,
    random_seed=seed,
    num_replicates=replicates,
    **GutenkunstThreePop_model.asdict())



model_dict = {"GutenkunstThreePop": GutenkunstThreePop_ts}
#             "TennessenTwoPop": TennessenTwoPop_ts}


# ## Add MNMs
# Loop through models/replicates and randomly generate MNMs.

run_archie = False

if __name__ == "__main__":

    
    replicates = sys.argv[1]

    # MNM simulation parameters
    try:
        print("Simulating with MNMs on " + model_label)
        mnm_dist = sys.argv[2]
        mnm_frac = sys.argv[3]
    except IndexError:
        print("Simulating without MNMs on " + model_label)
        mnm_dist = 100
        mnm_frac = -100


    for model_label, model in model_dict.items():
        
        
        for j, ts in enumerate(model):

            prefix = model_label + "_mnm" + str(mnm_dist) + "-" + str(mnm_frac) + "_"

            # add MNMs and output results
            add_mnms(ts, 
                     output_dir=msprime_dir, 
                     model_label=model_label, 
                     mnm_dist=mnm_dist, 
                     mnm_frac=mnm_frac, 
                     rep_label=str(j))


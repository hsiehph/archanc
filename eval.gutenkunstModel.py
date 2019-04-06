
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
import math, time

def add_mnms(ts, model_label, output_dir="./", rep_label=0, mnm_dist=100, mnm_frac=0.015):

    if mnm_frac != -100:
        prefix = model_label + "_mnm" + str(mnm_dist) + "-" + str(mnm_frac) + "_" + str(rep_label)
    else:
        prefix = model_label + "_womnm_%s" % str(rep_label)

    mnm_dict = {}

    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

#    dict_IDs = {}
    l_IDs = []
    for popID in ["AFR", "EUR", "EA"]:
#        dict_IDs[popID] = []
        for indID in range(50):
#            dict_IDs[popID].append("%s_ind%s" % (popID, indID))
            l_IDs.append("%s_ind%s" % (popID, indID))

    header.extend(l_IDs)

#
#    if rep_label == str(1):
#        for popID in ["AFR", "EUR", "EA"]:
#            with open("%s%s.%s.indID" % (output_dir, model_label, popID), "w") as text_file:
#                for indID in dict_IDs[popID]:
#                    text_file.write(indID + "\n")
#    else:
#        time.sleep(5)

#    with open(output_dir + prefix + ".vcf", "w") as text_file:
    with sys.stdout as text_file:
        text_file.write("\t".join(header) + "\n") 

        for variant in ts.variants():
            pos = round(variant.site.position)
            # output current variant to VCF
            snv = ["1", str(pos), "1_%s" % pos, "A", "G", "50", "PASS", "VT=SNP", "GT"] 
            snv.extend([("|").join([str(g) for g in variant.genotypes[i: i+2]]) for i in range(0, len(variant.genotypes), 2)])
            text_file.write("\t".join(snv) + "\n") 

            random.seed(variant.site.position)
            if random.random() < mnm_frac:
                # make sure a mnm is at least 10bp from its counterpart.
                # this is to avoid the internal filter of Sprime.
                dist = random.randint(10, mnm_dist)
                mnm_cand = variant.site.position + dist 
                mnm_cand_r = str(round(mnm_cand))
                snv = ["1", str(mnm_cand_r), "1_%s" % mnm_cand_r, "A", "G", "50", "PASS", "VT=MNM", "GT"] 
                snv.extend([("|").join([str(g) for g in variant.genotypes[i: i+2]]) for i in range(0, len(variant.genotypes), 2)])
                text_file.write("\t".join(snv) + "\n") 



if __name__ == "__main__":

    replicateID = int(sys.argv[1])

    # ## Simulate models
    # Simulate 150 samples (50 for each of African, European, and EA ancestry) under each of the specified models (with 1 sample)
    # coalescent simulation parameters
    replicates = 1
    sample_size = 100 #haploid, each pop
    length = 50000
    mu = 1.15e-8
    rr = 1e-8
    seed = replicateID

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


    for model_label, model in model_dict.items():
        # MNM simulation parameters
        try:
            mnm_dist = float(sys.argv[2])
            mnm_frac = float(sys.argv[3])
            sys.stderr("Simulating with MNMs on " + model_label)
        except IndexError:
            mnm_dist = 100
            mnm_frac = -100
            sys.stderr("Simulating without MNMs on " + model_label)

        for j, ts in enumerate(model):

            # add MNMs and output results
            add_mnms(ts, 
                     output_dir=msprime_dir, 
                     model_label=model_label, 
                     mnm_dist=mnm_dist, 
                     mnm_frac=mnm_frac, 
                     rep_label=str(replicateID))


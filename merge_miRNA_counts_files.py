#!/usr/bin/env python
# coding: utf-8

# In[7]:

import pandas as pd
import numpy as np
import os
import subprocess
import sys
import random
from io import StringIO


# In[12]:

def mergeMiRnaCounts(sample_list, mature_path, genome_path, output_dir):

    print ()
    # Reading sample list file
    samples=np.genfromtxt(sample_list, encoding='utf-8', dtype=str)
    print ("samples :" + str(samples))
    nFiles=len(samples)
    print("Number of samples:" + str(nFiles))
    print ()

    #create empty array to merge all samples count
    all_miRNAs=np.empty((2656,nFiles+1), dtype=object)
    print ("miRNA counts will be merged in: \"all_miRNAs\"")
    print ("all_miRNAs shape :" + str(all_miRNAs.shape))
    print()

    #For each sample:
    i=0
    for sample in samples:

        if (i==0):
            print ("Reading miRNA counts for sample: " + sample)

            #read mature miRNA counts
            f1_mature_file=mature_path + "/" + sample + "-maturemiRNA-counts.txt"
            os.system("grep -v '*' " + f1_mature_file + " > tmpfile && mv tmpfile " + f1_mature_file)
            f1_mature=np.loadtxt(f1_mature_file, dtype=str, delimiter='\t', skiprows=1, encoding='utf-8')
            print ("f1_mature shape: " + str(f1_mature.shape))
            f1_mature=f1_mature[f1_mature[:, 0].argsort()]

            #read genome mi RNA counts
            f1_genome=np.loadtxt(genome_path + "/" + sample + "-taggedBAMcounts.txt", dtype=str, delimiter='\t', skiprows=1, encoding='utf-8')
            print ("f1_genome shape: " + str(f1_genome.shape))
            f1_genome=f1_genome[f1_genome[:, 0].argsort()]
            print ()

            #testing that f1_mature ans f1_genome have same order
            l = [random.randint(0,len(all_miRNAs)) for i in range(10)]
            print("checking miRNA order in counts files")
            for n in l:
                if f1_mature[n,0] != f1_genome[n,0]:
                    print("mature miRNA and genome miRNA don't have same order")
                    return
            print("Mature and genome counts files are ordered")
            print ()

            #write miRNA counts to all_miRNAs array
            print("Writing sample counts to all_miRNAs table")
            all_miRNAs[:,0]=f1_mature[:,0]
            all_miRNAs[:,1]=f1_mature[:,1].astype(int)+f1_genome[:,1].astype(int)
            
            print ()

        else:
            # For the rest of samples
            print ("counter =" + str(i))
            print ("sample " + str(i+1) + ":" + sample) 
            #reading sample files
            print ("Reading miRNA counts for sample: " + sample)
            print ()
            f2_mature_file=mature_path + "/" + sample + "-maturemiRNA-counts.txt"
            os.system("grep -v '*' " + f2_mature_file + " > tmpfile && mv tmpfile " + f2_mature_file)
            f2_mature=np.loadtxt(f2_mature_file, dtype=str, delimiter='\t', skiprows=1, encoding='utf-8')
            f2_mature=f2_mature[f2_mature[:, 0].argsort()]
            print ("f2_mature shape :" + str(f2_mature.shape))
            f2_genome=np.loadtxt(mature_path + "/" + sample + "-maturemiRNA-counts.txt", dtype=str, delimiter='\t', skiprows=1, encoding='utf-8')
            f2_genome=f2_genome[f2_genome[:, 0].argsort()]
            print ("f2_genome shape :" + str(f2_genome.shape))

            #Checking that f2_mature and f2_genome have the same order
            print("checking miRNA order counts files")
            for n in l:
                if f2_mature[n,0] != f2_genome[n,0]:
                    print ("mature miRNA and genome miRNA don't have same order")
                    return
            print ("f2_mature and f2_genome are ordered")
            print ()

            #checking that miRNA order in next sample is the same than in first sample"
            print("checking miRNA order in sample : " + sample)
            for n in l:
                if f1_mature[n,0] != f2_mature[n,0]:
                    print("miRNA are odered in sample: " + sample)
                    return
            print("miRNA in sample: " + sample + " are ordered")
            print ()

            # merging conts of next sample to all_miRNAs
            print("Merging counts to all_miRNAs")
            all_miRNAs[:,i+1]=f2_mature[:,1].astype(int)+f2_genome[:,1].astype(int)
            print(" First two rows in all_miRNAs")
            print (all_miRNAs[0:2])
            print ()
        i+=1

    # Save all_miRNAs to a csv file
    print ( "Saving all_miRNAs to a file")
    header= np.insert(samples, 0, "miRNA")
    print ("File header: " + str(header))

    #np.savetxt(output_dir +"/all_miRNAs.csv", all_miRNAs, delimiter=',', fmt=fmt, header=header, comments='')
    np.savetxt(output_dir +"/all_miRNAs.csv", all_miRNAs, fmt='%s', delimiter=',', header=np.array2string(header, separator=','))
    print("miRNAs counts from MIRbase and genome alignments have been merged in file all_miRNAs.csv")

# In[ ]:

def main():

    print ()

    sample_list=sys.argv[1]
    mature_path=sys.argv[2]
    genome_path=sys.argv[3]
    output_dir=sys.argv[4]
    print ("sample list file: " + sample_list)
    print ("mature miRNA counts path: " + mature_path)
    print ("genome miRNA counts path: " + genome_path)

    mergeMiRnaCounts(sample_list, mature_path, genome_path, output_dir)

# Execute `main()` function
if __name__ == '__main__':
    main()


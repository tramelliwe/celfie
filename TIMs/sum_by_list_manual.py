# imports
import csv
import numpy as np
import pandas as pd
import sys
import collections
import os


def get_region_dict(file, tissues):
    """
    gets TIM regions to combine CpGs within
    """

    regions_dict = collections.defaultdict(list)  # intialize dict

    with open(file, "r") as input_file:
        list_file = csv.reader(input_file, delimiter="\t")

        for line in list_file:
            # for each line in the tim file get position info
            chrom, start, end = line[0], line[1], line[2]

            # add to a dictionary for easy lookup
            regions_dict[chrom].append(
                {
                    "start": int(start),
                    "end": int(end),
                    "meth": np.zeros(tissues),
                    "depth": np.zeros(tissues),
                }
            )

    return regions_dict




def get_methylation_counts(file, regions_dict):
    """
    add together the methylation counts for all CpGs in the selected region
    """
    c=0
    # file of CpGs
    with open(file, "r") as input_file:
        cpg_file = csv.reader(input_file, delimiter="\t")
        
        # get methylation read counts for each position
        for line in cpg_file:
            c+=1
            print(c)
            line = [element.replace(',', '.') for element in line]
            chrom, start, end = line[0], line[1], line[2]
            meth = np.array(line[3::2], dtype=np.float64)
            depth = np.array(line[4::2], dtype=np.float64)

            if chrom in regions_dict:
                for region in regions_dict[chrom]:
                    # check if the CpG is in any region
                    # if so, add its methylation counts
                    if int(start) >= region["start"] and int(end) <= region["end"]:
                        region["meth"] += meth
                        region["depth"] += depth

    return regions_dict


def write_bed_file(output_file, regions_dict):
    """
    write bed file of summed counts for all tissues
    """
    with open(output_file, "w") as output:
        bed_file = csv.writer(output, delimiter="\t",lineterminator="\n")
        
        for chrom in regions_dict:
            for region in regions_dict[chrom]:
                values = []

                for m, d in zip(region["meth"], region["depth"]):

                    # if data is missing output 0 count
                    if np.isnan(m):
                        m = 0.0
                    if np.isnan(d):
                        d = 0.0

                    values.append(m)
                    values.append(d)

                bed_file.writerow(
                    [chrom] + [region["start"]] + [region["end"]] + values
                )



list_file = "temp/celfie_original_atlas/top_500_markers_merged.bed"
=======
list_file = "test_celfie/markers_500_windowed.txt"
tissue_cpg_file = "temp/summed_tissues.txt"
output_file_name = "temp/summed_tissues_markers.bed"

tissues = 39

folder = ("C:/Users/cwillemart/Downloads/bed/bed/")
output_folder = os.path.join(folder,"summed_over_marker/")
files = os.listdir(folder)
file_paths = [os.path.join(folder, file) for file in files if file.endswith("bed")]

for file_path in file_paths:
    tissue_cpg_file = file_path
    sample = os.path.basename(file_path)
    output_file_name = output_folder + sample
    regions = get_region_dict(list_file, tissues)
    get_methylation_counts(tissue_cpg_file, regions)  # get methylation read counts of Cpgs within region
    write_bed_file(output_file_name, regions)

=======
list_file = "top_500_markers_merged.bed"
directory_files = "20240417_output/bed"
tissues = 1

for file in os.listdir(directory_files):
    tissue_cpg_file = os.path.join(directory_files,file)
    output_file_name = os.path.join(directory_files,os.path.splitext(file)[0]+"_summed.bed")
    regions = get_region_dict(list_file, tissues)
    get_methylation_counts(tissue_cpg_file, regions)  # get methylation read counts of Cpgs within region
    write_bed_file(output_file_name, regions)
   
>>>>>>> Stashed changes

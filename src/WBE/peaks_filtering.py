import csv
import numpy as np
import cvxpy as cp
import glpk
import copy
import time
import math
import sys

def read_barcode_csv_file(file_path):
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        mutations = next(csv_reader)  # Store headers separately
        mutations = mutations[1:]
        haplotypes = []
        hap_mut_matrix = []
        for row in csv_reader:
            hap_mut_matrix.append(row[1:])
            haplotypes.append(row[0])   
        #Convert to NumPy array of integers
        hap_mut_matrix = np.array([[int(elem) for elem in row] for row in hap_mut_matrix])
    
    return mutations, haplotypes, hap_mut_matrix

def read_vcf_file(file_path):
    with open(vcf_file, 'r') as file:
        vcf_reader = csv.reader(file, delimiter='\t')
        af_values = []
        depth_values = []
        for row in vcf_reader:
            if not row[0].startswith('#'):  # Skip header rows
                af_field = row[-2]
                af_value = af_field.split(';')[0].split('=')[1]
                af_values.append(float(af_value))
                depth_values.append(int(row[-1])+1)

    depth = np.log2(depth_values)
    depth = depth / np.max(depth)
    return np.array(af_values), depth

def read_condensed_csv_file(file_path):
    # Create an empty dictionary to store the data
    condensed_nodes_map = {}
    with open(file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)   
        # Iterate through each row in the CSV file
        for row in csv_reader:
            # Extract the key (CONDENSED-*) and values from the row
            key, *values = row
            # Store the values in a list in the dictionary
            condensed_nodes_map[key] = values

    return condensed_nodes_map

def read_loss_mutations_file(file_path):
    # Create an empty dictionary to store the data
    loss_mutations = []
    with open(file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)   
        # Iterate through each row in the CSV file
        for row in csv_reader:
            loss_mutations.append(row[0])

    return loss_mutations

def write_vcf_file(file_path, mut_hap, haplotypes, mutations):
    with open(file_path, 'w') as file:
        file.write("##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(haplotypes)))
        
        i = 0
        while i < len(mutations):
            same_pos_list = []
            hap_presence = copy.deepcopy(mut_hap[i]) 
            mut_string = mutations[i][-1]
            vcf_string = "NC_045512v2\t" + mutations[i][1:-1] + "\t" + mutations[i]
            for j in range(len(mutations[i+1:])):
                if mutations[i+1+j][1:-1] == mutations[i][1:-1]:
                    same_pos_list.append(i+1+j)
            for idx in same_pos_list:
                vcf_string += "," + mutations[idx]
                mut_string += "," + mutations[idx][-1]
                for k in range(len(mut_hap[i])):
                    if mut_hap[idx][k] == 1:
                        hap_presence[k] = idx - i + 1
            vcf_string += "\t" + mutations[i][0] + "\t" + mut_string + "\t.\t.\t.\t.\t"
            file.write(vcf_string)
            file.write("\t".join(str(present) for present in hap_presence))
            file.write("\n")
            i += 1 + len(same_pos_list)

def write_csv_file(data, file_path):
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)

def cp_solve(A, b, d):
    x = cp.Variable(A.shape[1])
    # Perform element-wise multiplication with broadcasting
    weighted_A =  A * d[:, np.newaxis]
    weighted_b = b * d
    cost = cp.norm(weighted_A @ x - weighted_b, 1)
    constraints = [sum(x) == 1, x >= 0]
    prob = cp.Problem(cp.Minimize(cost), constraints)
    prob.solve(verbose=False)
    return x.value, cost.value

def solve_abundance(hap_mut_matrix, read_af, depth_values, haplotypes, mutations):
    #Setting the variables for the regression problem
    A = hap_mut_matrix.T
    #Running the loop twice to ensure abundances sum up to 1
    for i in range(2):
        sol, ref_cost = cp_solve(A, read_af, depth_values)
        #Make the proportion to be 0 for anything < eps
        sol[sol < eps] = 0
        #Removing haplotypes that are < 0 in abundance
        hap_idx_remove = np.where(sol == 0)[0]
        sol = [val for i, val in enumerate(sol) if i not in hap_idx_remove]
        A = np.delete(A, hap_idx_remove, axis = 1)
        haplotypes = [val for i, val in enumerate(haplotypes) if i not in hap_idx_remove]
    
    #Remove rows (mutations) with all zeros (due to removal of haplotypes)
    zero_rows = np.all(A == 0, axis=1)
    mut_idx_remove = np.where(zero_rows)[0]
    A = np.delete(A, mut_idx_remove, axis=0)
    mutations = [val for i, val in enumerate(mutations) if i not in mut_idx_remove]
    return A, haplotypes, sol, mutations

#Check arguments
if len(sys.argv) != 3:
    printf("USAGE: python lineage_abundance.py <file_prefix> <directory>")
    sys.exit(1)

# Start time
start_time = time.time()
eps = 1e-3
file_prefix = sys.argv[1]
directory = sys.argv[2]

# Reading File
barcode_file_path = directory + "/" + file_prefix + "_barcode.csv"
mutations, haplotypes, hap_mut_matrix = read_barcode_csv_file(barcode_file_path)
print(f'\nPeaks given for filtering: {len(haplotypes)}')

vcf_file = directory + "/" + file_prefix + "_read_data.vcf"
af_values, depth_values = read_vcf_file(vcf_file)

condensed_file_path = directory + "/" + file_prefix + "_condensed_nodes.csv"
condensed_nodes_map = read_condensed_csv_file(condensed_file_path)

loss_mutations_file_path = directory + "/" + file_prefix + "_loss_mutations.csv"

#Solving abundance
mut_hap_matrix, haplotypes, abundances, mutations = solve_abundance(hap_mut_matrix, af_values, depth_values, haplotypes, mutations)

#Abundance of lineages
lineages = {}
uncertain_lineages = {}
for i, hap in enumerate(haplotypes):
    #Check if it is a condensed node
    if "CONDENSED" in hap:
        condensed_lineages = set()
        #Get all the lineages from condensed nodes
        for nodes in condensed_nodes_map[hap]:
            split_parts = nodes.split('_')
            if split_parts[-1] not in condensed_lineages:
                condensed_lineages.add(split_parts[-1])
        #Add the abundance if single lineage is present
        if len(condensed_lineages) == 1:
            lin = condensed_lineages.pop()
            if lin not in lineages:
                lineages[lin] = abundances[i]
            else:
                lineages[lin] += abundances[i]
        #Multiple lineages
        else:
            abun = abundances[i]
            for lin in condensed_lineages:
                if abun not in uncertain_lineages:
                    uncertain_lineages[abun] = [lin]
                else:
                    uncertain_lineages[abun].append(lin)
    else:
        #Extract the part after the last '_'
        split_parts = hap.split('_')
        if split_parts[-1] not in lineages:
            lineages[split_parts[-1]] = abundances[i]
        else:
            lineages[split_parts[-1]] += abundances[i]

print("\nLINEAGE ABUNDANCE (ORIG):")
for lin, abun in lineages.items():
    print(lin, abun)

for abun, lin in uncertain_lineages.items():
    lin_str = ', '.join(lin)
    print(lin_str, abun) 

#Write abundance of haplotypes in csv
csv_write_header = ['Haplotype', 'Abundance']
csv_write_data = np.vstack((haplotypes, abundances)).T
csv_write_data = np.vstack((csv_write_header, csv_write_data))
write_csv_file(csv_write_data, directory + "/" + file_prefix + "_haplotype_abundance.csv")

#Write VCF
write_vcf_file(directory + "/" + file_prefix + "_haplotypes.vcf", mut_hap_matrix, haplotypes, mutations)

# End time
print(f"\nElapsed time: {time.time() - start_time} seconds")

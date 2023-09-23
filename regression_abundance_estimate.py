import csv
import numpy as np
import cvxpy as cp
import glpk
import copy
import time

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
    
    return np.array(af_values), np.log2(depth_values)


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
    sol, ref_cost = cp_solve(A, read_af, depth_values)
    #Make the proportion to be 0 for anything < eps
    sol[sol < eps] = 0

    #Removing haplotypes that are < 0 in abundance
    hap_idx_remove = np.where(sol == 0)[0]
    sol = [val for i, val in enumerate(sol) if i not in hap_idx_remove]
    A = np.delete(A, hap_idx_remove, axis = 1)
    haplotypes = [val for i, val in enumerate(haplotypes) if i not in hap_idx_remove]
    
    #To ensure abundance sums up to 1
    sol, _ = cp_solve(A, read_af, depth_values)
    sol[sol < eps] = 0
    #Removing haplotypes that are < 0 in abundance
    hap_idx_remove = np.where(sol == 0)[0]
    sol = [val for i, val in enumerate(sol) if i not in hap_idx_remove]
    A = np.delete(A, hap_idx_remove, axis = 1)
    haplotypes = [val for i, val in enumerate(haplotypes) if i not in hap_idx_remove]
    
    #Remove rows with all zeros (due to removal of haplotypes)
    zero_rows = np.all(A == 0, axis=1)
    mut_idx_remove = np.where(zero_rows)[0]
    A = np.delete(A, mut_idx_remove, axis=0)
    mutations = [val for i, val in enumerate(mutations) if i not in mut_idx_remove]

    #Write VCF
    haplotypes = [hap + "_READ_1_29903" for hap in haplotypes]
    vcf_file = "my_vcf_haplotypes.vcf"
    write_vcf_file(vcf_file, A, haplotypes, mutations)
    
    return haplotypes, sol


# Start time
start_time = time.time()
eps = 1e-2

# Reading File
barcode_file_path = 'my_vcf_barcode.csv'
mutations, haplotypes, hap_mut_matrix = read_barcode_csv_file(barcode_file_path)

vcf_file = 'my_vcf_abundance.vcf'
af_values, depth_values = read_vcf_file(vcf_file)

#Solving abundance
haplotypes, abundances = solve_abundance(hap_mut_matrix, af_values, depth_values, haplotypes, mutations)

#Abundance of lineages
lineages = {}
for i, hap in enumerate(haplotypes):
    split_parts = hap.split('_')
    #Extract the part after the last '_'
    if split_parts[-4] not in lineages:
        lineages[split_parts[-4]] = abundances[i]
    else:
        lineages[split_parts[-4]] += abundances[i]

print("\nLINEAGE ABUNDANCE (ORIG):")
for lin, abun in lineages.items():
    if abun:
        print(lin, abun)

#Write abundance of haplotypes in csv
csv_write_header = ['Haplotype', 'Abundance']
csv_write_data = np.vstack((haplotypes, abundances)).T
csv_write_data = np.vstack((csv_write_header, csv_write_data))
write_csv_file(csv_write_data, "my_vcf_hap_abundance.csv")

# End time
print(f"\nElapsed time: {time.time() - start_time} seconds")
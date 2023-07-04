import csv
import numpy as np
import cvxpy as cp
import copy
import time

def read_csv_file(file_path):
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        mutations = next(csv_reader)  # Store headers separately
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

        for row in vcf_reader:
            if not row[0].startswith('#'):  # Skip header rows
                info_field = row[-1]
                af_value = info_field.split(';')[0].split('=')[1]
                af_values.append(float(af_value))
    
    return np.array(af_values)

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

def cp_solve(A, b):
    x = cp.Variable(A.shape[1])
    cost = cp.norm(A @ x - b, 1)
    constraints = [sum(x) == 1, x >= 0]
    prob = cp.Problem(cp.Minimize(cost), constraints)
    prob.solve(verbose=False)
    return x.value

def add_mutations_to_haplotypes(A, x_val, b, haplotypes): 
    A_copy = np.copy(A)
    ref_cost = np.linalg.norm(A @ x_val - b, 1)
    sol = x_val
    hap_new = copy.deepcopy(haplotypes)
    hap_added = 0
    
    #Removing haplotypes that are < 0 in abundance
    sol[sol < eps] = 0
    hap_idx_remove = np.where(sol == 0)[0]
    sol = [val for i, val in enumerate(sol) if i not in hap_idx_remove]
    A_copy = np.delete(A_copy, hap_idx_remove, axis = 1)
    hap_new = [val for i, val in enumerate(hap_new) if i not in hap_idx_remove]

    while True:
        hap_idx = -1
        residual = abs(A_copy @ sol - b)
        site = np.argmax(residual)
        for i, val in enumerate(A_copy[site]):
            #Copy the column at the end of the array and add the mutation there
            A_copy = np.concatenate((A_copy, A_copy[:, i][:, np.newaxis]), axis=1)
            # Invert the mutation in current haplotype
            A_copy[site, -1] = 1 - val
            #Find the new value of x
            curr_sol = cp_solve(A_copy, b)
            curr_cost = np.linalg.norm(A_copy @ curr_sol - b, 1)
            #Removing the last column that was added
            A_copy = A_copy[:, :-1]
            #Update the ref_cost on passing convergence condition
            if ((ref_cost - curr_cost) > (thresh * ref_cost)):
                ref_cost = copy.deepcopy(curr_cost)
                hap_idx = i
                sol = curr_sol
        
        #Check if cost could be reduced 
        if hap_idx >= 0:
            #Add the hap_idx column at the end of A_copy
            A_copy = np.concatenate((A_copy, A_copy[:, hap_idx][:, np.newaxis]), axis=1)
            # Invert the mutation in current haplotype
            A_copy[site, -1] = 1 - A_copy[site, hap_idx]
            hap_added += 1
            hap_new = np.append(hap_new, f"NEW_{hap_added}")
            
            #Removing haplotypes that are < 0 in abundance
            sol[sol < eps] = 0
            hap_idx_remove = np.where(sol == 0)[0]
            sol = [val for i, val in enumerate(sol) if i not in hap_idx_remove]
            A_copy = np.delete(A_copy, hap_idx_remove, axis = 1)
            hap_new = [val for i, val in enumerate(hap_new) if i not in hap_idx_remove]
            #print(f'Site: {site}, Hap: {hap_idx}, cost: {ref_cost}')
        else:
            if not hap_added:
                print("No Haplotype added!!!\n")
            break

    return A_copy, np.abs(sol), hap_new

def solve_abundance(hap_mut_matrix, read_af, haplotypes, mutations):
    #Setting the variables for the regression problem
    A = hap_mut_matrix.T
    orig_sol = cp_solve(A, read_af)
    #Make the proportion to be 0 for anything < eps
    orig_sol[orig_sol < eps] = 0
    
    #Add new haplotype names
    A_new, sol, hap_new = add_mutations_to_haplotypes(A, np.abs(orig_sol), read_af, haplotypes)
    
    #Remove rows with all zeros (due to removal of haplotypes)
    zero_rows = np.all(A_new == 0, axis=1)
    mut_idx_remove = np.where(zero_rows)[0]
    A_new = np.delete(A_new, mut_idx_remove, axis=0)
    mutations_new = [val for i, val in enumerate(mutations) if i not in mut_idx_remove]
    hap_new = [hap + "_READ_1_29903" for hap in hap_new]

    #Write VCF
    vcf_file = "my_vcf_haplotypes.vcf"
    write_vcf_file(vcf_file, A_new, hap_new, mutations_new)
    
    return hap_new, sol, orig_sol


# Start time
start_time = time.time()
eps = 1e-2
thresh = 0.04

# Reading File
barcode_file_path = 'my_vcf_barcode.csv'
mutations, haplotypes, hap_mut_matrix = read_csv_file(barcode_file_path)
mutations = mutations[1:]

vcf_file = 'my_vcf_abundance.vcf'
af_values = read_vcf_file(vcf_file)

#Solving abundance
new_haplotypes, new_abundances, abundances = solve_abundance(hap_mut_matrix, af_values, haplotypes, mutations)

#Abundance of lineages
lineages = {}
for i, hap in enumerate(haplotypes):
    split_parts = hap.split('_')
    #Extract the part after the last '_'
    if split_parts[-1] not in lineages:
        lineages[split_parts[-1]] = abundances[i]
    else:
        lineages[split_parts[-1]] += abundances[i]

print("\nLINEAGE ABUNDANCE (REGRESSION):")
for lin, abun in lineages.items():
    if abun:
        print(lin, abun)

#Abundance of haplotypes
print("\nHAPLOTYPE ABUNDANCE(Initial + New):")
for i in range(len(new_haplotypes)):
    print(new_haplotypes[i], new_abundances[i])

# End time
print(f"\nElapsed time: {time.time() - start_time} seconds")
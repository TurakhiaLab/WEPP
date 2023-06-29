import csv
import numpy as np
import cvxpy as cp
import copy

def read_csv_file(file_path):
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        mutations = next(csv_reader)  # Store headers separately
        haplotypes = []
        peak_mut_matrix = []
        
        for row in csv_reader:
            haplotypes.append(row[0])
            peak_mut_matrix.append(row[1:])
    
    return mutations, haplotypes, peak_mut_matrix


def read_vcf_file(file_path):
    with open(vcf_file, 'r') as file:
        vcf_reader = csv.reader(file, delimiter='\t')
        af_values = []
        for row in vcf_reader:
            if not row[0].startswith('#'):  # Skip header rows
                info_field = row[-1]
                af_value = info_field.split(';')[0].split('=')[1]
                af_values.append(float(af_value))
    return af_values


def write_vcf_file(file_path, mut_peak, haplotypes, mutations):
    with open(file_path, 'w') as file:
        file.write("##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(haplotypes)))
        i = 0
        while i < len(mutations):
            same_pos_list = []
            peak_presence = copy.deepcopy(mut_peak[i]) 
            mut_string = mutations[i][-1]
            vcf_string = "NC_045512v2\t" + mutations[i][1:-1] + "\t" + mutations[i]
            for j in range(len(mutations[i+1:])):
                if mutations[i+1+j][1:-1] == mutations[i][1:-1]:
                    same_pos_list.append(i+1+j)
            for idx in same_pos_list:
                vcf_string += "," + mutations[idx]
                mut_string += "," + mutations[idx][-1]
                for k in range(len(mut_peak[i])):
                    if mut_peak[idx][k] == 1:
                        peak_presence[k] = idx - i + 1
            vcf_string += "\t" + mutations[i][0] + "\t" + mut_string + "\t.\t.\t.\t.\t"
            file.write(vcf_string)
            file.write("\t".join(str(int(present)) for present in peak_presence))
            file.write("\n")
            i += 1 + len(same_pos_list)


def add_mutations_to_haplotypes(A, x_val, b): 
    hap_indices = []
    A_copy = np.copy(A)
    ref_cost = np.linalg.norm(A.astype('float64') @ x_val - b, 1)
    sol = x_val
    while True:
        hap_idx = -1
        residual = abs(A_copy.astype('float64') @ sol - b)
        site = np.argmax(residual)
        for i in range(len(A_copy[site])):
            #Only add mutation to haplotype where it is absent
            if not int(A_copy[site, i]):
                #Copy the column at the end of the array and add the mutation there
                A_copy = np.concatenate((A_copy, A_copy[:, i][:, np.newaxis]), axis=1)
                A_copy[site, -1] = 1
                #Find the new value of x
                x = cp.Variable(A_copy.shape[1])
                cost = cp.norm(A_copy @ x - b, 1)
                constraints = [sum(x) == 1, x >= 0]
                prob = cp.Problem(cp.Minimize(cost), constraints)
                prob.solve(verbose=False)
                curr_cost = np.linalg.norm(A_copy.astype('float64') @ x.value - b, 1)
                #Removing the last column that was added
                A_copy = A_copy[:, :-1]
                #Update the ref_cost
                if (abs(curr_cost - ref_cost) > 1e-9) and (curr_cost < ref_cost):
                    ref_cost = copy.deepcopy(curr_cost)
                    hap_idx = i
                    sol = x.value
        #Check if cost could be reduced 
        if hap_idx >= 0:
            #Add the hap_idx column at the end of A_copy
            A_copy = np.concatenate((A_copy, A_copy[:, hap_idx][:, np.newaxis]), axis=1)
            A_copy[site, -1] = 1
            hap_indices.append(hap_idx)
        else:
            if A_copy.shape == A.shape:
                print("No Haplotype added!!!")
            break
    return A_copy, hap_indices, np.abs(sol)


def solve_abundance(peak_mut_matrix, read_values, haplotypes, mutations):
    #Setting the variables for the regression problem
    eps = 1e-3
    A = np.array(np.array(peak_mut_matrix).T)
    b = np.array(read_values)
    x = cp.Variable(A.shape[1])
    cost = cp.norm(A @ x - b, 1)
    constraints = [sum(x) == 1, x >= 0]
    prob = cp.Problem(cp.Minimize(cost), constraints)
    prob.solve(verbose=False)
    sol = x.value
    
    #Add new haplotype names
    A_new, hap_idx_add, sol = add_mutations_to_haplotypes(A, np.abs(sol), b)
    hap_new = copy.deepcopy(haplotypes)
    count_dict = {}
    for idx in hap_idx_add:
        if idx in count_dict:
            count_dict[idx] += 1
        else:
            count_dict[idx] = 1
        hap_new += f'NEW_{count_dict[idx]}_{haplotypes[idx]}'
    hap_new = [hap + "_READ_1_29902" for hap in hap_new]

    #Make the proportion to be 0 for anything < eps
    sol[sol < eps] = 0
    #Removing haplotypes that are < 0 in abundance
    hap_idx_remove = np.where(sol == 0)[0]
    A_new = np.delete(A_new, hap_idx_remove, axis = 1)
    hap_new = [val for i, val in enumerate(hap_new) if i not in hap_idx_remove]
    sol = [val for i, val in enumerate(sol) if i not in hap_idx_remove]
    
    #Remove rows with all zeros (due to removal of haplotypes)
    A_new = A_new.astype(int)
    zero_rows = np.all(A_new == 0, axis=1)
    mut_idx_remove = np.where(zero_rows)[0]
    A_new = np.delete(A_new, mut_idx_remove, axis=0)
    mutations_new = [val for i, val in enumerate(mutations) if i not in mut_idx_remove]

    #Write VCF
    vcf_file = "my_vcf_haplotypes.vcf"
    write_vcf_file(vcf_file, A_new, hap_new, mutations_new)
    
    return hap_new, sol


# Reading File
barcode_file_path = 'my_vcf_barcode.csv'
mutations, haplotypes, peak_mut_matrix = read_csv_file(barcode_file_path)
mutations = mutations[1:]

vcf_file = 'my_vcf_abundance.vcf'
af_values = read_vcf_file(vcf_file)

#Solving abundance
haplotypes, abundances = solve_abundance(peak_mut_matrix, af_values, haplotypes, mutations)

for i in range(len(haplotypes)):
    print(haplotypes[i], abundances[i])
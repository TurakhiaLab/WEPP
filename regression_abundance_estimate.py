import csv
import numpy as np
import cvxpy as cp
import copy

def read_csv_file(file_path):
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        mutations = next(csv_reader)  # Store headers separately
        peak_nodes = []
        peak_mut_matrix = []
        
        for row in csv_reader:
            peak_nodes.append(row[0])
            peak_mut_matrix.append(row[1:])
    
    return mutations, peak_nodes, peak_mut_matrix


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


def write_vcf_file(file_path, mut_peak, peak_nodes, mutations):
    with open(file_path, 'w') as file:
        file.write("##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(peak_nodes)))
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


def add_mutations_to_haplotypes(A, sol, b): 
    ref_cost = np.linalg.norm(A.astype('float64') @ sol - b, 1)
    A_copy = np.copy(A)
    residual = abs(A_copy.astype('float64') @ sol - b)
    sorted_sites = np.argsort(residual)[::-1] 
    peak_indices = []
    for site in sorted_sites:
        peak_idx = -1
        for i in range(len(A_copy[site])):
            #Only add mutation to haplotype where it is absent
            if not int(A_copy[site, i]):
                A_copy[site, i] = 1
                curr_cost = np.linalg.norm(A_copy.astype('float64') @ sol - b, 1)
                if (abs(curr_cost - ref_cost) > 1e-9) and (curr_cost < ref_cost):
                    ref_cost = copy.deepcopy(curr_cost)
                    peak_idx = i
                #Bring the original mutation back
                A_copy[site, i] = 0
        if peak_idx >= 0:
            #Change the mutation in haplotype bringing the biggest change
            A_copy[site, peak_idx] = 1
            #Add unique peak indexes in the list
            if peak_idx not in peak_indices:
                peak_indices.append(peak_idx)
    
    return A_copy, peak_indices


def solve_abundance(peak_mut_matrix, read_values, peak_nodes, mutations):
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
    #Make the proportion to be 0 for anything < eps
    sol[sol < eps] = 0
    #Removing peaks that are < 0 in abundance
    peak_idx_remove = np.where(sol == 0)[0]
    A_new = np.delete(A, peak_idx_remove, axis = 1)
    peak_nodes_new = [val for i, val in enumerate(peak_nodes) if i not in peak_idx_remove]

    #Adding peaks that were modified to reduce cost
    A_hap_modified, peak_idx_add = add_mutations_to_haplotypes(A, sol, b)
    A_new = np.concatenate((A_new, A_hap_modified[:, peak_idx_add]), axis = 1)
    peak_nodes_new += ["NEW_" + peak_nodes[idx] for idx in peak_idx_add]
    peak_nodes_new = [peak + "_READ_0_29902" for peak in peak_nodes_new]

    #Remove rows with all zeros
    A_new = A_new.astype(int)
    zero_rows = np.all(A_new == 0, axis=1)
    mut_idx_remove = np.where(zero_rows)[0]
    A_new = np.delete(A_new, mut_idx_remove, axis=0)
    mutations_new = [val for i, val in enumerate(mutations) if i not in mut_idx_remove]

    #Write VCF
    vcf_file = "my_vcf_haplotypes.vcf"
    write_vcf_file(vcf_file, A_new, peak_nodes_new, mutations_new)

    #if np.array_equal(A, A_hap_modified):
    #    print("No New Haplotype found!!!")
    return sol


# Reading File
barcode_file_path = 'my_vcf_barcode.csv'
mutations, peak_nodes, peak_mut_matrix = read_csv_file(barcode_file_path)
mutations = mutations[1:]

vcf_file = 'my_vcf_abundance.vcf'
af_values = read_vcf_file(vcf_file)

#Solving abundance
abundances = solve_abundance(peak_mut_matrix, af_values, peak_nodes, mutations)

#for i in range(len(peak_nodes)):
#    print(peak_nodes[i], abundances[i])

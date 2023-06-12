import csv
import numpy as np
import cvxpy as cp

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


def solve_abundance(peak_mut_matrix, peak_nodes, read_values, eps):
    # set up and solve demixing problem
    A = np.array(np.array(peak_mut_matrix).T)
    b = np.array(read_values)
    x = cp.Variable(A.shape[1])
    cost = cp.norm(A @ x - b, 1)
    constraints = [sum(x) == 1, x >= 0]
    prob = cp.Problem(cp.Minimize(cost), constraints)
    prob.solve(verbose=False)
    sol = x.value
    # extract lineages with non-negligible abundance
    sol[sol < eps] = 0
    return sol

# Reading File
barcode_file_path = 'my_vcf_barcode.csv'
mutations, peak_nodes, peak_mut_matrix = read_csv_file(barcode_file_path)
mutations = mutations[1:]

vcf_file = 'Freyja_test.vcf'
af_values = read_vcf_file(vcf_file)


#Solving abundance
abundances = solve_abundance(peak_mut_matrix, peak_nodes, af_values, 1e-3)

for i in range(len(peak_nodes)):
    print(peak_nodes[i], abundances[i])

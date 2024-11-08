import csv

# Open the file for reading
with open('my_output_latest.txt', 'r') as file:
    # Read the lines from the file
    lines = file.readlines()

# Extract lineages and abundances
for l in lines:
    if l.split()[0] == "lineages":
        lineages = l.split()[1:]
    elif l.split()[0] == "abundances":
        abundances = l.split()[1:]

# Combine lineages and abundances
result = ['Lineage,Abundance']
result += [f"{lineage},{abundance}" for lineage, abundance in zip(lineages, abundances)]

#Write CSV
csv_write_header = ['Lineage', 'Abundance']

with open('my_vcf_freyja_results.csv', 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for item in result:
        csv_writer.writerow(item.split(','))

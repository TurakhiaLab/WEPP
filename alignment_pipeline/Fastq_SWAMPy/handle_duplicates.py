import sys

if len(sys.argv) != 3:
    print("Usage: python handle_duplicates.py <input.vcf> <output.vcf>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

read_count = {}

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
            continue
        
        fields = line.split('\t')
        read_name = fields[2]  
        
        if read_name in read_count:
            read_count[read_name] += 1
        else:
            read_count[read_name] = 1
        
        modified_read_name = f"{read_name}_{read_count[read_name]}"
        fields[2] = modified_read_name
        modified_line = '\t'.join(fields)
        outfile.write(modified_line)

print(f"Processed {input_file} and wrote modified VCF to {output_file}")

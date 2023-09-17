import sys

if len(sys.argv) != 2:
    print("Usage: python3 process_rc.py <fastq>")
    exit()

fastq = sys.argv[1]
out_fastq = fastq[:-6] + "_processed.fastq"

f = open('./intermediate_files/alignment.sam', 'r')
sam_lines = f.readlines()
sam_lines = sam_lines[3:]
f.close()

f = open(fastq, 'r')
fastq_lines = f.readlines()
f.close()

for i in range(len(sam_lines)):
    lines = sam_lines[i].strip().split('\t')
    read_name = lines[0]
    flag = int(lines[1])
    
    # Check for reverse complement
    if flag & 16:
        # Reverse complement read
        read = fastq_lines[4*i + 1].strip()
        read = read[::-1]
        read = read.translate(str.maketrans('ATCG', 'TAGC'))
        fastq_lines[4*i + 1] = read + '\n'
        
        # Reverse complement quality
        qual = fastq_lines[4*i + 3].strip()
        qual = qual[::-1]
        fastq_lines[4*i + 3] = qual + '\n'

f = open(out_fastq, 'w')
f.write(''.join(fastq_lines))
f.close()

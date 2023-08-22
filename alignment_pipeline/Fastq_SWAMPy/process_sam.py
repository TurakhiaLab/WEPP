import sys

def reverse(s):
    """Reverse a string."""
    return s[::-1]

def revcom(seq):
    """Get the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement[base] for base in reverse(seq))

def main():
    for line in sys.stdin:
        # If the line starts with '@', it's a header and should be printed as-is
        if line.startswith('@'):
            print(line, end='')
            continue

        # Split the SAM line into its fields
        fields = line.strip().split("\t")

        # Check if the 16th bit in the FLAG (fields[1]) is set
        if int(fields[1]) & 16:
            # Modify the SEQ and QUAL fields for reverse complement alignment
            fields[9] = revcom(fields[9])
            fields[10] = reverse(fields[10])
            fields[1] = str(int(fields[1]) - 16)  # Update the FLAG to indicate forward strand alignment

        # Print the modified SAM line
        print("\t".join(fields))

if __name__ == "__main__":
    main()

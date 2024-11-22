from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqIO import write


#### shows the identifier ID only 
fasta_file = '/Users/labuser3/Downloads/All_PR_ALigned.fasta'
for record in SeqIO.parse(fasta_file, 'fasta'): 
    print(f"ID: {record.id}")
  
## counting every line which was 3321   
with open(fasta_file, 'r') as file:
    line_count =sum(1 for line in file )
print(f"{line_count}")


## counted the numbers of sequences which was true including the reference : 1661
identifier = sum( 1 for _ in SeqIO.parse(fasta_file, 'fasta'))
print({identifier})


def get_amino_acid_at_position(fasta_file, position):
    # Read the first sequence from the FASTA file
    first_seq = next(SeqIO.parse(fasta_file, "fasta")).seq

    # Ensure the position is within the length of the sequence
    if position < len(first_seq):
        amino_acid = first_seq[position]
        print(f"Amino acid at position {position} is: {amino_acid}")
    else:
        print(f"Position {position} is out of range for the first sequence.")

#position = 605  # 0-index based, so position 91 is index 90

def filter_sequences(input_file, output_file, position_amino_acid_dict):

    filtered_sequences = []

# Iterate through sequences in the input file
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        print(f"Processing sequence {record.id} of length {len(sequence)}")  # Debug print
        
        keep = True
        
        # Check if the sequence meets the criteria
        for pos, amino_acid in position_amino_acid_dict.items():
            if pos >= len(sequence):
                print(f"Position {pos + 1} is out of range for sequence {record.id} with length {len(sequence)}")
                keep = False
                break
            if sequence[pos] != amino_acid:
                print(f"Sequence {record.id} does not meet the criteria at position {pos + 1}")
                print(f"Expected {amino_acid} but found {sequence[pos]} at position {pos + 1}")
                keep = False
                break
        
        if keep:
            print(f"Sequence {record.id} meets all criteria")
            filtered_sequences.append(record)
    
    # Write the filtered sequences to the output file
    SeqIO.write(filtered_sequences, output_file, "fasta")
    print(f"Filtered sequences have been saved to {output_file}")
    print(f"Total filtered sequences: {len(filtered_sequences)}")

# Example usage
input_file = "All_PR_ALigned.fasta"
output_file = "filtered_sequences.fasta"

# Convert given positions to zero-based and create the dictionary
position_amino_acid_dict = {
   459: 'D',  
    467: 'T',  
    504: 'E',  
    477: 'M',  
}

filter_sequences(input_file, output_file, position_amino_acid_dict)

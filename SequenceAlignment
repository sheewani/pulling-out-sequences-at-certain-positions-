from Bio import SeqIO

def filter_sequences_by_reference(input_file, output_file, reference_seq, position_amino_acid_dict):
    """
    Filters sequences based on specified positions and amino acids relative to a reference sequence.

    """
    filtered_sequences = []
    
   
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        print(f"Processing sequence {record.id} of length {len(sequence)}")  # Debug print
        
        keep = True
        
        # Check if the sequence meets the criteria
        for pos, amino_acid in position_amino_acid_dict.items():
            if pos-1 >= len(sequence):
                print(f"Position {pos} is out of range for sequence {record.id} with length {len(sequence)}")
                keep = False
                break
            if sequence[pos-1] != amino_acid:
                print(f"Sequence {record.id} does not meet the criteria at position {pos}")
                print(f"Expected {amino_acid} but found {sequence[pos-1]} at position {pos}")
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
reference_seq = "MKFLLQLLLLADPTKLDPSDYVGFTFFVGAMAMMAASAFFFLSLNQFNKKWRTSVLVSGLITFIAAVHYWYMRDYWFAIQESPTFFRYVDWVLTVPLMCVEFYLILKVAGAKPALMWKLIVFSVIMLVTGYFGEAVFQDQAALWGAISGAAYFYIVYEIWLGSAKKIAVAAGGDILKAHKILCWFVLVGWAIYPLGYMLGTDGWYTSILGKGSVDIAYNIADAINKIGFGLVIYALAVKKNEVEVV"  # Replace with your actual reference sequence

# Given positions are 1-based, convert them accordingly
position_amino_acid_dict = {
    90: 'D',  # 91
    94: 'T',  # 95
    101: 'E',  # 102
    98: 'M',  # 99
    226: 'K'  # 227
}

filter_sequences_by_reference(input_file, output_file, reference_seq, position_amino_acid_dict)

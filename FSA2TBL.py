def fasta_to_tbl_with_translation(nucleotide_fasta, amino_fasta, output_file, gene_name, product, note, codon_start=1):

    # Read nucleotide FASTA
    with open(nucleotide_fasta, 'r') as nuc_fasta, open(amino_fasta, 'r') as aa_fasta, open(output_file, 'w') as tbl:
        nuc_sequence = ""
        aa_sequence = ""
        nuc_sequence_name = ""
        aa_sequence_name = ""
        
        aa_lines = aa_fasta.readlines()  # Load the amino acid FASTA file
        aa_idx = 0  # Track the line index for the AA file
        
        for line in nuc_fasta:
            if line.startswith(">"):  # Find the FASTA header in nucleotide file
                if nuc_sequence_name:  # Write the previous sequence feature 
                    bp_start = 1
                    bp_end = len(nuc_sequence)
                    
                    
                    # gene feature
                    tbl.write(f"<1    >{bp_end}    gene\n")
                    tbl.write(f"                        gene           {gene_name}\n")
                    
                    # CDS feature
                    tbl.write(f"<1    >{bp_end}    CDS\n")
                    tbl.write(f"                        product        {product}\n")
                    tbl.write(f"                        note           {note}\n")
                    tbl.write(f"                        codon_start    {codon_start}\n")
                    
                    # Use the corresponding amino acid sequence as translation
                    aa_translation = aa_sequence.strip()  # Strip any trailing newlines
                    tbl.write(f"                        translation    {aa_translation}\n")
                    

                
                # Start a new nucleotide sequence
                nuc_sequence_name = line[1:].strip()
                tbl.write(f">Feature {nuc_sequence_name}\n")
                nuc_sequence = ""  # Reset sequence for the next header
                
                # Move to the corresponding amino acid header and sequence
                while not aa_lines[aa_idx].startswith(">"):  # Find the corresponding AA header
                    aa_idx += 1
                aa_sequence_name = aa_lines[aa_idx][1:].strip()  # Capture AA header
                aa_idx += 1  # Move to the corresponding sequence
                aa_sequence = aa_lines[aa_idx].strip()  # Capture the AA sequence
                aa_idx += 1
            else:
                # Append the current line (sequence) to nucleotide sequence string
                nuc_sequence += line.strip()

        # Handle the last sequence in the file
        if nuc_sequence_name:
            bp_start = 1
            bp_end = len(nuc_sequence)
            
            
            # Add the gene feature
            tbl.write(f"<1    >{bp_end}    gene\n")
            tbl.write(f"                        gene           {gene_name}\n")
            
            # Add the CDS feature
            tbl.write(f"<1    >{bp_end}    CDS\n")
            tbl.write(f"                        product        {product}\n")
            tbl.write(f"                        note           {note}\n")
            tbl.write(f"                        codon_start    {codon_start}\n")
            tbl.write(f"                        translation    {aa_sequence.strip()}\n")
            


#####(nucleotide_fasta, amino_fasta, output_file, gene_name, product, note, codon_start)
fasta_to_tbl_with_translation("AQP1_E0002.fsa", "AQP1_E0002_AA_Codon1.fasta", "AQP1_E0002.tbl", "AQP1", "Aquaporin 1", "RmAQP1 Exon 1", 1)

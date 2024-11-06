import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline

def translate_sequences(input_fasta, output_fasta):
    records = SeqIO.parse(input_fasta, "fasta")
    translated_records = []
    for record in records:
        translated_seq = record.seq.translate(to_stop=True)
        translated_record = SeqRecord(translated_seq, id=record.id, description="translated")
        translated_records.append(translated_record)
    SeqIO.write(translated_records, output_fasta, "fasta")

def align_sequences(input_fasta, output_fasta):
    mafft_cline = MafftCommandline(input=input_fasta)
    stdout, stderr = mafft_cline()
    with open(output_fasta, "w") as handle:
        handle.write(stdout)

def back_translate(aligned_protein_fasta, original_nucleotide_fasta, output_fasta):
    aligned_proteins = SeqIO.parse(aligned_protein_fasta, "fasta")
    original_nucleotides = SeqIO.to_dict(SeqIO.parse(original_nucleotide_fasta, "fasta"))

    with open(output_fasta, "w") as out_f:
        for protein_record in aligned_proteins:
            nucleotide_seq = original_nucleotides[protein_record.id].seq
            back_translated_seq = ""
            nucleotide_index = 0
            for aa in protein_record.seq:
                if aa == "-":
                    back_translated_seq += "---"
                else:
                    back_translated_seq += str(nucleotide_seq[nucleotide_index:nucleotide_index + 3])
                    nucleotide_index += 3
            out_f.write(f">{protein_record.id}\n{back_translated_seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Align protein-coding nucleotide sequences using amino acid translations.")
    parser.add_argument('-i', '--input', required=True, help='Input nucleotide FASTA file')
    parser.add_argument('-p', '--program', required=True, choices=['F'], help='Alignment program (F for MAFFT)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for the aligned files')
    args = parser.parse_args()

    input_fasta = args.input
    output_prefix = args.output
    translated_fasta = f"{output_prefix}_translated.fasta"
    aligned_protein_fasta = f"{output_prefix}_aligned_protein.fasta"
    back_translated_fasta = f"{output_prefix}_aligned_nucleotide.fasta"

    # Translate nucleotide sequences to protein sequences
    translate_sequences(input_fasta, translated_fasta)

    # Align protein sequences
    align_sequences(translated_fasta, aligned_protein_fasta)

    # Back-translate aligned protein sequences to nucleotide sequences
    back_translate(aligned_protein_fasta, input_fasta, back_translated_fasta)

    print(f"Alignment completed. Aligned nucleotide sequences written to {back_translated_fasta}")

if __name__ == "__main__":
    main()
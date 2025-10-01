# --- Anwesha Awasthi ---
# Biopython Project Pipeline (Nucleotide Only)

# Import required modules
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import gc_fraction
import urllib.error

# Always set email for NCBI access
Entrez.email = "anweshaawasthi@gmail.com"

# Example nucleotide accession (change to your assigned one)
nuc_accession = "NM_001301717"

# --- 1. Fetch Nucleotide Sequence from NCBI ---
print("\n--- Fetching Nucleotide Sequence ---")
handle = Entrez.efetch(db="nucleotide", id=nuc_accession, rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

print("Accession:", record.id)
print("Description:", record.description)
print("Sequence length:", len(record.seq))

# --- 2. Basic Sequence Analysis ---
seq = record.seq
print("\n--- Sequence Analysis ---")
print("Original Sequence:", seq)
print("Complement:", seq.complement())
print("Reverse Complement:", seq.reverse_complement())
print("Transcription (DNA → RNA):", seq.transcribe())
print("Translation (RNA → Protein):", seq.translate(to_stop=True))

# GC content
gc_content = gc_fraction(seq) * 100
print("GC Content:", round(gc_content, 2), "%")

# --- 3. Global and Local Alignments ---
aligner = PairwiseAligner()

# Global Alignment
aligner.mode = "global"
global_aln = aligner.align(seq, seq)[0]
print("\n--- Global Alignment ---")
print("Score:", global_aln.score)
print(global_aln)

# Local Alignment
aligner.mode = "local"
local_aln = aligner.align(seq, seq)[0]
print("\n--- Local Alignment ---")
print("Score:", local_aln.score)
print(local_aln)

# --- 4. BLASTn for Nucleotide Sequence ---
try:
    print("\n--- Running BLASTn on database nt ---")
    result_handle_n = NCBIWWW.qblast("blastn", "nt", seq)
    blast_record_n = NCBIXML.read(result_handle_n)

    # Display top 3 alignments
    print("\nTop 3 Alignments for BLASTn:")
    for alignment in blast_record_n.alignments[:3]:
        print(f"\nSequence: {alignment.title}")
        print(f"Length: {alignment.length}")
        for hsp in alignment.hsps:
            print(f"Score: {hsp.score}, E-value: {hsp.expect}")
            print(f"Query: {hsp.query}")
            print(f"Match: {hsp.match}")
            print(f"Subject: {hsp.sbjct}")

except urllib.error.URLError as e:
    print(f"URL Error: {e.reason}")
except Exception as e:
    print(f"An error occurred: {e}")

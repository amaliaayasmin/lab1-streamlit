import streamlit as st
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Set up your Entrez email
Entrez.email = "amaliayasmin@graduate.utm.my"  # Replace with your own email address

# 1. Custom Function: Retrieve Protein Sequence Data
def retrieve_data(uniprot_id):
    try:
        handle = Entrez.efetch(db='protein', id=uniprot_id, rettype='fasta', retmode='text')
        record = SeqIO.read(handle, 'fasta')
        handle.close()
        return record
    except Exception as e:
        st.error(f"Error retrieving data: {e}")
        return None

# 2. Custom Function: Perform Basic Analysis on Protein Sequence
def get_basic_analysis(sequence):
    seq_analysis = ProteinAnalysis(str(sequence))

    sequence_length = len(sequence)
    amino_acid_composition = seq_analysis.count_amino_acids()
    molecular_weight = seq_analysis.molecular_weight()
    isoelectric_point = seq_analysis.isoelectric_point()
    
    return {
        "Sequence Length": sequence_length,
        "Amino Acid Composition": amino_acid_composition,
        "Molecular Weight": molecular_weight,
        "Isoelectric Point": isoelectric_point
    }

# Streamlit App Interface
st.title('Lab 1 - AMALIA YASMIN')  

# Input field for UniProt ID
protein_id = st.text_input('Enter UniProt ID')
retrieve = st.button('Retrieve')

# When 'Retrieve' button is clicked
if retrieve:
    if protein_id:
        # Call the function to retrieve protein data
        record = retrieve_data(protein_id)
        
        if record:
            
            st.subheader("Protein Information")
            st.text(f"Description: {record.description}")
            st.text(f"Sequence Length: {len(record.seq)}")  
            st.text(f"Sequence: {record.seq}") 

            analysis_outcome = get_basic_analysis(record.seq)

            st.subheader("Protein Analysis")
            st.write("Sequence Length:", analysis_outcome["Sequence Length"])  
            st.write("Amino Acid Composition:", analysis_outcome["Amino Acid Composition"])
            st.write("Molecular Weight:", analysis_outcome["Molecular Weight"])
            st.write("Isoelectric Point:", analysis_outcome["Isoelectric Point"])
    else:
        st.warning('Please enter a UniProt ID')

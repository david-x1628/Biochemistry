"""Solves tasks for Katharina's biochemistry exercises."""

import io
import Bio
from Bio import SeqIO
from Bio import SeqRecord
#from Bio import Align
import requests  # for 'API calls', i. e. for retrieving information from an internet endpoint
import numpy  # for matrix
import pandas as pd  # for matrix as 'DataFrame', i. e. as tabular data type


print("Testdurchlauf:\n")


# TASK 1

def fetch_uniprot_fasta(accession: str, timeout: int = 10) -> Bio.SeqRecord | None:
    """Fetche and parse FASTA from UniProt.
    
    Parameters:
        accession (str): Accession of wanted protein or DNA.
        timeout (int): Timeout for HTTP request.
        
    Returns:
        record (SeqRecord): SeqRecord object of specified protein or DNA or
        None if fetching or parsing fails.
    """
    
    # Set endpoint for API call
    # -> API call = communication between programs
    url = f'https://rest.uniprot.org/uniprotkb/{accession}.fasta'
    

    try:  # try to run the code below and if an error occurs, run the code below 'except'  
        
        # API call:
        # requests.get(url) gets the information provided at the url;
        # 'timeout' determines that if there is no response, the call will
        # be canceled;
        # response from the url endpoint is stored in 'response'
        response = requests.get(url, timeout=timeout).text

        # convert response to FASTA format amd store in 'record'
        record = SeqIO.read(io.StringIO(response), 'fasta')

        return record

    except:  # run if an exception (= error at runtime) occurs
        return None


# example: HUMAN Hemoglobin subunit alpha
P69905 = fetch_uniprot_fasta('P69905')  # quasi wird f(P69905) der Name P69905 gegeben ()

# example: HUMAN Alpha-1-globin
D1MGQ2 = fetch_uniprot_fasta('D1MGQ2')

# example: HUMAN Relaxin-3
Q8WXF3 = fetch_uniprot_fasta('Q8WXF3')

# example: HUMAN Alpha-lactalbumin
P00709 = fetch_uniprot_fasta('P00709')

# example: HUMAN Baculoviral IAP repeat-containing protein 5
O15392 = fetch_uniprot_fasta('O15392')

print("Wir nutzen als Beispiel:\n\n", P69905, "\n\n", D1MGQ2)


# TASK 2

def summarize_record(record: Bio.SeqRecord) -> dict[str, str | int]: 
    """Summarise record in dictionary.
    
    Parameters:
        record (SeqRecord): Record to be summarised.
    
    Returns:
        summary (dict): Summarised record containing record's
            id, length, description and seq_start (first 60
            residues) as keys.
    """

    summary = {
        'id': record.id,  # .id accesses 'id' attribute of record
        'length': len(record.seq),  # len() determines length of 'sequence' attribute of record
        'description': record.description,  # .description accesses attribute 'description' of record
        'seq_start': record.seq[:60]  #  [:59] accesses the first 60 positions of the sequence (process is called 'slicing')  
        }
    
    return summary


#example: HUMAN Hemoglobin subunit alpha
P = summarize_record(P69905)

#example: HUMAN Alpha-1-globin
D = summarize_record(D1MGQ2)

print("\n\nDie Zusammenfassungen als Dictionaries sind:\n")
print(P)
print(D)


# TASK 3.1

def pairwise_identity(seq1: Bio.Seq, seq2: Bio.Seq):
    """Calculate percent identity between two sequences.
    
    Parameters: 
        seq1 (Seq): Sequence to compare to seq2.
        seq2 (Seq): Sequence to compare to seq1.
    
    Returns:
        percent_identity (int): Percent identity between seq1 and seq2.
    """

    # check if the length of the two given sequences is the same
    if len(seq1) != len(seq2):
        raise ValueError("The given sequences do not have the same length.")

    aligner = Bio.Align.PairwiseAligner()  # initialise PairwiseAligner object
    score = aligner.score(seq1, seq2)  # calculate score of the given sequences
    percent_identity = (score / len(seq1)) * 100  # calculate identity in percent

    return percent_identity


#print("\nDie prozentuale Identität der beiden Beispiele beträgt: ", pairwise_identity(P69905.seq, D1MGQ2.seq))


# TASK 3.2  -  unlabeled DataFrame matrix

def build_identity_matrix(records: list[Bio.SeqRecord]) -> pd.DataFrame:
    """Build matrix of percentage identities of all records.
    
    Parameters:
        records (list[SeqRecord]): All records to be included in the matrix.
        
    Returns:
        identity_matrix (DataFrame): Matrix containing all identity alignment
            comparisons of all given records with each other.
    """

    data = []
    for r in records:
        series = pd.Series([pairwise_identity(r.seq, s.seq) for s in records])
        data.append(series)
    
    identity_matrix = pd.concat(data, axis=1)

    return identity_matrix


# TASK 3.2  -  with labelled DataFrame matrix

def build_identity_matrix_labelled(records: list[Bio.SeqRecord]) -> pd.DataFrame:
    """Build matrix of percentage identities of all records.
    
    Parameters:
        records (list[SeqRecord]): All records to be included in the matrix.
        
    Returns:
        identity_matrix (DataFrame): Matrix containing all identity alignment
            comparisons of all given records with each other.
    """

    data = {'': [r.id for r in records]}
    for r in records:
        identities = []
        for s in records:
            identity = pairwise_identity(r.seq, s.seq)
            identities.append(identity)
        
        data.update({str(r.id): identities})
    
    identity_matrix = pd.DataFrame(data=data)

    return identity_matrix


records = records = [P69905, Q8WXF3, P00709, O15392]

print("\n\n\n")
print("Gewichtematrix mit einigen Beispielen:\n\n")
#print(build_identity_matrix(records))
print("\n\n\n")


# TASK 4

def find_conserved_columns(records: Bio.SeqRecord) -> list[int]:
    """Finds all positions that are conserved accross all provided records.
    
    Conserved refers to the same residue in every sequence (ignoring gaps).

    Parameters:
        records (SeqRecord): Records of which the conserved columns
            are to be searched.
    
    Returns:
        conserved_columns (list[int]): List containing all conserved
            positions of the specified records.
    """

    i = 0
    conserved_columns = []

    for position in records[0].seq:  # loop through the first sequence
        letters = [position]  # initialise list and store one position to compare the others to
        
        for record in records:  # loop through all records
           
            if record.id == records[0].id:  # skip first record because first serves as comparison
                continue
            
            if record.seq[i] in letters:  # check if the i. letter is equal to 'position'
                
                if i not in conserved_columns:  # if already appended, do not append twice

                    conserved_columns.append(i)  # append position that is equal among all compared sequences

            else:  # for every position that is not equal among all sequences

                try: 

                    conserved_columns.pop(conserved_columns.index(i))  # delete this position from the list

                except ValueError:  # filter for those positions that never made it to the list

                    pass

                break  # end loop for this position since it is not equal among all sequences
            
        i += 1
    
    return conserved_columns

print("Und welche Positionen sind in mehreren gegebenen Sequenzen gleich?\n\n")
print(find_conserved_columns([P69905, P00709]))
        
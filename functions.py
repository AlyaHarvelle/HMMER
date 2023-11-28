# Import the libraries
import os
import re
import json
import ipywidgets as widgets
from collections import Counter
import scipy.stats
from scipy.cluster.hierarchy import linkage, dendrogram
import numpy as np
import pandas as pd
import subprocess
from Bio.Blast import NCBIXML
from Bio import AlignIO, SeqIO
from Bio.AlignIO import ClustalIO
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
import subprocess

directory = os.getcwd()

# Set input parameters
name = 'alina'
server = 'ecate'

# Parses the json file and returns a pandas dataframe
def json_parser(line):
    rows = list()
    for key, value in line.items():
        if key == 'acc':
            acc = value
        else:
            evidence, feature, source = key.split('-')
            for k, v in value.items():
                if k == 'regions':
                    for item in range(len(v)):
                        reg = v[item]
                        start, end = reg[0], reg[1]
                        rows.append({
                            'acc': acc,
                            'evidence': evidence,
                            'feature': feature,
                            'source': source,
                            'start': start,
                            'end': end
                        })

    return rows

# Parses BLAST alignment in XML format and returns a pandas dataframe
def blast_parser(input_file):
    max_query_length = 0
    data = []
    with open(input_file) as f:
        blast_records = NCBIXML.parse(f)

        # Iterate PSIBLAST rounds
        for blast_record in blast_records:
            query_id = blast_record.query

            # Get the length of the query sequence
            query_length = int(blast_record.query_letters)
            if query_length > max_query_length:
                max_query_length = query_length

            # Iterate alignments (here just one)
            for i, alignment in enumerate(blast_record.alignments):
                subject_id = alignment.title

                # Iterate pairwise alignments
                for hsp in alignment.hsps:
                    data.append((query_id.split('|')[1].split(' ')[0],
                                 subject_id.split('|')[1].split(' ')[0],
                                 query_length,
                                 hsp.align_length,
                                 hsp.query,
                                 hsp.sbjct,
                                 hsp.query_start,
                                 hsp.query_end,
                                 hsp.sbjct_start,
                                 hsp.sbjct_end,
                                 hsp.expect,
                                 hsp.score,
                                ))

                    # Skip duplicated subjects
                    break

        df = pd.DataFrame(data, columns=['query_id', 'subject_id', 'query_len', 'hsp_len', 
                                         'query_seq','subject_seq', 'query_start', 'query_end',
                                         'subject_start', 'subject_end',
                                         'eval', 'bit_score'])
        
        df = df.sort_values('eval', ascending=True)
        grouped = df.groupby('query_id')['subject_id'].nunique().reset_index(name='count')
        df = pd.merge(df, grouped, on='query_id')

        return df

# Calculate the statistics
# Returns the values of occupancy and entropy for each alignment
def stats_calculation(seqs, q_id):
    
    data = []
    aa = 'ACDEFGHIKLMNPQRSTVWY'

    for i, column in enumerate(seqs.T):

        count = Counter(column)
        try:
            count.pop('-')
        except KeyError:
            pass
        count_sorted = sorted(count.items(), key=lambda x:x[1], reverse=True)

        non_gap = np.count_nonzero(column != '-')
        occupancy = non_gap / column.size

        probabilities = [count.get(k, 0.0) / column.size for k in aa]

        entropy = scipy.stats.entropy(probabilities, base=20)
        data.append([i, q_id, occupancy, entropy, count_sorted])

    df_calc = pd.DataFrame(data, columns=['pos', 'query_id', 'occupancy', 'entropy', 'counts'])
    return df_calc

# Calculate statistics of the MSAs - returns a list of dataframes for each MSA
def process_files(folder_path, method):
    statistics_list = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith('.fasta'):
            protein_id = os.path.splitext(filename)[0].split('_')[0]
            seq_data = print_dis_seqs(folder_path, method, protein_id)
            stats = stats_calculation(seq_data, protein_id)
            statistics_list.append(stats)
            
    statistics_df = pd.concat(statistics_list, ignore_index=True)
            
    return statistics_df

# Function retrieving a FASTA sequence from ID
def get_fasta(upac): 
    result = subprocess.run(['ssh', 'ecate', '/software/packages/ncbi-blast/latest/bin/blastdbcmd', '-db', '/db/blastdb/uniprot/uniprot.fasta', '-entry', upac], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return result.stdout

# ClustalOmega MSA generator
def clustalo_generator(input_folder, output_folder):
    # Iterate over all input files in the input folder
    for input_file in os.listdir(input_folder):
        if input_file.endswith('.fasta'):
            # Extract ID from the file name
            id_split = os.path.splitext(input_file)[0].split('_')[0]
            print(f'MSA ClustalOmega is generated for the {id_split} protein')
            output_file = os.path.join(output_folder, f'{id_split}_clustal.fasta') # name of the output file

            # Define the ClustalOmega command
            clustalomega_cline = ClustalOmegaCommandline(
                infile=os.path.join(input_folder, input_file),
                outfile=output_file,
                outputorder='input-order',
                verbose=False,
                auto=True,
                force=True)

            # Run the ClustalOmega command
            subprocess.run(str(clustalomega_cline), shell=True)
            
# Remove gaps from the query sequence of the alignment
def remove_extra_gaps(alignment):
    # Get the first sequence
    first_seq = alignment[0].seq

    # Find the positions of non-gap characters in the first sequence
    non_gap_positions = [i for i, base in enumerate(first_seq) if base != '-']

    # Create a list to hold SeqRecord objects for the filtered alignment
    filtered_seqs = []

    for seq_record in alignment:
        # Extract and join the non-gap characters
        filtered_seq = Seq(''.join(seq_record.seq[i] for i in non_gap_positions), seq_record.seq.alphabet)
        
        # Create a new SeqRecord with the filtered sequence
        filtered_record = SeqRecord(
            seq=filtered_seq,
            id=seq_record.id,
            description=seq_record.description)
        
        filtered_seqs.append(filtered_record)

    # Create a new alignment with filtered SeqRecord objects
    new_alignment = MultipleSeqAlignment(filtered_seqs, alphabet=alignment._alphabet)

    return new_alignment

# Apply gaps removing to all files
def process_folder(input_folder):
    # Iterate over all input files in the input folder
    for input_file in os.listdir(input_folder):
        if input_file.endswith('.fasta'):
            # Set file paths
            file_path = os.path.join(input_folder, input_file)

            # Read the alignment and remove extra gaps
            alignment = AlignIO.read(file_path, 'fasta')
            filtered_alignment = remove_extra_gaps(alignment)

            # Write the filtered alignment to a new file
            with open(file_path, "w") as output_handle:
                AlignIO.write(filtered_alignment, output_handle, 'fasta')


# Saves disordered regions from MSA output in separate FASTA files
def select_dis_regions(msa_file, query_id, start_positions, end_positions, region, output_file):
    sep_sequences = []

    with open(msa_file, 'r') as msa_handle:
        msa_records = list(SeqIO.parse(msa_handle, 'fasta'))

    for i, (start, end) in enumerate(zip(start_positions, end_positions), start=1):
        records = []
        query_record_id = str(query_id)
        subject_record_ids = []

        for j, record in enumerate(msa_records):
            sequence = record.seq
            seq_length = len(sequence)
            if seq_length >= start > 0 and end <= seq_length:
                trimmed_sequence = sequence[start - 1: end]  # Adjust the slicing here

                if j == 0:
                    record_id = query_record_id
                    description = ''
                else:
                    if record.id not in subject_record_ids:
                        subject_record_ids.append(record.id)
                    record_id = subject_record_ids[-1]
                    description = record.description

                disordered_record = SeqIO.SeqRecord(trimmed_sequence, id=record_id, description=description)
                records.append(disordered_record)
            else:
                print(f'Invalid region: start={start}, end={end}, Seq Length: {seq_length}')

        if records:
            SeqIO.write(records, output_file, 'fasta')

        sep_sequences.extend([record.seq for record in records])  # Extend the collected sequences

    return sep_sequences

# Creates a numpy array from the .fasta file
def get_seqs(aligns):
    seqs = []
    if isinstance(aligns, str):  # If the input is a file path
        with open(aligns) as f:
            for record in AlignIO.read(f, 'fasta'):
                seqs.append(list(record.seq))
        seqs = np.array(seqs, dtype='str')
    else:  # If the input is a list of sequences
        seqs = np.array(aligns, dtype='str')
    return seqs

# Prints the sequences of disordered regions of an alignment
def print_dis_seqs(directory, align_type, query_id):
    all_seqs = []

    for file_name in os.listdir(directory):
        if f'{query_id}_' in file_name and f'_disordered' not in file_name:
            alignment_file = os.path.join(directory, file_name)
            seqs = []

            with open(alignment_file) as f:
                for record in AlignIO.read(f, 'fasta'):
                    seqs.append(np.array(list(record.seq), dtype='str'))

            all_seqs.append(np.array(seqs))

    if len(all_seqs) == 1:
        return np.array(all_seqs[0])
    else:
        return all_seqs

# Function to plot occupancy and entropy distribution for each MSA
def plot_occupancy_entropy_distribution(stats, num_rows=4):
    unique_query_ids = np.sort(stats['query_id'].unique())

    # Calculate the number of subplots in each row
    num_plots_in_row = np.ceil(len(unique_query_ids) / num_rows).astype(int)

    # Create subplots for each unique query_id
    fig, axes = plt.subplots(num_rows, num_plots_in_row, figsize=(20, 8))

    for i, query_id in enumerate(unique_query_ids):
        subset = stats[stats['query_id'] == query_id]
        row_index = i // num_plots_in_row  # Determine the row for the current instance
        col_index = i % num_plots_in_row   # Determine the column for the current instance

        sns.kdeplot(subset['occupancy'], shade=True, ax=axes[row_index, col_index], label='Occupancy', warn_singular=False)
        sns.kdeplot(subset['entropy'], shade=True, ax=axes[row_index, col_index], label='Entropy', warn_singular=False)

        axes[row_index, col_index].set_ylabel('')
        axes[row_index, col_index].set_xlabel('')

        axes[row_index, col_index].set_title('ID: {}'.format(query_id))

    fig.legend(loc='upper left', labels=['Occupancy', 'Entropy'])

    plt.suptitle('Occupancy and entropy distribution for each Uniprot ID')
    plt.tight_layout()
    plt.show()

# Preprocess the Intepro outputs
def pfam_processing(filename):
    data = []
    with open(filename) as file:
        for line in file:
            row = line.strip().split("\t")
            uniprot_id = row[0]
            pfam_id = row[1]
            ipr_id = row[2]
            start_pfam = row[3]
            end_pfam = row[4]
                
            # Choose a Pfam ID
            if pfam_id.startswith("PF"):
                data.append([uniprot_id, pfam_id, ipr_id, start_pfam, end_pfam])

    columns = ["uniprot_id", "pfam_id", "ipr_id", "start_pfam", "end_pfam"]
    pfam = pd.DataFrame(data, columns=columns)
    pfam["start_pfam"] = pd.to_numeric(pfam["start_pfam"])
    pfam["end_pfam"] = pd.to_numeric(pfam["end_pfam"])
    pfam["length_pfam"] = pfam["end_pfam"] - pfam["start_pfam"] + 1

    return pfam

# Calculate redundancy 
def calculate_red(msa_file, output_file, threshold, word_size, id_split):

    cd_hit_path = '/Users/alina/cd-hit/cd-hit'

    # Run CD-Hit to cluster the sequences and remove redundancy
    cmd = f'{cd_hit_path} -i {msa_file} -o {output_file} -c {threshold} -n {word_size} > /dev/null'
    subprocess.call(cmd, shell=True)

    # Store the non-redundant sequences in a list
    non_redundant_sequences = []
    with open(output_file, 'r') as output_handle:
        for record in SeqIO.parse(output_handle, 'fasta'):
            non_redundant_sequences.append(record)

    # Write the sequences to the output file
    with open(output_file, 'w') as final_handle:
        SeqIO.write(non_redundant_sequences, final_handle, 'fasta')

    # Count the number of sequences in the original MSA and the non-redundant MSA
    total_sequences = sum(1 for record in SeqIO.parse(msa_file, 'fasta'))
    non_redundant_sequences_count = len(non_redundant_sequences)

    # Count the effective sequences (Nf)
    Nf = non_redundant_sequences_count / total_sequences
    print(f'Non-redundant seqs for {id_split}:', non_redundant_sequences_count, 
          'Total no. of seqs:', total_sequences,
          'Ratio:', '{:.2f}'.format(Nf))
    
    return

# Function to plot occupancy and entropy
def plot_disordered_regions(id, l_cl, l_nr, dis_calc_num):
    # Define the data frame names
    dis_calc_name = f"dis_calc{dis_calc_num}"
    dis_calc_nr_name = f"dis_calc_nr{dis_calc_num}"

    # Define the figure and axes
    fig, ax = plt.subplots(3, 2, figsize=(10, 10))

    # KDE plot of occupancy/entropy of disordered region
    sns.kdeplot(eval(dis_calc_name)['occupancy'], shade=True, ax=ax[0, 0], warn_singular=False)
    sns.kdeplot(eval(dis_calc_nr_name)['occupancy'], shade=True, ax=ax[0, 0], warn_singular=False)
    sns.kdeplot(eval(dis_calc_name)['entropy'], shade=True, ax=ax[0, 1], warn_singular=False)
    sns.kdeplot(eval(dis_calc_nr_name)['entropy'], shade=True, ax=ax[0, 1], warn_singular=False)

    # Histogram of occupancy/entropy redundant vs non-redundant
    ax[1, 0].hist([eval(dis_calc_name)['occupancy'], eval(dis_calc_nr_name)['occupancy']])
    ax[1, 1].hist([eval(dis_calc_name)['entropy'], eval(dis_calc_nr_name)['entropy']])

    # Scatterplot of occupancy/entropy redundant vs non-redundant
    ax[2, 0].scatter(np.sort(eval(dis_calc_name)['occupancy']), eval(dis_calc_name)['pos'])
    ax[2, 0].scatter(np.sort(eval(dis_calc_nr_name)['occupancy']), eval(dis_calc_nr_name)['pos'])

    ax[2, 1].scatter(np.sort(eval(dis_calc_name)['entropy']), eval(dis_calc_name)['pos'][::-1])
    ax[2, 1].scatter(np.sort(eval(dis_calc_nr_name)['entropy']), eval(dis_calc_nr_name)['pos'][::-1])

    # Add x-axis labels
    x_axis_labels = ['Occupancy', 'Entropy']
    for i in range(3):
        for j in range(2):
            ax[i, j].set_xlabel(x_axis_labels[j])

    # Add the legend
    legend_labels = [f'Initial MSA ({l_cl} sequences)', f'Non-redundant MSA ({l_nr} sequences)']
    for ax in ax.flat:
        ax.legend(legend_labels)

    plt.suptitle(f'KDE plot, histogram, and scatterplot of occupancy and entropy distribution for {id}_{dis_calc_num+1} protein')
    plt.tight_layout()
    plt.show()
    
# Extract the table information from the hmmsearch output
def process_hmmsearch_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Extract the column names
    column_names = lines[12].split()
    inclusion_threshold_index = lines.index('  ------ inclusion threshold ------\n')

    # Extract the data rows
    data_rows = [line.split() for line in lines[14:inclusion_threshold_index]]
    data_rows = [row[:9] + [' '.join(row[9:])] for row in data_rows]

    # Create the DataFrame
    stats = pd.DataFrame(data_rows, columns=column_names)

    # Print the total number of hits and unique sequences
    unique_sequences = stats.Sequence.nunique()

    return stats

# Extract the information regarding hmm and alignment regions
def extract_table_from_output(hmm_result_file):
    with open(hmm_result_file, 'r') as f:
        lines = f.readlines()

        # Extract a query ID
        query_id = None
        extract_info = False
        result_data = []

        for line in lines:
            if line.startswith(">>"):
                query_id = line.strip()[2:].split(" ")[1]
                extract_info = True  # Set the flag to True when query ID is found
            elif extract_info:
                line_data = line.split()
                if len(line_data) >= 2 and line_data[1] == "!":
                    try:
                        hmm_from = int(line_data[6])
                        hmm_to = int(line_data[7])
                        hmm_length = (hmm_to - hmm_from + 1)
                        ali_from = int(line_data[9])
                        ali_to = int(line_data[10])
                        ali_length = (ali_to - ali_from + 1)
                        env_from = int(line_data[12])
                        env_to = int(line_data[13])
                        env_length = (env_to - env_from + 1)
                        result_data.append([query_id, hmm_from, hmm_to, hmm_length,
                                            ali_from, ali_to, ali_length, 
                                            env_from, env_to, env_length])
                    except ValueError:
                        pass
                    extract_info = False

        # Create and return the DataFrame
        result_df = pd.DataFrame(result_data, columns=['id', 'hmm_from', 'hmm_to', 'hmm_length',
                                                       'ali_from', 'ali_to', 'ali_length',
                                                       'env_from', 'env_to', 'env_length'])
        return result_df

# Plot overlaps between HMM and alignment regions
def plot_overlaps(hmm_results, curated_query, id_dis):
    unique_uniprot_ids = hmm_results['Sequence'].unique()

    if len(unique_uniprot_ids) == 0:
        print("No overlapping regions to plot.")
        return

    # Plot overlapping regions
    fig, ax = plt.subplots(figsize=(10, 0.25 * len(unique_uniprot_ids)))

    # Plot the regions using 'hmm_from' and 'hmm_to' for Pfam, and 'ali_from' and 'ali_to' for Disprot
    ax.hlines(hmm_results['Sequence'], hmm_results['hmm_from'], hmm_results['hmm_to'], linewidth=10, color='blue', label='HMM region', linestyle='-')
    ax.hlines(hmm_results['Sequence'], hmm_results['ali_from'], hmm_results['ali_to'], linewidth=10, color='green', label='Alignment region', linestyle='-')

    ax.set_yticks(hmm_results['Sequence'])
    ax.set_yticklabels(hmm_results['Sequence'])

    # Create a custom range for the y-axis based on the unique 'uniprot_id' values
    y_axis_range = range(len(unique_uniprot_ids))
    ax.set_ylim(min(y_axis_range) - 1, max(y_axis_range) + 1)
    
    # Add annotations at specific positions
    for idx, row in hmm_results.iterrows():
        x_pos = (row['hmm_from'] + row['hmm_to']) / 2  # Calculate the x position for annotation using 'hmm' coordinates
        y_pos = row['Sequence']  # Use the 'uniprot_id' as the y position

    plt.title(f'Overlaps between the HMM and alignment regions of the {id_dis} protein')
    plt.xlabel('Position')
    plt.ylabel('UniProt ID')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.show()
    
# Calculate the overlaps percentage
def calculate_overlap(row_pfam, row_disprot):
    start_pfam = row_pfam['start_pfam']
    end_pfam = row_pfam['end_pfam']
    start_hmm = row_pfam['ali_from']
    end_hmm = row_pfam['ali_to']
    start_disprot = row_disprot['start'] 
    end_disprot = row_disprot['end']
    disprot_length = end_disprot - start_disprot + 1

    overlap_pfam_disprot = min(end_pfam, end_disprot) - max(start_pfam, start_disprot) + 1
    overlap_hmm_disprot = min(end_hmm, end_disprot) - max(start_hmm, start_disprot) + 1

    if overlap_pfam_disprot > 0:
        percentage_overlap_pfam = (overlap_pfam_disprot / disprot_length) * 100
    else:
        percentage_overlap_pfam = 0
    
    if overlap_hmm_disprot > 0:
        percentage_overlap_hmm = (overlap_hmm_disprot / disprot_length) * 100
    else:
        percentage_overlap_hmm = 0

    return percentage_overlap_pfam, percentage_overlap_hmm, overlap_pfam_disprot, overlap_hmm_disprot
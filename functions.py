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
                    data.append((query_id.split("|")[1].split(" ")[0],
                                 subject_id.split("|")[1].split(" ")[0],
                                 query_length,
                                 hsp.align_length,
                                 hsp.query,
                                 hsp.match,
                                 hsp.sbjct,
                                 hsp.query_start,
                                 hsp.query_end,
                                 hsp.sbjct_start,
                                 hsp.sbjct_end,
                                 hsp.identities,
                                 hsp.positives,
                                 hsp.gaps,
                                 hsp.expect,
                                 hsp.score,
                                ))

                    # Skip duplicated subjects
                    break

        df = pd.DataFrame(data, columns=["query_id", "subject_id", "query_len", "hsp_len", "query_seq",
                                               "match_seq", "subject_seq", "query_start", "query_end",
                                               "subject_start", "subject_end", "identity", "positive",
                                               "gaps", "eval", "bit_score"])
        df = df.sort_values('eval', ascending=True)
        grouped = df.groupby('query_id')['subject_id'].nunique().reset_index(name='count')
        df = pd.merge(df, grouped, on='query_id')

        return df

# Returns the values of occupancy and entropy for each alignment
def stats_calculation(seqs):
    
    data = []
    aa = "ACDEFGHIKLMNPQRSTVWY"

    for i, column in enumerate(seqs.T):

        count = Counter(column)
        try:
            count.pop('-')
        except KeyError:
            pass
        count_sorted = sorted(count.items(), key=lambda x:x[1], reverse=True)

        non_gap = np.count_nonzero(column != "-")
        occupancy = non_gap / column.size

        probabilities = [count.get(k, 0.0) / column.size for k in aa]

        entropy = scipy.stats.entropy(probabilities, base=20)
        data.append([i, occupancy, entropy, count_sorted])

    df_calc = pd.DataFrame(data, columns=['pos', 'occupancy', 'entropy', 'counts'])
    return df_calc

# Function retrieving a FASTA sequence from ID
def get_fasta(upac): 
    result = subprocess.run(['ssh', 'ecate', '/software/packages/ncbi-blast/latest/bin/blastdbcmd', '-db', '/db/blastdb/uniprot/uniprot.fasta', '-entry', upac], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return result.stdout

# Returns a FASTA file containing disordered regions only.
# def trim_dis_regions(file, query_id, start_positions, end_positions):
# 
#     records = []
#     for record in SeqIO.parse(file, "fasta"):
#         sequence = record.seq
#         trimmed_sequence = ""
#         previous_end = 0  # Track the end position of the previous region
#         for start, end in zip(start_positions, end_positions):
#             if len(sequence) >= start > previous_end and end <= len(sequence):
#                 trimmed_sequence += '-' * (start - previous_end - 1) + sequence[start - 1: end]
#                 previous_end = end  # Update the previous end position
#             else:
#                 print("Invalid region: start =", start, "end =", end)
#         trimmed_sequence += '-' * (len(sequence) - previous_end)  # Add trailing gaps
#         record.seq = trimmed_sequence
#         records.append(record)
#     SeqIO.write(records, f"{directory}/results/alignments/output_files/disordered/{query_id}_disordered.fasta", "fasta")

#     trimmed_sequences = [record.seq for record in records]  # Collect trimmed sequences
#     return trimmed_sequences

# Saves disordered regions from MSA output in separate FASTA files
def select_dis_regions(msa_file, query_id, start_positions, end_positions, output_directory):

    sep_sequences = []  # Collect trimmed sequences

    with open(msa_file, "r") as msa_handle:
        msa_records = list(SeqIO.parse(msa_handle, "fasta"))

    for i, (start, end) in enumerate(zip(start_positions, end_positions), start=1):
        records = []
        query_record_id = str(query_id)
        subject_record_ids = []

        for j, record in enumerate(msa_records):
            sequence = record.seq
            if len(sequence) >= start > 0 and end <= len(sequence):
                trimmed_sequence = sequence[start - 1: end]

                if j == 0:
                    record_id = query_record_id
                    description = ""
                else:
                    if record.id not in subject_record_ids:
                        subject_record_ids.append(record.id)
                    record_id = subject_record_ids[-1]
                    description = record.description

                disordered_record = SeqIO.SeqRecord(trimmed_sequence, id=record_id, description=description)
                records.append(disordered_record)
            else:
                print(f"Invalid region: start={start}, end={end}")

        if records:
            output_file_separate = os.path.join(output_directory, f"{query_id}_{i}.fasta")
            SeqIO.write(records, output_file_separate, "fasta")

        sep_sequences.extend([record.seq for record in records])  # Extend the collected sequences

    return sep_sequences

# Creates a numpy array from the .fasta file
def get_seqs(aligns):
    seqs = []
    if isinstance(aligns, str):  # If the input is a file path
        with open(aligns) as f:
            for record in AlignIO.read(f, "fasta"):
                seqs.append(list(record.seq))
        seqs = np.array(seqs, dtype="str")
    else:  # If the input is a list of sequences
        seqs = np.array(aligns, dtype="str")
    return seqs

# Prints the sequences of disordered regions of an alignment
def print_dis_seqs(directory, align_type, query_id):
    all_seqs = []

    for file_name in os.listdir(directory):
        if f'{query_id}_' in file_name and f'_disordered' not in file_name:
            alignment_file = os.path.join(directory, file_name)
            seqs = []

            with open(alignment_file) as f:
                for record in AlignIO.read(f, "fasta"):
                    seqs.append(np.array(list(record.seq), dtype="str"))

            all_seqs.append(np.array(seqs))

    if len(all_seqs) == 1:
        return np.array(all_seqs[0])
    else:
        return all_seqs

# Processes the output of the protein2ipr results with Pfam instances
def read_and_filter_pfam_data(filename, filtered_curated):
    data = []
    with open(filename) as file:
        for line in file:
            row = line.strip().split("\t")
            if len(row) >= 6:
                uniprot_id = row[0]
                ipr_id = row[1]
                description = row[2]
                pfam_id = row[3]
                start_pos = row[4]
                end_pos = row[5]
                
                # Choose a Pfam ID
                if pfam_id.startswith("PF"):
                    data.append([uniprot_id, ipr_id, description, pfam_id, start_pos, end_pos])

    columns = ["uniprot_id", "ipr_id", "description", "pfam_id", "start_pos", "end_pos"]
    pfam = pd.DataFrame(data, columns=columns)
    pfam["start_pos"] = pd.to_numeric(pfam["start_pos"])
    pfam["end_pos"] = pd.to_numeric(pfam["end_pos"])
    pfam['length'] = pfam['end_pos'] - pfam['start_pos'] + 1

    return pfam

# Calculate redundancy 
def calculate_Nf(msa_file, threshold, id_dis):

    output_file = f"/Users/alina/HMM/results/alignments/input_files/non-redundant/Nf_{id_dis}.fasta"
    cd_hit_path = "/Users/alina/cd-hit/cd-hit"

    # Run CD-HIT to cluster the sequences (excluding the first line) and remove redundancy
    cmd = f"{cd_hit_path} -i {msa_file} -o {output_file} -c {threshold} -n 4 > /dev/null"
    subprocess.call(cmd, shell=True)

    # Read the first line from the original MSA file
    with open(msa_file, "r") as msa_handle:
        first_record = next(SeqIO.parse(msa_handle, "fasta"))

    # Temporarily store the non-redundant sequences in a list
    non_redundant_sequences = []
    with open(output_file, "r") as output_handle:
        for record in SeqIO.parse(output_handle, "fasta"):
            non_redundant_sequences.append(record)

    # Write the non-redundant sequences to the output file
    with open(output_file, "w") as final_handle:
        SeqIO.write([first_record] + non_redundant_sequences, final_handle, "fasta")

    # Count the number of sequences in the MSA and the non-redundant MSA
    total_sequences = sum(1 for record in SeqIO.parse(msa_file, "fasta"))
    non_redundant_sequences_count = len(non_redundant_sequences)

    # Calculate the effective sequences (Nf)
    Nf = non_redundant_sequences_count / total_sequences
    print("The number of non-redundant sequences:", non_redundant_sequences_count)
    print("The total number of sequences:", total_sequences)
    print("The ratio of non-redundant sequences (Nf):", "{:.2f}".format(Nf))

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
#     print(f"The number of unique sequences: {unique_sequences}")

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
                        result_data.append([query_id, hmm_from, hmm_to, hmm_length])
                    except ValueError:
                        pass
                    extract_info = False

        # Create and return the DataFrame
        result_df = pd.DataFrame(result_data, columns=['id', 'hmm_from', 'hmm_to', 'hmm_length'])
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
#     ax.hlines(hmm_results['Sequence'], hmm_results['env_from'], hmm_results['env_to'], linewidth=10, color='pink', label='Envelope region', linestyle='-')

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
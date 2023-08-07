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

# directory = "/Users/alina/HMM/"
directory = os.getcwd()

def json_parser(line):
    """
    Parses .json files and returns a pandas dataframe
    """
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

def blast_parser(input_file):
    """
    Parses BLAST alignment in XML format and returns a pandas dataframe
    """
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


# def blast_iterator(out_file, query_sequence, selected):
#     with open(out_file, "w") as fout:
#         max_row = len(query_sequence)
#         mapped_seq = ["-"] * max_row
#         for index, row in selected.iterrows():
#             c = 0
#             query_start = row["query_start"]
#             for l_q, l_s in zip(row['query_seq'], row['subject_seq']):
#                 if l_q != " " and l_q != '-': # if the initial aa from query is not empty or gapped
#                     mapped_seq[query_start + c - 1] = l_s if l_s != " " else "-" # assign aa to subject
#                     c += 1
#             fout.write(">{}\n{}\n".format(row["subject_id"], "".join(mapped_seq)))
#             return fout

# def build_msa_from_blast(path, selected_id, query_sequence, selected_init):
#     out_file = f'{directory}results/alignments/output_files/blast/{selected_id}_blast.fasta'
#
#     with open(out_file, "w") as fout:
#         mapped_seq = ["-"] * len(query_sequence)
#
#         # Write the header line for the query sequence
#         fout.write(">{}\n".format(selected_id))
#
#         # Map the query sequence to the mapped_seq list
#         c = 0
#         for l_q in query_sequence:
#             if l_q != " " and l_q != '-':
#                 mapped_seq[c] = l_q
#                 c += 1
#
#         # Write the query_mapped_seq sequence to the output file
#         fout.write("{}\n".format("".join(mapped_seq)))
#
#         # Map the subject sequences to the mapped_seq list and write to the output file
#         for index, row in selected_init.iterrows():
#             c = 0
#             query_start = row["query_start"]
#             for l_q, l_s in zip(row['query_seq'], row['subject_seq']):
#                 if l_q != " " and l_q != '-':
#                     mapped_seq[query_start + c - 1] = l_s if l_s != " " else "-"
#                     c += 1
#             fout.write(">{}\n{}\n".format(row["subject_id"], "".join(mapped_seq)))

def stats_calculation(seqs):
    """
    Returns the values of occupancy and entropy for each alignment
    :param seqs: alignment in a dataframe
    """
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

def get_fasta(upac): # a function retrieving a FASTA sequence from ID
    """
    Returns a sequence for a particular query ID from the remote cluster
    :param upac: query ID
    """
    result = subprocess.run(['ssh', 'echidna', '/software/packages/ncbi-blast/latest/bin/blastdbcmd', '-db', '/db/blastdb/uniprot/uniprot.fasta', '-entry', upac], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return result.stdout

def trim_dis_regions(file, query_id, start_positions, end_positions):
    """
    Returns a FASTA file containing disordered regions only.
    :param file: an initial FASTA alignment
    :param query_id: query ID
    :param start_positions: the number of column/columns representing the start of the disordered region
    :param end_positions: the number of column/columns representing the end of the disordered region
    """
    records = []
    for record in SeqIO.parse(file, "fasta"):
        sequence = record.seq
        trimmed_sequence = ""
        previous_end = 0  # Track the end position of the previous region
        for start, end in zip(start_positions, end_positions):
            if len(sequence) >= start > previous_end and end <= len(sequence):
                trimmed_sequence += '-' * (start - previous_end - 1) + sequence[start - 1: end]
                previous_end = end  # Update the previous end position
            else:
                print("Invalid region: start =", start, "end =", end)
        trimmed_sequence += '-' * (len(sequence) - previous_end)  # Add trailing gaps
        record.seq = trimmed_sequence
        records.append(record)
    SeqIO.write(records, f"{directory}/results/alignments/output_files/disordered/{query_id}_disordered.fasta", "fasta")

    trimmed_sequences = [record.seq for record in records]  # Collect trimmed sequences
    return trimmed_sequences

def select_dis_regions(msa_file, query_id, start_positions, end_positions, output_directory):
    """
    Saves disordered regions from MSA output in separate FASTA files
    """

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
    print(f"The number of unique sequences: {unique_sequences}")

    return stats

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

                if pfam_id.startswith("PF"):
                    data.append([uniprot_id, ipr_id, description, pfam_id, start_pos, end_pos])

    columns = ["uniprot_id", "ipr_id", "description", "pfam_id", "start_pos", "end_pos"]
    pfam = pd.DataFrame(data, columns=columns)
    pfam["start_pos"] = pd.to_numeric(pfam["start_pos"])
    pfam["end_pos"] = pd.to_numeric(pfam["end_pos"])
    pfam['length'] = pfam['end_pos'] - pfam['start_pos'] + 1

    pfam_filtered = pfam[~pfam.apply(lambda x: any(
        (x['start_pos'] > curated_region['end']) or (x['end_pos'] < curated_region['start'])
        for _, curated_region in filtered_curated.iterrows()), axis=1)]

    return pfam_filtered

def calculate_Nf(msa_file, threshold, id_dis):
    # Read the first line from the original MSA file
    with open(msa_file, "r") as msa_handle:
        first_record = next(SeqIO.parse(msa_handle, "fasta"))

    output_file = f"/Users/alina/HMM/results/alignments/input_files/non-redundant/Nf_{id_dis}.fasta"
    cd_hit_path = "/Users/alina/cd-hit/cd-hit"

    # Run CD-HIT to cluster the sequences (excluding the first line) and remove redundancy
    cmd = f"{cd_hit_path} -i {msa_file} -o {output_file} -c {threshold}"
    subprocess.call(cmd, shell=True)

    # Temporarily store the non-redundant sequences in a list
    non_redundant_sequences = []
    with open(output_file, "r") as output_handle:
        for record in SeqIO.parse(output_handle, "fasta"):
            non_redundant_sequences.append(record)

    # Write the non-redundant sequences (including the first line) to the output file
    with open(output_file, "w") as final_handle:
        SeqIO.write([first_record] + non_redundant_sequences, final_handle, "fasta")

    # Count the number of sequences in the MSA and the non-redundant MSA
    total_sequences = sum(1 for record in SeqIO.parse(msa_file, "fasta"))
    non_redundant_sequences_count = len(non_redundant_sequences)

    # Calculate the effective sequences (Nf)
    Nf = non_redundant_sequences_count / total_sequences
    print("The number of non-redundant sequences:", non_redundant_sequences_count)
    print("The total number of sequences:", total_sequences)
    print("Number of effective sequences (Nf):", "{:.2f}".format(Nf))

    return
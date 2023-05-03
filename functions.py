from collections import Counter
import scipy.stats
import numpy as np
import pandas as pd
import subprocess
from Bio.Blast import NCBIXML

def parser(input_file):
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
                    data.append((query_id.split(" ")[0],
                                    subject_id.split(" ")[0],
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


def blast_iterator(selected, out_file):
    with open(out_file, "w") as fout:
        max_row = len(query_sequence)
        mapped_seq = ["-"] * max_row
        for index, row in selected.iterrows():
            c = 0
            query_start = row["query_start"]
            for l_q, l_s in zip(row['query_seq'], row['subject_seq']):
                if l_q != " " and l_q != '-': # if the initial aa from query is not empty or gapped
                    mapped_seq[query_start + c - 1] = l_s if l_s != " " else "-" # assign aa to subject
                    c += 1
            fout.write(">{}\n{}\n".format(row["subject_id"], "".join(mapped_seq)))


def calculate(seqs):
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
    result = subprocess.run(['ssh', 'echidna', '/software/packages/ncbi-blast/latest/bin/blastdbcmd', '-db', '/db/blastdb/uniprot/uniprot.fasta', '-entry', upac], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return result.stdout


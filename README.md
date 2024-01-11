# Hidden Markov Models (HMM) for disordered protein domains

Author: Skrylnik Alina.

This repository contains the following files: 

- `msa_analysis.ipynb`: analysis of the multiple sequence alignments.
- `hmm_analysis.ipynb`: analysis of the generated Hidden Markov Models.
- `functions.py`: contains the functions used in the previously listed notebooks.

## Data acquisition and processing

We open the folder `proj`. All scripts are run from the directory `proj/scripts`.

- `proj/`
  - `databases/`
    - `curated_disprot.csv`: A data frame containing information about Uniprot IDs, start and end positions of the disordered regions.
    - `filtered.tsv.gz`: An extract from protein2ipr database, containing only Pfam domains.
    - `uniprot_filtered/`: Curated Uniprot instances containing the disordered regions data frames.
      - `curated_uniprot_disprot_0.csv`
      - `curated_uniprot_disprot_1.csv`
      - ...
      - `parts/`: Data frames of Uniprot instances split for each query ID.
         - `curated_uniprot_disprot_(Uniprot_ID_1).csv`
         - `curated_uniprot_disprot_(Uniprot_ID_2).csv`
         - ...
  - `msa/`
    - `unali_seqs/`: Unaligned sequences, serving as input for MSA build.
      - `(Uniprot_ID_1)_input.fasta`
      - `(Uniprot_ID_2)_input.fasta`
      - ...
    - `clustal_msa/`: ClustalOmega MSA results.
      - `(Uniprot_ID_1)_clustal.fasta`
      - `(Uniprot_ID_2)_clustal.fasta`
      - ...
      - `disordered/`: Trimmed MSA, containing only disordered regions.
         - `(Uniprot_ID_1)_(region).fasta`
         - `(Uniprot_ID_2)_(region).fasta`
         - ...
  - `hmm/`
    - `hmmbuild/`: Hidden Markov Models.
      - `(Uniprot_ID_1)_(region).hmm`
      - `(Uniprot_ID_2)_(region).hmm`
      - ...
    - `hmmsearch/`: HMM search results against RP75%, in text format.
      - `hmmsearch_clustal_(Uniprot_ID_1)_(region).txt`
      - `hmmsearch_clustal_(Uniprot_ID_2)_(region).txt`
      - ...
         - `preprocessed/`: Processed hmmsearch output in CSV format.
            - `hmmsearch_clustal_df_(Uniprot_ID_1)_(region).csv`
            - `hmmsearch_clustal_df_(Uniprot_ID_2)_(region).csv`
            - ...
         - `combined/`: Merged files from the `preprocessed` folder.
            - `combined_hmmsearch_results.csv`
    - `pfam/`
      - `protein2ipr_clustal.tsv`: hmmsearch instances found in InterPro.
      - `pfam_disprot.csv`: Merged hmmsearch and InterPro results.
  - `scripts/`
    - `blast_parser.py`: Converts XML files to data frames.
    - `df_filter.sh`: Selects rows from `curated_disprot.csv`.
    - `split_csv.py`: Splits CSV files for each MSA.
    - `start_array_jobs_msa.sh`, `array_job_msa.sh`: Retrieves unaligned sequences for MSA construction.
    - `clustalo_generator.py`: Builds ClustalOmega MSA.
    - `cut_msa.py`: Selects disordered regions from MSA and saves them in .fasta format.
    - `hmmbuild.py`: Builds HMM from ClustalOmega MSA.
    - `start_array_jobs_hmm.sh`, `array_job_hmm.sh`: Runs hmmsearch for multiple HMM files simultaneously.
    - `hmmsearch_prep.py`: Converts hmmsearch results to CSV format.
    - `protein2ipr.py`: Iterates over a subset of the Pfam database (only with domains).
    - `pfam_hmm_merge.py`: Combines hmmsearch and Pfam files.

We run the scripts one by one in the directory `root/scripts`.

1 `python blast_parser.py`
2 `bash df_filter.sh ../databases/curated_disprot.csv ../databases/uniprot_filtered`
3 `python split_csv.py ../databases/uniprot_filtered/ ../databases/uniprot_filtered/parts/`
4 `sbatch start_array_jobs_msa.sh ../databases/uniprot_filtered/parts ../msa/unali_seq/` - TO FIX
5 `python clustalo_generator.py ../msa/unali_seq ../msa/clustal_msa`
6 `python cut_msa.py`
7 `python hmmbuild.py`
8 `sbatch start_array_jobs_hmm.sh clustal`
9 `python hmmsearch_prep.py`
10 `./protein2ipr.py ../hmm/hmmsearch/combined/ ../databases/filtered.tsv.gz ../hmm/pfam/protein2ipr_clustal.tsv` - search only in domains
11 `python pfam_hmm_merge.py`


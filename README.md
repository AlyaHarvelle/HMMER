# Hidden Markov Models (HMM) for disordered protein domains

Author: Skrylnik Alina.

This repository contains the following files: 

- `msa_preprocessing.ipynb`, `hmm_preprocessing.ipynb`: extract the files to work with in `msa_analysis` and `hmm_analysis` notebooks.
- `msa_analysis.ipynb`: analysis of the multiple sequence alignments.
- `hmm_analysis.ipynb`: analysis of the generated Hidden Markov Models.
- `functions.py`: contains the functions used in the previously listed notebooks.

To prepare the files for the following analysis, one needs to start with the `msa_preprocessing.ipynb`, then `hmm_preprocessing.ipynb`. 
The corresponding analysis is conducted in `msa_analysis.ipynb` and `hmm_analysis.ipynb`.

## Data acquisition and processing

# Project Structure

We open the folder **proj**. All scripts are run there.

- `proj/`
  - `databases/`
    - `curated_disprot.csv`: Uniprot IDs, start and end positions of the disordered regions.
    - `filtered.tsv.gz`: Extract from protein2ipr db, containing only domains.
    - `uniprot_filtered/`: Curated Uniprot instances containing the disordered regions.
      - `curated_uniprot_disprot_0.csv`: Data frame for Uniprot instances with disordered regions.
      - `curated_uniprot_disprot_1.csv`: ...
      - ...
  - `msa/`
    - `unali_seqs/`: Unaligned sequences, serving as input for MSA build.
      - `*_input.fasta`: Unaligned sequences file for a specific Uniprot ID.
      - ...
    - `clustal_msa/`: ClustalOmega MSA results.
      - `*_clustal.fasta`: ClustalOmega MSA file for a specific Uniprot ID.
      - ...
    - `disordered/`: Trimmed MSA, containing only disordered regions.
      - `*_*.fasta`: Trimmed MSA file for a specific Uniprot ID and region.
      - ...
  - `hmm/`
    - `hmmbuild/`: Hidden Markov Models.
      - `*.hmm`: HMM file for a specific Uniprot ID.
      - ...
    - `hmmsearch/`: HMM search results against RP75%, in text format.
      - `hmmsearch_clustal_*_**.txt`: HMM search result file for a specific Uniprot ID and region.
      - ...
    - `preprocessed/`: Processed hmmsearch output in CSV format.
      - `hmmsearch_clustal_df_*_**.csv`: CSV file for a specific Uniprot ID and region.
      - ...
    - `combined/`: Merged files from the `preprocessed` folder.
      - `combined_hmmsearch_results.csv`: Merged CSV file.
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


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

The data is presented with the following structure:

proj/
│
├── databases/
│   ├── curated_disprot.csv
│   ├── filtered.tsv.gz
│   └── uniprot_filtered/
│       ├── curated_uniprot_disprot_0.csv
│       ├── curated_uniprot_disprot_1.csv
│       ├── ...
│       └── parts/
│           ├── curated_uniprot_disprot_(Uniprot_ID_1).csv
│           ├── curated_uniprot_disprot_(Uniprot_ID_2).csv
│           └── ...
│       
│
├── msa/
│   ├── unali_seqs/
│   │   ├── *_input.fasta
│   │   └── ...
│   ├── 
│   │   ├── *_clustal.fasta
│   │   └── ...
│   └── clustal_msa/
│       └── *_*.fasta
│
├── hmm/
│   ├── hmmbuild/
│   │   ├── *.hmm
│   │   └── ...
│   ├── hmmsearch/
│   │   ├── hmmsearch_clustal_*_**.txt
│   │   └── ...
│   ├── preprocessed/
│   │   ├── hmmsearch_clustal_df_*_**.csv
│   │   └── ...
│   ├── combined/
│   │   └── combined_hmmsearch_results.csv
│   └── pfam/
│       ├── protein2ipr_clustal.tsv
│       └── pfam_disprot.csv
│
└── scripts/
    ├── blast_parser.py
    ├── df_filter.sh
    ├── split_csv.py
    ├── start_array_jobs_msa.sh
    ├── array_job_msa.sh
    ├── clustalo_generator.py
    ├── cut_msa.py
    ├── hmmbuild.py
    ├── start_array_jobs_hmm.sh
    ├── array_job_hmm.sh
    ├── hmmsearch_prep.py
    ├── protein2ipr.py
    └── pfam_hmm_merge.py

We run the scripts one by one in the directory `root/scripts`.


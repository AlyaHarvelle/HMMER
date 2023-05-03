import subprocess

def get_fasta(upac): # a function retrieving a FASTA sequence from ID
    result = subprocess.run(['ssh', 'echidna', '/software/packages/ncbi-blast/latest/bin/blastdbcmd', '-db', '/db/blastdb/uniprot/uniprot.fasta', '-entry', upac], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return result.stdout
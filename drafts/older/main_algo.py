import os
import time
import argparse
from ahocorasick import Automaton
from quickdna import DnaSequence
from Bio import SeqIO

DB_FILE = '/app/Sarit_CDR3_database.csv'

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process a FASTQ file to find CDR3 sequences.")
    parser.add_argument("fastq_file", help="Path to the FASTQ file")
    return parser.parse_args()

def fastq_to_aa(input_fastq_file):
    """Generates amino acid sequences from FASTQ."""
    with open(input_fastq_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fastq"):
            seq = DnaSequence(str(record.seq))
            frames = seq.translate_all_frames()
            for frame in frames:
                yield str(frame)

def read_cdr3_1M_into_sorted_list():
    with open(DB_FILE, 'r') as source_file:
        cdr3_1M_list = []
        source_file.readline()  # delete the first line of the file
        for line in source_file:
            line = line.strip()
            line = line.split(',')[2]  # remove the second column from the line
            cdr3_1M_list.append(line)
    return sorted(set(cdr3_1M_list))

def find_multiple_aho(strings, substrings):
    """Search for multiple substrings in a list of strings using Aho-Corasick."""
    A = Automaton()
    for idx, sub in enumerate(substrings):
        A.add_word(sub, (idx, sub))
    A.make_automaton()

    found_substrings = []
    for string in strings:
        found_substrings.extend([sub for pos, (_, sub) in A.iter(string) if len(sub) > 7])
    
    return found_substrings

if __name__ == "__main__":
    args = parse_arguments()
    FASTQ_FILE = args.fastq_file
    FASTQ_NAME = os.path.basename(FASTQ_FILE)
    FASTQ_NAME_WITHOUT_EXT = os.path.splitext(FASTQ_NAME)[0]
    RESULT_DIR = f'/app/algo_results'
    os.makedirs(RESULT_DIR, exist_ok=True)
    OUTPUT_CDR3_FILE = f'{RESULT_DIR}/CDR3_results_{FASTQ_NAME_WITHOUT_EXT}.txt'

    if not os.path.exists(RESULT_DIR):
        os.makedirs(RESULT_DIR)

    print(f'## Start algo on file: {FASTQ_FILE} ##')
    start = time.time()

    print(f'Converting FASTQ to amino acids...')
    aa_sequences = fastq_to_aa(FASTQ_FILE)

    # Load CDR3 sequences from the database
    cdr3MegaList = read_cdr3_1M_into_sorted_list()

    print(f'Start search CDR3 from DB in amino acid results')
    results = find_multiple_aho(aa_sequences, cdr3MegaList)

    # Write the results to output file
    with open(OUTPUT_CDR3_FILE, 'w') as output_file:
        for result in set(results):
            output_file.write(f"{result}\n")

    elapsed_time = time.time() - start
    minutes = elapsed_time // 60
    print(f'Time taken for algo: {minutes} minutes')

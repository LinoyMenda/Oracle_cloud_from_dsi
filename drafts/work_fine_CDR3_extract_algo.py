import os
import subprocess
import logging
import time
import argparse
from ahocorasick import Automaton
from quickdna import DnaSequence
from Bio import SeqIO

###
DB_FILE = '/mnt/cdr_fs_us/Sarit_CDR3_database.csv'
WORKDIR_PATH = '/app'
RESULT_DIR = '/app/algo_results'
SRA_TOOLKIT_PATH =f'{WORKDIR_PATH}/sratoolkit.3.1.1-centos_linux64/bin'

### Create a logger
logger = logging.getLogger('my_logger')
logger.setLevel(logging.DEBUG)  # Set the log level to DEBUG
file_handler = logging.FileHandler(f'{WORKDIR_PATH}/log_file.txt')
file_handler.setLevel(logging.DEBUG)  # Set the handler's log level to DEBUG
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
# Add the handler to the logger
logger.addHandler(file_handler)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download SRA file as a FASTQ file and then find CDR3 sequences.")
    parser.add_argument("sra_id", help="SRA acc id number")
    return parser.parse_args()

def download_sra_file(acc_id):
    try:
        logger.debug(f"Start downloaded: {acc_id} file")
        fastq_download_command = f'{SRA_TOOLKIT_PATH}/fasterq-dump'
        # Run the prefetch command
        download_sra_file = subprocess.run([fastq_download_command, acc_id], check=True)
        logger.debug(f"{download_sra_file}")
        logger.debug(f" ~ Successfully downloaded: {acc_id} file ~")
        fastq_file_path = f'{WORKDIR_PATH}/{acc_id}.fastq'
        return fastq_file_path
    except subprocess.CalledProcessError as e:
        logger.debug(f"Error downloading {acc_id} file: {e}")
        raise  # This will stop the program by re-raising the exception

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
    sra_id = args.sra_id
    FASTQ_FILE = download_sra_file(sra_id)

    os.makedirs(RESULT_DIR, exist_ok=True)
    output_cdr3_file = f'{RESULT_DIR}/CDR3_results_{sra_id}.txt'

    logger.debug(f'## Start algo on file: {FASTQ_FILE} ##')
    start = time.time()

    logger.debug(f'Converting FASTQ to amino acids...')
    aa_sequences = fastq_to_aa(FASTQ_FILE)

    # Load CDR3 sequences from the database
    cdr3MegaList = read_cdr3_1M_into_sorted_list()

    logger.debug(f'Start search CDR3 from DB in amino acid results')
    results = find_multiple_aho(aa_sequences, cdr3MegaList)

    # Write the results to output file
    with open(output_cdr3_file, 'w') as output_file:
        for result in set(results):
            output_file.write(f"{result}\n")
    
    logger.debug(f'Remove FASTQ file')
    os.remove(FASTQ_FILE)

    elapsed_time = time.time() - start
    minutes = elapsed_time // 60
    logger.debug(f'~~~ Time taken for algo: {minutes} minutes ~~~')


import subprocess
import logging
import os
### This Pipleine divide into 2 section: 
### First one : download sra file from dataset.
### For the above section need to used sra toolkit, download from here:
### wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-centos_linux64.tar.gz
### tar -xzvf sratoolkit.3.1.1-centos_linux64.tar.gz
### cd sratoolkit.3.1.1-centos_linux64/bin
### ./vdb-config -i
### ./prefetch /dsi/sbm/linoym/oracle_cloud/top10_SRR.txt
### Second: run "Tmer" algo in order to extract cdr3 sequnces
### remove file and repeat the process

WORKDIR_PATH = '/app'
SRA_TOOLKIT_PATH =f'{WORKDIR_PATH}/sratoolkit.3.1.1-centos_linux64/bin'
SRA_ACC_FILE_PATH = F'{WORKDIR_PATH}/top10_SRR.txt'
FASTQֹ_DOWNLOAD_COMMAND = f'{SRA_TOOLKIT_PATH}/fasterq-dump'

if __name__ == "__main__":

    # Create a logger
    logger = logging.getLogger('my_logger')
    logger.setLevel(logging.DEBUG)  # Set the log level to DEBUG
    file_handler = logging.FileHandler(f'{WORKDIR_PATH}/log_file.txt')
    file_handler.setLevel(logging.DEBUG)  # Set the handler's log level to DEBUG
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    # Add the handler to the logger
    logger.addHandler(file_handler)

    # Read the accession numbers from 'ars_file.txt'
    with open(SRA_ACC_FILE_PATH, 'r') as file:
        accession_numbers = file.read().splitlines()

    # Loop over the accession numbers and download each one using fasterq-dump
    for acc in accession_numbers:
        try:
            logger.debug(f"Start downloaded: {acc} file")
            # Run the prefetch command
            download_sra_file = subprocess.run([FASTQֹ_DOWNLOAD_COMMAND, acc], check=True)
            logger.debug(f"{download_sra_file}")
            logger.debug(f" ~ Successfully downloaded: {acc} file ~")
        except subprocess.CalledProcessError as e:
            logger.debug(f"Error downloading {acc} file: {e}")
        file_path = f'{WORKDIR_PATH}/{acc}.fastq'
        try:
            algo_command = ['python', f'{WORKDIR_PATH}/main_algo.py', file_path]
            logger.debug(f"Start algo on {acc} file")
            result = subprocess.run(algo_command, check=True, text=True, capture_output=True)
            # Print the output and errors (if any)
            logger.debug(result.stdout)
        except subprocess.CalledProcessError as e:
            logger.debug(f"An error occurred while executing the algo on fle{acc}: {e}")
        os.remove(file_path)
        logger.debug(f'Delete {acc} fastq file')
    logger.debug("~End script~")


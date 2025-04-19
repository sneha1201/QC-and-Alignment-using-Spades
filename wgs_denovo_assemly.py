import os
import subprocess
import argparse
import pandas as pd
from multiprocessing import Pool

def run_command(command, log_file):
    """Run a shell command and log output."""
    print(f"Running: {command}")
    with open(log_file, "a") as log:
        log.write(f"\n[COMMAND] {command}\n")
        subprocess.run(command, shell=True, stdout=log, stderr=subprocess.STDOUT)

def run_parallel_command(cmd_log_pair):
    """Wrapper for multiprocessing."""
    command, log_file = cmd_log_pair
    run_command(command, log_file)

def create_dir(directory):
    """Create output directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def main():
    parser = argparse.ArgumentParser(description="WGS de novo Assembly Pipeline")
    parser.add_argument("--sample_info", required=True, help="TSV file with columns: Sample, File Name, Pair")
    parser.add_argument("--fastq_folder", required=True, help="Folder with input FASTQ files")
    parser.add_argument("--output_folder", required=True, help="Folder to store output results")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads per tool")
    parser.add_argument("--jobs", type=int, default=2, help="Number of samples to run in parallel")
    args = parser.parse_args()

    sample_info = pd.read_csv(args.sample_info, sep="\t")
    fastq_folder = os.path.abspath(args.fastq_folder)
    output_folder = os.path.abspath(args.output_folder)
    threads = args.threads
    parallel_jobs = args.jobs

    log_file = os.path.join(output_folder, "pipeline.log")
    with open(log_file, "w") as log:
        log.write("WGS De Novo Assembly Pipeline Log\n")

    # Define output directories
    qc_folder = os.path.join(output_folder, "1_QC")
    qc_before = os.path.join(qc_folder, "qc_before")
    qc_after = os.path.join(qc_folder, "qc_after")
    trimmed_folder = os.path.join(output_folder, "2_Trimmed")
    spades_folder = os.path.join(output_folder, "3_Spades_Assembly")

    for folder in [qc_before, qc_after, trimmed_folder, spades_folder]:
        create_dir(folder)

    # Pre-trim QC
    run_command(f"fastqc {fastq_folder}/*.fastq.gz -o {qc_before} -t {threads}", log_file)
    run_command(f"seqkit stats -a {fastq_folder}/*.fastq.gz -o {qc_before}/stats.tsv", log_file)
    run_command(f"multiqc {qc_before} -o {qc_before}", log_file)

    # Trimming
    trim_commands = []
    for sample in sample_info["Sample"].unique():
        files = sample_info[sample_info["Sample"] == sample]
        r1 = os.path.join(fastq_folder, files[files["Pair"] == "Pair1"]["File Name"].values[0])
        r2 = os.path.join(fastq_folder, files[files["Pair"] == "Pair2"]["File Name"].values[0])
        out1 = os.path.join(trimmed_folder, f"{sample}_1_trimmed.fastq.gz")
        out2 = os.path.join(trimmed_folder, f"{sample}_2_trimmed.fastq.gz")
        html = os.path.join(trimmed_folder, f"{sample}_fastp.html")
        json = os.path.join(trimmed_folder, f"{sample}_fastp.json")

        cmd = (
            f"fastp -i {r1} -I {r2} -o {out1} -O {out2} -x "
            f"--detect_adapter_for_pe -w {threads} -h {html} -j {json}"
        )
        trim_commands.append((cmd, log_file))

    with Pool(parallel_jobs) as pool:
        pool.map(run_parallel_command, trim_commands)

    # Post-trim QC
    run_command(f"fastqc {trimmed_folder}/*.fastq.gz -o {qc_after} -t {threads}", log_file)
    run_command(f"seqkit stats -a {trimmed_folder}/*.fastq.gz -o {qc_after}/stats.tsv", log_file)
    run_command(f"multiqc {qc_after} -o {qc_after}", log_file)

    # SPAdes Assembly
    spades_commands = []
    for sample in sample_info["Sample"].unique():
        left = os.path.join(trimmed_folder, f"{sample}_1_trimmed.fastq.gz")
        right = os.path.join(trimmed_folder, f"{sample}_2_trimmed.fastq.gz")
        sample_spades_out = os.path.join(spades_folder, sample)  # Sample-specific subdirectory
        create_dir(sample_spades_out)

        cmd = (
            f"spades.py -1 {left} -2 {right} -o {sample_spades_out} -t {threads} --careful"
        )
        spades_commands.append((cmd, log_file))

    with Pool(parallel_jobs) as pool:
        pool.map(run_parallel_command, spades_commands)

    # Assembly Stats (QUAST in same assembly folders)
    quast_commands = []
    for sample in sample_info["Sample"].unique():
        contigs = os.path.join(spades_folder, sample, "scaffolds.fasta")  # Reference the sample-specific contigs file
        quast_out = os.path.join(spades_folder, sample, "assembly_stats")
        create_dir(quast_out)

        cmd = f"quast {contigs} -o {quast_out} -t {threads} "
        quast_commands.append((cmd, log_file))

    with Pool(parallel_jobs) as pool:
        pool.map(run_parallel_command, quast_commands)

    print("✅ Pipeline finished successfully!")

if __name__ == "__main__":
    main()


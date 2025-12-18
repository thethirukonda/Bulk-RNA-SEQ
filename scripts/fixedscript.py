import subprocess
import time
import os

sra_numbers = [
    "SRR7179504", "SRR7179505", "SRR7179506", "SRR7179507",
    "SRR7179508", "SRR7179509", "SRR7179510", "SRR7179511",
    "SRR7179520", "SRR7179521", "SRR7179522", "SRR7179523",
    "SRR7179524", "SRR7179525", "SRR7179526", "SRR7179527",
    "SRR7179536", "SRR7179537", "SRR7179540", "SRR7179541"
]

os.makedirs("fastq", exist_ok=True)

# Convert all to FASTQ
for sra_id in sra_numbers:
    # CORRECT PATH for your setup
    sra_path = f"/Users/aakash/Downloads/rnaseq/{sra_id}/{sra_id}.sra"
    
    # Check if file exists first
    if not os.path.exists(sra_path):
        print(f"  Skipping {sra_id} - file not found at {sra_path}")
        continue
    
    print(f"\n=== Generating FASTQ for: {sra_id} ===")
    fastq_dump_cmd = (
        f"fastq-dump --outdir fastq --gzip --skip-technical "
        f"--readids --read-filter pass --dumpbase --split-3 --clip {sra_path}"
    )
    print("Command:", fastq_dump_cmd)
    start_time = time.time()
    subprocess.call(fastq_dump_cmd, shell=True)
    end_time = time.time()
    elapsed_min = (end_time - start_time) / 60
    print(f" FASTQ generation time for {sra_id}: {elapsed_min:.2f} minutes")

print("\n All done! Check the 'fastq' folder for your files.")


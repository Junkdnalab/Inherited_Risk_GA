# library
import subprocess
# read in sra id
with open("/drive-pool/data/peter_data/sc_data/brca/GSE168837/gse168837_sra.txt") as sraFile:
    sra = sraFile.read().splitlines()
# download .sra files
for sra_id in sra:
    print("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    print("The command used was: " + prefetch)
    subprocess.call(prefetch, shell = True)
# extract .sra files from above 
for sra_id in sra:
    print("Generating fastq for: " + sra_id)
    fastq_dump = "fastq-dump --split-files --outdir fastq --gzip " + sra_id
#    fastq_dump = "fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip /drive-pool/data/peter_data/sc_data/brca/GSE168837/" + sra_id + "/" + sra_id + ".sra"
    print("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell = True)

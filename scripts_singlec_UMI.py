# Workflow for the pyrho analyses of the wolves data

from gwf import Workflow
import glob
import os
import os.path
import itertools

# import pandas as pd

gwf = Workflow()


def download_data():
    inputs = []
    outputs = [
        "/u/home/m/mica20/project-collaboratory/running_star/hgmm_100_fastqs.tar"
    ]
    options = {
        "memory": "8g",
        "cores": 1,
        "walltime": "23:59:59",
    }
    # Step 1: get data
    spec = """
	wget http://cf.10xgenomics.com/samples/cell-exp/1.3.0/hgmm_100/hgmm_100_fastqs.tar;
	tar -xf hgmm_100_fastqs.tar;
	cat fastqs/hgmm_100_S1_L00?_R1_001.fastq.gz > hgmm_100_R1.fastq.gz;
	cat fastqs/hgmm_100_S1_L00?_R2_001.fastq.gz > hgmm_100_R2.fastq.gz;
	"""
    return inputs, outputs, options, spec


def barcodes(fastq):
    inputs = []
    outputs = ["/u/home/m/mica20/project-collaboratory/whitelist.txt"]
    options = {
        "memory": "8g",
        "cores": 1,
        "walltime": "23:59:59",
    }
    # Step 2: Identify correct cell barcodes
    spec = """
	umi_tools whitelist --stdin {fastq} \
	                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	                    --set-cell-number=100 \
	                    --log2stderr > whitelist.txt;
	""".format(
        fastq=fastq
    )
    return inputs, outputs, options, spec


def extract_barcodes_umis(fastq1, fastq2, whitelist):
    inputs = []
    outputs = []
    options = {
        "memory": "8g",
        "cores": 1,
        "walltime": "23:59:59",
    }
    # Step 3: Extract barcdoes and UMIs and add to read names
    spec = """
	umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	                  --stdin {fastq1} \
	                  --stdout hgmm_100_R1_extracted.fastq.gz \
	                  --read2-in {fastq2} \
	                  --read2-out=hgmm_100_R2_extracted.fastq.gz \
	                  --whitelist={whitelist};
	""".format(
        fastq1=fastq1, fastq2=fastq2, whitelist=whitelist
    )
    return inputs, outputs, options, spec


directory_analyses = "/u/home/m/mica20/project-collaboratory/running_star"
print("Downloading samples")
gwf.target_from_template(f"Download_data", download_data())

print("Extracting barcode from samples")
gwf.target_from_template(
    f"Barcode_data", barcodes(fastq=f"{directory_analyses}/hgmm_100_R1.fastq.gz")
)

print("Extract barcdoes and UMIs and add to read names")
gwf.target_from_template(
    f"Barcode_data_umis",
    extract_barcodes_umis(
        fastq1=f"{directory_analyses}/hgmm_100_R1.fastq.gz",
        fastq2=f"{directory_analyses}/hgmm_100_R2.fastq.gz",
        whitelist=f"{directory_analyses}/whitelist.txt",
    ),
)

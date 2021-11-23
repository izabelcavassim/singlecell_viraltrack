# Workflow for the pyrho analyses of the wolves data

from gwf import Workflow
import glob
import os
import os.path
import itertools

# import pandas as pd

gwf = Workflow()


def concatenate_data(filename):
    inputs = []
    outputs = [
        f"/u/home/m/mica20/project-collaboratory/running_star/{filename}_R1.fastq.gz",
        f"/u/home/m/mica20/project-collaboratory/running_star/{filename}_R2.fastq.gz"
    ]
    options = {
        "memory": "8g",
        "cores": 1,
        "walltime": "23:59:59",
    }
    # Step 1: get the data
    spec = """
	cat {filename}/*R1_001.fastq.gz > {filename}_R1.fastq.gz;
	cat {filename}/*_R2_001.fastq.gz > {filename}_R2.fastq.gz;
	""".format(filename=filename)
    print(spec)
    return inputs, outputs, options, spec


def barcodes(fastq, filename, working_dir):
    inputs = [f'{working_dir}/{filename}_R1.fastq.gz', f'{working_dir}/{filename}_R1.fastq.gz']
    outputs = [f'{working_dir}/whitelist_{filename}.txt']
    options = {
        "memory": "8g",
        "cores": 1,
        "walltime": "23:59:59",
    }
    # Step 2: Identify correct cell barcodes
    spec = """
	umi_tools whitelist --stdin {fastq} \
	--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--set-cell-number=5500 \
    --plot-prefix \
	--log2stderr > whitelist_{filename}.txt;
	""".format(
        fastq=fastq,
        filename=filename
    )
    print(spec)
    return inputs, outputs, options, spec


def extract_barcodes_umis(filename, fastq1, fastq2, whitelist, working_dir):
    inputs = [f'{working_dir}/{filename}_R1.fastq.gz', f'{working_dir}/{filename}_R1.fastq.gz', f'{working_dir}/whitelist_{filename}.txt']
    outputs = [f'{working_dir}/{filename}_5500_R1_extracted.fastq.gz', f'{working_dir}/{filename}_5500_R2_extracted.fastq.gz']
    options = {
        "memory": "8g",
        "cores": 1,
        "walltime": "23:59:59",
    }
    # Step 3: Extract barcdoes and UMIs and add to read names
    spec = """
	umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
	--stdin {fastq1} \
	--stdout {filename}_5500_R1_extracted.fastq.gz \
	--read2-in {fastq2} \
	--read2-out={filename}_5500_R2_extracted.fastq.gz \
	--whitelist={whitelist};
	""".format(
        fastq1=fastq1, fastq2=fastq2, whitelist=whitelist, filename=filename
    )
    print(spec)
    return inputs, outputs, options, spec


# Two folders with single cell data
sc_folder = ["595_BILs_hTCR", "595_BILs_5GEX"]


# Submitting jobs
directory_analyses = "/u/home/m/mica20/project-collaboratory/running_star"
for fol in sc_folder:
    print("Concatenating samples")
    gwf.target_from_template(f"Concatenate_data_{fol}", concatenate_data(filename=fol))

    print("Extracting barcode from samples")
    gwf.target_from_template(
        f"Barcode_data_{fol}", 
        barcodes(
            fastq=f"{directory_analyses}/{fol}_R1.fastq.gz",
            filename=fol,
            working_dir=directory_analyses
        )
    )

    print("Extract barcdoes and UMIs and add to read names")
    gwf.target_from_template(
        f"Barcode_data_umis_{fol}",
        extract_barcodes_umis(
            filename = fol,
            fastq1=f"{directory_analyses}/{fol}_R1.fastq.gz",
            fastq2=f"{directory_analyses}/{fol}_R2.fastq.gz",
            whitelist=f"{directory_analyses}/whitelist_{fol}.txt",
            working_dir=f"{directory_analyses}"
        ),
    )

import os
import pandas as pd
import subprocess

# Input GWAS paths table
INPUT_FILE = "all_traits_gwaspath.txt"

# METAL templates
TEMPLATE = """SCHEME STDERR
WEIGHT N
AVERAGEFREQ ON
MINMAXFREQ ON
VERBOSE OFF
GENOMICCONTROL ON
CUSTOMVARIABLE TotalSampleSize

{marker_sections}

OUTFILE STDERR_{trait} .tbl
ANALYZE HETEROGENEITY

QUIT
"""

MARKER_TEMPLATE = """MARKER SNP
ALLELE A1 A2
FREQ Freq
PVAL p
EFFECT b
STDERR se
WEIGHT N
LABEL TotalSampleSize as N
PROCESS {path}"""

# Load GWAS table
df = pd.read_csv(INPUT_FILE, sep="\t")
all_traits = sorted(df["trait"].unique())

# Generate METAL config and run for each trait
for trait in all_traits:
    trait_data = df[df["trait"] == trait]
    
    marker_sections = "\n\n".join([
        MARKER_TEMPLATE.format(path=row["path"])
        for _, row in trait_data.iterrows()
    ])
    
    content = TEMPLATE.format(marker_sections=marker_sections, trait=trait)
    
    os.makedirs(trait, exist_ok=True)
    config_file = os.path.join(trait, f"{trait}_meta.txt")
    with open(config_file, "w") as f:
        f.write(content)
    
    log_file = os.path.join(trait, "metal.log")
    with open(log_file, "w") as log:
        subprocess.run(["metal", config_file], stdout=log, stderr=log)

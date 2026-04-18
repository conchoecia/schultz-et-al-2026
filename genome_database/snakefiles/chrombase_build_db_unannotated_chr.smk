"""
Program  : GenDB_build_db_unannotated.snakefile
Language : snakemake
Date     : 2023-10-01
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.

Description:
  This program:
    - takes the list of unannotated genomes,
    - annotates the genomes with common animal linkage groups,
    - and prepares a database of those linkage groups for ODP.
  The program requires that either a directory of the annotated and unannotated genome lists be provided.
  Otherwise the user has to specify specific paths to those tsv files.

Usage instructions:
  - There are currently no usage instructions. This is a work in progress.

Updates:
  - August 2025 - updates for chrombase
"""

import pandas as pd
from datetime import datetime
import yaml

# This block imports fasta-parser as fasta
import os
import sys
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
dependencies_path = os.path.join(snakefile_path, "dependencies")
sys.path.insert(1, dependencies_path)
import fasta

src_path = os.path.join(snakefile_path, "src")
sys.path.insert(1, src_path)
import GenDB

import math
from pathlib import Path
GiB = 1024**3
def size_gib(p: str) -> float:
    return max(0.001, Path(p).stat().st_size / GiB)  # float GiB, tiny floor

# figure out where bin is because we need to use some outside tools
bin_path = os.path.join(snakefile_path, "bin")

configfile: "config.yaml"
config["tool"] = "odp_ncbi_genome_db"

# Do some logic to see if the user has procided enough information for us to analyse the genomes
if ("directory" not in config) and ("accession_tsvs" not in config):
    raise IOError("You must provide either a directory of the annotated and unannotated genome lists, or a list of the paths to those tsv files. Read the config file.")

config["tempdir"] = "/tmp"
# check that the tempdir exists
if "tempdir" not in config:
    raise IOError("You must provide a temporary directory to store temporary files. Read the config file for instructions.")

# Strip all trailing slashes from the tempdir
#  20250901 NOTE - I'm not sure why this is here, maybe in case the tempdir has space characters?
#                  This seems like rule written for some specific change originally.
config["temp"] = config["tempdir"].rstrip("/").rstrip("\\")
if not os.path.isdir(config["tempdir"]):
    raise IOError("The temporary directory you provided does not exist. {}".format(config["tempdir"]))

config = GenDB.opening_logic_GenDB_build_db(config, chr_scale = True, annotated = False,
                                            tempignore_path = os.path.join(snakefile_path, "data/temporary_ignore_list.txt"))
# Print some info about the files that we found.
printed = False
if not printed:
    GenDB.print_gendb_config_summary(config, chr_scale=True, annotated=True)
    printed = True

# One key feature of this script is that we will map proteins from
#  LG databases to annotate those genomes with the LG identities.
#  Before we accept the user's input, we need to parse the supplied
#  directories to make sure that they are valid.
LG_to_db_directory_dict = {}
LG_to_rbh_dfs           = {}
LG_outfiles             = []
if len(LG_to_db_directory_dict) == 0: # only do this once
    if "LG_db_directories" not in config:
        raise ValueError("You have not specified the directories of LG databases. Please add the key 'LG_db_directories' to your config file.")
    for thisdirectory in config["LG_db_directories"]:
        # The directory must exist
        if not os.path.isdir(thisdirectory):
            raise IOError("The directory of LG databases you provided does not exist. {}".format(thisdirectory))
        # There must be a directory named aligned
        if not os.path.isdir(os.path.join(thisdirectory, "aligned")):
            raise IOError("The directory of LG databases you provided does not contain a directory named 'aligned'. {}".format(thisdirectory))
        # The directory name is the name of the LG, just get the last part of the path
        LG_name = os.path.basename(os.path.normpath(thisdirectory))
        # The directory must contain a .rbh file
        if not any([x.endswith(".rbh") for x in os.listdir(thisdirectory)]):
            raise IOError("The directory of LG databases you provided does not contain a .rbh file. {}".format(thisdirectory))
        # Now we read in the .rbh file as a pandas df to get the LG names from the `.rbh` column.
        rbhfilepath = [x for x in os.listdir(thisdirectory) if x.endswith(".rbh")][0]
        df = pd.read_csv(os.path.join(thisdirectory, rbhfilepath), sep="\t")
        # For each rbh entry in the dataframe, there should be a fasta file in the aligned directory
        fasta_list = [x for x in os.listdir(os.path.join(thisdirectory, "aligned")) if x.endswith(".fasta")]
        for rbh in df["rbh"].unique():
            if rbh + ".fasta" not in fasta_list:
                raise IOError("The directory of LG databases you provided does not contain a fasta file for the rbh entry {}. {}".format(rbh, thisdirectory))
        # This appears to be a legitimate LG database directory. Add it to the dictionary.
        # Also add a pandas dataframe to the dictionary.
        LG_to_db_directory_dict[LG_name] = thisdirectory
        LG_to_rbh_dfs[LG_name] = df
        outfile = config["tool"] + "/input/LG_proteins/{}.fasta".format(LG_name)
        LG_outfiles.append(outfile)

# I am not sure if this is actually used
resource_scopes:
    mem_mb         = "local",
    time           = "local",
    runtime        = "local",
    download_slots = "global"

wildcard_constraints:
    taxid="[0-9]+",

# first we must load in all of the files. Only do it once
rule all:
    input:
        expand("NCBI_odp_db.unannotated.{LG_name}.yaml",
               LG_name=LG_to_db_directory_dict.keys()),
        expand("NCBI_odp_sp_list.unannotated.{LG_name}.txt",
               LG_name=LG_to_db_directory_dict.keys()),
        # the creation of the yaml file implies the other files were created
        #expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz", assemAnn=config["assemAnn"]),
        ##LG_outfiles,
        #expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom.gz",
        #       assemAnn=config["assemAnn"], LG_name=LG_to_db_directory_dict.keys()),
        #expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep.gz",
        #       assemAnn=config["assemAnn"], LG_name=LG_to_db_directory_dict.keys()),

rule dlChrs:
    """
    We have selected the unannotated genomes to download.
    For these genomes we need to find a way to annotate them.

    To specifically download the chromosome-scale scaffolds, there is a series of commands with the NCBI datasets tool.
      I found this set of instructions after opening a ticket on NCBI's github page: https://github.com/ncbi/datasets/issues/298

    Notes:
      - 12-31-2023 - Using the command line had too many edge cases that didn't work, so I resorted to using python to do a more careful job.
        This verifies that the files are downloaded and unzipped correctly, and contain all of the expected sequences.
        Therefore, we do not need additional verification steps for the assembly file.
    """
    output:
       fasta   = temp(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta"),
       allscaf = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.scaffold_df.all.tsv",
       chrscaf = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.scaffold_df.chr.tsv"
    retries: 3
    params:
        datasets = os.path.join(bin_path, "datasets"),
        outdir   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/",
    threads: 1
    group: "dlgz"
    resources:
        mem_mb = GenDB.dlChrs_get_mem_mb, # 1 GB of RAM
        runtime = 20,
        download_slots = 1
    run:
        result = GenDB.download_unzip_genome(wildcards.assemAnn, params.outdir,
                                             params.datasets, chrscale = True)
        if result != 0:
            raise ValueError("The download of the genome {} failed.".format(wildcards.assemAnn))

rule gzip_fasta_file:
    """
    In this step zip the fasta file to conserve space.
    """
    input:
        genome = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta",
    output:
        genome = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz"
    threads: 1
    group: "dlgz"
    resources:
        mem_mb  = 1000, # 1 GB of RAM
        runtime = lambda wildcards: GenDB.gzip_get_time(config["assemAnn_to_scaflen"][wildcards.assemAnn])
    shell:
        """
        echo "Gzipping the fasta file."
        gzip < {input.genome} > {output.genome}
        """

rule generate_LG_fasta_sequence:
    """
    Currently, in the LG database, the sequences are only available as alignments.
    Here, we just concatenate all of the alignments, then strip the gaps.
    """
    input:
        LG_dir = lambda wildcards: LG_to_db_directory_dict[wildcards.LG_name]
    output:
        fasta = config["tool"] + "/input/LG_proteins/{LG_name}.fasta"
    resources:
        mem_mb  = 1000, # 1 GB of RAM
        runtime = 5
    threads: 1
    run:
        with open(output.fasta, "w") as o:
            for fastafile in os.listdir(input.LG_dir + "/aligned"):
                for record in fasta.parse(input.LG_dir + "/aligned/" + fastafile):
                    o.write(">{}\n{}\n".format(record.id, record.seq.replace("-", "")))

rule miniprot:
    """
    This handles all of the miniprot steps. Both the indexing and the mapping.
    Doing it this way prevents keeping a ton of temporary data on the hard drive.
    It is kind of wonky in that it may crash a few times before it works, but it eventually
      works as the amount of ram increases with each run.
    Subsequent jobs are usually scheduled on the same nodes, so the index files are often still there.
    """
    input:
        pep    = config["tool"] + "/input/LG_proteins/{LG_name}.fasta",
        genome = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz",
    output:
        paf  = config["tool"] + "/output/mapped_reads/{assemAnn}/{LG_name}_to_{assemAnn}.paf"
    threads: 8
    retries: 8
    params:
        mpi_suffix = "{LG_name}_{assemAnn}.filt.fasta.gz.mpi" # this is the temporary index file
    resources:
        tmpdir  = config["tempdir"], # the place where the temporary index file will be stored
        mem_mb  = GenDB.miniprot_get_mem_mb, # The RAM usage can blow up during indexing. Often > 10GB. 6Gbp genomes need more than 20GB of RAM.
        runtime  = lambda wildcards, input: int(math.ceil( # makes this run 30 minutes for every 1/2 GiB of genome size
            5 + 30 * max(0.001, Path(input.genome).stat().st_size / GiB)
            ))
    shell:
        """
        INDEXFILE=$TMPDIR/{params.mpi_suffix}
        # don't index if the file already exists
        if [ ! -f ${{INDEXFILE}} ]; then
            # INDEXING step
            echo "Indexing the genome: Genome {wildcards.assemAnn} for LG {wildcards.LG_name}"
            miniprot -t {threads} -d ${{INDEXFILE}} {input.genome}
        fi

        # MAPPING step
        echo ""
        echo "Mapping the proteins to the genome. {wildcards.LG_name} to {wildcards.assemAnn}"
        miniprot -t {threads} ${{INDEXFILE}} {input.pep} > {output.paf}

        # REMOVE the index
        echo ""
        echo "Removing the index file: ${{INDEXFILE}}"
        rm ${{INDEXFILE}}
        """

rule filter_paf_for_longer_scaffold:
    """
    Many times, a single protein will map equally well to multiple scaffolds.
    In this case, we just pick the longest of the scaffolds.

    To do this, we load of a file of all the chrom sizes first to do the comparisons.
    """
    input:
        paf = config["tool"] + "/output/mapped_reads/{assemAnn}/{LG_name}_to_{assemAnn}.paf",
    output:
        paf = config["tool"] + "/output/mapped_reads/{assemAnn}/{LG_name}_to_{assemAnn}.filt.paf"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 1
    run:
        paf_colnames = ["query",  "qlen", "qstart", "qend", "strand",
                         "target", "tlen", "tstart", "tend", "matches",
                         "alnlen", "mapq", "AS", "ms"]

        entries = []
        with open(input.paf, "r") as f:
            for line in f:
                # skip the header lines
                if line.startswith("#"):
                    continue
                # split the line
                line = line.strip().split("\t")[:len(paf_colnames)]
                line[-2] = line[-2].replace("AS:i:", "")
                line[-1] = line[-1].replace("ms:i:", "")
                for colindex in [1,2,3,6,7,8,9,10,11,12,13]:
                    line[colindex] = int(line[colindex])
                entries.append(line)
        # load in the paf file, only read in the first 14 columns
        pafdf = pd.DataFrame(entries, columns=paf_colnames)
        # Sort by query, then ms, then tlen. This will put the best hits at the top.
        # Then, drop duplicates of the query column. This will keep only the best hit.
        pafdf = pafdf.sort_values(by=["query", "ms", "tlen"], ascending=False
                 ).drop_duplicates(subset=["query"])
        # save the filtered paf file to the outout, but don't use the index or header
        pafdf.to_csv(output.paf, sep="\t", index=False, header=False)

rule paf_to_chrom_and_pep:
    """
    This script converts the filtered paf file to the chrom file.
    There are some decisions to make here about which proteins are best

    The chrom format is: protein_name, scaf_name, strand, scaf_start, scaf_end
    BFGX8636T1      sca1    +       1       1246
    BFGX0001T1      sca1    -       2059    2719
    BFGX0002T1      sca2    +       6491    12359
    BFGX0003T1      sca2    -       12899   18848

    Also outputs a protein fasta file of the entries.
    """
    input:
        paf = config["tool"] + "/output/mapped_reads/{assemAnn}/{LG_name}_to_{assemAnn}.filt.paf",
        pep  = config["tool"] + "/input/LG_proteins/{LG_name}.fasta"
    output:
        chrom = temp(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom"),
        pep   = temp(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep")
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 5
    run:
        paf_colnames = ["query",  "qlen", "qstart", "qend", "strand",
                         "target", "tlen", "tstart", "tend", "matches",
                         "alnlen", "mapq", "AS", "ms"]

        # read in the paf file as a pandas dataframe
        df = pd.read_csv(input.paf, sep="\t", header=None, names=paf_colnames)
        # Split the query column into two parts. The first part is up to the last
        # _ character, the second part is everything else. Save the first part as
        # the rbh column.
        df["rbh"]     = df["query"].str.rsplit("_", n=1, expand=True)[0]
        # the second part is the species name
        df["species"] = df["query"].str.rsplit("_", n=1, expand=True)[1]
        # sort based on the rbh column, then the ms column, then the tlen column
        df = df.sort_values(by=["rbh", "ms", "tlen"], ascending=False)
        # for now just keep the first entry for each rbh
        df = df.drop_duplicates(subset=["rbh"])
        chromdf = df[["query", "target", "strand", "tstart", "tend"]]
        # sort by target, then tstart, then tend
        chromdf = chromdf.sort_values(by=["target", "tstart", "tend"])
        # save the chromdf as a tab delimited file
        chromdf.to_csv(output.chrom, sep="\t", index=False, header=False)

        # now we generate the protein fasta file
        with open(output.pep, "w") as o:
            for record in fasta.parse(input.pep):
                if record.id in chromdf["query"].values:
                    o.write(">{}\n{}\n".format(record.id, record.seq))

rule gzChrom:
    """
    just gzips the chrom file
    """
    input:
        chrom = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom",
        pep   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep"
    output:
        chrom = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom.gz",
        pep   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep.gz"
    threads: 1
    resources:
        mem_mb  = 1000, # shouldn't take a lot of RAM.
        runtime = 5
    shell:
        """
        # first gzip the chrom file
        gzip < {input.chrom} > {output.chrom}
        # second gzip the pep file. This will be bigger.
        gzip < {input.pep} > {output.pep}
        """

rule generate_assembled_config_entry:
    """
    Print out a small piece of a yaml file specifically for ODP.
    These will be gathered and concatenated later.
    """
    input:
        unannot_genomes = ancient(config["unannotated_genome_chr_tsv"]),
        genome          = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz",
        protein = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep.gz",
        chrom   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom.gz"
    output:
        yaml   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.yaml.part",
    threads: 1
    resources:
        runtime = 5,
        mem_mb  = 1000
    params:
    run:
        # load in the dataframe of the annotated genomes
        df = pd.read_csv(input.unannot_genomes, sep="\t")
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()

        row = df.loc[df["Assembly Accession"] == wildcards.assemAnn]
        taxid = row["Organism Taxonomic ID"].values[0]

        # s is the output string.
        # h is headspace. Currently 2 spaces because the odp yaml files are compoased like so:
        #species:
        #  Sp_name:
        #    taxid:
        #    genus:
        #    species:
        #    assembly_accession:
        #    proteins:
        #    chrom:
        #    genome:
        #    minscaflen:
        #    ....

        # genus
        try:
            genus = row["Organism Name"].values[0].split(" ")[0]
        except:
            genus = row["Organism Name"].values[0]
        # species
        try:
            species = row["Organism Name"].values[0].split(" ")[1]
        except:
            species = "None"
        # cleanup the species name
        species = species.replace("sp.", "sp")
        assemAnn = wildcards.assemAnn

        orgnamvalue = row["Organism Name"].values[0]

        # now we strip unwanted characters from the genus, species name, and assembly accession
        strip_these_chars = [" ", "_", "-"]
        for char in strip_these_chars:
            genus = genus.replace(char, "")
            species = species.replace(char, "")
            assemAnn = assemAnn.replace(char, "")

        minscaflen = 100000
        if "|" in [genus, species, taxid, wildcards.assemAnn]:
            raise ValueError("You cannot have a | character in the genus, species, taxid, or assembly accession.")
        for thischar in strip_these_chars:
            if thischar in [genus, species, taxid, wildcards.assemAnn]:
                raise ValueError("You cannot have a {} character in the genus, species, taxid, or assembly accession.".format(thischar))
        # We currently cannot have "_" characters in the species name because of problems it causes downstream with snakemake.
        spstring = "{}{}-{}-{}".format(genus, species, taxid, assemAnn)
        # If there are more than 3 '-' characters in the string, then we have a problem.
        if spstring.count("-") > 3:
            raise ValueError("There are too many '-' characters in the string: {}".format(spstring))
        h = "  "
        s = ""
        s += h + "{}:\n".format(spstring)
        s += h + h + "assembly_accession: {}\n".format(wildcards.assemAnn)
        s += h + h + "taxid:              {}\n".format(taxid)
        s += h + h + "genus:              {}\n".format(genus)
        s += h + h + "species:            {}\n".format(species)
        s += h + h + "proteins:           {}\n".format(os.path.abspath(input.protein))
        s += h + h + "chrom:              {}\n".format(os.path.abspath(input.chrom))
        s += h + h + "genome:             {}\n".format(os.path.abspath(input.genome))
        s += h + h + "minscaflen:         {}\n".format(minscaflen)

        # write the string to the output
        with open(output.yaml, "w") as f:
            f.write(s)

rule collate_assembled_config_entries:
    """
    Concatenate all of the yaml files together into one big file.
    """
    input:
        yaml_parts = expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{{LG_name}}.yaml.part",
                            assemAnn=config["assemAnn"])
    output:
        yaml = "NCBI_odp_db.unannotated.{LG_name}.yaml"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 5
    shell:
        """
        echo "species:" > {output.yaml}
        cat {input.yaml_parts} >> {output.yaml}
        """

rule generate_species_list_for_timetree:
    """
    Generates a list of species that are in this database.
    One species per line.
    """
    input:
        yaml    = "NCBI_odp_db.unannotated.{LG_name}.yaml"
    output:
        sp_list = "NCBI_odp_sp_list.unannotated.{LG_name}.txt"
    threads: 1
    resources:
        runtime = 5,
        mem_mb  = 1000
    run:
        # open the yaml file into a dictionary
        with open(input.yaml, "r") as f:
            yaml_dict = yaml.load(f, Loader=yaml.FullLoader)
        # open the output file for writing
        with open(output.sp_list, "w") as f:
            # loop through the dictionary and write the species names to the file
            for sp in yaml_dict["species"]:
                f.write("{}\t{}\n".format(yaml_dict["species"][sp]["genus"], yaml_dict["species"][sp]["species"]))

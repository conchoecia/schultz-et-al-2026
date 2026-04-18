"""
Date: 2023-12-19 This is an updated version of the original script that focuses on genomes identified as chromosome-scale by NCBI.

Build a database of chromosome-scale genome assemblies using the NCBI website

PREREQUISITES:
  - This script requires the ete4 toolkit to be installed. This can be done with conda:
    https://etetoolkit.org/download/
  - You must also preload the taxonomy toolkit with the following commands:
    https://etetoolkit.github.io/ete/tutorial/tutorial_taxonomy.html#setting-up-local-copies-of-the-ncbi-and-gtdb-taxonomy-databases
    ```
    from ete4 import NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    ```

20250901 - TODO - there is a problem where sometimes assemblies appear to be annotated, but there are no clear peptide files:
         - This needs to be edited in this file. The whole process of deciding what is annotated would be helpful.
         - GCA_964026615.1 is one of those
         - GCA_919967415.2 is another
         - GCA_964023275.1 Norana najaformis is another that lacks a protein file.
"""

# Some specific NCBI taxids cause problems with the NCBI datasets tool.
# This one, GCA_900186335.3, causes a parsing error: https://github.com/ncbi/datasets/issues/300
hardcoded_ignore_accessions = ["GCA_900186335.3",
                               "GCA_000002165.1", # This is the 2009 Celera Genomics mouse genome. It was found to be contaminated, and too large.
                               "GCF_000002265.2" # This is the Celera Genomics rat genome. It is currently suppressed.
                               ]
import datetime
from datetime import timedelta
from ete4 import NCBITaxa
import os
import pandas as pd
import sys
from datetime import datetime

# import fasta parser from dependencies
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
bin_path          = os.path.join(snakefile_path, "bin")
dependencies_path = os.path.join(snakefile_path, "dependencies")

src_path = os.path.join(snakefile_path, "src")
sys.path.insert(1, src_path)
from InferDB import append_to_list, fields_to_print
from InferDB import return_stats_string
from InferDB import prefer_SRA_over_others, prefer_refseq, prefer_assembly_with_higher_N50
from InferDB import remove_specific_GCAs, get_best_contig_L50_assembly
from InferDB import legal_True_final_group_df, filter_raw_genome_df
from InferDB import load_and_cleanup_NCBI_datasets_tsv_df, dataset_summary_table

configfile: "chrombase.config.yaml"

if "datetime" not in config:
    # previously We recorded down to the hour and minute, but this felt a bit excessive
    #config["datetime"] = datetime.now().strftime('%Y%m%d%H%M')
    config["datetime"] = datetime.now().strftime('%Y%m%d')


config["tool"] = "odp_ncbi_genome_scraper"

wildcard_constraints:
    taxid="[0-9]+",
    datetime="[0-9]+"

rule all:
    input:
        # The files that will be used for downloading the database.
        expand(config["tool"] + "/output/for_downloading_humanfilt/annotated_chr_notembargoed_{datetime}.tsv",
                datetime = config["datetime"]),
        expand(config["tool"] + "/output/for_downloading_humanfilt/unannotated_chr_notembargoed_humanfilt_{datetime}.tsv",
                datetime = config["datetime"]),
        ## report of the final dataset
        #expand(config["tool"] + "/input/report_{taxid}_{datetime}.tsv",
        #        taxid = config["taxids"], datetime = config["datetime"]),
        ## report and plot of the filtered dataset
        #expand(config["tool"] + "/input/report_history_filtered_{taxid}_{datetime}.tsv",
        #        taxid = config["taxids"], datetime = config["datetime"]),
        #expand(config["tool"] + "/input/report_history_filtered_{taxid}_{datetime}.pdf",
        #        taxid = config["taxids"], datetime = config["datetime"]),
        ## report and plot of the unfiltered dataset. This is more appropriate for looking at the growth of NCBI genomes.
        #expand(config["tool"] + "/input/report_history_raw_{taxid}_{datetime}.tsv",
        #        taxid = config["taxids"], datetime = config["datetime"]),
        #expand(config["tool"] + "/input/report_history_raw_{taxid}_{datetime}.pdf",
        #        taxid = config["taxids"], datetime = config["datetime"]),

rule download_json:
    output:
        genome_report = config["tool"] + "/input/{taxid}_{datetime}.json"
    params:
        datasets = os.path.join(bin_path, "datasets")
    threads: 1
    resources:
        time    = 5, # 5 minutes
        runtime = 5,
        mem_mb = 1000
    shell:
        """
        {params.datasets} summary genome taxon --as-json-lines {wildcards.taxid} > {output.genome_report}
        """

rule format_json_to_tsv:
    input:
        genome_report = config["tool"] + "/input/{taxid}_{datetime}.json",
    output:
        report_tsv = config["tool"] + "/input/{taxid}_{datetime}.tsv" # NOTE - long-term should be temp()
    params:
        fields = ",".join(fields_to_print),
        dataformat = os.path.join(bin_path, "dataformat")
    threads: 1
    resources:
        time    = 5, # 5 minutes
        runtime = 5,
        mem_mb = 1000
    shell:
        """
        {params.dataformat} tsv genome --inputfile {input.genome_report} --fields {params.fields} > {output.report_tsv}
        """

rule accession_list_raw:
    """
    This gets a list of all the accessions from the genome tsv.
    The intention with this is to download to do a dry run.

    Format: One accession per line in the file.
    """
    input:
        report_tsv = config["tool"] + "/input/{taxid}_{datetime}.tsv"
    output:
        accession_list = config["tool"] + "/input/accession_list_{taxid}_{datetime}.txt"
    threads: 1
    resources:
        time = 5,
        runtime = 5,
        mem_mb = 1000
    run:
        df = pd.read_csv(input.report_tsv, sep="\t", dtype=str)
        # output the unique entries of the Assembly Accession column
        accession_list = sorted(df["Assembly Accession"].unique())
        with open(output.accession_list, "w") as f:
            for accession in accession_list:
                f.write(accession + "\n")

rule dryrun_download_all_annotations:
    """
    Uses the dryrun mode of the datasets tool to download the json file for those assemblies.
    This json will have all the info about whether the necessary annotations exist.

    Use the datasets download --dehydrated tool.
    datasets download genome accession GCA_964023275.1,GCA_964026615.1  --include gff3,protein --dehydrated --filename temp2.zip
    """
    input:
        accession_list = config["tool"] + "/input/accession_list_{taxid}_{datetime}.txt"
    output:
        zipped = config["tool"] + "/input/dehydrated_annot_{taxid}_{datetime}.zip"
    params:
        datasets = os.path.join(bin_path, "datasets")
    threads: 1
    resources:
        time    = 10,
        runtime = 10,
        mem_mb  = 1000
    shell:
        """
        {params.datasets} download genome accession --inputfile {input.accession_list} \
            --include gff3,protein --dehydrated --filename {output.zipped}
        """

rule unzip_curation_json:
    """
    We only need one file from this download - the file that has the info about whether
    both the gff and protein fasta files are present.
    """
    input:
        zipped = config["tool"] + "/input/dehydrated_annot_{taxid}_{datetime}.zip"
    output:
        json = config["tool"] + "/input/dehydrated_annot_{taxid}_{datetime}_dataset_catalog.json"
    threads: 1
    resources:
        time    = 2,
        runtime = 2,
        mem_mb  = 1000
    shell:
        """
        unzip -p {input.zipped} ncbi_dataset/data/dataset_catalog.json > {output.json}
        """

rule dataset_catalog_json_to_tsv:
    """
    - What this does is convert the dataset catalog JSON file to a TSV file.
    - The TSV contains three columns: accession, contains_GFF3, and contains_PROTEIN_FASTA.
    - This is used to note whether something is properly annotated, since the normal json from
      NCBI doesn't contain information on whether both the GFF3 and protein fasta file are present.
    """
    input:
        json = config["tool"] + "/input/dehydrated_annot_{taxid}_{datetime}_dataset_catalog.json"
    output:
        tsv = config["tool"] + "/input/dehydrated_annot_{taxid}_{datetime}_dataset_catalog.tsv"
    threads: 1
    resources:
        time    = 2,
        runtime = 2,
        mem_mb  = 1000
    run:
        import json
        import sys
        import csv
        with open(input.json) as f:
            data = json.load(f)
        with open(output.tsv, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["accession", "contains_GFF3", "contains_PROTEIN_FASTA"])
            for asm in data["assemblies"]:
                accession = asm.get("accession")
                if accession is None:
                    continue
                filetypes = {f["fileType"] for f in asm.get("files", [])}
                contains_gff = "GFF3" in filetypes
                contains_protein = "PROTEIN_FASTA" in filetypes
                writer.writerow([accession, int(contains_gff), int(contains_protein)])

rule mark_annotations:
    """
    We now mark which assemblies are properly annotated - which ones have both gff3 and protein fasta
      files. We will mark the ones that have both of these files as being annotated. The presence of
      an annotation in the NCBI json doesn't record whether both of these essential files are present.
    """
    input:
        files_tsv = config["tool"] + "/input/dehydrated_annot_{taxid}_{datetime}_dataset_catalog.tsv",
        report_tsv = config["tool"] + "/input/{taxid}_{datetime}.tsv" # NOTE - long-term should be temp()
    output:
        report_tsv = config["tool"] + "/input/{taxid}_{datetime}_annotmarked.tsv"
    threads: 1
    resources:
        time    = 2,
        runtime = 2,
        mem_mb  = 4000
    run:
        files_df  = pd.read_csv(input.files_tsv, sep="\t")
        files_df["contains_GFF3"]          = files_df["contains_GFF3"].astype(int)
        files_df["contains_PROTEIN_FASTA"] = files_df["contains_PROTEIN_FASTA"].astype(int)

        report_df = pd.read_csv(input.report_tsv, sep="\t", low_memory=False)

        # map from files_df → report_df, fill missing with 0, and cast to int
        report_df["has_GFF"] = (
            report_df["Assembly Accession"]
            .map(files_df.set_index("accession")["contains_GFF3"])
            .fillna(0)
            .astype(int)
        )

        report_df["has_protein"] = (
            report_df["Assembly Accession"]
            .map(files_df.set_index("accession")["contains_PROTEIN_FASTA"])
            .fillna(0)
            .astype(int)
        )

        # AND logic (already int 0/1, so just multiply)
        report_df["is_annotated"] = (report_df["has_GFF"] & report_df["has_protein"]).astype(int)

        # reorder columns
        new_order = ["Assembly Accession", "has_GFF", "has_protein", "is_annotated"]
        report_df = report_df[new_order + [c for c in report_df.columns if c not in new_order]]

        # save
        report_df.to_csv(output.report_tsv, sep="\t", index=False)

rule get_representative_genomes:
    """
    This rule is responsible for parsing the genome report and selecting the representative genomes
     for each species.

    Currently this does not support multiple genomes for one species.
    """
    input:
        report_tsv = config["tool"] + "/input/{taxid}_{datetime}_annotmarked.tsv",
        assembly_ignore_list = os.path.join(snakefile_path, "data/assembly_ignore_list.txt")
    output:
        report                 = config["tool"] + "/input/report_{taxid}_{datetime}.tsv",
        representative_genomes = config["tool"] + "/input/selected_genomes_{taxid}_{datetime}.tsv"
    threads: 1
    resources:
        time    = 5, # 5 minutes
        runtime = 5,
        mem_mb = 1000
    run:
        # first we get the list of things to ignore
        # remove things that are in assembly_ignore_list.txt
        ignore_list = []
        with open(input.assembly_ignore_list, "r") as f:
            for line in f:
                # remove leading and trailing whitespace
                line = line.strip()
                # ignore lines that start with a comment character
                if line.startswith("#"):
                    continue
                # ignore empty lines
                if len(line) == 0:
                    continue
                # Add this line to the set of assemblies to ignore
                # Get the entry as the string until the first whitespace
                # Should be tab or a space character.
                entry = line.split()[0]
                ignore_list.append(entry)

        # add hardcoded_ignore_accessions to the ignore_list
        ignore_list = set(ignore_list + hardcoded_ignore_accessions)

        # Load in and clean up the dataframe.
        # This mostly changes certain column types to ints for consistent processing later on.
        df = load_and_cleanup_NCBI_datasets_tsv_df(input.report_tsv, ignore_list)

        # Perform filtering of the dataset to make sure that we have the best possible genomes.
        # This does the actual removal of rows of assemblies
        df = filter_raw_genome_df(df, hardcoded_ignore_accessions)

        # print out the dataset summary table
        summarydf = dataset_summary_table(df)
        summarydf.to_csv(output.report, sep="\t", index=False)

        # make a new column called Lineage. Get the NCBI Taxa lineage from ete4 NCBITaxa
        ncbi = NCBITaxa()
        taxid_dict = {taxid: ";".join([str(x) for x in ncbi.get_lineage(taxid)]) for taxid in list(df["Organism Taxonomic ID"].unique())}
        df["Lineage"] = df["Organism Taxonomic ID"].map(taxid_dict)

        # save the dataframe to the output
        df.to_csv(output.representative_genomes, sep="\t", index=False)

# this checkpoint triggers re-evaluation of the DAG
checkpoint split_into_annotated_and_unannotated_and_chr_nonchr:
    """
    This rule is responsible for parsing the genome report and selecting the representative genomes

     for each species.

    Currently this does not support multiple genomes for one species.
    """
    input:
        representative_genomes = expand(config["tool"] + "/input/selected_genomes_{taxid}_{datetime}.tsv", taxid = config["taxids"], datetime = config["datetime"])
    output:
        # DO NOT CHANGE THE ORDER OF THESE FILES. THE FUNCTION get_assemblies(wildcards) DEPENDS ON IT
        annotated_genomes_chr       =       config["tool"] + "/output/annotated_genomes_chr_{datetime}.tsv",
        annotated_genomes_nonchr    =       config["tool"] + "/output/annotated_genomes_nonchr_{datetime}.tsv",
        unannotated_genomes_chr     =       config["tool"] + "/output/unannotated_genomes_chr_{datetime}.tsv",
        unannotated_genomes_nonchr  =       config["tool"] + "/output/unannotated_genomes_nonchr_{datetime}.tsv"
    threads: 1
    resources:
        time   = 5, # 5 minutes
        mem_mb = 1000
    run:
        list_of_annotated_chr_dfs      = []
        list_of_annotated_nonchr_dfs   = []
        list_of_unannotated_chr_dfs    = []
        list_of_unannotated_nonchr_dfs = []
        # load in the dataframe
        for thisfile in input.representative_genomes:
            df = pd.read_csv(thisfile, sep="\t")
            # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
            df.columns = df.columns.str.strip()

            # annotated genomes.
            # use the "annotated" and "chrscale" columns to determine what goes into what
            annot_chrom_df    = df.loc[df["annotated"] == True].loc[df["chrscale"] == True]
            annot_notchrom_df = df.loc[df["annotated"] == True].loc[df["chrscale"] == False]
            list_of_annotated_chr_dfs.append(      annot_chrom_df)
            list_of_annotated_nonchr_dfs.append(annot_notchrom_df)

            # unannotated genomes.
            # use the "annotated" and "chrscale" columns to determine what goes into what
            unannot_chrom_df    = df.loc[df["annotated"] == False].loc[df["chrscale"] == True]
            unannot_notchrom_df = df.loc[df["annotated"] == False].loc[df["chrscale"] == False]
            list_of_unannotated_chr_dfs.append(      unannot_chrom_df)
            list_of_unannotated_nonchr_dfs.append(unannot_notchrom_df)

        # put all the annotated dataframes together
        annot_chr_df    = pd.concat(list_of_annotated_chr_dfs)
        annot_nonchr_df = pd.concat(list_of_annotated_nonchr_dfs)
        unann_chr_df    = pd.concat(list_of_unannotated_chr_dfs)
        unann_nonchr_df = pd.concat(list_of_unannotated_nonchr_dfs)

        # Remove duplicate rows that may have been picked up by nested taxids.
        #   For example, if we pick "Metazoa" and "Arthropoda", there will be many duplicate rows.
        annot_chr_df    = annot_chr_df.drop_duplicates()
        annot_nonchr_df = annot_nonchr_df.drop_duplicates()
        unann_chr_df    = unann_chr_df.drop_duplicates()
        unann_nonchr_df = unann_nonchr_df.drop_duplicates()

        # save the dataframe to the output
        annot_chr_df.to_csv(     output.annotated_genomes_chr,      sep="\t", index=False)
        annot_nonchr_df.to_csv(  output.annotated_genomes_nonchr,   sep="\t", index=False)
        unann_chr_df.to_csv(     output.unannotated_genomes_chr,    sep="\t", index=False)
        unann_nonchr_df.to_csv(  output.unannotated_genomes_nonchr, sep="\t", index=False)

rule filter_embargoed_genomes:
    """
    This goes through the annotated and unannotated genome files and filters these based on
      embargo status. It uses a conservative option to future-proof databases s.t. even if
      annotations are released at the end of the embargo period.

    outfiles for one:
      annotated_chr_20250829_all_20250903.tsv
      annotated_chr_20250829_all_20250903_fewerColumns.tsv
      annotated_chr_20250829_embargoed_20250903.tsv
      annotated_chr_20250829_embargoed_20250903_fewerColumns.tsv
      annotated_chr_20250829_notembargoed_20250903.tsv
      annotated_chr_20250829_notembargoed_20250903_fewerColumns.tsv
      annotated_chr_20250829_report.txt
    """
    input:
        annotat = config["tool"] + "/output/annotated_genomes_chr_{datetime}.tsv",
        unannot = config["tool"] + "/output/unannotated_genomes_chr_{datetime}.tsv"
    output:
        annotated_all   = config["tool"] + "/output/embargo_filtered/annotated_chr_all_{datetime}.tsv",
        annotated_emb   = config["tool"] + "/output/embargo_filtered/annotated_chr_embargoed_{datetime}.tsv",
        annotated_non   = config["tool"] + "/output/embargo_filtered/annotated_chr_notembargoed_{datetime}.tsv",
        annotated_rep   = config["tool"] + "/output/embargo_filtered/annotated_chr_report_{datetime}.txt",
        unannotat_all   = config["tool"] + "/output/embargo_filtered/unannotated_chr_all_{datetime}.tsv",
        unannotat_emb   = config["tool"] + "/output/embargo_filtered/unannotated_chr_embargoed_{datetime}.tsv",
        unannotat_non   = config["tool"] + "/output/embargo_filtered/unannotated_chr_notembargoed_{datetime}.tsv",
        unannotat_rep   = config["tool"] + "/output/embargo_filtered/unannotated_chr_report_{datetime}.txt"
    threads: 1
    params:
        genbargo = os.path.join(dependencies_path, "genbargo/filter_assemblies.py"),
        out_dir = config["tool"] + "/output/embargo_filtered/",
        annotated_prefix   = "annotated_chr",
        unannotated_prefix = "unannotated_chr"
    resources:
        time    = 5, # 5 minutes
        runtime = 5,
        mem_mb = 1000
    shell:
        """
        # check the annotated genomes first
        python {params.genbargo} -t {input.annotat} -c \
          -p {params.annotated_prefix} -d {params.out_dir}

        # Now check the unannotated genomes
        python {params.genbargo} -t {input.unannot} -c \
          -p {params.unannotated_prefix} -d {params.out_dir}
        """

rule db_for_downloading_limited_human:
    """
    This filters the humans from the nonannotated, non-embargoed genomes.
    Takes that filtered file and the annotated, chr-scale genomes and saves them to a new directory.
    """
    input:
        annotated_non   = config["tool"] + "/output/embargo_filtered/annotated_chr_notembargoed_{datetime}.tsv",
        unannotat_non   = config["tool"] + "/output/embargo_filtered/unannotated_chr_notembargoed_{datetime}.tsv"
    output:
        annotated_non   = config["tool"] + "/output/for_downloading_humanfilt/annotated_chr_notembargoed_{datetime}.tsv",
        unannotat_non   = config["tool"] + "/output/for_downloading_humanfilt/unannotated_chr_notembargoed_humanfilt_{datetime}.tsv"
    threads: 1
    resources:
        time    = 5, # 5 minutes
        runtime = 5,
        mem_mb = 1000
    run:
        # copy the annotated file (python)
        import shutil
        shutil.copyfile(input.annotated_non, output.annotated_non)

        # For the non-annotated, read in with pandas and filter "Homo sapiens"
        #  out of the column: "Assembly BioSample Description Organism Name"
        df = pd.read_csv(input.unannotat_non, sep="\t")
        df = df.loc[df["Organism Name"] != "Homo sapiens"]
        df = df.loc[df["Organism Taxonomic ID"] != 9606]
        df.to_csv(output.unannotat_non, sep="\t", index=False)

rule history_of_assemblies_filtered:
    """
    Takes in the genome report and outputs the stats of the dataframe for different points in time.
      - This rule is for the filtered version of the dataset, in which each species can occur in multiple categories.
      - This is better for looking at the datasets that will be used for whole-genome comparisons.
    Step backward 7 days in time until we run out of assemblies to consider.
    """
    input:
        report_tsv = config["tool"] + "/input/{taxid}_{datetime}.tsv",
    output:
        report     = config["tool"] + "/input/report_history_filtered_{taxid}_{datetime}.tsv",
    threads: 1
    resources:
        time    = 5, # 5 minutes
        runtime = 5,
        mem_mb = 1000
    params:
        day_step = 7
    run:
        historical_view_dfs = []
        # load in and clean up the dataframe
        df = load_and_cleanup_NCBI_datasets_tsv_df(input.report_tsv, hardcoded_ignore_accessions)

        # go back in time 7 days at a time until we run out of assemblies to consider
        # go back until January 2000
        jan2000 = datetime.strptime("2000-01-01", '%Y-%m-%d')
        current_date = datetime.today().strftime('%Y-%m-%d')
        while current_date >= jan2000.strftime('%Y-%m-%d'):
            # print to sys.stderr in a progress-bar type configuration that just prints out the date
            print("   Filtering on or before date: {}".format(current_date), file=sys.stderr, end="\r")
            # filter the dataframe to only include assemblies that were released before seven_days_ago
            dftemp = df.loc[df["Assembly Release Date"] <= current_date]
            # annotations may have come at a later date
            dftemp.loc[dftemp["Annotation Release Date"] > current_date, "Annotation Release Date"] = np.nan
            # if there are no assemblies left, then we are done
            # print out the dataset summary table
            dftemp = filter_raw_genome_df(dftemp, [], suppress_text = True)
            summarydf = dataset_summary_table(dftemp, assembly_release_date = current_date)
            historical_view_dfs.append(summarydf)
            current_date = datetime.strptime(current_date, '%Y-%m-%d')
            current_date = (current_date - timedelta(days=params.day_step)).strftime('%Y-%m-%d')
        # print a newline to keep the progress bar on the screen/in the log
        print("   Filtering on or before date: {}".format(current_date), file=sys.stderr)
        # concatenate all the dataframes together
        concatdf = pd.concat(historical_view_dfs)
        # save to the output
        concatdf.to_csv(output.report, sep="\t", index=False)

rule assembly_report_plot_filtered:
    """
    Make a pdf of the filtered dataset assembly report.
    """
    input:
        report          = config["tool"] + "/input/report_history_filtered_{taxid}_{datetime}.tsv",
        plotting_script = os.path.join(snakefile_path, "scripts/plot_NCBI_genomes_history.py")
    output:
        pdf             = config["tool"] + "/input/report_history_filtered_{taxid}_{datetime}.pdf",
    threads: 1
    resources:
        time    = 1, # 5 minutes
        runtime = 1,
        mem_mb = 500
    shell:
        """
        python {input.plotting_script} -i {input.report} -o {output.pdf}
        """

rule history_of_assemblies_raw:
    """
    Takes in the genome report and outputs the stats of the dataframe for different points in time.
      - This rule is for the raw version of the dataset, in which each species can only occur in one category.
      - This is more accurate for looking at the absolute growth of the NCBI dataset over time.
    Step backward 7 days in time until we run out of assemblies to consider.
    """
    input:
        report_tsv = config["tool"] + "/input/{taxid}_{datetime}.tsv",
    output:
        report     = config["tool"] + "/input/report_history_raw_{taxid}_{datetime}.tsv",
    threads: 1
    resources:
        time    = 5, # 5 minutes
        runtime = 5,
        mem_mb = 1000
    params:
        day_step = 7
    run:
        historical_view_dfs = []
        # load in and clean up the dataframe
        df = load_and_cleanup_NCBI_datasets_tsv_df(input.report_tsv, hardcoded_ignore_accessions)

        # go back in time 7 days at a time until we run out of assemblies to consider
        # go back until January 2000
        jan2000 = datetime.strptime("2000-01-01", '%Y-%m-%d')
        current_date = datetime.today().strftime('%Y-%m-%d')
        while current_date >= jan2000.strftime('%Y-%m-%d'):
            # print to sys.stderr in a progress-bar type configuration that just prints out the date
            print("   Filtering on or before date: {}".format(current_date), file=sys.stderr, end="\r")
            # filter the dataframe to only include assemblies that were released before seven_days_ago
            dftemp = df.loc[df["Assembly Release Date"] <= current_date]
            # if there are no assemblies left, then we are done
            # print out the dataset summary table
            # sort by Assembly Accession, then Annotation release date, preferring the ones with annotations first
            dftemp = dftemp.sort_values(by=["Assembly Accession", "Annotation Release Date"], ascending=[True, False])
            # drop duplicates, keeping the first one
            dftemp = dftemp.drop_duplicates(subset=["Assembly Accession"], keep="first")
            # add the columns called "chrscale" and "annotated", set them to False
            dftemp = dftemp.assign(chrscale=False, annotated=False)
            # for "Annotation Release Date", make the cells after current date NaN
            dftemp.loc[dftemp["Annotation Release Date"] > current_date, "Annotation Release Date"] = np.nan
            # if Annotation Release Date is NaN, then set "annotated" to False, otherwise True. We already handled False, just do True
            dftemp.loc[~dftemp["Annotation Release Date"].isna(), "annotated"] = True
            # if Assembly level == "Chromosome", then set "chrscale" to True
            dftemp.loc[dftemp["Assembly Level"] == "Chromosome", "chrscale"] = True
            summarydf = dataset_summary_table(dftemp, assembly_release_date = current_date)
            historical_view_dfs.append(summarydf)
            current_date = datetime.strptime(current_date, '%Y-%m-%d')
            current_date = (current_date - timedelta(days=params.day_step)).strftime('%Y-%m-%d')
        # print a newline to keep the progress bar on the screen/in the log
        print("   Filtering on or before date: {}".format(current_date), file=sys.stderr)
        # concatenate all the dataframes together
        concatdf = pd.concat(historical_view_dfs)
        # save to the output
        concatdf.to_csv(output.report, sep="\t", index=False)

rule assembly_report_plot_raw:
    """
    Make a pdf of the raw dataset assembly report.
    """
    input:
        report          = config["tool"] + "/input/report_history_raw_{taxid}_{datetime}.tsv",
        plotting_script = os.path.join(snakefile_path, "scripts/plot_NCBI_genomes_history.py")
    output:
        pdf             = config["tool"] + "/input/report_history_raw_{taxid}_{datetime}.pdf",
    threads: 1
    resources:
        time  = 1, # 5 minutes
        runtime = 1,
        mem_mb = 500
    shell:
        """
        python {input.plotting_script} -i {input.report} -o {output.pdf}
        """

## TODO -implement this rule to summarize all the chromosome-scale genomes
#rule generate_report_allchr:
#    """
#    Generate a report that summarizes the genome provenance of all the chromosome-scale genomes.
#    """
#    input:
#        annotated_non   = config["tool"] + "/output/for_downloading_humanfilt/annotated_chr_notembargoed_{datetime}.tsv",
#        unannotat_non   = config["tool"] + "/output/for_downloading_humanfilt/unannotated_chr_notembargoed_humanfilt_{datetime}.tsv"
#    output:
#        report          = config["tool"] + "/output/for_downloading_humanfilt/report_allchr_{datetime}.txt"
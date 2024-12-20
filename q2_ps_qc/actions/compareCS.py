#!/usr/bin/env python
import os
import pandas as pd
import fastatools as ft
import time
from collections import defaultdict

def compareCS(
        ctx,
        metadata_file,
        codename_column,
        parent_codename_column,
        zscores_file,
        fasta_file,
        fullname_column=None,
        generate_epitope_data=False,
        min_zscore=8,
        min_zscore_diff=0.5,
        max_zscore_diff=0.1,
        pep_seq_len=30,
        histogram_bins=6,
        min_epitope_size=7,
        data_output_dir="./output_data"
):
    start_time = time.time()

    assert not os.path.exists(data_output_dir), f"{data_output_dir} already exists! Please move or delete it and try again."
    os.mkdir(data_output_dir)

    scatter_plot = ctx.get_action("ps-plot", "compareCS_scatter")
    histogram = ctx.get_action("ps-plot", "compareCS_histogram")

    # read in fasta file
    fasta_dict = ft.read_fasta_dict(fasta_file)

    # read in metadata file
    metadata = pd.read_csv(metadata_file, sep="\t")

    # validate input columns
    if codename_column not in metadata.columns:
        raise ValueError(f"{codename_column} column not found in metadata file.")
    elif parent_codename_column not in metadata.columns:
        raise ValueError(f"{parent_codename_column} column not found in metadata file.")
    elif fullname_column not in metadata.columns:
        if generate_epitope_data:
            raise ValueError(f"{fullname_column} column not found in metadata file. Turn off the generation of epitope data if fullname is not included.")
        else:   
            # extract columns
            metadata = metadata[[codename_column, parent_codename_column]]
    else:
        metadata = metadata[[fullname_column, codename_column, parent_codename_column]]

    # validate zscores (doing this to avoid using qiime2 format types honestly)
    with open(zscores_file) as fh:
        for _, line in zip(range(1), fh):
            if not line.startswith("Sequence name\t"):
                raise ValueError('TSV does not start with "Sequence name"')
        
    # read in zscores file
    zscores = pd.read_csv(zscores_file, sep="\t").set_index('Sequence name')

    metadata_has_s_versions = metadata[metadata[parent_codename_column].notna()]

    scatterplot, chart_data = generate_scatterplot(scatter_plot, zscores, metadata_has_s_versions, parent_codename_column, codename_column, min_zscore, fasta_dict)
    chart_data.to_csv(os.path.join(data_output_dir, "reactivity_plot_data.tsv"), sep="\t", index=False)

    
    if generate_epitope_data:
        codename_2_fullname, ordered_fullnames = get_fullnames(metadata_has_s_versions, fullname_column, parent_codename_column)

        # map epitopes that are more reactive with the C versions and those that are more reactive with the S versions
        max_distance = pep_seq_len-min_epitope_size
        raw_epitope_df, formatted_epitope_df = map_epitopes(chart_data, codename_2_fullname, ordered_fullnames, fullname_column, min_zscore_diff, max_distance, fasta_dict)

        raw_epitope_df.to_csv(os.path.join(data_output_dir, "raw_epitope_data.tsv"), sep="\t")
        formatted_epitope_df.to_csv(os.path.join(data_output_dir, "formatted_epitope_data.tsv"), sep="\t", index=False)


    # sample_name: category_name: peptide_sequences
    sample_category_peptides = seperate_peptides_into_categories(chart_data, min_zscore_diff, max_zscore_diff, fasta_dict)

    c_count_summary_df = get_c_count_summaries(sample_category_peptides)
    c_count_summary_df.to_csv(os.path.join(data_output_dir, "c_count_summary.tsv"), sep="\t", index=False)

    histogram = generate_c_count_histogram(histogram, sample_category_peptides, histogram_bins, pep_seq_len)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {round(elapsed_time, 2)} seconds")

    return scatterplot, histogram


def generate_c_count_histogram(
    histogram,
    sample_category_peptides: dict, 
    histogram_bins: int,
    pep_seq_len: int
):
    # automatically populate dict stucture
    c_count_data = defaultdict(lambda: histogram_dict(pep_seq_len))

    # create an option for all samples
    c_count_data["all"]

    for sample_name, category_2_parent_peptides in sample_category_peptides.items():
        for category_name, parent_peptides in category_2_parent_peptides.items():
            for parent_pep_seq in parent_peptides:
                for pos, aa in enumerate(parent_pep_seq):
                    if pos>=pep_seq_len:
                        raise IndexError(f"A peptide sequence length exceeds pep-seq-len: {pep_seq_len}")
                    if aa.lower() == 'c':
                        c_count_data[sample_name][category_name]["C count"][pos] += 1
                        c_count_data["all"][category_name]["C count"][pos] += 1

    histogram_data = {
        "Sample Name": list(),
        "Category Name": list(),
        "C count": list(),
        "Position": list()
    }

    for sample_name, category_dicts in c_count_data.items():
        for category_name, count_dict in category_dicts.items():
            for i in range(pep_seq_len):
                histogram_data["Sample Name"].append(sample_name)
                histogram_data["Category Name"].append(category_name)
                histogram_data["C count"].append(count_dict["C count"][i])
                histogram_data["Position"].append(count_dict["Position"][i])
    
    chart, = histogram(
        sample_names = histogram_data["Sample Name"],
        category_names = histogram_data["Category Name"],
        c_counts = histogram_data["C count"],
        positions = histogram_data["Position"],
        num_bins = histogram_bins
    )

    return chart


def create_count_structure(pep_seq_len):
    return {
        "C count": [0] * pep_seq_len,
        "Position": list(range(pep_seq_len))
    }


def histogram_dict(pep_seq_len):
    return defaultdict(lambda: create_count_structure(pep_seq_len))


# TODO: make this generalizable
def get_fullnames(metadata, fullname_column, parent_codename_column):
    # map each base full name to codename, just use C columns
    # note: need to remove " CtoS" from S version peptides' fullnames
    codename_2_fullname = defaultdict()
    ordered_fullnames = list()
    for i, row in metadata.iterrows():
        fullname = row[fullname_column][0:-len(" CtoS")]
        ordered_fullnames.append(fullname)
        codename_2_fullname[row[parent_codename_column]] = fullname
    
    return codename_2_fullname, ordered_fullnames


def generate_scatterplot(scatter_plot, zscores, metadata, parent_codename_column, codename_column, min_zscore, fasta_dict):
    # get variants
    variants = zscores.columns.to_list()

    # format data to have {c_version: s_version}
    c_to_s = {row[parent_codename_column]: row[codename_column] for i, row in metadata.iterrows()}

    # filter zscore matrix to only have scores for peptides with a c and s version
    filtered_zscores = zscores.loc[zscores.index.isin(list(c_to_s.keys())+list(c_to_s.values()))]

    # create data dict
    data_dict = {
        "C Z score": list(),
        "S Z score": list(),
        "C codename": list(),
        "S codename": list(),
        "Parent Sequence": list(),
        "Sample Name": list()
    }

    for variant in variants:
        for c_version, s_version in c_to_s.items():
            c_zscore = float(filtered_zscores.loc[c_version, variant])
            s_zscore = float(filtered_zscores.loc[s_version, variant])
            # check if either meets the min zscore threshold
            if c_zscore >= min_zscore or s_zscore >= min_zscore:
                data_dict["C Z score"].append(c_zscore)
                data_dict["S Z score"].append(s_zscore)
                data_dict["C codename"].append(c_version)
                data_dict["S codename"].append(s_version)
                data_dict["Parent Sequence"].append(fasta_dict[c_version])
                data_dict["Sample Name"].append(variant)

    # generate graph
    chart, = scatter_plot(
        x = data_dict["C Z score"],
        y = data_dict["S Z score"],
        c_codenames = data_dict["C codename"],
        s_codenames = data_dict["S codename"],
        parent_sequences = data_dict["Parent Sequence"],
        sample_names = data_dict["Sample Name"]
    )

    # output chart data
    return chart, pd.DataFrame.from_dict(data_dict)

def get_c_count_summaries(
    sample_category_peptides: dict,
    c_higher_category_name: str = "C_Higher",
    s_higher_category_name: str = "S_Higher",
    similar_category_name: str = "Similar"
):
    summary_data = list()

    for sample_name, category_2_parent_peptides in sample_category_peptides.items():

        # keep track of counts in this format = category: num_cysteines: num_peptides
        peptide_counts = {
            c_higher_category_name: defaultdict(int),
            s_higher_category_name: defaultdict(int),
            similar_category_name: defaultdict(int)
        }

        for category_name, parent_peptides in category_2_parent_peptides.items():
            for parent_pep_seq in parent_peptides:
                # count number of peptides that have that number of c's in the parent peptide
                peptide_counts[category_name][get_c_count(parent_pep_seq)] += 1
        
        # find the max number of cysteines
        max_num_c = get_max_num_c(peptide_counts)

        # write to output data
        for num_c in range(1, max_num_c+1):
            for category in peptide_counts.keys():
                number_peptides = peptide_counts[category][num_c]
                number_total_peptides = len(category_2_parent_peptides[category])
                proportion_peptides = number_peptides / number_total_peptides if number_total_peptides > 0 else 0
                summary_data.append((sample_name, category, num_c, number_peptides, round(proportion_peptides, 2)))
    
    c_count_summary_df = pd.DataFrame(summary_data, columns=["Sample", "Category", "Number_Cysteines", "Number_Peptides", "Proportion_Peptides"])

    return c_count_summary_df


def seperate_peptides_into_categories(
    data_df: pd.DataFrame,
    min_zscore_diff: float,
    max_zscore_diff: float,
    fasta_dict: dict,
    c_higher_category_name: str = "C_Higher",
    s_higher_category_name: str = "S_Higher",
    similar_category_name: str = "Similar"
):
    # group by sample
    grouped_data = data_df.groupby("Sample Name")

    # sample_name: category_name: peptide_sequences
    sample_category_peptides = dict()

    for sample_name, sample_df in grouped_data:
        sample_category_peptides[sample_name] = {
            c_higher_category_name: list(),
            s_higher_category_name: list(),
            similar_category_name: list()
        }

        # loop through each peptide
        for i, row in sample_df.iterrows():
            parent_pep = row["C codename"]
            c_zscore = row["C Z score"]
            s_zscore = row["S Z score"]

            category_name = None

            # check that the difference ratio meets the min difference threshold
            if (abs(c_zscore-s_zscore)/max([c_zscore, s_zscore])) >= min_zscore_diff:
                # check if c is more reactive
                if c_zscore > s_zscore:
                    category_name = c_higher_category_name

                # check if s is more reactive
                elif c_zscore < s_zscore:
                    category_name = s_higher_category_name

            # check if both versions are "equally" reactive
            elif (abs(c_zscore-s_zscore)/max([c_zscore, s_zscore])) <= max_zscore_diff:
                category_name = similar_category_name
            
            # check if a category name was set
            if category_name:
                parent_pep_seq = fasta_dict[parent_pep]
                # seperate the raw peptide sequences into categories
                sample_category_peptides[sample_name][category_name].append(parent_pep_seq)

    return sample_category_peptides 


# find the max number of cysteines in summary peptide counts dictionary
def get_max_num_c(peptide_counts):
    max_num = 0
    for category in peptide_counts.keys():
        all_num_c = list(peptide_counts[category].keys())
        if all_num_c:
            max_c_for_this_category = max(all_num_c)
            if max_c_for_this_category > max_num:
                max_num = max_c_for_this_category
    
    return max_num


# get the count of cysteine amino acids in a peptide sequence
def get_c_count(peptide_seq):
    return peptide_seq.lower().count('c')


 # map epitopes that are more reactive with the C versions and those that are more reactive with the S versions
def map_epitopes(
        data_df: pd.DataFrame, 
        codename_2_fullname: dict, 
        ordered_fullnames: list, 
        fullname_column: str,
        min_zscore_diff: float,
        max_distance: int,
        fasta_dict: dict
    ):
    raw_epitope_df = pd.DataFrame(ordered_fullnames, columns=[fullname_column]).set_index(fullname_column)

    # group by sample
    grouped_data = data_df.groupby("Sample Name")

    epitopes_data = list()

    for sample_name, sample_df in grouped_data:
        raw_sample_reactivity = list()

        # get reactivity data
        c_reactivity_data = list()
        s_reactivity_data = list()

        # loop through each peptide
        for i, row in sample_df.iterrows():
            # note: codename_2_fullname uses C codenme
            fullname = codename_2_fullname[row["C codename"]]

            sequence_name, start_pos, end_pos = extract_fullname_data(fullname)

            c_zscore = row["C Z score"]
            s_zscore = row["S Z score"]

            # check that the difference ratio meets the min difference threshold
            if (abs(c_zscore-s_zscore)/max([c_zscore, s_zscore])) >= min_zscore_diff:
                # check if c is more reactive
                if c_zscore > s_zscore:
                    reactive = "C"
                    c_reactivity_data.append([sequence_name, start_pos, end_pos, row["C codename"], c_zscore])
                # check if s is more reactive
                elif c_zscore < s_zscore:
                    reactive = "S"
                    s_reactivity_data.append([sequence_name, start_pos, end_pos, row["S codename"], s_zscore])
                # otherwise, tie
                else:
                    reactive = "Tie"
                    
                raw_sample_reactivity.append([fullname, reactive])

        # add reativity data to raw_epitope_df
        new_raw_sample_df = pd.DataFrame(raw_sample_reactivity, columns=[fullname_column, sample_name]).set_index(fullname_column)
        raw_epitope_df = raw_epitope_df.merge(new_raw_sample_df, how='outer', left_index=True, right_index=True)

        # create data frame for c and s reactivity
        c_epitopes = None
        s_epitopes = None
        c_reactivity_df = pd.DataFrame(c_reactivity_data, columns=["SequenceName", "StartPos", "EndPos", "Peptide", "Zscore"])
        if not c_reactivity_df.empty:
            c_epitopes, opposite_epitopes = get_epitopes(c_reactivity_df, max_distance, sample_df, "C", "S")
            for i in range(len(c_epitopes)):
                epitope = c_epitopes[i]
                sequence_name = epitope[0]
                peps = epitope[1]
                zscores = epitope[2]
                inferred_epitope = get_inferred_epitope(fasta_dict, peps.split(","))

                opposite_epitope = opposite_epitopes[i]
                opposite_peps = opposite_epitope[1]
                opposite_zscores = opposite_epitope[2]
                opposite_inferred_epitope = get_inferred_epitope(fasta_dict, opposite_peps.split(","))

                epitopes_data.append((sample_name, sequence_name, "C", inferred_epitope, peps, zscores, "S", opposite_inferred_epitope, opposite_peps, opposite_zscores))

        s_reactivity_df = pd.DataFrame(s_reactivity_data, columns=["SequenceName", "StartPos", "EndPos", "Peptide", "Zscore"])
        if not s_reactivity_df.empty:
            s_epitopes, opposite_epitopes = get_epitopes(s_reactivity_df, max_distance, sample_df, "S", "C")
            for i in range(len(s_epitopes)):
                epitope = s_epitopes[i]
                sequence_name = epitope[0]
                peps = epitope[1]
                zscores = epitope[2]
                inferred_epitope = get_inferred_epitope(fasta_dict, peps.split(","))

                opposite_epitope = opposite_epitopes[i]
                opposite_peps = opposite_epitope[1]
                opposite_zscores = opposite_epitope[2]
                opposite_inferred_epitope = get_inferred_epitope(fasta_dict, opposite_peps.split(","))

                epitopes_data.append((sample_name, sequence_name, "S", inferred_epitope, peps, zscores, "C", opposite_inferred_epitope, opposite_peps, opposite_zscores))

    formatted_epitope_df = pd.DataFrame(epitopes_data, columns=["Sample", "SequenceName", "MoreReactiveVersion", "MoreReactiveInferredEpitope", "MoreReactiveSupportingPeptides", "MoreReactiveZscores", "LessReactiveVersion", "LessReactiveInferredEpitope", "LessReactiveSupportingPeptides", "LessReactiveZscores"])
        
    return raw_epitope_df, formatted_epitope_df


# makes assumptions about column names
def get_epitopes(reactivity_df: pd.DataFrame, max_distance: int, sample_df: pd.DataFrame, index: str, opposite: str):
    data_df = sample_df.set_index(f"{index} codename")

    # group by sequence name
    sequence_groups = reactivity_df.groupby('SequenceName')
    all_epitopes = list()
    opposite_epitopes = list()

    for sequence_name, sequence_df in sequence_groups:
        # order by position
        sequence_df.sort_values(by='StartPos', ascending=True, inplace=True)

        epitope_groups = list()
        epitope_group = list()

        # create epitopes 
        for i, row in sequence_df.iterrows():
            # add first one to an epitope group
            if i != 0:
                if row["StartPos"]-prev_pep_start_pos > max_distance:
                    epitope_groups.append(tuple(epitope_group))
                    epitope_group.clear()

            epitope_group.append(i)
            
            prev_pep_start_pos = row["StartPos"]
        
        for group in epitope_groups:
            # add peptides and zscores for found epitope
            peptides = [sequence_df.loc[i]["Peptide"] for i in group]
            zscores = [str(sequence_df.loc[i]["Zscore"]) for i in group]

            all_epitopes.append((
                sequence_name,
                ",".join(peptides), 
                ",".join(zscores)
            ))

            # add peptides and zscores for opposite of found epitope
            opposite_peptides = [data_df.loc[pep][f"{opposite} codename"] for pep in peptides]
            opposite_zscores = [str(data_df.loc[pep][f"{opposite} Z score"]) for pep in peptides]

            opposite_epitopes.append((
                sequence_name,
                ",".join(opposite_peptides), 
                ",".join(opposite_zscores)
            ))
    
    return all_epitopes, opposite_epitopes


# TODO: make this generalizable
# extracts sequence name, start position, and end position from fullname
def extract_fullname_data(fullname: str):
    parts = fullname.split("_")
    return "_".join(parts[0:-2]), int(parts[-2]), int(parts[-1])


# Assumes that peptides are in order by start position
def get_inferred_epitope(fasta_dict, peptides):
    # use the first and last peptide sequences
    first_seq = fasta_dict[peptides[0]]
    last_seq = fasta_dict[peptides[-1]]

    start_pos = 0

    while not first_seq[start_pos:] == last_seq[0:len(last_seq) - start_pos]:
        start_pos += 1
    
    return first_seq[start_pos:]

#!/usr/bin/env python

import pandas as pd
import optparse
import argparse
import csv
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats
import numpy as np
import qiime2
# from q2_pepsirf.format_types import Visualization
# from q2_ps_qc.format_types import CorrelationMatrix

plt.style.use('ggplot')
plt.close("all")


def find_replicate_pair(matrix, index, base_name):
    while index < len(matrix):
        current_replicate = matrix[index].split('\t')[0]

        if base_name in current_replicate:
            return current_replicate

        index += 1

    return None


def generate_bad_corr_tsv(data, bad_corr_file_name, bad_corr_replicates):
    score_fh = open(data, "r")
    bad_corr_fh = open(bad_corr_file_name, "w")
    scores = score_fh.readlines()
    replicates = scores[0].replace("\n", "").split("\t")
    replicates.pop(0)

    bad_corr_fh.write("Sequence name\t")

    for replicate in bad_corr_replicates:
        bad_corr_fh.write(replicate)
        bad_corr_fh.write("\t")

    row_index = 1
    try:
        while True:
            # Create a list of all the scores in the current row
            # Replace all 'nan' values with '0' and remove first element as it is not a score
            row = scores[row_index].replace("\n", "").replace('nan', '0').split("\t")

            # Write the sequence name followed by a tab, then pop the sequence name from the row list
            bad_corr_fh.write(row[0])
            bad_corr_fh.write("\t")
            row.pop(0)

            # Loop through each score in the row list
            score_index = 0
            while score_index < len(row):
                # Get the name of the replicate associated with the current score
                replicate = replicates[score_index]

                # Check if the replicate matches the replicate in the predicted correlation file
                if replicate in bad_corr_replicates:
                    bad_corr_fh.write(row[score_index])
                    bad_corr_fh.write("\t")

                score_index += 1

            bad_corr_fh.write("\n")
            row_index += 1
    except EOFError and IndexError:
        pass


def generate_corr_matrix(ctx, source, data, log_normalization = False, correlation_threshold = 0.8):
    LN_CONSTANT = 11

    # Open the file with replicate scores
    score_fh = open(data, "r")
    scores = score_fh.readlines()

    # Get the names of all the replicates in the input file
    replicates = scores[0].replace("\n", "").split("\t")
    replicates.pop(0)

    repScatters_tsv = ctx.get_action('ps-plot', 'repScatters_tsv')

    first_pair_index = 0
    second_pair_index = 0
    looking_for_second_pair = False
    index = 0
    temp_index = 0
    distance_between_indices = 0
    replicate_pair_dict = {}
    bad_corr_replicates = []
    good_corr_replicates = []
    try:
        while True:
            current_replicate = replicates[index]

            if not looking_for_second_pair:
                replicate_pair_dict = {}
                base_sequence_name = current_replicate.split('_A')[0].split('_B')[0].split('_P3')[0].split('_P4')[0]

                temp_index = index
                first_pair_index = index

                replicate_pair_dict[current_replicate] = []

                looking_for_second_pair = True
                index += 1
                continue

            # Check if current replicate has matching base sequence name
            if looking_for_second_pair and current_replicate.split('_A')[0].split('_B')[0].split('_P3')[0].split('_P4')[0] \
                    == base_sequence_name:
                second_pair_index = index

                replicate_pair_dict[current_replicate] = []

                # Get row

                # Loop through each of the rows in the scores file
                row_index = 1
                try:
                    while True:
                        row = scores[row_index].replace("\n", "").replace('nan', '0').split("\t")
                        row.pop(0)

                        # Get the score from the index associated with the first pair and add it to the first pair entry
                        first_pair_score = float(row[first_pair_index])
                        replicate_pair_dict[replicates[first_pair_index]].append(first_pair_score)

                        # Get the score from the index associated with the second pair and add it to the second pair entry
                        second_pair_score = float(row[second_pair_index])
                        replicate_pair_dict[replicates[second_pair_index]].append(second_pair_score)

                        row_index += 1

                except EOFError and IndexError:
                    pass

                # Create a data frame & convert it into a correlation matrix
                data_frame = pd.DataFrame(data=replicate_pair_dict)
                corr_matrix = [data_frame.corr(method='pearson')]

                print(corr_matrix)

                for matrix in corr_matrix:
                    score_found = False
                    temp_score = '2.0'
                    for replicate in matrix:
                        if score_found and temp_score != '2.0':
                            if float(temp_score) < 0.8:
                                bad_corr_replicates.append(replicate)
                        for score in matrix.get(replicate):
                            # TODO: create a check for if all four entries in matrix are 1.0
                            if score != 1.0 and not score_found:
                                temp_score = str(score)
                                score_found = True

                                if score < 0.8:
                                    bad_corr_replicates.append(replicate)
                                elif score >= 0.8:
                                    good_corr_replicates.append(replicate)

                looking_for_second_pair = False

                if distance_between_indices > 0:
                    distance_between_indices = 0
                    index = temp_index
            else:
                distance_between_indices += 1
                if index == len(replicates) - 1:
                    print("Second pair not found; either it's already been found or there is no second pair")
                    looking_for_second_pair = False
                    distance_between_indices = 0
                    index = temp_index

            index += 1

            # Create correlation matrix

    except EOFError and IndexError:
        pass

    print(bad_corr_replicates)

    score_fh.close()

    generate_bad_corr_tsv(data, "bad_corr.tsv", bad_corr_replicates)

    correlation_vis, = repScatters_tsv(
		source = source,
		pn_filepath = None,
		plot_log = False,
		zscore_filepath = "bad_corr.tsv",
		col_sum_filepath = None,
		facet_charts = False,
		xy_threshold = None
	)

    return correlation_vis
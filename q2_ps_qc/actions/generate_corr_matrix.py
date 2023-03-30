#!/usr/bin/env python

import pandas as pd
import optparse
import argparse
import csv
import os
import inout as io

import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats
import numpy as np
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
            bad_corr_fh.write(row[0])
            bad_corr_fh.write("\t")
            row.pop(0)

            # Loop through each score in the current row
            # Skip the first element as it does not contain a score
            score_index = 0
            while score_index < len(row):
                # Get the base name of the peptide associated with the current score
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


def pep_correlation(data, log_normalization, correlation_threshold, output):
    LN_CONSTANT = 11
    with open("rep_corr_score_matrix.txt", "r") as heather_fh:
        # Open the file with replicate scores
        score_fh = open(data, "r")
        output_fh = open(output, "w")
        scores = score_fh.readlines()
        predicted_corr = heather_fh.readlines()

        # Write the column titles for the output file
        output_fh.write("Sequence name\tPredicted Correlation\tPearson\tSpearman\tKendall\n")

        # Get the names of all the peptides in the input file
        replicates = scores[0].replace("\n", "").split("\t")
        replicates.pop(0)

        index = 1
        end_of_replicate_pair = False
        found_replicates = []
        good_corr_replicates = []
        bad_corr_replicates = []
        while index < len(predicted_corr):
            # Get the current replicate
            predicted_replicate = predicted_corr[index].split('\t')[0]
            if predicted_replicate in found_replicates:
                break

            found_replicates.append(predicted_replicate)

            base_replicate_name = predicted_replicate.split('_A')[0].split('_B')[0].split('_P3')[0]
            predicted_score = predicted_corr[index].split('\t')[1].replace('\n', '\t')

            other_predicted_replicate = find_replicate_pair(predicted_corr, index + 1, base_replicate_name)
            if other_predicted_replicate is None:
                print("Warning: replicate " + predicted_replicate + " does not have a matching pair; skipping over it")
                index += 1
                continue
            found_replicates.append(other_predicted_replicate)

            # Write the current replicate & its associated predicted correlation score
            # To the output file
            output_fh.write(base_replicate_name)
            output_fh.write('\t')
            output_fh.write(predicted_score)

            # Loop through the file until it reaches its end
            row_index = 1
            corr_dict = {}
            try:
                while True:
                    # Create a list of all the scores in the current row
                    # Replace all 'nan' values with '0' and remove first element as it is not a score
                    row = scores[row_index].replace("\n", "").replace('nan', '0').split("\t")
                    row.pop(0)

                    # Loop through each score in the current row
                    # Skip the first element as it does not contain a score
                    score_index = 0
                    while score_index < len(row):
                        # Get the base name of the peptide associated with the current score
                        replicate = replicates[score_index]

                        # Check if the replicate matches the replicate in the predicted correlation file
                        if replicate == predicted_replicate or replicate == other_predicted_replicate:
                            # Add to dictionary if it hasn't already been added
                            if replicate not in corr_dict:
                                corr_dict[replicate] = []

                            if log_normalization:
                                if float(row[score_index]) >= LN_CONSTANT:
                                    LN_CONSTANT = float(row[score_index]) + 1

                                    corr_dict = {}
                                    score_index = 0
                                    print("RESETTING SCORE INDEX...")
                                    continue

                                # Add the current score to the dictionary using its associated replicate as a key
                                corr_dict[replicate].append(np.log10(float(row[score_index]) + LN_CONSTANT)
                                                            - np.log10(LN_CONSTANT))
                            else:
                                # Add the current score to the dictionary using its associated replicate as a key
                                corr_dict[replicate].append(float(row[score_index]))

                        score_index += 1

                    row_index += 1
            except EOFError and IndexError:
                pass

            # Create a data frame & convert it into a correlation matrix
            data_frame = pd.DataFrame(data=corr_dict)
            corr_matrix = [data_frame.corr(method='pearson')]

            for matrix in corr_matrix:
                score_found = False
                temp_score = '2.0'
                for replicate in matrix:
                    if score_found and temp_score != '2.0':
                        if float(temp_score) < 0.8:
                            print(temp_score, "\t", replicate)
                            bad_corr_replicates.append(replicate)
                        elif float(temp_score) >= 0.8:
                            print("This is a good score")
                    for score in matrix.get(replicate):
                        # TODO: create a check for if all four entries in matrix are 1.0
                        if score != 1.0 and not score_found:
                            temp_score = str(score)
                            output_fh.write(str(score))
                            output_fh.write("\t")
                            score_found = True

                            if score < 0.8:
                                bad_corr_replicates.append(replicate)
                            elif score >= 0.8:
                                good_corr_replicates.append(replicate)

            output_fh.write("\n")
            index += 2
            if end_of_replicate_pair:
                end_of_replicate_pair = False
            else:
                end_of_replicate_pair = True

    score_fh.close()
    output_fh.close()

    generate_bad_corr_tsv(data, "bad_corr.tsv", bad_corr_replicates)


if __name__ == "__main__":
    p = optparse.OptionParser()

    # Input controls
    p.add_option('-d', '--data',
                 help='Matrix containing data that will be used to generate correlations [none, REQ]')

    # Correlation controls
    p.add_option('-l', '--log_normalization', default=False,
                 help='Run a log normalization on each of the sets of scores before running a correlation test on them'
                      '[False]')
    p.add_option('-t', '--correlation_threshold', default=0.8,
                 help='Set a threshold value; anything below the value will be considered a bad correlation score, and'
                      'anything above will be considered a good correlation score [0.8]')

    # Output Controls
    p.add_option('-o', '--output', help='File name for output file [none, REQ]')

    opts, args = p.parse_args()

    pep_correlation(opts.data, opts.log_normalization, opts.correlation_threshold, opts.output)

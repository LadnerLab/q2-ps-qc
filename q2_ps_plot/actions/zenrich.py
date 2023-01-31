#!/usr/bin/env python
import pandas as pd
import numpy as np
import altair as alt
import glob, os, subprocess
import tempfile
import csv, itertools
from collections import defaultdict

from pandas.core.dtypes.missing import isnull

from q2_pepsirf.format_types import PepsirfContingencyTSVFormat, Zscore, PepsirfInfoSNPNFormat
import qiime2
from q2_types.feature_table import BIOMV210Format

# Name: _make_pairs_file
# Process: makes a pairs file for the purpose of inputting it
# into enrich
# Method Inputs/Parameters: column and outpath
# Method Outputs/Returned: list of pairs
# Dependencies: itertools
def _make_pairs_file(column, outpath):
    series = column.to_series()
    pairs = {k: v.index for k,v in series.groupby(series)}
    result = []
    for _, ids in pairs.items():
        result.append(list(itertools.combinations(ids, 2)))
    with open(outpath, 'w') as fh:
        for pair in result:
            for a, b in pair:
                fh.write(a + '\t' + b + '\n')
    return result

# Name: zenrich
# Process: runs the enrich module to collect enriched peptides at
# different z score thresholds, then creates an interactive graph with
# the enriched peptides highlighted and the sample available in the 
# dropdown menu
# Method inputs/parameters: output_dir, data, zscores, negative_controls,
# negative_id, highlight_probes, source, pn_filepath, peptide_metadata
# tooltip, color_by, negative_data, step_z_thresh, upper_z_thresh
# lower_z_thresh, exact_z_thresh, exact_cs_thresh, pepsirf_binary
# Dependencies: os, pandas, csv, glob, numpy, altair, subprocess, tempfile,
# and defaultdict
def zenrich(output_dir: str,
            data: PepsirfContingencyTSVFormat,
            zscores: PepsirfContingencyTSVFormat,
            negative_controls: list = None,
            negative_id: str = None,
            highlight_probes: PepsirfInfoSNPNFormat = None,
            source: qiime2.CategoricalMetadataColumn = None,
            pn_filepath: str = None,
            peptide_metadata: qiime2.Metadata = None,
            tooltip: list=['Species', 'SpeciesID'],
            color_by: str = 'z_score_threshold',
            negative_data: pd.DataFrame = None,
            step_z_thresh: int=5,
            upper_z_thresh: int=30,
            lower_z_thresh: int=5,
            exact_z_thresh: list=None,
            exact_cs_thresh: str = '20',
            pepsirf_binary: str="pepsirf") -> None:

    old = os.getcwd() # TODO: bug in framework, remove this when fixed

    #collect absolute filepath of pairs file
    if pn_filepath:
        pairsFile = os.path.abspath(pn_filepath)

    #collect absolute filepath of
    if os.path.isfile(pepsirf_binary):
        pepsirf_binary = os.path.abspath(pepsirf_binary)
        pepsirf_binary = "'%s'" % (pepsirf_binary)

    # collect probe names in a list if higlighted is provided
    peptides = []
    if highlight_probes:
        with open(str(highlight_probes), 'r') as fin:
                    for row in fin:
                        pep = row.rstrip("\n")
                        peptides.append(pep)

    #create temporary directory to work in
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)

        #create pairs file
        if source:
            pairsFile = os.path.join(tempdir, 'pairs.tsv')
            _make_pairs_file(source, pairsFile)

        #flip data frame
        data_file = str(data)
        data = data.view(pd.DataFrame)
        data = data.transpose()

        #put zscores data into dataframe and flip
        zBiom = zscores.view(BIOMV210Format)
        zData = zBiom.view(pd.DataFrame)
        zData = zData.transpose()

        #set negative matrix
        if negative_data is None:
            negative_data = data
        else:
            negative_data = negative_data.transpose()

        if negative_id and not negative_controls:
                negative_controls = []
                neg_cont = negative_data.columns
                for samp in neg_cont:
                    if samp.startswith(negative_id):
                        negative_controls.append(samp)

        #set max rows for altair to none
        alt.data_transformers.enable('default', max_rows=None)
        
        #values created by user input
        threshFile = "tempThreshFile.tsv"
        
        #set values for pepsirf
        outSuffix = "_tempEnriched.txt"
        failOut = "enrichFail.txt"
        outputDir = []
        
        #convert colsum data to pandas data frame
        peptideNames = negative_data.index
        
        #get list of thresholds
        if not exact_z_thresh:
            threshRange = list(reversed(range(lower_z_thresh, upper_z_thresh, step_z_thresh)))
        else:
            threshRange = exact_z_thresh

        #iterate through thresholds and generate threshold file to be used to get enriched peptide
        for score in threshRange:
            outputDir.append(str(score))
            #write the threshold file for specified zscore
            with open(threshFile, 'w', newline='') as out_file:
                    tsv_writer = csv.writer(out_file, delimiter='\t')
                    tsv_writer.writerow([str(zscores), score])
                    tsv_writer.writerow([str(data_file), str(exact_cs_thresh)])
            
            #run p enrich module
            cmd = "%s enrich -t %s -s '%s' -x %s -o %s -f %s >> enrich.out" % (pepsirf_binary, threshFile, pairsFile, outSuffix, str(score), failOut)
            subprocess.run(cmd, shell=True)

        #create empty dicionaries to collect all sample names and sample/peptide combos
        dropNames = {}
        pepSamp = {}

        #create empty dictionary of lists to collect enriched peptides
        enrDict = defaultdict(list)
        columnNms = ['Peptide','sample','sample_value', 'negative_control', 'z_score_threshold', 'Zscores']
        for nm in columnNms:
            enrDict[nm]

        #iterate through created directories
        for oD in outputDir:
            enrFiles = glob.glob("%s/*%s" % (oD, outSuffix))
            
            #iterate through enriched files to gather the peptides
            for file in enrFiles:
                
                #get individual sample names
                sNames = os.path.basename(file).split(outSuffix)[0].split("~")
                
                #get joined sample name
                sample = "~".join(sNames)

                #generate a list of sample names
                if sample not in dropNames:
                    dropNames[sample] = ""
                
                #collect information for the pandas data frame
                with open(file, 'r') as fin:
                    for row in fin:
                        pep = row.rstrip()
                        if pep:
                            if (pep,sample) not in pepSamp:
                                x = np.mean([float(negative_data[sn][pep]) for sn in negative_controls])
                                y = np.mean([float(data[sn][pep]) for sn in sNames])
                                z = [zData[sn][pep] for sn in sNames]
                                zToStr = ', '.join([str(elem) for elem in z])

                                #add all enriched peptide info to enriched dictionary
                                enrDict['Peptide'].append(pep)
                                enrDict['sample'].append(sample)
                                enrDict['sample_value'].append(np.log10(y+1))
                                enrDict['negative_control'].append(np.log10(x+1))
                                enrDict['z_score_threshold'].append(oD)
                                enrDict['Zscores'].append(zToStr)
                                pepSamp[(pep,sample)] = ""
            
            #check if there were any samples that had no enriched peptides
            if os.path.exists(os.path.join(oD, failOut)): 

                #add sample names that have no enriched peptides to dropNames dict
                with open(os.path.join(oD, failOut), 'r') as fh:
                    lines = fh.readlines()
                    c = 0
                    for ln in lines:
                        fail = ln.rstrip("\n").split("\t")
                        samp = "~".join(fail[0].split(", "))
                        if c > 0 and fail[1] == 'No enriched peptides' and samp not in dropNames:
                            dropNames[samp] = ''
                        c += 1

        #convert eriched dictionary to enriched data frame
        enrichedDf = pd.DataFrame(enrDict)

        #collect peptide metadata information and add it to the enriched data frame
        if peptide_metadata:
            metaDf = peptide_metadata.to_dataframe()
            enrichedDf = metaDf.merge(enrichedDf, how="right", right_on="Peptide", left_index=True)

            #add pepsirf to tooltip list and add data type to each item
            tooltip.insert(0, 'Peptide')
            tooltip.insert(1, 'Zscores')
            i = 0
            for title in tooltip:
                tooltip[i] += ':N'
                i += 1
        else:
            tooltip = ['Peptide:N', 'Zscores:N']

        # if highlighted probes provided, collect all of the peptides and put into new df
        if highlight_probes:
            highlightDf = enrichedDf.loc[enrichedDf["Peptide"].isin(peptides)]

        #create empty directory to collect all heatmap information
        hmDic = defaultdict(list)
        columnNms = ['sample', 'bin_x_start', 'bin_x_end', 'bin_y_start', 'bin_y_end', 'count']
        for nm in columnNms:
            hmDic[nm]

        #collect x axis data
        x = np.array([np.mean([float(negative_data[sn][pn]) for sn in negative_controls if not pd.isnull(negative_data[sn][pn])]) for pn in peptideNames])
        xLog = np.log10(x+1)

        #iterate through samples to fill Data Frame
        for sample in dropNames.keys():
            sNames = sample.split('~')
            y = np.array([np.mean([float(data[sn][pn]) for sn in sNames if not pd.isnull(data[sn][pn])]) for pn in peptideNames])
            yLog = np.log10(y+1)
            heatmap, xedges, yedges = np.histogram2d(xLog, yLog, bins=(70,70))

            #fill binned data frame for heatmap generation by row
            for x in range(0, heatmap.shape[0]):
                for y in range(0, heatmap.shape[1]):
                    count = heatmap[x,y]
                    #do not add data with count of 0
                    if count == 0.0:
                        continue
                    bin_x_start = xedges[x]
                    bin_x_end = xedges[x+1]
                    bin_y_start = yedges[y]
                    bin_y_end = yedges[y+1]

                    #add all x and y bins to heatmap dictionary
                    hmDic['sample'].append(sample)
                    hmDic['bin_x_start'].append(bin_x_start)
                    hmDic['bin_x_end'].append(bin_x_end)
                    hmDic['bin_y_start'].append(bin_y_start)
                    hmDic['bin_y_end'].append(bin_y_end)
                    hmDic['count'].append(count)

        #convert heatmap dictionary to data frame
        heatmapDf = pd.DataFrame(hmDic)

        #find axis ratio for chart width an height
        xy_max = heatmapDf[['bin_x_end', 'bin_y_end']].max()
        ratio = (xy_max[0] / xy_max[1]) - 1 
        chartHeight = 500
        chartWidth =  chartHeight + (20 * ratio)

        #create a sorted list for dropNames
        dropList = list(dropNames.keys())
        dropList.sort()

        #create dropdown specs
        sample_dropdown = alt.binding_select(options=dropList, name='Sample Select')
        sample_select = alt.selection_single(fields=['sample'], bind=sample_dropdown, name="sample", init={'sample': dropList[0]})
        
        # set color by as nominal
        color_by = color_by + ":N"

        #create scatterplot of enriched peptides
        scatter = alt.Chart(enrichedDf).mark_circle(size=50).encode(
            x = alt.X('negative_control:Q', title = "Negative Control log10(value+1.0)"),
            y = alt.Y('sample_value:Q', title = "Sample log10(value+1.0)"),
            color = alt.Color(color_by,
                scale=alt.Scale(range=['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']),
                sort=threshRange,
                legend=alt.Legend(title='Z Score Thresholds')),
            tooltip = tooltip
        ).add_selection(
            sample_select
        ).transform_filter(
            sample_select
        )

        #create heatmap of all peptides
        heatmapChart = alt.Chart(heatmapDf, width=chartWidth, height=chartHeight).mark_rect().encode(
            alt.X('bin_x_start:Q', scale=alt.Scale(domain=(0,xy_max[0]))),
            alt.X2('bin_x_end:Q'),
            alt.Y('bin_y_start:Q', scale=alt.Scale(domain=(0,xy_max[1]))),
            alt.Y2('bin_y_end:Q'),
            alt.Color('count:Q', scale = alt.Scale(scheme='greys'), legend=alt.Legend(title='Point Frequency'))
        ).transform_filter(
            sample_select
        )

        #layer scatterplot and heatmap
        finalChart = alt.layer(heatmapChart, scatter).properties(title="Z Score Threshold Variance")

        # if highlighted probes provided, create square scatterplot
        if highlight_probes:
            highlightedChart = alt.Chart(highlightDf).mark_square(size=60).encode(
                x = alt.X('negative_control:Q', title = "Negative Control log10(value+1.0)"),
                y = alt.Y('sample_value:Q', title = "Sample log10(value+1.0)"),
                color = alt.Color(color_by,
                    scale=alt.Scale(range=['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']),
                    sort=threshRange,
                    legend=alt.Legend(title='Z Score Thresholds')),
                tooltip = tooltip
            ).transform_filter(
                sample_select
            )

            #layer the square scatterplot over the final chart
            finalCharts = alt.layer(finalChart, highlightedChart)

            #save the plot
            finalCharts.save(os.path.join(output_dir, "index.html"))

        #otherwise, save the final chart without the higlighted probes
        else:
            finalChart.save(os.path.join(output_dir, "index.html"))

    os.chdir(old)  # TODO: found a bug in the framework, remove when fixed
    
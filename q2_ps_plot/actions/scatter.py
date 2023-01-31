from altair.vegalite.v4.schema.channels import Tooltip
import pandas as pd
import altair as alt
import numpy as np
from collections import defaultdict
import qiime2, itertools, os, csv
import random
import math
from q2_pepsirf.format_types import MutantReferenceFileFmt

# Name: _make_pairs_list
# Process: creates a list of sample replicates
# Method inputs/parameters: column
# Method outputs/Returned: list of pairs
# Dependencies: itertools
# function to collect pairwise pairs of samples
def _make_pairs_list(column, type):
    if type == "sourceCol":
        series = column.to_series()
        pairs = {k: v.index for k,v in series.groupby(series)}
        result = []
        for _, ids in pairs.items():
            combination = list(itertools.combinations(ids, 2))
            if len(combination) > 0:
                result.append(tuple(combination))
            else:
                result.append(ids)
    elif type == "pnFile":
        result = []
        with open(column) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter="\t")
            for row in csv_reader:
                result.append(tuple([tuple(row)]))
    return result

# Name: repScatters
# Process: creates an interactive scatterplot/heatmap for reps of 
# Col-sum data or reps of z score data
# Method inputs/parameters: output_dir, source, plot_log, zscore, col_sum,
# facet_charts
# Dependencies: os, pandas, numpy, altair, and defaultdict
def repScatters(
    output_dir: str,
    source: qiime2.CategoricalMetadataColumn = None,
    pn_filepath: str = None,
    plot_log: bool = False,
    zscore: pd.DataFrame = None,
    col_sum: pd.DataFrame = None,
    facet_charts: bool = False,
    xy_threshold: int = None)->None:

    # check wether zscore matrix or colsum matrix was provided
    if zscore is not None:
        data = zscore.transpose()
    if col_sum is not None:
        data = col_sum.transpose()

    # start the default dict and add the following column names or keys
    hmDic = defaultdict(list)
    columnNms = ['sample', 'bin_x_start', 'bin_x_end', 'bin_y_start', 'bin_y_end', 'count']
    for nm in columnNms:
        hmDic[nm]

    # start the default dict for collecting scatterplot values and add the following keys
    scatterDict = defaultdict(list)
    scatter_columnNms = ['sample', 'x_value', 'y_value']
    for nm in scatter_columnNms:
        scatterDict[nm]

    # call function to collect samples pairs
    if source:
        pairs = _make_pairs_list(source, "sourceCol")
    elif pn_filepath:
        pairs = _make_pairs_list(pn_filepath, "pnFile")

    # initialize empty list to collet sample names for dropdown menu
    samples = []

    # start loop to collect heatmap values
    for lst in pairs:
        for tpl in lst:

            #set x and y values from the dataframe
            if type(tpl) is tuple:
                samp = '~'.join(tpl)
                samples.append(samp)
                # x = [ x for x in data[tpl[0]] if not np.isnan(x)]
                # y = [ y for y in data[tpl[1]] if not np.isnan(y)]
                x = []
                y = []
                xData = list(data[tpl[0]])
                yData = list(data[tpl[1]])
                peptides = list(data.index)
                for i in range(len(xData)):
                    if not np.isnan(xData[i]) and not np.isnan(yData[i]):
                        x.append(xData[i])
                        y.append(yData[i])

                        #if xy-threshold provided collect values matching threshold
                        if xy_threshold and xData[i] > xy_threshold and yData[i] > xy_threshold:
                            scatterDict["x_value"].append(xData[i])
                            scatterDict["y_value"].append(yData[i])
                            scatterDict['peptide'].append(peptides[i])
                            scatterDict['sample'].append(samp)
            else:
                samp = tpl
                samples.append(samp)
                x = [ x for x in data[tpl] if not np.isnan(x)]
                y = [ y for y in data[tpl] if not np.isnan(y)]


            if plot_log:
                x = np.log10(np.array(x)+1)
                y = np.log10(np.array(y)+1)
                xTitle = "Sample1 log10(score + 1)"
                yTitle = "Sample2 log10(score + 1)"
            else:
                xTitle = "Sample1"
                yTitle = "Sample2"

            #generate heatmap values
            heatmap, xedges, yedges = np.histogram2d(x, y, bins=(70,70))
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

                    #place heatmap values in the dictionary created
                    hmDic['sample'].append(samp)
                    hmDic['bin_x_start'].append(bin_x_start)
                    hmDic['bin_x_end'].append(bin_x_end)
                    hmDic['bin_y_start'].append(bin_y_start)
                    hmDic['bin_y_end'].append(bin_y_end)
                    hmDic['count'].append(count)
    
    #convert the heatmap to a dataframe
    hmDf = pd.DataFrame(hmDic)

    #convert the scatterplot dict to a dataframe
    scatterDf = pd.DataFrame(scatterDict)

    # set the dropdown values
    sample_dropdown = alt.binding_select(options=samples, name='Sample Select')
    sample_select = alt.selection_single(fields=['sample'], bind=sample_dropdown, name="sample", init={'sample': samples[0]})

    if not facet_charts:
        # create scatterplot chart with attached dropdown menu
        heatmapChart = alt.Chart(hmDf).mark_rect().encode(
                alt.X('bin_x_start:Q', title=xTitle),
                alt.X2('bin_x_end:Q'),
                alt.Y('bin_y_start:Q', title=yTitle),
                alt.Y2('bin_y_end:Q'),
                alt.Color('count:Q', scale = alt.Scale(scheme='plasma'))
        ).add_selection(
            sample_select
        ).transform_filter(
            sample_select
        )
    else:
        # create scatterplot chart with attached dropdown menu
        heatmapChart = alt.Chart(hmDf).mark_rect().encode(
                alt.X('bin_x_start:Q', title=xTitle),
                alt.X2('bin_x_end:Q'),
                alt.Y('bin_y_start:Q', title=yTitle),
                alt.Y2('bin_y_end:Q'),
                alt.Color('count:Q', scale = alt.Scale(scheme='plasma'))
        ).facet(
            facet='sample:N',
            columns=5
        )
    
    #if xy-threshold then create the scatterplot
    if xy_threshold:
        Scatter = alt.Chart(scatterDf).mark_circle(opacity=1, size=100, color="#009E73").encode(
            x = alt.X("x_value:Q", title = xTitle),
            y = alt.Y('y_value:Q', title = yTitle),
            tooltip = ['sample', 'peptide', 'x_value', 'y_value']
        ).transform_filter(
            sample_select
        )

        chart = alt.layer(heatmapChart, Scatter)
        chart.save(os.path.join(output_dir, "index.html"))

    if not xy_threshold:
        # save the chart to index.html to be saved to a .qzv file
        heatmapChart.save(os.path.join(output_dir, "index.html"))

# Name: mutantScatters
# Process: creates an interactive scatterplot/boxplot for mutant peptides
# Method inputs/parameters: output_dir, source, metadata, zscore, reference_file,
# peptide_header, reference_header, x_axis_header, category_header, label_header
# x_axis_label, y_axis_label, min_wobble, max_wobble, wobble, scatter_only,
# scatter_boxplot, boxplot_only
# Dependencies: math, pandas, altair, and random
def mutantScatters(
    output_dir: str,
    source: qiime2.CategoricalMetadataColumn,
    metadata: pd.DataFrame,
    zscore: pd.DataFrame,
    reference_file: MutantReferenceFileFmt = None,
    peptide_header: str = "FeatureID",
    reference_header: str = "Reference",
    x_axis_header: str = 'Position',
    category_header: str = 'Category',
    label_header: str = 'Label',
    x_axis_label: str = 'Position',
    y_axis_label: str = 'Mean Zscore',
    min_wobble: float = -0.4,
    max_wobble: float = 0,
    wobble: bool = False,
    scatter_only: bool = True,
    scatter_boxplot: bool = False,
    boxplot_only: bool = False
)->None:

    # convert metadata to a dataframe and reset the index
    metadata = metadata.to_dataframe().reset_index()

    # reverse the zscore dataframe
    zscore = zscore.transpose()

    #Collect replicate names and source names and collect mean values
    series = source.to_series()
    pairs = {k: v.index for k,v in series.groupby(series)}
    delete = []
    for key, lst in pairs.items():
        zscore[key] = zscore[list(lst.values)].mean(axis=1)
        for item in lst:
            delete.append(item)
    
    #delete unused columns
    zscore = zscore.drop(columns = delete)

    #reorganized the dataframe
    zscore = zscore.stack().to_frame().reset_index()
    zscore = zscore.rename(columns={'Sequence name': peptide_header, 'level_1': 'Samples', 0 : 'value'})

    #merge zscore and metadata and create a copy of the dataframe
    df = pd.merge(zscore, metadata, how="right", on=peptide_header)
    dfCopy = df.copy()

    if reference_file:
        #collect information for labeling
        References = []
        labeling = []
        CodeName = []
        Position = []
        with open(str(reference_file)) as file:
            tsv_file = csv.reader(file, delimiter="\t")
            c = 0
            for line in tsv_file:
                c+=1
                if c != 1:
                    References.append(line[0])
                    ct = 0
                    for item in line[1]:
                        ct+=1
                        Position.append(ct)
                        labeling.append(item)
                        CodeName.append(line[0])
                
        #create pandas dataframe with collected info
        labelDict = {reference_header: CodeName, x_axis_header: Position, 'labeling': labeling}
        label = pd.DataFrame(labelDict)

    else:
        References = list(df['Reference'].unique())
    #create reference dataframe
    refernceDict = {reference_header: References, peptide_header: References}
    ReferenceDf = pd.DataFrame(refernceDict)
    ReferenceDf = pd.merge(zscore, ReferenceDf, how="right", on=peptide_header)

    #if wobble selected wobble positions
    if wobble:

        df[x_axis_header] = df[x_axis_header].astype(float)
        for index, row in df.iterrows():
            if not math.isnan(df.at[index, x_axis_header]):
                df.at[index, x_axis_header + "_value"] = float(df.at[index, x_axis_header]) + random.uniform(min_wobble, max_wobble)

    # collect sample and reference names
    samples = list(df[reference_header].unique())
    samp = list(df['Samples'].unique())

    # set the dropdown values
    sample_dropdown = alt.binding_select(options=samples, name='Reference Select')
    sample_select = alt.selection_single(fields=[reference_header], bind=sample_dropdown, name="Reference", init={reference_header: samples[0]})

    # set the dropdown values
    samp_dropdown = alt.binding_select(options=samp, name='Sample Select')
    samp_select = alt.selection_single(fields=['Samples'], bind=samp_dropdown, name="Samples", init={'Samples': samp[0]})

    #set max rows for altair to none
    alt.data_transformers.enable('default', max_rows=None)

    # if not boxplot only selected create scatterplot
    if not boxplot_only:
        if wobble:
            chart = alt.Chart(df).mark_circle(size=50).encode(
                x = alt.X(x_axis_header + "_value:Q", title = x_axis_label),
                y = alt.Y('value:Q', title = y_axis_label),
                color = alt.Color(category_header + ":N",
                scale=alt.Scale(range=['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7'])),
                tooltip = [label_header, peptide_header, category_header]
            ).add_selection(
                sample_select,
                samp_select
            ).transform_filter(
                sample_select &
                samp_select
            ).properties(
                width=1000,
                height=300
            )
        else:
            chart = alt.Chart(df).mark_circle(size=50).encode(
            x = alt.X(x_axis_header + ":Q", title = x_axis_label),
            y = alt.Y('value:Q', title = y_axis_label),
            color = alt.Color(category_header + ":N",
            scale=alt.Scale(range=['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7'])),
            tooltip = [label_header, peptide_header, category_header]
        ).add_selection(
            sample_select,
            samp_select
        ).transform_filter(
            sample_select &
            samp_select
        ).properties(
            width=1000,
            height=300
        )

        #create dashed reference line
        line = alt.Chart(ReferenceDf).mark_rule(
            strokeDash=[4,5]
        ).encode(y='value:Q').transform_filter(
            sample_select &
            samp_select
        )

        #layer the line and scatterplot
        chart1 = alt.layer(line, chart)

        if reference_file:
            #add the text chart underneath the scatter plot
            text = alt.Chart(label).mark_text(
                align='center'
            ).encode(
                x = alt.X(x_axis_header + ':Q', title=x_axis_label),
                text='labeling:N'
            ).properties(
                width=1000
            ).add_selection(
                sample_select
            ).transform_filter(
                sample_select
            )

            # layer the text chart underneath the scatterplot
            chart2 = alt.vconcat(chart1, text)

        # if only scatter, save the plot
        if scatter_only and reference_file:
            chart2.save(os.path.join(output_dir, "index.html"))
        elif scatter_only:
            chart1.save(os.path.join(output_dir, "index.html"))

    # if scatter_boxplot is chosen create the boxplot with the dropdown menu
    if scatter_boxplot:
        boxplot = alt.Chart(dfCopy).mark_boxplot().encode(
            x=alt.X(x_axis_header + ':O', title=x_axis_label),
            y=alt.Y('value:Q', title=y_axis_label)
        ).properties(
            width=1000,
            height=300
        ).add_selection(
            sample_select,
            samp_select
        ).transform_filter(
            sample_select &
            samp_select
        )

        if reference_file:
            #layer the chart and save
            layeredChart = alt.vconcat(chart2, boxplot)
            layeredChart.save(os.path.join(output_dir, "index.html"))
        else:
            #layer the chart and save
            layeredChart = alt.vconcat(chart1, boxplot)
            layeredChart.save(os.path.join(output_dir, "index.html"))

    # if boxplot only create a boxplot chart
    if boxplot_only:
        boxplot = alt.Chart(dfCopy).mark_boxplot().encode(
            x=alt.X('Position:O', title=x_axis_label),
            y=alt.Y('value:Q', title=y_axis_label),
            column = 'Reference:N',
            row = 'Samples:N'
        )

        # save the boxplot chart
        boxplot.save(os.path.join(output_dir, "index.html"))

---
description: >-
  Open Science Data Repository (OSDR) “Data Visualization” portal that allows
  users to interact with the data from space-related studies within GeneLab's
  processed database.
---

# OSDR Data Visualization Portal

[Link to Data Visualization Portal](https://visualization.genelab.nasa.gov/data/) ([https://visualization.genelab.nasa.gov/data/](https://visualization.genelab.nasa.gov/data/))

```{button-link} https://visualization.genelab.nasa.gov/data/
:color: primary
:outline:
Open the Data Visualization Portal
```

The portal encompasses various visualization types, including Gene Expression query tables, Dendrograms, Heatmaps, Ideogram, Gene Set Enrichment Analysis (GSEA) and interactive plots including PCA, Pair, and Volcano plots. Each tool offers the flexibility to adjust parameters and explore specific aspects of the data.

* The OSDR has many data visualization tools suitable for multiple applications and stakeholders.&#x20;
  * **Data & Tools: Use the dropdown menu to access the** “[Data Visualization](https://visualization.genelab.nasa.gov/data/)” hyperlink to open the OSDR multi-study data visualization application.

[Link to detailed Multi-study visualization tutorial (PDF)](https://web.archive.org/web/20240615223406/https://genelab.nasa.gov/sites/default/files/2024-03/Data%20Visualization%20Portal%20JT_AU.pdf) — *Data Visualization Portal Tutorial*, J. Davis and A. Uriarte-Acuna. The original genelab.nasa.gov copy is no longer available (it now redirects to a NASA landing page), so this links to the Internet Archive's copy from June 2024.

### **Metadata Dashboard** <a href="#id-18i9g98tw37j" id="id-18i9g98tw37j"></a>

[Link to Metadata Dashboard](https://visualization.genelab.nasa.gov/data/)

```{button-link} https://visualization.genelab.nasa.gov/data/
:color: primary
:outline:
Open the Metadata Dashboard
```

**The metadata dashboard is designed to help users narrow search results for experimental data. It provides various tools for filtering and displaying results.**

*   The main tools for filtering the study table's results are the Pie charts and the filters on the left side of the dashboard. Each section of the Pie chart acts as a separate section of filters, and when a filter from the Pie chart is selected the results containing that factor will automatically populate in the studies table below. A user can make one selection on each Pi chart to narrow the results in the studies table. On the left side of the dashboard are a series of filtering options that can be used to identify similar studies. At the bottom of the page is a table summarizing the main metadata that describes how these studies were conducted.

    * At the bottom left of the table, there is a blue “Visualize Study” button that will launch any combination of studies that can be selected into a new tab in your browser for further analysis

### &#x20;<a href="#uk50q5k2hiyu" id="uk50q5k2hiyu"></a>

### Selecting “processed studies” to find higher-order data. <a href="#uk50q5k2hiyu" id="uk50q5k2hiyu"></a>

* Users can select factors as filters on the Pie charts or the left side of the dashboard, both sections will be updated to show the selected filters and the studies table will be updated to show the relevant studies.
  * At the top of the filter panel on the left of the screen, you can see the “Show only studies with processed data tab”, selecting this option and then pressing the “Apply” filter button that excludes studies without processed data.

### &#x20;<a href="#e5g8py1a7trn" id="e5g8py1a7trn"></a>

### Bar charts of factors <a href="#e5g8py1a7trn" id="e5g8py1a7trn"></a>

**The bar graph displays the results per factor.** When a factor is selected, it is broken down into types. After selecting a filter, the user can click "apply" to activate the filter and update the results on the graph. A red "X" appears next to the crosshair of the Pi chart to clear filters.

* Each pie chart comes with a crosshair tool located at the bottom left. Selecting the crosshair displays a bar graph showing the factors listed in the pie chart.
  * **Example factors bar plot:** In this example, the tissue factor has been viewed as a bar plot, and then the “Root” is selected, this then opens a new table that includes the 2 tissue ontology terms that contributed to the original “Root” ontology term. In this example “Root” is made up of both “plant root” and “root”

### &#x20;<a href="#baeehoyiz78k" id="baeehoyiz78k"></a>

### Selecting a study or group of studies for visualization <a href="#baeehoyiz78k" id="baeehoyiz78k"></a>

**Once a user has selected a study or multiple studies they can press the "Visualize Study" button to be directed to the data visualization tools.**

* You may also select multiple studies to visualize simultaneously in which case a user will be directed to a Multi-Study preview page before being directed to the data visualization tools.
  * **Filtered down to one study example:** A series of filters applied narrowed the selection to 1 study that fulfils all the user's parameters. The red dashed lines highlight the box that allows the user to select a study, the screen image also shows how a selected study appears when it is highlighted and also identifies where the blue “Visualize Study” button can be found.

### &#x20;<a href="#s7aljci4w0hb" id="s7aljci4w0hb"></a>

### Studies Metadata Menu as a Table <a href="#s7aljci4w0hb" id="s7aljci4w0hb"></a>

**At the bottom of the page below the pie charts is a table that lists the studies resulting from the selected filters from above.**

*   The table includes the following information for each study: OSD, Title, Assay, Organism, Tissue, and Factor. By default, the studies will be listed in order of OSD-# from smallest to largest, but the order can be reversed based on each information category by double-clicking on the study title.

    * **First 10 in chronological order:** If no list is applied then the table default is to show these in numerical order based on the OSD accession number.

### Selecting studies for Multi-study analysis: <a href="#adcagqrqpdky" id="adcagqrqpdky"></a>

The multi-study page is used to initialize the parameters for data visualization of the multiple studies. Researchers can uncover intricate patterns of gene expression associated with specific conditions or treatments across a variety of experiments. Below are detailed instructions on how to effectively navigate and utilize the Multi-Study Page.

* For your initial test, let's use Arabidopsis thaliana (plant) studies as an example.
* Start by selecting "Arabidopsis thaliana" as the organism of interest. Since combining DNA microarray assays is not supported, ensure filtering by both "Arabidopsis thaliana" and "RNA sequencing" in the assay technology type.
* Choose two different Arabidopsis studies that encompass various organ/tissue types. For instance, select "OSD-120" and "OSD-217."
* Mark the checkboxes beside the selected studies in the studies table.
* Click the "Visualize Study" button to proceed.
* In this example OSD-120 & OSD-217 have been selected, normalized with DESeq2 and their first 2 principle components plotted in a 2D scatter plot. To the right the is a factor selection tab where you can add new factors that will appear as columns beside the OSD-###. Clicking on the “Expand table” button reveals a table showing all the replicates and the metadata the user loads.

### &#x20;<a href="#n1ennvxv3ysp" id="n1ennvxv3ysp"></a>

### Multi-study Data Normalization: <a href="#n1ennvxv3ysp" id="n1ennvxv3ysp"></a>

**A dialogue box will appear to prompt you for data normalization options.** The default selection is often "DESeq2" for normalization, but you can also choose "No Normalization."

* **E-mail notification:** If desired, you can enter your email address to receive a notification when the studies have been combined and normalized. Alternatively, proceed without entering an email address.
* **Normalization and PCA Insights:** Understand normalization details by clicking the Information button next to the normalization method. View the PCA chart, which provides insights into data distribution after normalization.
* **Factor Selection and Differential Gene Analysis:** Under "Factor Selection," choose variables to generate a factors table for differential gene analysis. Select parameters, characteristics, or factors from the dropdown list to add to the table.
* **Exploring the Multi-Study Page:** A PCA chart for data visualization will be included on the multi-study page. Utilize the PCA chart options to tailor your visualization based on specific criteria.
* **Modifying Normalization Method:** If you wish to change the normalization method (e.g., from "DESeq2" to "No Normalization"), click "Change Normalization Method."
* **Sample Selection for Gene Expression Analysis:** Select specific samples by clicking the "Select labelled, expand table" button. Choose samples based on factors added during factor selection.
* **Visualizing and Downloading Results:** Click "Visualize Studies" to proceed to visualization plots or download the accounts table. Depending on your selection, enter your email address for a notification upon completion.
* **Exploring Visualization Plots:** Upon completion, the page will direct you to a range of visualization plots and graphs for your data analysis.

### &#x20;<a href="#oek6t1v60xxg" id="oek6t1v60xxg"></a>

### Multi-study Sidebar Functions <a href="#oek6t1v60xxg" id="oek6t1v60xxg"></a>

**When a user has selected the study/studies to visualize, they will be directed to the data visualization tools, where a sidebar of helpful tools is provided on the left side of the screen.**

The "Study Details" button is located at the top of the sidebar. This button pulls up a display with the study information including a small description. The display also includes a tab labelled “samples” that a user can press to see the individual samples and additional information for the study.

Below the study details button is a label for each plot provided for a user within the data visualization tool. Clicking these labels will automatically direct the user to the plot associated with the label.

At the bottom of the sidebar is the default Group selection is utilized for each plot. A user can modify the groups are selected by pressing the "Modify Groups" button. This button will prompt the user to select the individual groups that a user would like to see displayed on each plot.

A feature exclusive to multi-study visualization is the option to download the combined Differential Gene Expression (DGE) table.

Users can then select a threshold based on quantitative factors of gene expression such as fold change, P-value and Adjusted P-value.

When a user accesses the multi-study visualization the DGE table will have several options to export the information at the top of the table.

### Plotly <a href="#id-6klo26n9so4x" id="id-6klo26n9so4x"></a>

**Plotly is a third-party software that uses data provided by OSDR to create the interactive visualizations displayed.**

* At the top-right corner of each plot will be options to help a user better visualize the data. The house logo within the options will reset the axes of the plots back to default.
  * Users also have the option to zoom in/out on each plot as well as auto-scale the graphic. There are two tools provided for data point selection, the lasso tool and the box select tool. Each of these tools provides a shape that will select any data points that fall within them. Lastly, there is a download button in the shape of a camera that will let you download the plot as a PNG file.

#### &#x20;<a href="#pf6ggq2bhe0b" id="pf6ggq2bhe0b"></a>

### PCA Plots <a href="#zgf6azkbfb8w" id="zgf6azkbfb8w"></a>

**PCA stands for Principal Component Analysis, and this type of plot is used to reduce the dimensionality of large experiments to simplify the process of analyzing the data points.**

[Link to read about PCA plots](https://builtin.com/data-science/step-step-explanation-principal-component-analysis)

* Each PCA plot will include options for a 2D and 3D representation of the data. The default selection is a 3D representation on an "X", "Y", and "Z" axis. In the upper left corner of the plot area select the "2D" button and then press "Update" The graph will update to display the data on an "X", and "Y" axis only.
  * The "Color by Factor" feature allows users to select a specific factor from the study for representation on the graph to allow for an easier comparison between differences in the data. Select the "Color by Factor" drop-down menu, within the drop-down menu select one factor, then press the "Update" button.

### **Selecting groups in PCA plots**

**PCR plots contain a lot of information, on some occasions selecting a fraction of the samples and plotting their PCA variations in 2D can often be a clear and simple method to highlight clustering.**

* In this example, the "Cell Line" factor was selected from the drop-down. The results will now be represented by colors matching the factor selected. In this example, (OSD-154) the colors represent the different cell lines from the experiment and clearly show how the cell lines could be a factor in the differences between the data points.
  * Another feature within the PCA plot tool allows users to hide factors by selecting the label located on the right side of the plot. The two labels provided are the cell lines "GM15036" which is represented by blue, and "GM15510" represented by the color orange. Click on the label "GM15036" and the data points will be hidden as shown below. Click on the label "GM15036" a second time and the data will return to the graph.

### Pair Plots <a href="#a8pdjzjwz1gr" id="a8pdjzjwz1gr"></a>

[Link to read about Pair Plots](https://medium.com/analytics-vidhya/pairplot-visualization-16325cd725e6)

Pair plots are used for Exploratory Data Analysis, where the plot visualizes the data in order to find a relationship between variables that can be continuous or categorical.

* A “Pair Plot” is used to understand the best set of features to explain a relationship between two variables. It also helps to form simple classification models by making linear separations in a dataset.
  * **Enhancing difference visualization:** The default display for the pair plot will be the comparison between two sets of data with a % difference color threshold of 20%. Two plots will be displayed on the dashboard for the ability to compare multiple sets of data simultaneously. In this example on the right, we can see an increase in the threshold reduced the number loci annotated as red. Clicking each of the drop-down menus will allow a user to change which axis the sample data is displayed on.

**Refining Pair plot correlation visualization**

Upon selecting the green "Samples" button located on the plot, one is presented with the opportunity to adjust several correlation coefficient values. These values serve to quantify both the strength and direction of the correlation that exists between data points.

* Coefficients closer to +1 indicate a strong positive correlation, while those near -1 suggest a strong negative correlation. A coefficient near 0 signifies a weak or no correlation. Selecting a different coefficient dynamically alters the spread of the data, with higher coefficients leading to tighter clustering and lower ones resulting in broader scattering.
* Users also can view different data correlations by clicking the green "Samples" button at the top of the plot and comparing it to other study providing a correlation coefficient summary. Clicking this button will change the dropdown to show multiple correlation coefficients for a pair of samples from the OSDR.

### Volcano Plots <a href="#v6dx63l0d9d" id="v6dx63l0d9d"></a>

[Link to read about Volcano Plots](https://www.htgmolecular.com/blog/2022-08-25/understanding-volcano-plots)

**The name volcano plot comes from its resemblance to a volcanic eruption with the most significant points at the top, like spewed pieces of molten lava.**

* A volcano plot is useful for identifying events that differ significantly between two groups of experimental subjects. The default display for Volcano Plots will have the -Log10(Adj P Value) with an Adj P Value threshold of 0.05 and a Log2 FC threshold of 1.00 as shown below.
  * **Each point on the graph represents a gene:** The log2-fold differences between the groups are plotted on the x-axis and the -log10 p-value differences are plotted on the y-axis. The horizontal dashed line represents the significance threshold specified in the analysis, usually derived using a multiple testing correction.

#### Customizing volcano plots <a href="#xqfjx7fvgehs" id="xqfjx7fvgehs"></a>

**Adjusting the cutoff threshold for coloring the volcano plot and adjusting the distribution of data can increase the ease of data integration.**

* Users can change the type of data displayed on the Y axis, and the options from the dropdown menu include "P Value, Adjusted P Value, and -Log10(P Value)".
  * Below is an example of the "P-value", adjusted P value and log10 (P-values) display for a volcano plot. In addition, the ability to change the P value threshold is available and the image below shows a P value and Log2 FC thresholds are also provided. In this example the threshold is increased to P values of 0.25 and a Log2 Fold Change threshold of 3 is applied.

### Heatmaps <a href="#id-90l4f2mpkzz5" id="id-90l4f2mpkzz5"></a>

[Link to read about Heatmaps](https://www.htgmolecular.com/blog/2023-05-03/understanding-heat-maps-in-gene-expression-profiling)

Heatmaps allow researchers to quickly and easily identify patterns of gene expression that are associated with specific conditions or treatments and use color coding to indicate the magnitude of values.

* By measuring the number of RNA molecules produced by genes in a particular sample, researchers can determine the level of gene expression.
  * The default settings for the heat map are shown in the image below. The heatmap links genes depending on how alike they are based on the conditions set in the experiment.
  * **Log2 Transformation:** The Log2 transformation is available to enhance the display of genes with more pronounced differences. Applying Log2 transformation can reveal subtler variations in gene expression between conditions.

**Filtering Heatmaps**

Users can filter the genes displayed on the heatmap based on their significance or fold change.

* Filtering by Significance (Adjusted p values) and log2 foldchange thresholds:
  * This filtering helps highlight genes with statistically significant expression changes. Analyze the heatmap to identify gene expression patterns associated with specific conditions or treatments.
  * Users can switch off the clustering of rows or columns of genes causing the dendrograms to be removed. This also affects how the heatmap links genes based on similarity. In some cases, toggling off these types of clustering can help focus on specific groups of genes and their expression patterns. Pay attention to color intensity and clustering.

**Heat map clustering options**

Users can select a clustering method to display results.

* Choosing Clustering Method: The default is often set to UPGMA (Unweighted Pair Group Method with Arithmetic Mean). Clustering helps group genes with similar expression profiles, making patterns more apparent so it can be useful to compare multiple options.
  * In this example we compare UPGMA to the Nearest Point Algorithm, Furthest Point Algorithm and the Ward variance minimization

#### &#x20;<a href="#gehyxb9h3uh2" id="gehyxb9h3uh2"></a>

### Ideogram <a href="#k0o83ew6adb9" id="k0o83ew6adb9"></a>

**Ideograms provide a schematic representation of chromosomes, and they are used to show the relative size of the chromosomes and their characteristic banding patterns.**

[Link to read about Ideograms](https://en.wikipedia.org/wiki/Ideogram)

* The ideogram available through OSDR offers three options for customization.
  * Two options are related to how the user would like to filter significant genes, the first choice being genes with an adjusted P value less than the value set in the text box, and the second choice would be filtered by a Log2 FC value greater than the value set in the text box. The third option to customize the plot allows the user to change the layout of annotations for the Ideogram from a drop-down menu.

### &#x20;<a href="#id-4ngxvslairnl" id="id-4ngxvslairnl"></a>

### DGE Table <a href="#id-2ub25g3cd45c" id="id-2ub25g3cd45c"></a>

**Each study will have an associated Differential Gene Expression (DGE) table available that includes information on each sample from the study.**

To export a Differential Gene Expression (DGE) table from the study visualization page, follow these steps:

1. **Locate the DGE Table:** Scroll down to the bottom of the study visualization page. There, you'll find the Differential Gene Expression table containing valuable data.
2. **Copy the Table to Clipboard:** Identify the "Copy" button within the DGE table. It should be prominently displayed. Click the "Copy" button. This action will copy the entire table, including all data, headers, and values, to your device's clipboard which can then be placed into Word documents and/or spreadsheets.
3. **Or save and export:** As CSV, Excel, PDF, or send to Printer for a physical copy:

* To save the data in various file formats, look for the corresponding buttons.
  * **For a CSV file:** Locate and click the "CSV" button. This will prompt a download of the DGE table data in CSV format to your device.
  * **For an Excel file:** Look for the "Excel" button. Click it to initiate the download of the DGE table data in Excel format (XLSX) to your device.
  * **For a PDF file:** Find and select the "PDF" button. This action will convert the DGE table into a PDF file that you can save to your device.
  * **For Printing:** Spot the "Print" button. Clicking this will open a new window displaying a printer-friendly version of the DGE table. You can then use your browser's print functionality to print the table directly.

**Note:** Choose the method that best suits your needs to access and analyze the DGE data efficiently.

#### &#x20;<a href="#id-8vsniiswkqfg" id="id-8vsniiswkqfg"></a>

#### Multiple study GSEA <a href="#d7vfdcantbds" id="d7vfdcantbds"></a>

**GSEA stands for gene set enrichment analysis, a method to identify gene groups that are overrepresented in a large gene set. It uses statistics to pinpoint significantly enriched or depleted gene classes.**

[Link to read about GSEA](https://www.pnas.org/doi/10.1073/pnas.0506580102)

On the Gene Lab Visualization Portal, you'll find a dedicated GSEA section for each study. Within this GSEA section, there are various parameters you can customize:

1. **Choose genesets:** Select the database geneset enrichment analysis. The default is "KEGG 2019," which is recommended.
2. **Permutations:** Decide the number of permutations you desire and whether they're based on phenotypes or gene sets.
3. **Gene Number Range:** Adjust the minimum and maximum gene sizes. Increasing the minimum size omits genes with fewer than 15 data points, same for the maximum size.
4. **Weighted Score Type:** Defaults to one, representing the t-test. Alternatively, choose signal-to-noise, fold change, or log2 fold change.
5. **Statistical Method:** Select your preferred statistical method. The default is the t-test.

To update the plot with your changes, simply click "Update." A range of plot types is available:

* **Normalized Enrichment Score (NES) Table:** View different gene sets in a table format. Export this table using the options at the top.
* **NES Plot**: The default plot displays normalized enrichment scores based on gene sets.
* **Dot Plot:** Similar to NES Plot, it showcases the top six gene sets based on false discovery rate (FDR). FDR indicates the likelihood that a result is valid, e.g., FDR of 0.25 means a 25% chance of validity.
* **Enrichment Plot**: This reveals the fold change distribution of the top three gene sets with an FDR of under 0.25.
* **Network Plot:** Visualize relationships between gene sets using a network plot.
* **GSEA Info:** For in-depth details about GSEA creation, statistics, and plot documentation.

With these steps, you can effectively navigate and utilize the GSEA section, gaining insights into gene set enrichment analysis for your study. A more detailed explanation can be found in the earlier section on the OSDR single study data visualization application.

### &#x20; <a href="#l8m8wfdlxlsd" id="l8m8wfdlxlsd"></a>

### &#x20;<a href="#id-9fsa75cgexjv" id="id-9fsa75cgexjv"></a>

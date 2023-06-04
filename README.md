# Challenge3_Figures_Code

Code to generate Challenge 3 figures for LRGASP paper.  The data is available from Sage Synapse. 

To generate the figures, obtain code from GitHub:

```
git clone git@github.com:LRGASP/Challenge3_Figures_Code.git
cd Challenge2_Figures_Code
```

The data file is `Challenge3_Figures_Data.zip`, Synapse id  `XXX`.
Install synapseclient if necessary: `pip install synapseclient`

Download, extract, and run the R programs to build figures into the `output` directory:

```
synapse get xxx
unzip -q Challenge3_Figures_Data.zip
Rscript Code_figures_Challenge3_LRGASP_paper_figure.R
Rscript Code_supplementary_figures_Challenge3.R
```




# XomicsEnrich
Integration omics summary data and enrichment analysis
`Requirement Dependencies`: MATLAB 2018b; SPM 12;

#Table of content
1. Use instructions and parameter selection.
2. Explanation of the output files.
3. Posthoc enrichment analysis with R code.


## 1. Use instructions and parameter selection:
Run the "XomicsEnrich/main/XomicsEnrich.m" in the matlab, and you can refer to a GUI like the below figure:
![Input](https://github.com/CaoLuolong/XomicsEnrich/blob/main/pics/GUI.png)

### 1.1 Select the input file.
A text file could be read by matlab(.txt, .xls, ...), and contains at least 2 column, the first column is a list of brain regions, the following column(s) is corresponding region's eigenvalue (trait).For example:
![Input](https://github.com/CaoLuolong/XomicsEnrich/blob/main/pics/input_example.png)

Here we collected 5 types of brain atlas totally. All brain region IDs of different atlas were united coded and recorded in "./rawdata/brain_atlas/".

### 1.2 Select brain atlas and parameters.
`Brain Atlas`: DK68_aparcaseg; 500.aparc; power264PNC; AAL1_ROI_MNI_V4; AAL3v1_1mm  
`Sphere`: left; map2left; whole; map2whole  
`Interpolate`: no; nearest; gauss  

### 1.3 Input PPI parameters.
`PPI database`: String full version11.5; String full version12.0; String physical version11.5; String physical version12.0;  
`PPI neighber`: a integer k (default=1), the k-th neighbour of a gene in the PPI.  
`PPI quantile`: a decimal number j between 0 and 1 (default=0.5), the threshold of connection edge in the PPI.  

### 1.4 Input partial least squares regression parameters.
`PLS dimension`: the dimention of PLS analysis. Typically the first 5 dim could explain enough of the trait variance.  
`PLS permutation`: the number of permutation times in analysis. Typically 1000 times permutation enables the PLS analysis have a constant dim-signaficance result.  
`PLS bootstrap`: the number of bootstrap times in analysis. Typically 3000 times bootstrap enables the PLS analysis have a constant factor-weight result.  
### 1.5 The output path.
`Output path`: the result of analysis will be stored at this path automately.

### 1.6 Load and Analysis.

## 2. Explanation of the output files.
1) a "Run.log" file records he parameters and time of this analysis program.  
2) a .mat file contains gene expression matrix and PET map of biomoleculars.  
3) a .fig file depicts the variance proportation of trait explained by genes/PETs on PLS analysis.  
4) a .mat file contains gene/PET weights on trait. This depict a gene could impact the trait through which biomolecules(PETs).  
5) a .txt file contains: all genes' association with trait by totaly 10 integration methods (rows * columns).  
6) a .txt file contains: all genes' association with trait by 6 of the 10 methods.  
7) a .txt file contains: all genes' association with trait by addPPI method only.  

## 3. Posthoc enrichment analysis with R code.
The enrichment analysis can refer to "./utils/pathway_enrichment.R" (run in R session), and you can do GO term, KEGG pathway, or celltype enrichent analysis and plot results of interest. 

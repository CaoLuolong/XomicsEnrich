# XomicsEnrich
Integration omics summary data and enrichment analysis

## Useage:
Run the "XomicsEnrich/main/XomicsEnrich.m" in the matlab, and youcan refer to a GUI like the below figure:
![Input](https://github.com/CaoLuolong/XomicsEnrich/blob/main/pics/GUI.png)

### 1. Select the input file.
A text file could be read by matlab(.txt, .xls, ...), and contains at least 2 column, the first column is a list of brain regions, the following column(s) is corresponding region's eigenvalue (trait).For example:
![Input](https://github.com/CaoLuolong/XomicsEnrich/blob/main/pics/input_example.png)

### 2. Select brain atlas and parameters.
`Brain Atlas`: DK68_aparcaseg; 500.aparc; power264PNC; AAL1_ROI_MNI_V4; AAL3v1_1mm
`Sphere`: left; ap2left; whole; map2whole
`Interpolate`: no; nearest; gauss

### 3. Input PPI parameters.
`PPI neighber`: a integer k (default=1), the k-th neighbour of a gene in the PPI.
`PPI quantile`: a decimal number j between 0 and 1 (default=0.75), the threshold of connection edge in the PPI.

### 4. Input partial least squares regression parameters.

### 5. Select the output path.

### 6. Load AHBA gene expression data.

### 7. Load [STRING](https://cn.string-db.org/) PPI data.

### 8. Run analysis.

## Output results.
A .mat file contains gene expression matrix and PET map of biomoleculars; a .txt file contains: all genes' association with trait by 14 integration methods (rows * columns).

## Posthoc enrichment analysis with R code.

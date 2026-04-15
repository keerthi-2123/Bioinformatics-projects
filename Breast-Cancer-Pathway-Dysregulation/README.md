## Breast Cancer Pathway Dysregulation

**Biological Question:** Are DNA repair genes more dysregulated than TP53 checkpoint genes in breast cancer?
**Dataset:** TCGA Breast Cancer gene expression (Normal vs Tumor samples)

**Methods:** Fold change analysis + independent t-test (scipy) + 
heatmap visualization (seaborn)

**Workflow:** 1. Loaded normal and tumor datasets
              2. Extracted selected genes
              3. Calculated mean expression values
              4. Computed fold change (Tumor - Normal)
              5. Performed Welch's t-test
              6. Generated plots

**Key Finding:** 1. Several DNA repair genes showed elevated dysregulation in tumor tissue.
                 2. TP53 checkpoint genes displayed mixed regulation patterns.
                 3. These findings suggest altered DNA damage response signaling in breast cancer.
                 
**Future improvements:** 1. Add pathway enrichment analysis
                         2. Include additional cancer datasets
                         3. Perform clustering and PCA
                         4. Add survival analysis
**Tools:** Python, pandas, matplotlib, seaborn, scipy

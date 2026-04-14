import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

# Load datasets
normal = pd.read_csv("BC-TCGA-Normal.txt", sep="\t")
tumor = pd.read_csv("BC-TCGA-Tumor.txt", sep="\t")

# Gene groups
tp53_genes = ["TP53", "MDM2", "CDKN1A", "BAX", "CHEK2"]
repair_genes = ["BRCA1", "BRCA2", "RAD51", "PALB2", "ATM", "FANCD2"]

all_genes = tp53_genes + repair_genes

# Filter rows
normal_f = normal[normal["Hybridization REF"].isin(all_genes)].copy()
tumor_f = tumor[tumor["Hybridization REF"].isin(all_genes)].copy()

# Mean expression
normal_f["Normal_Mean"] = normal_f.iloc[:,1:].mean(axis=1)
tumor_f["Tumor_Mean"] = tumor_f.iloc[:,1:].mean(axis=1)

# Merge
result = pd.DataFrame({
    "Gene": normal_f["Hybridization REF"].values,
    "Normal": normal_f["Normal_Mean"].values,
    "Tumor": tumor_f["Tumor_Mean"].values
})

# Fold change
result["FoldChange"] = result["Tumor"] - result["Normal"]

# Pathway labels
result["Pathway"] = result["Gene"].apply(
    lambda x: "TP53_Checkpoint" if x in tp53_genes else "DNA_Repair"
)

print(result)

# Save CSV
result.to_csv("gene_expression_summary.csv", index=False)

# T-test
tp53_vals = result[result["Pathway"]=="TP53_Checkpoint"]["FoldChange"]
repair_vals = result[result["Pathway"]=="DNA_Repair"]["FoldChange"]

t_stat, p_val = ttest_ind(tp53_vals, repair_vals, equal_var=False)

print("\nT-test Results")
print("T-statistic:", t_stat)
print("P-value:", p_val)

# Save stats
pd.DataFrame({
    "T_statistic":[t_stat],
    "P_value":[p_val]
}).to_csv("pathway_statistics.csv", index=False)

# BAR GRAPH
plt.figure(figsize=(10,6))
sns.barplot(data=result, x="Gene", y="FoldChange", hue="Pathway")
plt.axhline(0, linestyle="--")
plt.title("Tumor vs Normal Fold Change")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("foldchange_plot.png", dpi=300)
plt.show()

# HEATMAP
heat = result.set_index("Gene")[["Normal","Tumor"]]

plt.figure(figsize=(8,6))
sns.heatmap(heat, annot=True, cmap="coolwarm", center=0)
plt.title("Expression Heatmap")
plt.tight_layout()
plt.savefig("heatmap.png", dpi=300)
plt.show()
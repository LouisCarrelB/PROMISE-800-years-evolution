import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

gene_name = "ENSG00000109339/"
GENE = "../DATA/" + gene_name

file_path_original = GENE + "inter/scores_ali_can_updated.csv"
file_path_alternative = GENE + "inter/scores_ali_alt_updated.csv"
df_original = pd.read_csv(file_path_original)
df_alternative = pd.read_csv(file_path_alternative)

counts_alt = df_alternative['Regne'].value_counts()
counts_o = df_original['Regne'].value_counts()

plt.figure(figsize=(10, 6))
colors = plt.cm.Paired(range(len(counts_alt)))
wedges, texts, autotexts = plt.pie( 
    counts_alt, 
    autopct=lambda pct: '{:.1f}%'.format(pct) if pct > 5 else '',
    textprops=dict(color="w"),
    colors=colors,
    startangle=140,
)
plt.legend(wedges, counts_alt.index, title="Règnes", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
plt.title('Proportion de règnes dans alt')
plt.show()

plt.figure(figsize=(10, 6))
colors = plt.cm.Paired(range(len(counts_o)))
wedges, texts, autotexts = plt.pie(
    counts_o, 
    autopct=lambda pct: '{:.1f}%'.format(pct) if pct > 5 else '',
    textprops=dict(color="w"),
    colors=colors,
    startangle=140,
)
plt.legend(wedges, counts_o.index, title="Règnes", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
plt.title('Proportion de règnes dans can')
plt.show()




common_classes = set(df_original['Class']).intersection(set(df_alternative['Class']))
common_regnes = set(df_original['Regne']).intersection(set(df_alternative['Regne']))


venn2(subsets=(len(set(df_original['Class'])) - len(common_classes),
               len(set(df_alternative['Class'])) - len(common_classes),
               len(common_classes)),
      set_labels=('can Embranchements', 'alt Embranchements'))

plt.title('Venn Diagram for Embranchement')
plt.show()

venn2(subsets=(len(set(df_original['Regne'])) - len(common_regnes),
               len(set(df_alternative['Regne'])) - len(common_regnes),
               len(common_regnes)),
      set_labels=('can Regnes', 'alt Regnes'))

plt.title('Venn Diagram for Regnes')
plt.show()




#for the notebook later

import pandas as pd

df = pd.read_table("../data/gtotree/Genomes_summary_info.tsv")
df["is_complete"] = ["Yes","No","Yes","Yes","Yes","Yes","No","Yes","Yes","Yes","Yes","No","Yes","No","No","No","No"]
df["is_endosymb"] = ["No","Yes","Yes","Yes","Yes","No","Yes","Yes","No","Yes","No","Yes","Yes","Yes","Yes","Yes","Yes"]

df["assembly_id"] = df["assembly_id"].apply(lambda name: '_'.join(name.split("_")[:-2]))

print(df)

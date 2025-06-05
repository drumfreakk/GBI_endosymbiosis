
import pandas as pd
import os
import matplotlib.pyplot as plt

#Import the species table & info
species_df = pd.read_table("../data/gtotree/Genomes_summary_info.tsv")
species_df["is_complete"] = ["Yes","No","Yes","Yes","Yes","Yes","No","Yes","Yes","Yes","Yes","No","Yes","No","No","No","No"]
species_df["is_endosymb"] = ["No","Yes","Yes","Yes","Yes","No","Yes","Yes","No","Yes","No","Yes","Yes","Yes","Yes","Yes","Yes"]

species_df["assembly_id"] = species_df["assembly_id"].apply(lambda name: '_'.join(name.split("_")[:-2]))


#Import the COG annotations into a dataframe
cog = pd.read_table("../data/functional_information/cog-20.def.tab", names=["ID", "CATEGORY", "NAME", "SYMBOL", "PATHWAYS", "PUBMED", "PDB"], encoding="latin")

cog_cats = pd.read_table("../data/functional_information/fun-20.tab", names=["CATEGORY", "??", "CATEGORY_DESCRIPTION"])

cog = cog.join(cog_cats.set_index("CATEGORY"), on="CATEGORY")

print(cog)

#Import the attributes from the gff files into a dataframe
data = {}

species = os.listdir("attrs")

for specy in species:
	i = 0
	with open("attrs/"+specy, "r") as file:
		for line in file:
	#		print("\n\n", i)
			attrs_linked = line.strip().split(";")
			attrs = [["Species",specy]] + [i.split("=") for i in attrs_linked]
			had_keys = []
			for keyval in attrs:
				if keyval[0] == "db_xref":
					keyval[0] = "cog_id"
					keyval[1] = keyval[1].split(":")[1]
				if keyval[0] in data:
	#				print("IN", keyval)
					data[keyval[0]].append(keyval[1])
				else:
	#				print("NOT", keyval)
					#Initialise an array with empty values
					data[keyval[0]] = [keyval[1] if j==i else "" for j in range(i+1)]
				had_keys.append(keyval[0])
			for key in list(set(data.keys()).difference(had_keys)):
				data[key].append("")
			i+=1

genes = pd.DataFrame.from_dict(data)

#pd.options.display.max_columns = None

#Add the COG annotations to the gff annotations
genes = genes.join(cog.set_index("ID"), on="cog_id")
#print(genes)

categories=['T','L','G','E','P','R','F','K','C','H','J','M','O','V','I','S','Q','N','U','W','X','D','A','Z']
columns = ["Species"] + categories
#columns=["Species",'T','L','G','E','P','R','F','K','C','H','J','M','O','V','I','S','Q','N','U','W','X','D','A','Z']
counts = {}
for i in columns:
	counts[i] = [0]
counts = pd.DataFrame.from_dict(counts)


def get_gene_category_counts(species):
	counts = {}
	for i in columns:
		counts[i] = [0]
#counts = {"Species": species}

	counts["Species"] = species
	for index, row in genes.loc[genes['Species'] == species].iterrows():
		for category in row["CATEGORY"]:
			counts[category][0] += 1
	return pd.DataFrame.from_dict(counts)

for specy in species:
	count = get_gene_category_counts(specy)
	counts = pd.concat([counts, count], ignore_index=True)

counts = counts.drop(0)

#Now counts is a dataframe with a row per species, with columns for the amount of genes in a given COG category, and wether they are complete and endosymbiotic
counts = counts.join(species_df[["assembly_id","is_complete","is_endosymb"]].set_index('assembly_id'), on='Species')

print(counts)


# Sum all of the gene counts in the different categories
free_living_mean = counts.loc[counts['is_endosymb'] == "No", categories].mean().to_frame(name="Free Living")
symbiont_complete_mean = counts.loc[(counts['is_endosymb'] == "Yes") & (counts['is_complete'] == "Yes"), categories].mean().to_frame(name="Symbiont (complete)")
symbiont_incomplete_mean = counts.loc[(counts['is_endosymb'] == "Yes") & (counts['is_complete'] == "No"), categories].mean().to_frame(name="Symbiont (incomplete)")

free_living_std = counts.loc[counts['is_endosymb'] == "No", categories].std().to_frame(name="Free Living")
symbiont_complete_std = counts.loc[(counts['is_endosymb'] == "Yes") & (counts['is_complete'] == "Yes"), categories].std().to_frame(name="Symbiont (complete)")
symbiont_incomplete_std = counts.loc[(counts['is_endosymb'] == "Yes") & (counts['is_complete'] == "No"), categories].std().to_frame(name="Symbiont (incomplete)")

mean_counts = free_living_mean.join(symbiont_complete_mean.join(symbiont_incomplete_mean))
std_counts = free_living_std.join(symbiont_complete_std.join(symbiont_incomplete_std))
print(mean_counts)
print(std_counts)
mean_counts.plot(kind='bar', yerr=std_counts, rot=0)
plt.show()




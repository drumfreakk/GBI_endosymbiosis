
import pandas as pd
import os


#Import the COG annotations into a dataframe
cog = pd.read_table("../data/functional_information/cog-20.def.tab", names=["ID", "CATEGORY", "NAME", "SYMBOL", "PATHWAYS", "PUBMED", "PDB"], encoding="latin")


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

df = pd.DataFrame.from_dict(data)

#pd.options.display.max_columns = None

#Add the COG annotations to the gff annotations
df = df.join(cog.set_index("ID"), on="cog_id")
print(df)




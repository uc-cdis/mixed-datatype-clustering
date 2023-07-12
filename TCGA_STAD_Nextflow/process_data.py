# The writefile command writes the code into the specified file rather than running it

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
import sys

# The value to use as a mask is specified as a command line argument
# since it is used in multiple scripts
mask_na = float(sys.argv[1])

with open("expression_data.json", "r") as f:
    decoded = json.load(f)
matrix = decoded["data"]["quantDataMatrix"]

ga = pd.DataFrame(matrix[1:], columns=matrix[0]).set_index("Gene/Aliquot")
oldnames = list(ga.columns)
newnames = [x.split(":")[1] for x in oldnames]
ga.rename(columns=dict(zip(oldnames, newnames)), inplace=True)
ga = ga.sort_index(axis=1)

for col in ga.keys():
    ga[col] = pd.to_numeric(ga[col], errors="coerce")

ga = ga.fillna(mask_na)

# The processed data is saved to be fed into later steps in the workflow
ga.to_csv("processed_data.csv")

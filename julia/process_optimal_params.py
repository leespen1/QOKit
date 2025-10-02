import pandas as pd
import glob
import os

# mapping CSV columns to desired names
csv_columns = ["h_tweak_sub", "hc_tweak_add", "l_tweak_mul", "r_tweak_mul", "min_mse", "mean_mse"]

def cast_value(val: str):
    """Try to cast val to int or float if possible, otherwise keep as string."""
    try:
        if "." in val or "e" in val.lower():  # float-like
            return float(val)
        return int(val)  # int-like
    except ValueError:
        return val

def parse_value(val: str):
    """Parse a value, splitting ranges into min/max if needed."""
    if ":" in val:
        min_val, max_val = val.split(":", 1)
        return cast_value(min_val), cast_value(max_val)
    else:
        return cast_value(val)

# directory containing CSV files
files = glob.glob("saved_results/optimal*.csv")

all_dfs = []
for file in files:
    # --- read CSV content ---
    df = pd.read_csv(file, sep="\t", header=None)
    df.columns = csv_columns

    # --- clean filename ---
    fname = os.path.basename(file)
    clean_name = fname
    if clean_name.startswith("optimalparams_"):
        clean_name = clean_name[len("optimalparams_"):]
    if clean_name.endswith(".csv"):
        clean_name = clean_name[:-4]

    # --- extract metadata ---
    for token in clean_name.split("_"):
        if "=" in token:
            key, val = token.split("=", 1)
            parsed = parse_value(val)
            if isinstance(parsed, tuple):
                df[f"{key}_min"] = parsed[0]
                df[f"{key}_max"] = parsed[1]
            else:
                df[key] = parsed

    # optional: keep original filename
    df["filename"] = fname

    all_dfs.append(df)

# combine all dataframes
big_df = pd.concat(all_dfs, ignore_index=True)

print(big_df.dtypes)
print(big_df.head())


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Example: big_df is the DataFrame we created earlier

# List of numeric fields to plot
fields = ["h_tweak_sub", "hc_tweak_add", "l_tweak_mul", "r_tweak_mul"]
#fields = ["hc_tweak_add", "l_tweak_mul", "r_tweak_mul"]

# Iterate over each graph type
for graphtype in big_df['graphtype'].unique():
    df_graph = big_df[big_df['graphtype'] == graphtype]

    # Iterate over each probability
    for prob in df_graph['probability'].unique():
        df_prob = df_graph[df_graph['probability'] == prob]

        plt.figure(figsize=(12, 6))
        for field in fields:
            # plot each field vs numnodes
            sns.lineplot(
                data=df_prob,
                x='numnodes',
                y=field,
                marker='o',
                label=field
            )

        plt.title(f"{graphtype} graphs, probability={prob}")
        plt.xlabel("Number of nodes")
        plt.ylabel("Field value")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

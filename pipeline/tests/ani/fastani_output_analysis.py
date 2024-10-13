import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


RAW_OUTPUTS_DIR = r"C:\Users\yairs\OneDrive\Documents\University\Master\stuff\pactobacter ANI analysis\outputs"

all_dfs = []
for file_name in os.listdir(RAW_OUTPUTS_DIR):
    file_path = os.path.join(RAW_OUTPUTS_DIR, file_name)
    df = pd.read_csv(file_path, delimiter='\t', names=['query', 'subject', 'ani_value', 'orthologous_segments', 'total_segments'])
    all_dfs.append(df)

combined_df = pd.concat(all_dfs, ignore_index=True)

query_short = []
subject_short = []
for index, row in combined_df.iterrows():
    query_short.append(os.path.splitext(os.path.basename(row['query']))[0])
    subject_short.append(os.path.splitext(os.path.basename(row['subject']))[0])

combined_df['query'] = query_short
combined_df['subject'] = subject_short

ani_values_df = combined_df.pivot_table(index='query', columns='subject', values='ani_value')

plt.subplots(figsize=(40, 40))
plot = sns.heatmap(ani_values_df, cmap='Blues')
fig = plot.get_figure()
fig.savefig('heatmap.png')

ani_values_df.replace({100: np.nan}, inplace=True)
max_values = ani_values_df.max(axis=1)
max_columns = ani_values_df.idxmax(axis=1)
ani_values_df = ani_values_df.assign(max_value=max_values.values)
ani_values_df = ani_values_df.assign(max_column=max_columns.values)
ani_values_df.to_csv('pactobacter.csv')

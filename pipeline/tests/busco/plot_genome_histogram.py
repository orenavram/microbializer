import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

FILE = r"C:\Users\TalPNB22\OneDrive\Documents\University\Master\Posters\73 ecoli outputs\genomes_completeness.csv"
OUTPUT_FILE = r'C:\temp\pic.png'
title = 'Genome BUSCO completeness score'
xlabel = 'Genome BUSCO completeness score per genome'

df = pd.read_csv(FILE)
df.set_index('Genome', inplace=True)
sns.histplot(df, x=title, kde=True)
plt.title(f'Distribution of {title}', fontsize=20, loc='center', wrap=True)
plt.xlabel(xlabel, fontsize=15)
plt.ylabel('Genomes count', fontsize=15)
plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))  # make y-axis integer
plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=600)

plt.clf()

import os
import pandas as pd

DF_1_PATH = r"C:\Users\TalPNB22\Downloads\10\putative_orthologs_table.txt"
DF_2_PATH = r"C:\Users\TalPNB22\Downloads\73\putative_orthologs_table.txt"


def sort_df(path):
    df = pd.read_csv(path)

    for index, row in df.iterrows():
        for column in df.columns:
            if pd.notnull(row[column]):
                df.at[index, column] = ';'.join(sorted(row[column].split(';')))

    df = df.sort_values(by=df.columns.tolist()).reset_index(drop=True)

    df.to_csv(path, index=False)


def main():
    sort_df(DF_1_PATH)
    sort_df(DF_2_PATH)


if __name__ == '__main__':
    main()

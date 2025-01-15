import pandas as pd

PATH = r"C:\Users\TalPNB22\Downloads\putative_orthologs_table_73_old.txt"
OUTPUT_PATH = r"C:\Users\TalPNB22\Downloads\putative_orthologs_table_73_old_sorted.csv"

# PATH = r"C:\temp\compare_salmonella_300\old\putative_orthologs_table_old.txt"
# OUTPUT_PATH = r"C:\temp\compare_salmonella_300\old\putative_orthologs_table_old_sorted.txt"


def main():
    df = pd.read_csv(PATH)
    df_sorted_cells = df.map(lambda x: ';'.join(sorted(x.split(';'))) if not pd.isna(x) else x)
    df_sorted_rows = df_sorted_cells.sort_values(by=list(df.columns[1:])).reset_index(drop=True)
    df_sorted_rows['OG_name'] = [f'OG_{i}' for i in range(len(df_sorted_rows.index))]
    df_sorted_rows.to_csv(OUTPUT_PATH, index=False)


if __name__ == '__main__':
    main()

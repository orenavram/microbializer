import os
import pandas as pd

# BASE_PATH = r"/groups/pupko/yairshimony/microbializer_prod/pipeline/data/kegg"
BASE_PATH = r"C:\repos\microbializer\pipeline\data\kegg"

RAW_DATA_PATH = os.path.join(BASE_PATH, 'raw_data')
CLEAN_DATA_PATH = os.path.join(BASE_PATH, 'cleaned_data')

RAW_KO_LIST_PATH = os.path.join(RAW_DATA_PATH, 'ko_list')
CLEANED_KO_LIST_PATH = os.path.join(CLEAN_DATA_PATH, 'ko_list')
PROKARYOTE_PROFILES_LIST_PATH = os.path.join(RAW_DATA_PATH, 'profiles', 'prokaryote.hal')


def main():
    with open(PROKARYOTE_PROFILES_LIST_PATH, 'r') as fp:
        prokaryote_profiles_list = fp.readlines()

    prokaryote_profiles_list = [os.path.splitext(profile.strip())[0] for profile in prokaryote_profiles_list if profile]

    ko_list_df = pd.read_csv(RAW_KO_LIST_PATH, sep='\t')
    ko_list_df = ko_list_df[ko_list_df['threshold'] != '-']
    ko_list_df = ko_list_df[ko_list_df['knum'].isin(prokaryote_profiles_list)]
    ko_list_df = ko_list_df[['knum', 'threshold', 'score_type', 'definition']]

    os.makedirs(CLEAN_DATA_PATH, exist_ok=True)
    ko_list_df.to_csv(CLEANED_KO_LIST_PATH, index=False)


if __name__ == "__main__":
    main()



PATH_1 = r"C:\Users\TalPNB22\Downloads\original_100\GCF_041463645.1_vs_GCF_041463755.1.m8"
PATH_2 = r"C:\Users\TalPNB22\Downloads\one_mmseqs_100\GCF_041463645.1_vs_GCF_041463755.1.m8"


def get_pairs(path):
    pairs = set()
    with open(path) as fp:
        for line in fp:
            tokens = line.split(',')
            if 'score' not in tokens:
                gene1, gene2 = tokens[:2]
                pairs.add(tuple(sorted((gene1, gene2))))

    return pairs


def main():
    pairs_1 = get_pairs(PATH_1)
    pairs_2 = get_pairs(PATH_2)

    diff_1_to_2 = pairs_1 - pairs_2
    diff_2_to_1 = pairs_2 - pairs_1

    print(f"Number of pairs in PATH_1 but not in PATH_2: {len(diff_1_to_2)}. pairs: {diff_1_to_2}")
    print(f"Number of pairs in PATH_2 but not in PATH_1: {len(diff_2_to_1)}. pairs: {diff_2_to_1}")


if __name__ == '__main__':
    main()

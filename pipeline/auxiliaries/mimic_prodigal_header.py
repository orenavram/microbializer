import os


def mimic_prodigal_output(orfs_dir, output_orf_file_extension):
    for file_name in os.listdir(orfs_dir):
        file_path = os.path.join(orfs_dir, file_name)

        # edit headers of ORFs to match the structure of prodigal output
        fixed_content = ''
        with open(file_path, 'r') as orfs_file:
            for line in orfs_file:
                if line.startswith('>'):
                    fixed_content += f'{line.strip()} # START # END # 1 # \n'
                else:
                    fixed_content += line

        # override the old file with the fixed content
        with open(file_path, 'w') as f:
            f.write(fixed_content)

        # change file name to match the output of step 2
        os.rename(file_path, f'{os.path.splitext(file_path)[0]}.{output_orf_file_extension}')


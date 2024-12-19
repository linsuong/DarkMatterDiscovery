from file_join import join_dat_files
from file_cleanup_v3 import clean_incomplete_rows

dates = [
        'Dec2',
        'Dec13',
        'Dec15',
        'Dec15_Linux',
        'Dec16_Linux'
        ]

for date in dates:
    folder_path = f"scans/5-D_scans/run_{date}/dat_files"
    output_file = f"scans/5-D_scans/run_{date}/combined_{date}.dat.gz"
    join_dat_files(folder_path, output_file)
    
for date in dates:
    input_file = f"scans/5-D_scans/run_{date}/combined_{date}.dat.gz"
    output_file = f'scans/5-D_scans/cleaned_dat_files/combined_{date}_clean.dat.gz'
    
    clean_incomplete_rows(input_file, output_file, expected_columns=12)

join_dat_files('scans/5-D_scans/cleaned_dat_files', 'scans/5-D_scans/combined.dat.gz')
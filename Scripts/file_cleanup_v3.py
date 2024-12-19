import pandas as pd

def clean_incomplete_rows(file_path, output_path, sep=r'\s+', expected_columns=12):
    """
    Removes rows with incomplete data (fewer than expected columns) from a .dat file.

    Parameters:
    - file_path: Path to the input .dat file.
    - output_path: Path to save the cleaned .dat file.
    - sep: Separator for columns (default is whitespace).
    - expected_columns: Number of expected columns in the file.
    """
    cleaned_rows = []

    # Open and read the file line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Split line into components
            split_line = line.split()
            # Keep lines that have the correct number of columns
            if len(split_line) == expected_columns:
                cleaned_rows.append(line)

    # Write the cleaned rows to the output file
    with open(output_path, 'w') as output_file:
        for row in cleaned_rows:
            output_file.write(row)

    print(f"Cleaned data saved to {output_path}")

if __name__ == "__main__":
    dates = [
            'Dec2',
            'Dec13',
            'Dec15',
            'Dec15_Linux',
            'Dec16_Linux' 
            ]

    for date in dates:
        input_file = f"scans/5-D_scans/run_{date}/combined_{date}.dat"
        output_file = f'scans/5-D_scans/cleaned_dat_files/combined_{date}_clean.dat'
        
        clean_incomplete_rows(input_file, output_file, expected_columns=12)
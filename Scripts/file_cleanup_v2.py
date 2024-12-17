import os

def clean_dat_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Split line into columns
            columns = line.split()
            # Check if the row has the exact number of columns
            if len(columns) == expected_columns:
                outfile.write(line)

# Set the expected number of columns (count based on your header line)
expected_columns = 12  # Update this based on your data's actual column count

# Process multiple files
input_directory = "/Users/linusong/Repositories/DarkMatterDiscovery/scans/5-D_scans/run_Dec16_Linux"  # Folder containing .dat files
output_directory = "/Users/linusong/Repositories/DarkMatterDiscovery/scans/5-D_scans/"  # Folder for cleaned files

# Create output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Iterate through all .dat files in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith(".dat"):
        input_path = os.path.join(input_directory, filename)
        output_path = os.path.join(output_directory, f"cleaned_{filename}")
        clean_dat_file(input_path, output_path)
        print(f"Cleaned: {filename} -> {output_path}")
        

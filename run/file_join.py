import os

def join_dat_files(folder_path, output_file):
    if not os.path.isdir(folder_path):
        print("Invalid folder path.")
        return

    with open(output_file, 'w', encoding='utf-8') as outfile:
        header_written = False  # To track if the header is already written

        for file_name in sorted(os.listdir(folder_path)):  # Process files in sorted order
            if file_name.endswith('.dat'):  # Only process .dat files
                file_path = os.path.join(folder_path, file_name)

                try:
                    with open(file_path, 'r', encoding='utf-8') as infile:
                        lines = infile.readlines()
                        if not lines:
                            print(f"File {file_name} is empty. Skipping.")
                            continue

                        if not header_written:
                            # Write the header from the first file
                            outfile.write(lines[0].rstrip() + '\n')  # Ensure newline
                            header_written = True

                        # Write all non-header lines
                        for line in lines[1:]:
                            outfile.write(line.rstrip() + '\n')  # Ensure newline

                        print(f"Appended data from {file_name}")

                except Exception as e:
                    print(f"Error processing file {file_name}: {e}")

    print(f"All files successfully joined into {output_file}")
    
# Replace with your folder path and output file name
folder_path = "run/5-D_scans/"
output_file = "run/5-D_scans/combined.dat"
join_dat_files(folder_path, output_file)

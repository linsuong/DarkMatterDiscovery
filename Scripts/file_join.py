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
    
if __name__ == "__main__":
    #unfiltered data:
    dates = [
            'Dec2',
            'Dec13',
            'Dec15',
            'Dec15_Linux',
            'Dec16_Linux' 
            ]
    for date in dates:
        folder_path = f"scans/5-D_scans/run_{date}/dat_files"
        output_file = f"scans/5-D_scans/run_{date}/combined_{date}.dat"
        join_dat_files(folder_path, output_file)

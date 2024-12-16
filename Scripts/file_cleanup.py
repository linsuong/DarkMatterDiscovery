import os

def remove_partially_empty_rows(folder_path):
    if not os.path.isdir(folder_path):
        print("Invalid folder path.")
        return

    print(f"Scanning folder: {os.path.abspath(folder_path)}\n")

    # List all files (recursive)
    files_processed = 0
    for root, _, files in os.walk(folder_path):
        print(f"Checking directory: {root}")
        for file_name in files:
            if not file_name.endswith('.dat'):
                print(f"Skipping non-dat file: {file_name}")
                continue

            file_path = os.path.join(root, file_name)

            try:
                with open(file_path, 'r', encoding='utf-8') as file:
                    lines = file.readlines()

                print(f"Processing file: {file_path} (Total lines: {len(lines)})")

                cleaned_lines = []
                for i, line in enumerate(lines):
                    # Debug: Print raw line content
                    #print(f"Line {i}: {repr(line)}")

                    # Strip leading/trailing whitespace
                    stripped_line = line.strip()

                    # Remove blank lines (lines with only whitespace)
                    if not stripped_line:
                        print(f"Removed blank row {i}: {repr(line)}")
                        continue

                    # Split into columns and check for empty values
                    columns = stripped_line.split()
                    if all(col.strip() for col in columns):
                        cleaned_lines.append(line)
                    else:
                        print(f"Removed row {i}: {repr(line)} (partially empty)")

                # Rewrite the file only if changes are needed
                if len(cleaned_lines) != len(lines):
                    with open(file_path, 'w', encoding='utf-8') as file:
                        file.writelines(cleaned_lines)
                    print(f"Cleaned file: {file_name}")
                else:
                    print(f"No changes needed for file: {file_name}")

                files_processed += 1

            except Exception as e:
                print(f"Error processing file '{file_name}': {e}")

    if files_processed == 0:
        print("\nNo valid .dat files found to process.")
    else:
        print(f"\nTotal files processed: {files_processed}")

# Replace with your folder path
folder_path = "run/5-D_scans/run_Dec13"
remove_partially_empty_rows(folder_path)
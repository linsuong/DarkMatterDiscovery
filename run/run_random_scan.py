import os

def remove_partially_empty_rows(folder_path):
    if not os.path.isdir(folder_path):
        print("Invalid folder path.")
        return

    for file_name in os.listdir(folder_path):
        if file_name.endswith('.dat'):  # Only process .dat files
            file_path = os.path.join(folder_path, file_name)

            try:
                with open(file_path, 'r', encoding='utf-8') as file:
                    lines = file.readlines()

                print(f"\nProcessing file: {file_name}")

                # Identify header (optional, assumes first row is header)
                header = lines[0] if lines else ""
                cleaned_lines = [header]  # Keep the header intact

                # Process the data rows
                for i, line in enumerate(lines[1:], start=2):  # Start index at 2 for clarity
                    # Split the line into columns
                    columns = line.strip().split()
                    
                    # Debug: Print line and its columns
                    print(f"Line {i}: {repr(line)} -> Columns: {columns}")
                    
                    # Check if any column is empty or missing
                    if len(columns) == len(header.strip().split()) and all(columns):
                        cleaned_lines.append(line)  # Valid row
                    else:
                        print(f"Removed row {i} (partially empty): {repr(line)}")

                # Write cleaned lines back to the file
                with open(file_path, 'w', encoding='utf-8') as file:
                    file.writelines(cleaned_lines)

                print(f"Finished cleaning file: {file_name}")

            except UnicodeDecodeError:
                print(f"Encoding error in file: {file_name}. Try a different encoding.")
            except Exception as e:
                print(f"Error processing file {file_name}: {e}")

# Replace this with your folder path
folder_path = "run"
remove_partially_empty_rows(folder_path)

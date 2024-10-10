import time
import shutil
import os
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

class NewFileHandler(FileSystemEventHandler):
    def __init__(self, destination_folder):
        self.destination_folder = destination_folder
    
    def on_created(self, event):
        # Ignore directories
        if event.is_directory:
            return

        # Only process .txt files
        if not event.src_path.endswith('.txt'):
            return

        # Get the source path and file name
        src_path = event.src_path
        file_name = os.path.basename(src_path)
        dst_path = os.path.join(self.destination_folder, file_name)

        # Check if file already exists, and if so, modify the name
        dst_path = self._get_non_conflicting_name(dst_path)

        try:
            # Move the .txt file to the destination folder
            shutil.move(src_path, dst_path)
            print(f'Moved file: {src_path} -> {dst_path}')
        except Exception as e:
            print(f'Error moving file {src_path}: {e}')

    def _get_non_conflicting_name(self, path):
        """
        If a file with the same name exists, append '_1', '_2', etc., to the filename.
        """
        base, extension = os.path.splitext(path)
        counter = 1

        # Keep modifying the name until there's no conflict
        while os.path.exists(path):
            path = f"{base}_{counter}{extension}"
            counter += 1

        return path

if __name__ == "__main__":
    # Set your source and destination folders here
    source_folder = '~/Documents/CalcHEP/Work/results'
    destination_folder = '~/Repositories/DarkMatterDiscovery/Results'

    # Create the event handler and observer
    event_handler = NewFileHandler(destination_folder)
    observer = Observer()
    observer.schedule(event_handler, path=source_folder, recursive=False)

    # Start the observer
    observer.start()
    print(f"Monitoring '{source_folder}' for new .txt files...")

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        # Stop the observer if the user presses Ctrl+C
        observer.stop()
        print("\nMonitoring stopped.")
    observer.join()

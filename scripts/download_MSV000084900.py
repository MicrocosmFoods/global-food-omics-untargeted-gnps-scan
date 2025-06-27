import ftplib
import os
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlparse
from tqdm import tqdm
import threading
from queue import Queue
import time


class FTPConnectionPool:
    def __init__(self, hostname, max_connections=1):
        self.hostname = hostname
        self.max_connections = max_connections
        self.connections = Queue()
        self.lock = threading.Lock()

        # Initialize pool with connections
        for _ in range(max_connections):
            ftp = ftplib.FTP(hostname)
            ftp.login()
            ftp.set_pasv(True)
            self.connections.put(ftp)

    def get_connection(self):
        return self.connections.get()

    def return_connection(self, ftp):
        try:
            # Test if connection is still alive
            ftp.voidcmd("NOOP")
            self.connections.put(ftp)
        except:
            # If connection is dead, create new one
            try:
                ftp = ftplib.FTP(self.hostname)
                ftp.login()
                ftp.set_pasv(True)
                self.connections.put(ftp)
            except:
                pass

    def close_all(self):
        while not self.connections.empty():
            try:
                ftp = self.connections.get_nowait()
                ftp.quit()
            except:
                pass


class FTPDownloader:
    def __init__(self, ftp_url, max_workers=1):
        self.url = urlparse(ftp_url)
        self.max_workers = max_workers
        self._tqdm_lock = threading.Lock()
        self.connection_pool = FTPConnectionPool(self.url.hostname, max_connections=1)

    def get_ftp_dir_contents(self, ftp):
        """Get directory contents using the FTP LIST command"""
        contents = []
        ftp.retrlines("LIST", contents.append)
        return contents

    def list_files_recursively(self, path="."):
        """Recursively list all files in the FTP directory"""
        files = []

        try:
            # Get raw list of files
            with ftplib.FTP(self.url.hostname) as ftp:
                ftp.login()
                ftp.set_pasv(True)
                if self.url.path:
                    ftp.cwd(self.url.path)
                if path != ".":
                    ftp.cwd(path)

                entries = self.get_ftp_dir_contents(ftp)

                for entry in entries:
                    # Parse the FTP LIST output
                    # Format typically: "type permissions links owner group size month day time/year name"
                    parts = entry.split(None, 8)  # Split on whitespace, max 8 splits
                    if len(parts) < 9:
                        continue

                    entry_type = parts[0][0]
                    name = parts[8]  # This preserves the exact name from FTP

                    if entry_type == "d":  # Directory
                        new_path = f"{path}/{name}" if path != "." else name
                        subfiles = self.list_files_recursively(new_path)
                        files.extend(subfiles)
                    elif entry_type == "-":  # Regular file
                        size = int(parts[4])
                        file_path = f"{path}/{name}" if path != "." else name
                        files.append({"path": file_path, "size": size})

        except ftplib.error_perm as e:
            print(f"Permission error listing {path}: {str(e)}")
        except Exception as e:
            print(f"Error listing directory {path}: {str(e)}")

        return files

    def download_file(self, remote_path, local_path, size):
        """Download a single file with progress bar"""
        os.makedirs(
            os.path.dirname(local_path) if os.path.dirname(local_path) else ".",
            exist_ok=True,
        )

        ftp = self.connection_pool.get_connection()
        try:
            if self.url.path:
                ftp.cwd(self.url.path)

            with open(local_path, "wb") as f:
                with self._tqdm_lock:
                    progress = tqdm(
                        total=size,
                        unit="B",
                        unit_scale=True,
                        desc=os.path.basename(remote_path),
                        leave=False,
                    )

                def callback(data):
                    f.write(data)
                    with self._tqdm_lock:
                        progress.update(len(data))

                try:
                    ftp.retrbinary(f"RETR {remote_path}", callback, blocksize=8192)
                except Exception as e:
                    print(f"Error downloading {remote_path}: {str(e)}")
                    if os.path.exists(local_path):
                        os.remove(local_path)
                finally:
                    with self._tqdm_lock:
                        progress.close()
        finally:
            self.connection_pool.return_connection(ftp)

    def download_all(self, output_dir="."):
        """Download all files using thread pool"""
        try:
            print("Scanning for files...")
            files = self.list_files_recursively()
            print(f"Found {len(files)} files to download")

            with tqdm(
                total=len(files), desc="Overall Progress", unit="file"
            ) as overall_progress:
                with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                    futures = []
                    for file_info in files:
                        remote_path = file_info["path"]
                        local_path = os.path.join(output_dir, remote_path)

                        if (
                            os.path.exists(local_path)
                            and os.path.getsize(local_path) == file_info["size"]
                        ):
                            overall_progress.update(1)
                            continue

                        future = executor.submit(
                            self.download_file,
                            remote_path,
                            local_path,
                            file_info["size"],
                        )
                        future.add_done_callback(lambda p: overall_progress.update(1))
                        futures.append(future)

                    for future in futures:
                        future.result()

        finally:
            self.connection_pool.close_all()


if __name__ == "__main__":
    FTP_URL = "ftp://massive-ftp.ucsd.edu/v02/MSV000084900/"
    OUTPUT_DIR = "MSV000084900"
    downloader = FTPDownloader(FTP_URL)
    downloader.download_all(OUTPUT_DIR)

import os
import sys
from datetime import datetime
from Helpers import Helpers

class FileTracker:
    def __init__(self, input_dir, link_dir, max_batch_size):
        self.processed_files = [] #basecalled, .fastq
        self.input_dir = input_dir
        self.link_dir = link_dir
        self.unprocessed = []
        self.max_batch_size = max_batch_size
        self.reads_processed = 0
    
    """
    def insert_file(self, fname):
        if not fname in self.processed_files:
            self.processed_files.append(fname)
        else:
            print(str(datetime.now())+" :: FileTracker :: [insert_file] File "+fname+" was already processed, skipping!", file=sys.stderr) 
    """
                    
    def get_processed(self):
        return self.processed_files
        
    def get_unprocessed(self):
        all_from_dir = [filename for filename in os.listdir(self.input_dir) if os.path.splitext(filename)[1]==".fastq"]
        self.unprocessed = Helpers.alphanumsort([fastq for fastq in all_from_dir if fastq not in self.processed_files])
        self.unprocessed = self.unprocessed[:self.max_batch_size]
        return self.unprocessed
        
    def get_dir(self):
        return self.input_dir
        
    def create_links_unprocessed(self):
        count = 0
        self.unprocessed = self.get_unprocessed()
        for f in self.unprocessed:
            count += 1
            link_from = self.input_dir+f
            link_to = self.link_dir+f
            os.symlink(link_from, link_to)
        return count # num of new files
            
    def set_processed(self):
        #add new reads to count
        new_reads = 0
        for file in self.unprocessed:
            count = Helpers.count_reads_in_fastq(self.input_dir+file)
            self.reads_processed += count
            new_reads += count
        num_new_files = len(self.unprocessed)
        #set files as processed
        self.processed_files = self.processed_files + self.unprocessed
        self.unprocessed = []
        num_processed_files = len(self.processed_files)
        files_present = len([filename for filename in os.listdir(self.input_dir) if os.path.splitext(filename)[1]==".fastq"])
        print(str(datetime.now())+" :: FileTracker :: [set_processed] processed:", str(self.processed_files))
        print(str(datetime.now())+" :: FileTracker :: [set_processed] "+str(num_processed_files)+" files processed in total ("+str(num_new_files)+" from this iteration) out of "+str(files_present)+" files in the input directory ("+str((num_processed_files/files_present)*100)+" % of files processed)")
        print(str(datetime.now())+" :: FileTracker :: [set_processed] "+str(self.reads_processed)+" reads processed in total ("+str(new_reads)+" form this iteration)")
        
    def remove_links(self):
        for f in os.listdir(self.link_dir):
            os.unlink(self.link_dir+f)
            

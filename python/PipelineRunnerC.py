#pipeline that was previously in master with map to assembly loop
from PipelineRunner import PipelineRunner
import argparse
import os
import time
import sys
from datetime import datetime
from FileTracker import FileTracker
from Minimap import Minimap
from Miniasm import Miniasm
from ReadsSubsetMapToAssembly import ReadsSubsetMapToAssembly
from Helpers import Helpers
from Stats import Stats
from Quast import Quast
from threading import Thread
from threading import Condition
import copy

class PipelineRunnerC(PipelineRunner):
    def __init__(self, args):
        super().__init__(args)

    def minimap_to_asm(self):
        super().minimap_to_asm()
        if self.it>0:
            self.minimap.run_minimap_old_subsample_to_new_assembly()
        sys.stdout.flush()

    def subsample_reads(self):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [subsample_reads] SUBSAMPLE READS ")
        #subsample - new reads_subset object created each time.. 
        self.reads_subset = ReadsSubsetMapToAssembly(self.coverage_threshold, self.subsample_strategy, self.miniasm_fasta, self.it)
        self.reads_subset.parse_paf(self.minimap_paf)
        if self.it>0:
            self.reads_subset.parse_prev_subsample_paf(self.minimap_prev_subsample_paf) #so that we know coverage of new assembly by reads already selected for previous subsample
        self.reads_subset.subsamples_for_contigs() #subsample reads for each contig to covergae threshold
        self.reads_subset.write_fastq(self.link_dir, self.subsample_fastq, self.unmapped_fastq, self.out_dir+"subsample_filtered/subsampling_filtered_"+str(self.it)+".txt", self.out_dir+"subsample_filtered/subsampling_filtered_"+str(max(self.it-1, 0))+".txt")

    def map_to_assembly(self):
        waiting = 0
        while(self.ftracker.create_links_unprocessed() == 0): # new reads only, wait until new file exists or waiting threshold reached
            print(str(datetime.now())+" :: PipelineRunner :: [map_to_assembly] waiting for data..")
            waiting+=1
            if waiting >=self.max_waiting_count:
                print(str(datetime.now())+" :: PipelineRunner :: [map_to_assembly] no new data for "+str(self.max_waiting_count*5)+" seconds, exiting", flush=True)
                exit(0)
            time.sleep(self.wait_time)
        t1 = time.time()
        self.minimap.run_minimap_to_assembly()
        t2 = time.time()
        self.new_reads_subset = ReadsSubsetMapToAssembly(self.coverage_threshold, self.subsample_strategy, self.miniasm_fasta, self.it) # for new reads
        t3 = time.time()
        #self.new_reads_subset.parse_paf(self.minimap_paf) #NEW reads mapped to assembly
        self.new_reads_subset.parse_paf_for_map_to_assembly(self.minimap_paf)
        t4 = time.time()
        #if self.it == 1: #possible speed up, approximate coverage (i.e. from prev. it. may be ok)
        self.minimap.run_minimap_old_subsample_to_new_assembly() #run it again, because subsample has changed - can possibly be improved by running on new sbsampled reads only
        t5 = time.time()
        self.new_reads_subset.parse_prev_subsample_paf(self.minimap_prev_subsample_paf) #so that we know coverage of new assembly by reads already selected for previous subsample
        t6 = time.time()
        self.new_reads_subset.find_contig_overhangs()
        t7 = time.time()
        self.new_reads_subset.write_contig_overhangs_and_unmapped(self.link_dir, self.unmapped_fastq)  #APPEND
        t8 = time.time()
        self.ftracker.set_processed() #only interesting reads chosen from the files and stored in subsample.fastq
        self.ftracker.remove_links()
        self.log_times(str(self.it)+"_map_to_asm",t1,t2,t3,t4,t5,t6,t7,t8)
        

    def loop_map_to_assembly(self):
        for i in range(self.num_map_to_assembly_loops): # treba domysliet.. 
            self.map_to_assembly() 
            print(str(datetime.now())+" :: PipelineRunner :: [loop_map_to_assembly] num of reads in subsample.fastq:"+str(Helpers.count_reads_in_fastq(self.subsample_fastq)), flush=True)
            print(str(datetime.now())+" :: PipelineRunner :: [loop_map_to_assembly] num of reads in unmapped.fastq:"+str(Helpers.count_reads_in_fastq(self.unmapped_fastq)), flush=True)
            #os.system("cat "+self.subsample_fastq+" | wc -l")

    def main_loop(self):
        super().main_loop()
        self.last_it_start = datetime.now()
        
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [main_loop] MAPPING NEW READS TO ASSEMBLY [pre-filter] (num. loops: "+str(self.num_map_to_assembly_loops)+")", flush=True)
        self.loop_map_to_assembly()
        #enough new reads, run miniasm again 
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [main_loop] NEW ASSEMBLY", flush=True)
        t1,t2,t3,t4 = self.new_assembly(first=False)
        self.minimap_to_asm()
        t5 = time.time()
        self.subsample_reads()
        t6 = time.time()

        self.last_it_end = datetime.now()

        self.iteration_stats(self.it)

        t7 = time.time()
        self.clean_links_after_new_assembly(first=False)
        t8 = time.time()
        self.log_times(self.it,t1,t2,t3,t4,t5,t6,t7,t8)
        print(str(datetime.now())+" :: PipelineRunner :: [main_loop] num of reads in subsample.fastq:"+str(Helpers.count_reads_in_fastq(self.subsample_fastq)), flush=True)
        print(str(datetime.now())+" :: PipelineRunner :: [main_loop] num of reads in unmapped.fastq:"+str(Helpers.count_reads_in_fastq(self.unmapped_fastq)), flush=True)
        #os.system("cat "+self.subsample_fastq+" | wc -l")
        if self.is_last_assembly:
            print("\n"+str(datetime.now())+" :: PipelineRunner :: [new_assembly] last assembly done, EXITING")
            exit(0)
        self.it+=1
        print("\n\n\n")

    def link_files_for_new_assembly(self, first):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [link_files_for_new_assembly] LINKING FILES ")

        self.is_last_assembly=False 
        if first or self.num_map_to_assembly_loops==0: #######################
            waiting = 0
            while(self.ftracker.create_links_unprocessed() == 0): # new reads only, wait until new file exists or waiting threshold reached
                print(str(datetime.now())+" :: PipelineRunner :: [link_files_for_new_assembly] waiting for data..")
                waiting+=1
                if waiting >=self.max_waiting_count:
                    print(str(datetime.now())+" :: PipelineRunner :: [link_files_for_new_assembly] no new data for "+str(self.max_waiting_count*5)+" seconds, running last assembly and exiting", flush=True)
                    #exit(0)
                    self.is_last_assembly = True
                    break
                time.sleep(self.wait_time)
        
        if not first:
            self.link_previous(self.subsample_fastq, self.unmapped_fastq) #link subsampled and interesting new reads 

    def clean_links_after_new_assembly(self, first = False):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [clean_links_after_new_assembly] CLEAN LINKS, SET PROCESSED FILES ")
        #set files as processed and remove links at the end of iteration
        if first or self.num_map_to_assembly_loops==0: ###################
            self.ftracker.set_processed()
        self.ftracker.remove_links()
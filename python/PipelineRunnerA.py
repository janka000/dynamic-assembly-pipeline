#pipeline that was previously in master
from PipelineRunner import PipelineRunner
import argparse
import os
import time
import sys
from datetime import datetime
from FileTracker import FileTracker
from Minimap import Minimap
from Miniasm import Miniasm
from ReadsSubset import ReadsSubset
from Helpers import Helpers
from Stats import Stats
from Quast import Quast
from threading import Thread
from threading import Condition
import copy

class PipelineRunnerA(PipelineRunner):
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
        self.reads_subset = ReadsSubset(self.coverage_threshold, self.subsample_strategy, self.miniasm_fasta, self.it)
        self.reads_subset.parse_paf(self.minimap_paf)
        if self.it>0:
            self.reads_subset.parse_prev_subsample_paf(self.minimap_prev_subsample_paf) #so that we know coverage of new assembly by reads already selected for previous subsample
        self.reads_subset.subsamples_for_contigs() #subsample reads for each contig to covergae threshold
        self.reads_subset.write_fastq(self.link_dir, self.subsample_fastq, self.unmapped_fastq, self.out_dir+"subsample_filtered/subsampling_filtered_"+str(self.it)+".txt", self.out_dir+"subsample_filtered/subsampling_filtered_"+str(max(self.it-1, 0))+".txt")

    def main_loop(self):
        super().main_loop()
        self.last_it_start = datetime.now()
        
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [main_loop] MAPPING NEW READS TO ASSEMBLY [pre-filter] (num. loops: "+str(self.num_map_to_assembly_loops)+")", flush=True)
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

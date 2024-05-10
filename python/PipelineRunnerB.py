#"unbiased subsample" with heuristics
from PipelineRunner import PipelineRunner
import argparse
import os
import time
import sys
from datetime import datetime
from FileTracker import FileTracker
from Minimap import Minimap
from Miniasm import Miniasm
from ReadsSubsetUnbiased import ReadsSubsetB, ReadsSubsetBs, ReadsSubsetLongReadB, ReadsSubsetLongReadBs
from Helpers import Helpers
from Stats import Stats
from Quast import Quast
from threading import Thread
from threading import Condition
import copy

class PipelineRunnerB(PipelineRunner):
    def __init__(self, args):
        super().__init__(args)
        self.reads_prob = self.out_dir+"reads_prob.tsv"

    def parse_args(self,args):
        super().parse_args(args)
        self.reads_prob = self.out_dir+"reads_prob.tsv"
        self.simple_subsample = args.unbiased_simple
        self.rlbias = args.prefer_longer_reads

    def makedirs(self):
        super().makedirs()
        if not os.path.exists(self.reads_prob):
            with open(self.reads_prob, "w+"):
                pass

    def minimap_to_asm(self):
        print("\n"+str(datetime.now())+" :: PipelineRunnerB :: [minimap_to_asm] MINIMAP - MAP READS TO ASSEMBLY ")
        self.minimap.run_minimap_to_assembly(remove_subsample=False)
        sys.stdout.flush()

    def subsample_reads(self):
        print("\n"+str(datetime.now())+" :: PipelineRunnerB :: [subsample_reads] SUBSAMPLE READS ")
        #subsample - new reads_subset object created each time.. 
        if not self.simple_subsample and not self.rlbias:
            self.reads_subset = ReadsSubsetB(self.coverage_threshold, self.subsample_strategy, self.miniasm_fasta, self.reads_prob)
        elif self.simple_subsample and not self.rlbias: #simple
            self.reads_subset = ReadsSubsetBs(self.coverage_threshold, self.subsample_strategy, self.miniasm_fasta, self.reads_prob)
        elif self.rlbias and not self.simple_subsample: #take loner reads first
            self.reads_subset = ReadsSubsetLongReadB(self.coverage_threshold, self.subsample_strategy, self.miniasm_fasta, self.reads_prob)
        else: # rlbias and simple
            self.reads_subset = ReadsSubsetLongReadBs(self.coverage_threshold, self.subsample_strategy, self.miniasm_fasta, self.reads_prob)
        self.reads_subset.subsample(self.minimap_paf,self.link_dir, self.subsample_fastq, self.unmapped_fastq, self.out_dir+"subsample_filtered/subsampling_filtered_"+str(self.it)+".txt")

    def main_loop(self):
        super().main_loop()
        print("\n"+str(datetime.now())+" :: PipelineRunnerB :: [main_loop] NEW ASSEMBLY", flush=True)
        
        self.last_it_start = datetime.now()
        
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
        print(str(datetime.now())+" :: PipelineRunnerB :: [main_loop] num of reads in subsample.fastq:"+str(Helpers.count_reads_in_fastq(self.subsample_fastq)), flush=True)
        print(str(datetime.now())+" :: PipelineRunnerB :: [main_loop] num of reads in unmapped.fastq:"+str(Helpers.count_reads_in_fastq(self.unmapped_fastq)), flush=True)
        #os.system("cat "+self.subsample_fastq+" | wc -l")
        if self.is_last_assembly:
            print("\n"+str(datetime.now())+" :: PipelineRunnerB :: [new_assembly] last assembly done, EXITING")
            exit(0)
        self.it+=1
        print("\n\n\n")
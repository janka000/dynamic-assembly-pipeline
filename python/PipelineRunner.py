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
from StatsThreadSafe import StatsThreadSafe
from FastStats import FastStats
#from Quast import Quast
from threading import Thread
from threading import Lock
import threading
import copy
import shutil
import time
import traceback

class PipelineRunner:

     def __init__(self,args):
         self.it = 0
         self.ftracker = None
         self.minimap = None
         self.miniasm = None
         #self.quast = None
         self.reads_subset = None
         self.new_reads_subset = None
         self.last_it_start = datetime.now()
         self.last_it_end = datetime.now()
         self.matplotlib_lock = Lock() #matplotlib is not thread-safe, mutal exclusion necessary
         self.prev_stats_thread = None
         self.parallel_stats = False
         self.max_waiting_stats_threads = 5
         self.i_want_bam_alignments_in_stats = False #have to be patient (+ ~10min per iteration)


         print(str(datetime.now())+" :: PipelineRunner :: [init] parse args")
         self.parse_args(args)
         self.print_params()
         print(str(datetime.now())+" :: PipelineRunner :: [run] makedirs")
         self.makedirs()
         print(str(datetime.now())+" :: PipelineRunner :: [run] setup")
         self.setup()
         sys.stdout.flush()
     
     def parse_args(self, args): 
         #dirs
         if args.in_dir!=None:
            self.input_dir = args.in_dir
         if args.out_dir!=None:    
            self.out_dir = args.out_dir
            self.link_dir = args.out_dir+"link/"      
            self.stats_dir = args.out_dir+"stats/"
         self.base_dir = args.base_dir

         #pipeline
         self.subsample_fastq = self.out_dir+"subsample.fastq" #subsampled reads for contigs - sholud not be subsampled again
         self.unmapped_fastq = self.out_dir+"unmapped.fastq" #reads that were not mapped yet - to be subsampled later
         self.coverage_threshold = int(args.coverage_threshold)
         self.num_map_to_assembly_loops = int(args.mapping_loops)
         self.num_main_loops = int(args.max_loops)
         self.max_waiting_count = int(args.max_waiting_count)
         self.wait_time = int(args.wait_time)
         self.max_batch_size = int(args.batch_size)
         self.subsample_strategy = args.subsample_strategy #probabilistic/take_all

         #minimap
         self.minimap_paf = self.out_dir+"minimap_assembly.paf.gz"
         self.minimap_threads = int(args.minimap_threads)
         self.minimap_K = args.minimap_K
         self.minimap_prev_subsample_paf = self.out_dir+"minimap_assembly_prev_it_subsample.paf.gz"

         #miniasm
         self.miniasm_min_unitig_reads = int(args.miniasm_min_unitig_reads)
         self.miniasm_min_coverage = int(args.miniasm_min_coverage)
         self.miniasm_gfa = self.out_dir+"miniasm.gfa"
         self.miniasm_fasta = self.out_dir+"assembly.fasta"

         #debug/track
         self.track_finalassembly =args.track_finalassembly
         self.have_assembly = args.have_assembly

         self.min_disk_space = args.min_disk_space
         self.time_file = self.out_dir+"times.tsv"
     
     def print_params(self):
         print("------------------------------------------------------------------------------")
         print("pipeline parameters:")
         print("------------------------------------------------------------------------------")
         print("in dir: ", self.input_dir)
         print("out dir:", self.out_dir)
         print("link dir: ", self.link_dir)
         print("stats dir: ", self.stats_dir)
         print("------------------------------------------------------------------------------")
         print("subsample fastq:", self.subsample_fastq)
         print("subsampling coverage threshold:", self.coverage_threshold)
         print("num map to assembly loops:", self.num_map_to_assembly_loops)
         print("max loops:", self.num_main_loops)
         print("wait time (check for new files interval):", self.wait_time)
         print("max waiting count (pipeline will terminate after no new files for [wait time interval]*[this number]):", self.max_waiting_count)
         print("max batch size:", self.max_batch_size)

         print("------------------------------------------------------------------------------")
         print("minimap")
         print("------------------------------------------------------------------------------")
         print("minimap threads:", self.minimap_threads)
         print("minimap K:", self.minimap_K)
         print("minimap paf:", self.minimap_paf)

         print("------------------------------------------------------------------------------")
         print("miniasm")
         print("------------------------------------------------------------------------------")
         print("min unitig reads:", self.miniasm_min_unitig_reads)
         print("min coverage:", self.miniasm_min_coverage)
         print("gfa path:", self.miniasm_gfa)
         print("fasta path:", self.miniasm_fasta)
         print("-------------------------------------------------------------------------------")

         print("debug")
         print("-------------------------------------------------------------------------------")
         #print("reference for quast:", self.reference)
         print("assembly provided? :", self.have_assembly)
         print("final assembly:", self.track_finalassembly)
         print("remove non-essential files at the end of run:", self.min_disk_space)
         print("-------------------------------------------------------------------------------", flush=True)

     def makedirs(self):
        if not os.path.exists(self.input_dir):
            os.makedirs(self.input_dir)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if not os.path.exists(self.out_dir+"subsample_filtered/"):
            os.makedirs(self.out_dir+"subsample_filtered/")
        if not os.path.exists(self.out_dir+"assemblies/"):
            os.makedirs(self.out_dir+"assemblies/")
        if not os.path.exists(self.out_dir+"miniasm_filtered/"):
            os.makedirs(self.out_dir+"miniasm_filtered/")
        if not os.path.exists(self.link_dir):
            os.makedirs(self.link_dir)
        if not os.path.exists(self.subsample_fastq):
            with open(self.subsample_fastq, "w+"):
                pass
        #stats
        if not os.path.exists(self.stats_dir):
            os.makedirs(self.stats_dir)
        if not os.path.exists(self.stats_dir+"with_reference/"):
            os.makedirs(self.stats_dir+"with_reference/")
        if not os.path.exists(self.stats_dir+"without_reference/"):
            os.makedirs(self.stats_dir+"without_reference/")
        #without reference
        if not os.path.exists(self.stats_dir+"without_reference/coverage_graph_per_it/"):
            os.makedirs(self.stats_dir+"without_reference/coverage_graph_per_it/")
        if not os.path.exists(self.stats_dir+"without_reference/contig_lengths_graph_per_it/"):
            os.makedirs(self.stats_dir+"without_reference/contig_lengths_graph_per_it/")
        if not os.path.exists(self.stats_dir+"without_reference/prev_asm_to_this/"):
            os.makedirs(self.stats_dir+"without_reference/prev_asm_to_this/")
        if not os.path.exists(self.stats_dir+"without_reference/json_summary_per_it/"):
            os.makedirs(self.stats_dir+"without_reference/json_summary_per_it/")
        if not os.path.exists(self.stats_dir+"html/"):
            os.makedirs(self.stats_dir+"html/")
        #with reference
        if self.have_assembly:
            if not os.path.exists(self.stats_dir+"with_reference/assembly_to_ref_dotplots/"):
                os.makedirs(self.stats_dir+"with_reference/assembly_to_ref_dotplots/")
            if not os.path.exists(self.stats_dir+"with_reference/bam_subsample_reads_to_ref/"):
                os.makedirs(self.stats_dir+"with_reference/bam_subsample_reads_to_ref/")
            if not os.path.exists(self.stats_dir+"with_reference/assembly_to_ref_paf/"):
                os.makedirs(self.stats_dir+"with_reference/assembly_to_ref_paf/")
            if not os.path.exists(self.stats_dir+"with_reference/reads_to_ref_coverage/"):
                os.makedirs(self.stats_dir+"with_reference/reads_to_ref_coverage/")
            #if not os.path.exists(self.stats_dir+"with_reference/assembly_to_ref_bam/"):
            #    os.makedirs(self.stats_dir+"with_reference/assembly_to_ref_bam/")
            
            #igv track
            if not os.path.exists(self.stats_dir+"with_reference/igv/"):
                os.makedirs(self.stats_dir+"with_reference/igv/")
            if not os.path.exists(self.stats_dir+"with_reference/igv/subsampling_bw_per_it/"):
                os.makedirs(self.stats_dir+"with_reference/igv/subsampling_bw_per_it/")
            if not os.path.exists(self.stats_dir+"with_reference/igv/asm_to_asm_map_bed_per_it/"):
                os.makedirs(self.stats_dir+"with_reference/igv/asm_to_asm_map_bed_per_it/")
            
            
     def setup(self):
        self.ftracker = FileTracker(self.input_dir, self.link_dir, self.max_batch_size)
        self.minimap = Minimap(self.out_dir, self.link_dir, self.minimap_threads, self.minimap_K, self.ftracker)
        self.miniasm = Miniasm(self.out_dir, self.link_dir, self.base_dir, self.miniasm_min_unitig_reads, self.miniasm_min_coverage)
            #self.quast = Quast(self.stats_dir, self.out_dir, reference=self.track_finalassembly)
        #else:
        #    self.quast = Quast(self.stats_dir, self.out_dir, reference=None)

     def minimap_ava(self):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [minimap_ava] MINIMAP ALL-VS-ALL ", flush=True)
        self.minimap.run_minimap_all_vs_all()
        sys.stdout.flush()

     def miniasm_new_assembly(self):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [miniasm_new_assembly] MINIASM - NEW ASSEMBLY ",flush=True)
        self.miniasm.run_miniasm(self.it)
        sys.stdout.flush()
        self.miniasm.gfa_to_fasta()
        self.miniasm.log_assembly(self.it)
        sys.stdout.flush()

     def minimap_to_asm(self):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [minimap_to_asm] MINIMAP - MAP READS TO ASSEMBLY ",flush=True)
        self.minimap.run_minimap_to_assembly()
        sys.stdout.flush()

     def iteration_stats(self, it):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [iteration_stats] CALCULATE STATS, DRAW COVERAGE GRAPHS (logging to "+self.stats_dir+"stats_"+str(it)+".log"+")",flush=True)
        stats = None
        if self.i_want_bam_alignments_in_stats:
            stats = Stats(self.input_dir, self.out_dir, self.stats_dir, self.subsample_fastq, self.coverage_threshold, self.minimap, self.base_dir+"miniasm/", reference=self.track_finalassembly)
            #stats = StatsThreadSafe(self.input_dir, self.out_dir, self.stats_dir, self.subsample_fastq, self.coverage_threshold, self.minimap, self.base_dir+"miniasm/", reference=self.track_finalassembly)
        else:
            stats = FastStats(self.input_dir, self.out_dir, self.stats_dir, self.subsample_fastq, self.coverage_threshold, self.minimap, self.base_dir+"miniasm/", reference=self.track_finalassembly)
        stats.set_iteration(it)
        stats.start_end_times(self.last_it_start, self.last_it_end)
        processed_files = copy.deepcopy(self.ftracker.processed_files+self.ftracker.unprocessed) #v processed su subory z predoslych iteracii, v unprocessed z tejto
        num_unmapped_reads = Helpers.count_reads_in_fastq(self.unmapped_fastq)
        num_threads = len(threading.enumerate())
        if self.parallel_stats:
            print(str(datetime.now())+" :: PipelineRunner :: [iteration_stats] running stats in parallel, BE AWARE, this may require a lot of disk space in case there are too threads (they will wait, after threshold for max number of running threads is exceeded in Stats class), therefore, there is a limit set to self.max_waiting_stats_threads == "+str(self.max_waiting_stats_threads)+" - if there are more threads, main thread will wait for the previous stats thread to finish", flush=True)
            if num_threads > self.max_waiting_stats_threads: #if there are too many stats threads, wait for the previos one to finish, otherwise there will be too many subsample copies
                print(str(datetime.now())+" :: PipelineRunner :: [iteration_stats] too many stats threads running ("+num_threads+"), waiting for the last one to finish before continuing", flush=True)
                self.prev_stats_thread.join()
            thread = Thread(target=stats.stats_thread_task, args=(it,processed_files, num_unmapped_reads, self.prev_stats_thread, num_threads, self.matplotlib_lock))
            thread.start()
            self.prev_stats_thread = thread
        #thread.join() #do not run stats in parallel, option 1 (separate thread, but main waiting)
        else:
            print(str(datetime.now())+" :: PipelineRunner :: [iteration_stats] running parellel stats sequentially with main - waiting until all stats are finished for this iteration", flush=True)
            stats.stats_thread_task(it,processed_files,num_unmapped_reads, None, 1, self.matplotlib_lock)
            stats_log = self.stats_dir+"stats_"+str(it)+".log"
            with open(stats_log,"r") as slog:
                print(slog.read()) #print log also to stdout

     def subsample_reads(self):
        raise NotImplementedError
        
     def link_previous(self,subsample_fastq, unmapped_fastq):
        link_from = subsample_fastq
        link_to = self.link_dir+"subsample.fastq"
        os.symlink(link_from, link_to)

        link_from = unmapped_fastq
        link_to = self.link_dir+"unmapped.fastq"
        os.symlink(link_from, link_to)

     def link_files_for_new_assembly(self, first):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [link_files_for_new_assembly] LINKING FILES ",flush=True)

        self.is_last_assembly=False 
        waiting = 0
        while(self.ftracker.create_links_unprocessed() == 0): # new reads only, wait until new file exists or waiting threshold reached
            print(str(datetime.now())+" :: PipelineRunner :: [link_files_for_new_assembly] waiting for data..",flush=True)
            waiting+=1
            if waiting >=self.max_waiting_count:
                print(str(datetime.now())+" :: PipelineRunner :: [link_files_for_new_assembly] no new data for "+str(self.max_waiting_count*5)+" seconds, running last assembly and exiting", flush=True)
                #exit(0)
                self.is_last_assembly = True
                break
            time.sleep(self.wait_time)
        
        if not first:
            self.link_previous(self.subsample_fastq, self.unmapped_fastq) #link subsampled and interesting new reads 

     def clean_links_after_new_assembly(self, first=False):
        print("\n"+str(datetime.now())+" :: PipelineRunner :: [clean_links_after_new_assembly] CLEAN LINKS, SET PROCESSED FILES ", flush=True)
        #set files as processed and remove links at the end of iteration
        self.ftracker.set_processed()
        self.ftracker.remove_links()

     def new_assembly(self, first=False):
        t1 = time.time()
        self.link_files_for_new_assembly(first)
        t2 = time.time()
        self.minimap_ava()
        t3 = time.time()
        self.miniasm_new_assembly()
        t4 = time.time()
        return t1,t2,t3,t4
         
     def main_loop(self):
        print("\n_______________________________________________________________________________")
        print("-------------------------------------------------------------------------------")
        print(str(datetime.now())+" :: PipelineRunner :: [main_loop]: ITERATION "+str(self.it))
        print("-------------------------------------------------------------------------------", flush=True)
        #each version 
        pass

     def start(self):
        t1,t2,t3,t4 = self.new_assembly(first=True)
        self.minimap_to_asm()
        t5 = time.time()
        self.subsample_reads()
        t6 = time.time()
        self.iteration_stats(self.it)
        t7 = time.time()
        self.clean_links_after_new_assembly(first=True)
        t8 = time.time()
        self.log_times(0,t1,t2,t3,t4,t5,t6,t7,t8)
    
     def log_times(self,it,t1,t2,t3,t4,t5,t6,t7,t8):
         with open(self.time_file,"a+") as tf:
             print(str(it)+"\t"+str(t1)+"\t"+str(t2)+"\t"+str(t3)+"\t"+str(t4)+"\t"+str(t5)+"\t"+str(t6)+"\t"+str(t7)+"\t"+str(t8), file=tf, flush=True)
        
     def run(self):
        #first round
        start = time.time()
        print("\n_______________________________________________________________________________")
        print("-------------------------------------------------------------------------------")
        print(str(datetime.now())+" :: PipelineRunner :: [run] first assembly - ITERATION 0", flush=True)
        sys.stdout.flush()
        self.start()
        print(str(datetime.now())+" :: PipelineRunner :: [run] num of reads in subsample.fastq:"+str(Helpers.count_reads_in_fastq(self.subsample_fastq)), flush=True)
        print(str(datetime.now())+" :: PipelineRunner :: [run] num of reads in unmapped.fastq:"+str(Helpers.count_reads_in_fastq(self.unmapped_fastq)), flush=True)
        self.it+=1
        end = time.time()
        with open(self.out_dir+"iteration_times.tsv","a+") as it_times:
            print(str(0)+"\t"+str(start)+"\t"+str(end), file = it_times) 
        

        for i in range(self.num_main_loops): # treba domysliet..
            start = time.time()
            self.main_loop()
            end = time.time()
            with open(self.out_dir+"iteration_times.tsv","a+") as it_times:
                print(str(self.it-1)+"\t"+str(start)+"\t"+str(end), file = it_times) 

     def clean_at_the_end(self):
        if self.min_disk_space:
            to_remove = []
            to_remove.append(self.out_dir+"minimap_assembly_prev_it_subsample.paf.gz")
            #to_remove.append(self.stats_dir+"with_reference/assembly_to_ref_paf/")
            to_remove.append(self.stats_dir+"with_reference/assembly_to_ref_bam/")
            to_remove.append(self.stats_dir+"with_reference/igv/miniasm_sub_bw_per_it/")
            #to_remove.append(self.stats_dir+"with_reference/igv/subsampling_bw_per_it/")
            to_remove.append(self.out_dir+"minimap_ava.paf.gz")
            to_remove.append(self.out_dir+"minimap_assembly.paf.gz")
            to_remove.append(self.stats_dir+"with_reference/bam_subsample_reads_to_ref/")
            for rem in to_remove:
                if os.path.exists(rem):
                    if os.path.isfile(rem):
                        os.unlink(rem)
                    elif os.path.isdir(rem):
                        shutil.rmtree(rem)

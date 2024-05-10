import pysam
import json
import os
import time
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.patches as patches
from Helpers import Helpers
from Minimap import Minimap
import traceback
import random
from threading import Thread

class Stats:
    """
    class for computing stats for completed assembly
    """
    def __init__(self, in_dir, out_dir, stats_dir, subsample_fastq, coverage_threshold, minimap, miniasm_dir, reference = None):
        self.subsample_fastq = subsample_fastq
        self.sum_coverage_per_pos = None
        self.assembly_length = None
        self.num_reads = None
        self.stats_dir = stats_dir
        self.out_dir = out_dir
        self.contigs = [] #contig names
        self.dict = {}
        self.c_lengths = [] #dlzky contigov
        self.num_reads_per_file = None
        self.in_dir = in_dir
        self.paf_figures = 0
        self.coverage_threshold = coverage_threshold
        self.minimap = minimap
        self.miniasm_dir = miniasm_dir
        self.lock_log = stats_dir+"lock_log.txt"
        if reference is not None:
            self.reference = reference
            self.have_ref = True
        else:
            self.have_ref = False

        self.max_stats_threads = 5
        self.debug_or_have_time = False #nice stats plots, but for patient people only
        
    def set_iteration(self, it):
        self.it = it

        #files that exist:
        self.assembly_fasta = self.out_dir+"assemblies/assembly_"+str(it)+".fasta" #created by Miniasm class
        self.prev_it_assembly = self.out_dir+"assemblies/assembly_"+str(it-1)+".fasta"

        #have to be created now:
        self.subsample_fastq_copy = self.stats_dir+"subsample_stats_copy"+str(it)+".fastq"
        os.system("cp "+self.subsample_fastq+" "+self.subsample_fastq_copy)
        
        #will be created later - target:
        self.stats_file = self.stats_dir+"without_reference/json_summary_per_it/stats_"+str(it)+".json"
        self.bigwig = self.stats_dir+"with_reference/igv/subsampling_bw_per_it/subsample_"+str(it)+".bw"
        self.assembly_to_ref_bed = self.stats_dir+"with_reference/igv/asm_to_asm_map_bed_per_it/asm_to_asm_"+str(it)+".bed"
        self.assembly_to_ref_minidot_dotplot = self.stats_dir+"with_reference/assembly_to_ref_dotplots/dotplot_"+str(it)+".ps"
        self.minidot_dotplot_pdf = self.stats_dir+"with_reference/assembly_to_ref_dotplots/dotplot_"+str(it)+".pdf"
        self.subsample_to_ref_coverage_png = self.stats_dir+"with_reference/reads_to_ref_coverage/ref_cov_"+str(it)+".png"
        self.prev_asm_to_this_dotplot_png = self.stats_dir+"without_reference/prev_asm_to_this/prev_asm_to_this_dotplot_"+str(it)+".png"
        self.noref_contig_coverage = self.stats_dir+"without_reference/coverage_graph_per_it/coverage_"+str(it)+".png" 
        self.noref_contig_coverage_bar = self.stats_dir+"without_reference/coverage_graph_per_it/coverage_bar_"+str(it)+".png"
        self.contig_lengths_graph = self.stats_dir+"without_reference/contig_lengths_graph_per_it/contig_lengths_"+str(it)+".png" 
        self.assembly_to_ref_dotplot = self.stats_dir+"with_reference/assembly_to_ref_dotplots/paf_bam_dotplot_"+str(it)

        #needed (may be deleted later)
        #noref
        self.subsample_to_assembly_bam = self.stats_dir+"stats_subsample_to_assembly_"+str(it)+".bam"
        self.paf_asm_to_asm = self.stats_dir+"without_reference/prev_asm_to_this/map_"+str(it)+".paf"
        #ref
        self.genome_sizes = self.out_dir+"final_assembly.sizes"
        self.subsample_to_ref_bam = self.stats_dir+"with_reference/bam_subsample_reads_to_ref/subsample_to_ref_"+str(it)+".bam"
        self.assembly_to_ref_paf = self.stats_dir+"with_reference/assembly_to_ref_paf/align_to_ref"+str(it)+".paf"
        
        #log file
        self.stats_log = self.stats_dir+"stats_"+str(it)+".log"
        self.log_file = open(self.stats_log, "w+")
        print("------------------------- stats for it. "+str(it)+"-----------------------------", flush=True, file=self.log_file)
        print("\n"+str(datetime.now())+" :: Stats :: [iteration_stats] CALCULATE STATS, DRAW COVERAGE GRAPHS (logging to "+self.stats_dir+"stats_"+str(it)+".log"+")", flush=True, file=self.log_file)

    
    def stats_thread_task_no_ref(self,it,processed_files, matplotlib_lock):
        try:
            stats_log = self.stats_dir+"noref_stats_"+str(it)+".log"
            stats_log_f = open(stats_log, "a+")
            
            #can run in parallel with the rest of pipeline        
            self.create_no_ref_necessary_files(stats_log_f)
            
            if self.debug_or_have_time:
                print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_no_ref] drawing prev. asm to this dotplot.", file=stats_log_f)
                self.prev_asm_to_this_dotplot(matplotlib_lock)
            
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_no_ref] building dict.", flush=True, file=stats_log_f)
            self.build_dict(stats_log_f)
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_no_ref] setting processed files.", flush=True, file=stats_log_f)
            self.files_processed(processed_files)
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_no_ref] drawing contig lengths.", flush=True, file=stats_log_f)
            self.draw_contig_lengths(matplotlib_lock)
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_no_ref] drawing coverage graphs for contigs.", flush=True, file=stats_log_f)
            self.draw_coverage(matplotlib_lock,stats_log_f)
        except Exception as e: 
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_no_ref] error in noref stats thread", file=stats_log_f)
            print(e, file=stats_log_f)
            traceback.print_exc()
        finally:
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_no_ref] stats for it. "+str(it)+" DONE.", flush=True, file=stats_log_f)
            print("\n ------------[stats_thread_task_no_ref] end of no ref. stats for it. "+str(it)+" ---------", flush=True, file=stats_log_f)
            stats_log_f.close()
            if matplotlib_lock.locked(): #if anything went wrong.. but it should be released earlier even in cas of failure
                matplotlib_lock.release()
            #end of parallel block
                
    def stats_thread_task_ref(self,it,matplotlib_lock):
        try:
            stats_log = self.stats_dir+"ref_stats_"+str(it)+".log"
            stats_log_f = open(stats_log, "a+")

            self.create_ref_necessary_files(stats_log_f)
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_ref] minidot assembly to ref.", flush=True, file=stats_log_f)
            self.minidot_asm_to_ref()
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_ref] drawing ref subsample coverage.", flush=True, file=stats_log_f)
            self.ref_subsample_coverage(matplotlib_lock,stats_log_f)
            if self.debug_or_have_time:
                print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_ref] drawing assembly to ref iteration dotplot.",flush=True, file=stats_log_f)
                self.assembly_to_ref_iteration_dotplot(it, matplotlib_lock, stats_log_f) 
                print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_ref] creating bigwig from bam.", flush=True, file=stats_log_f)
                #self.create_reads_to_ref_bigwig()
        except Exception as e:
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_ref] error in ref stats thread", file=stats_log_f)
            print(e, file=stats_log_f)
            traceback.print_exc()
        finally:
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task_ref] ref stats for it. "+str(it)+" DONE.", flush=True, file=stats_log_f)
            print("\n ------------ [stats_thread_task_ref] end of ref. stats for it. "+str(it)+" ---------", flush=True, file=stats_log_f)
            stats_log_f.close()
            if matplotlib_lock.locked(): #if anything went wrong.. but it should be released earlier even in cas of failure
                matplotlib_lock.release()
        
    def stats_thread_task(self,it,processed_files, num_unmapped_reads, prev_stats_thread, num_threads, matplotlib_lock):

            stats_log = self.stats_dir+"noref_stats_"+str(it)+".log"
            stats_log_f = open(stats_log, "a+")
            self.num_unmapped_reads = num_unmapped_reads
            #wait for prev stats thread to finish if there are too many threads (so that there are not too many stats threads and they are not executed in random order)
            if num_threads >= self.max_stats_threads:#prev_stats_thread != None:
                if prev_stats_thread.is_alive():
                    print(str(datetime.now())+" :: Stats :: [stats_thread_task] joining previous stats thread... ", flush=True, file=stats_log_f)
                    prev_stats_thread.join()
                    print(str(datetime.now())+" :: Stats :: [stats_thread_task] previous stats thread joined. ", flush=True, file=stats_log_f)
            stats_log_f.close()

            thread_no_ref = Thread(target=self.stats_thread_task_no_ref, args=(it,processed_files, matplotlib_lock))
            thread_no_ref.start()

            if self.have_ref: 
                thread_ref = Thread(target=self.stats_thread_task_ref, args=(it,matplotlib_lock))
                thread_ref.start()

            #wait for the threads to finish
            thread_no_ref.join()
            thread_ref.join() 

            log = self.stats_dir+"stats_"+str(it)+".log"
            with open(log,"a+") as l:
                print("STATS with and without reference for each iteartion run in PARALLEL, this is merged log", file=l)
            no_ref_log = self.stats_dir+"noref_stats_"+str(it)+".log"
            ref_log = self.stats_dir+"ref_stats_"+str(it)+".log"
            tmp_log = self.stats_dir+"tmp_stats_"+str(it)+".log"
            os.system("cat "+log+" "+no_ref_log+" "+ref_log+" > "+tmp_log)
            os.system("mv "+tmp_log+" "+log)
            os.unlink(no_ref_log)
            os.unlink(ref_log)

            stats_log_f = open(log, "a+")
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task]  ref and noref threads joined.", flush=True, file=stats_log_f)
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task]  printing stats.", flush=True, file=stats_log_f)
            self.print_stats()
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task]  printing HTML.", flush=True, file=stats_log_f)
            self.print_html() 
            print("\n"+str(datetime.now())+" :: Stats :: [stats_thread_task]  cleaning temporary files.", flush=True, file=stats_log_f)
            self.clean_files()
            print("\n ------------ [stats_thread_task] end of stats for it. "+str(it)+" ---------", flush=True, file=stats_log_f)

            stats_log_f.close()
    
    def create_necessary_files(self):
        self.create_no_ref_necessary_files(self.log_file)
        if self.have_ref:
            self.create_ref_necessary_files(self.log_file)
        print("\n"+str(datetime.now())+" :: Stats :: [create_necessary_files] all bam and paf files created ", flush=True, file=self.log_file)

    def create_no_ref_necessary_files(self,log):
        print("\n"+str(datetime.now())+" :: Stats :: [create_necessary_files] creating files for stats without reference ", flush=True, file=log)
        def create_subsample_to_asm_bam(assembly,subsample):
            self.minimap.map_reads_to_assembly(subsample,assembly,self.subsample_to_assembly_bam)
            os.system("samtools index "+self.subsample_to_assembly_bam)
        def create_paf_asm_to_asm(prev_asm,this_asm):
            if self.it > 0: #previous does not exist otherwise
                self.minimap.asm_to_asm_paf(this_asm, prev_asm, self.paf_asm_to_asm)
        print(str(datetime.now())+" :: Stats :: [create_no_ref_necessary_files] creating subsample to assembly bam ", flush=True, file=log)
        create_subsample_to_asm_bam(self.assembly_fasta,self.subsample_fastq_copy)
        print(str(datetime.now())+" :: Stats :: [create_no_ref_necessary_files] creating assembly to assembly paf ", flush=True, file=log)
        create_paf_asm_to_asm(self.prev_it_assembly, self.assembly_fasta)

    def create_ref_necessary_files(self,log):
        print("\n"+str(datetime.now())+" :: Stats :: [create_necessary_files] creating files for stats with reference ", flush=True, file=log)
        def create_subsample_to_ref_bam(ref,subsample):
            self.minimap.map_reads_to_assembly(subsample,ref,self.subsample_to_ref_bam)
            os.system("samtools index "+self.subsample_to_ref_bam)
        def create_assembly_to_ref_paf(assembly, ref):
            self.minimap.asm_to_asm_paf(ref,assembly,self.assembly_to_ref_paf)
        def create_assembly_to_ref_bed(paf,bed):
            awk='{print $6 "\t" $8 "\t" $9 "\t" $1 ":l:" $2 ":eoffset:" ($2-$4) "\t" $10 "\t" $5}'
            return os.system("cat "+paf+" | awk -F '\t' '"+awk+"' > "+bed)
        def reference_stats(ref):
            if not os.path.exists(self.genome_sizes):
                os.system("faSize -detailed "+ref+" > "+self.genome_sizes)
        print(str(datetime.now())+" :: Stats :: [create_ref_necessary_files] creating genome.sizes ", flush=True, file=log)
        reference_stats(self.reference)
        print(str(datetime.now())+" :: Stats :: [create_ref_necessary_files] creating subsample to ref bam ", flush=True, file=log)
        create_subsample_to_ref_bam(self.reference, self.subsample_fastq_copy)
        print(str(datetime.now())+" :: Stats :: [create_ref_necessary_files] creating assembly to ref paf ", flush=True, file=log)
        create_assembly_to_ref_paf(self.assembly_fasta, self.reference)
        print(str(datetime.now())+" :: Stats :: [create_ref_necessary_files] creating assembly to ref bed (from paf) ", flush=True, file=log)
        create_assembly_to_ref_bed(self.assembly_to_ref_paf, self.assembly_to_ref_bed)


###############################################################################
#####        without reference
###############################################################################
        
    def prev_asm_to_this_dotplot(self, lock):
        prev_it = self.it - 1
        if prev_it < 0:
            return 
        self.dotplot_with_coverage(self.paf_asm_to_asm, self.coverage_threshold, self.prev_asm_to_this_dotplot_png, lock, self.lock_log, self.it, self.log_file, bam_reads_to_asm=self.subsample_to_assembly_bam) 

    def alignment_stats(self):
        alignment = pysam.AlignmentFile(self.subsample_to_assembly_bam, "rb", check_sq=False)
        if self.contigs == None:
            self.assembly_stats()
        self.num_aligned_reads = {} #number of reads for each contig
        self.avg_min_max_coverages = {} #avg, min and max coverage for a position in contig
        self.coverage_per_base = {}
        self.sum_coverage_per_pos = {}
        with alignment as al:
            for ctg in self.contigs:
                self.num_aligned_reads[ctg] = al.count(contig=ctg)
                covA, covC, covG, covT = al.count_coverage(contig=ctg,quality_threshold=0)
                self.coverage_per_base[ctg] = {}
                self.coverage_per_base[ctg]["A"] = list(covA)
                self.coverage_per_base[ctg]["C"] = list(covC)
                self.coverage_per_base[ctg]["G"] = list(covG)
                self.coverage_per_base[ctg]["T"] = list(covT)
                len_cov = len(covA)
                sum_coverage_per_pos_ctg = []
                sum = 0
                min = covA[0]+covC[0]+covG[0]+covT[0]
                max = min
                for i in range(len_cov-1):
                    c=covA[i]+covC[i]+covG[i]+covT[i]
                    sum_coverage_per_pos_ctg.append(c)
                    sum+=c
                    if c > max:
                        max = c
                    if c < min:
                        min = c
                avg = sum/len_cov
                self.avg_min_max_coverages[ctg]=[avg,min,max]
                self.sum_coverage_per_pos[ctg]=sum_coverage_per_pos_ctg

    def assembly_stats(self):
        while not os.path.exists(self.assembly_fasta): #wait while assembly does not exist
            time.sleep(1)
            print(str(datetime.now())+" :: Stats :: [assembly_stats] waiting for assembly...", flush=True)
        if os.stat(self.assembly_fasta).st_size == 0: #fasta file is empty
            self.num_contigs = 0
            self.contigs = []
            self.contig_lengths = {}
            self.assembly_length = 0
            return
        #fasta file is not empty
        with pysam.FastaFile(self.assembly_fasta) as fa:
            self.num_contigs = fa.nreferences # # contigs
            #self.contig_lengths = list(fa.lengths) # contig lengths
            self.contigs = list(fa.references) # contig names
            self.c_lengths = list(fa.lengths)
            self.contig_lengths = {}
            for i in range(len(self.c_lengths)):
                self.contig_lengths[self.contigs[i]] = self.c_lengths[i]
            self.assembly_length = sum(self.c_lengths)

    def subsample_stats(self):
        self.num_reads = 0
        self.total_reads_len = 0
        with pysam.FastxFile(self.subsample_fastq_copy) as fq:
            for entry in fq:
                self.num_reads+=1
                self.total_reads_len+=len(entry.sequence)

    def run_noref_base(self):
        if self.assembly_length == None:
            self.assembly_stats()
        if self.num_reads == None:
            self.subsample_stats()
        if self.sum_coverage_per_pos == None:
            self.alignment_stats()

    def build_dict(self,log):
        print(str(datetime.now())+" :: Stats :: [build_dict] building dict", flush=True, file=log)
        self.run_noref_base()
        if len(self.dict)>0:
            print(str(datetime.now())+" :: Stats :: [build_dict] skipping, dict. exists", flush=True, file=log)
            return
        self.dict["summary"] = {}
        self.dict["contigs"] = {}
        self.dict["summary"]["time_start"] = str(self.start)
        self.dict["summary"]["time_end"] = str(self.end)
        self.dict["summary"]["num of contigs"] = self.num_contigs
        self.dict["summary"]["num of subsampled reads"] = self.num_reads
        self.dict["summary"]["num unmapped reads"] = self.num_unmapped_reads
        self.dict["summary"]["sum contig length"] = self.assembly_length
        self.dict["summary"]["sum read length"] = self.total_reads_len
        if self.have_ref:
            self.dict["summary"]["N50"] = self.n50(self.assembly_fasta,self.reference)
        if self.assembly_length>0:
            self.dict["summary"]["avg coverage"] = self.total_reads_len/self.assembly_length
        else:
            self.dict["summary"]["avg coverage"] = 0
        for ctg in self.contigs:
            self.dict["contigs"][ctg] = {}
            self.dict["contigs"][ctg]["length"] = self.contig_lengths[ctg]
            self.dict["contigs"][ctg]["num aligned reads"] = self.num_aligned_reads[ctg]
            self.dict["contigs"][ctg]["avg,min,max coverages"] = self.avg_min_max_coverages[ctg]
            #self.dict["contigs"][ctg]["coverage per base"] = self.coverage_per_base[ctg]
            #self.dict["contigs"][ctg]["sum coverage per pos"] = self.sum_coverage_per_pos[ctg]

    def files_processed(self, processed_files):
        if self.num_reads_per_file == None:
            self.num_reads_per_file = Helpers.count_reads_in_fastq(self.in_dir+processed_files[0])
        self.dict["processed_files"] = processed_files
        num_all_files = len(os.listdir(self.in_dir))
        self.dict["catch-up"] = str(len(processed_files)/num_all_files*100)+"%"
        self.dict["all_files"] = num_all_files
        self.dict["summary"]["num processed files"] = len(processed_files)
        self.dict["summary"]["processed reads"]=len(processed_files)*self.num_reads_per_file

    def draw_coverage(self, lock, log):
        try:
            self.run_noref_base()

            fname = self.noref_contig_coverage  
            bar = self.noref_contig_coverage_bar
            #fname_ylim = self.stats_dir+"without_reference/coverage_graph_per_it/coverage_ylim"+str(self.it)+".png"
            def moving_average(arr,n):
                ma = []
                window = 1000
                n = n-1
                for i in range(len(arr)-window-1):
                    ma.append(sum(arr[i:i+window])/window)

                for i in range(n):
                    ma = moving_average(ma,1)    
                return ma
            
            print(str(datetime.now())+" :: Stats :: [draw coverage] waiting for matplotlib lock.", flush=True, file = log)
            if lock.acquire():
                with open(self.lock_log, "a+") as llog:
                    print("thread for it. " +str(self.it)+ " acquired lock in [draw_coverage]", flush=True, file=llog)
            else:
                with open(self.lock_log, "a+") as llog:
                    print("thread for it. " +str(self.it)+ " lock.acquire() returned false in [draw_coverage]", flush=True, file=llog)
            print(str(datetime.now())+" :: Stats :: [draw coverage] acquired matplotlib lock.", flush=True, file = log)

            plt.clf()

            #--------------------------
            #draw per-base coverage
            #plt.figure(figsize=(20,10))
            if len(self.contigs) > 1:
                fig, ax = plt.subplots(nrows = len(self.contigs), figsize=(10, len(self.contigs)*2))
                #print(ax)
                i = 0
                labels = self.contigs
                lengths = self.c_lengths
                labels_lengths = zip(labels, lengths)
                labels_lengths = sorted(labels_lengths, key=lambda x: x[1], reverse = True)
                labels, lengths = zip(*labels_lengths)
                sorted_contigs = labels
                for ctg in sorted_contigs:
                    #col = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
                    col = "blue"
                    #ma = moving_average(self.sum_coverage_per_pos[ctg],3) #too slow
                    ax[i].plot(self.sum_coverage_per_pos[ctg], label=ctg, color=col, linewidth=0.5)
                    ax[i].legend()
                    i+=1
                #if len(self.contig_lengths)<=15: #draw legend only if there are not too many contigs
                #    plt.legend() 
                plt.savefig(fname, bbox_inches='tight',dpi=100)
                #plt.ylim(0, 50) 
                #plt.savefig(fname_ylim, bbox_inches='tight',dpi=100)
            elif len(self.contigs) == 1:
                fig, ax = plt.subplots(nrows = 1, figsize=(10, 2))
                col = "blue"
                ctg = self.contigs[0]
                ma = moving_average(self.sum_coverage_per_pos[ctg],3)
                ax.plot(ma, label=ctg, color=col, linewidth=0.5)
                ax.legend()
                plt.savefig(fname, bbox_inches='tight',dpi=100)
            else: # 0 contigs
                plt.savefig(fname, bbox_inches='tight',dpi=1, figsize=(10,10))
                #plt.savefig(fname_ylim, bbox_inches='tight',dpi=1)
            plt.close()

            #-----------------------------
            #draw coverage bar graph
            plt.clf()

            labels = self.contigs
            lengths = self.c_lengths
            coverages = []
            for i in range(len(labels)):
                cov = sum(self.sum_coverage_per_pos[labels[i]]) #sum of per-base covarage values for contig
                coverages.append(cov/lengths[i])
            fig, ax = plt.subplots(figsize=(10,len(labels)*0.5))
            ax.barh(labels, width = coverages, color="green")
            ax.invert_yaxis()
            
            plt.savefig(bar, bbox_inches='tight',dpi=100)
            plt.close()
            
        except Exception as e: 
            print("Stats :: draw_coverage :: iteration "+str(self.it)+" coverage image draw failed", flush=True, file=self.log_file)
            print(e, flush=True, file=self.log_file)
        finally:
            if lock.locked():
                with open(self.lock_log, "a+") as llog:
                    print("thread for it. " +str(self.it)+ " releases lock in [draw_coverage]", flush=True, file=llog)
                lock.release()


    def draw_contig_lengths(self, lock):
        try:
            self.run_noref_base()
            fname = self.contig_lengths_graph
            labels = self.contigs
            lengths = self.c_lengths
            labels_lengths = zip(labels,lengths)
            labels_lengths = sorted(labels_lengths, key=lambda x: x[1], reverse = True)
            labels, lengths = zip(*labels_lengths)

            lock.acquire()
            with open(self.lock_log, "a+") as llog:
                print("thread for it. " +str(self.it)+ " acquired lock in [draw_contig_lengths]",flush=True, file=llog)

            plt.clf()
            fig, ax = plt.subplots(figsize=(10,len(labels)*0.5))
            ax.barh(labels, width = lengths)
            ax.invert_yaxis()
            
            plt.savefig(fname, bbox_inches='tight',dpi=100)

        except Exception as e: 
            print("Stats :: draw_contig_lengths :: iteration "+str(self.it)+" contig lengths image draw failed", flush=True, file=self.log_file)
            print(e, flush=True, file=self.log_file)
        finally:
            if lock.locked():
                with open(self.lock_log, "a+") as llog:
                    print("thread for it. " +str(self.it)+ " releases lock in [draw_contig_lengths]", flush=True, file=llog)
                lock.release()


#################################################################
    ##### WITH REFERENCE ####
#################################################################

    def minidot_asm_to_ref(self):
        minidot_ps = self.miniasm_dir+"minidot -f3 "+self.assembly_to_ref_paf+" > "+self.assembly_to_ref_minidot_dotplot
        minidot_pdf = "ps2pdf -dEPSCrop "+self.assembly_to_ref_minidot_dotplot+" "+self.minidot_dotplot_pdf
        return os.system(minidot_ps) or os.system(minidot_pdf)

    def create_reads_to_ref_bigwig(self):
        self.bam_to_bigwig(self.bigwig, self.subsample_to_ref_bam, self.genome_sizes)
    
    def bam_to_bigwig(self,out_bigwig, bam,genome_sizes):
        print(str(datetime.now())+" :: Stats :: [bam_to_bigwig] creating bigwig from bam: "+bam, flush=True, file=self.log_file)
        out_bedgraph = "/tmp/dyn_assembly_"+str(self.it)+"_reads_list_to_bigwig.bedgraph"
        bedgraph_sorted = "/tmp/dyn_assembly_"+str(self.it)+"_reads_list_to_bigwig.sorted.bedgraph"

        def bedgraph(bam, out_bedgraph, bedgraph_sorted):
            res = os.system("bedtools genomecov -ibam "+bam+" -bga -split > "+out_bedgraph)
            res += os.system("cat "+out_bedgraph+" | sort -k1,1 -k2,2n > "+bedgraph_sorted)
            return res
        
        def bedgraph_to_bigwig(bedgraph_sorted,genome_sizes,out_bigwig):
            return os.system("bedGraphToBigWig "+bedgraph_sorted+" "+genome_sizes+" "+out_bigwig)

        def clean(out_bedgraph,bedgraph_sorted):
            res = os.system("rm "+out_bedgraph)
            res += os.system("rm "+bedgraph_sorted)
            if res > 0:
                raise Exception("clean after reads list to bigwig failed")
        
        print(str(datetime.now())+" :: Stats :: [bam_to_bigwig] running bedgraph",flush=True, file=self.log_file)
        if bedgraph(bam, out_bedgraph, bedgraph_sorted) > 0:
            print("!!! bedgraph in bam_to_bigwig retruned non-zero code", file=self.log_file)
            #raise Exception("bedgraph in reads_list_to_bigwig retruned non-zero code")
        print(str(datetime.now())+" :: Stats :: [bam_to_bigwig] running bedgraph to bigwig",flush=True, file=self.log_file)
        if bedgraph_to_bigwig(bedgraph_sorted,genome_sizes,out_bigwig) > 0:
            print("!!! bedgraph_to_bigwig in fastq_to_bigwig retruned non-zero code", flush=True, file=self.log_file)
            #raise Exception("bed graph to bigwig in reads_list_to_bigwig retruned non-zero code")
        print(str(datetime.now())+" :: Stats :: [bam_to_bigwig] cleaning bedgraph files",flush=True, file=self.log_file)
        clean(out_bedgraph,bedgraph_sorted)    
    
    def ref_subsample_coverage(self, lock,log):
        self.reference_coverage(self.assembly_to_ref_bed, self.subsample_to_ref_coverage_png, self.subsample_to_ref_bam, lock,log)

    def assembly_to_ref_iteration_dotplot(self, lock):
        self.dotplot_with_coverage(self.assembly_to_ref_paf, self.coverage_threshold, self.assembly_to_ref_dotplot, self.lock_log, lock, self.it, bam_reads_to_asm=self.subsample_to_ref_bam)

    def reference_coverage(self, bed, out_png, bamfile, lock, log):
        print(str(datetime.now())+" :: Stats :: [reference coverage] generating plotting data.", flush=True, file = log)
        
        bam = pysam.AlignmentFile(bamfile, 'rb')
        contigs = bam.references
        contig_lengths = bam.lengths
        sum_coverage_per_pos = {}
        for contig_id in contigs:
            A,C,G,T = bam.count_coverage(contig_id)
            coverage = [a+c+g+t for a,c,g,t in zip(A,C,G,T)]
            sum_coverage_per_pos[contig_id] = coverage


        assemblies = {}
        for contig_id in contigs:
            assemblies[contig_id] = []
        with open(bed, "r") as bedfile:
            for line in bedfile:
                l = line.split()
                assemblies[l[0]].append([int(l[1]),int(l[2])])
       
        print(str(datetime.now())+" :: Stats :: [reference coverage] waiting for matplotlib lock.", flush=True, file = log)
        if lock.acquire():
            with open(self.lock_log, "a+") as llog:
                print("thread for it. " +str(self.it)+ " acquired lock in [reference_coverage]", flush=True, file=llog)
        else:
            with open(self.lock_log, "a+") as llog:
                print("thread for it. " +str(self.it)+ " lock.acquire() returned false in [reference_coverage]", flush=True, file=llog)
        print(str(datetime.now())+" :: Stats :: [reference coverage] matplotlib lock acquired, drawing reference coverage plot.", flush=True, file = log)

        try:
            plt.clf()
            if len(contigs) > 1:
                fig, ax = plt.subplots(nrows = len(contigs), figsize=(10, len(contigs)*2))
                #print(ax)
                i = 0
                labels = contigs
                lengths = contig_lengths
                labels_lengths = zip(labels, lengths)
                labels_lengths = sorted(labels_lengths, key=lambda x: x[1], reverse = True)
                labels, lengths = zip(*labels_lengths)
                sorted_contigs = labels
                for ctg in sorted_contigs:
                    col = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]) #random for each contig in ref.
                    #ma = moving_average(self.sum_coverage_per_pos[ctg],3) #too slow
                    ax[i].plot(sum_coverage_per_pos[ctg], label=ctg, color="blue", linewidth=0.5)
                    for asm in assemblies[ctg]:
                        col = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]) #random for each contig in assembly
                        ax[i].plot([asm[0],asm[1]],[-1,-1], linewidth = 5, color=col)
                    ax[i].legend()
                    i+=1
                #if len(self.contig_lengths)<=15: #draw legend only if there are not too many contigs
                #    plt.legend() 
                plt.savefig(out_png, bbox_inches='tight',dpi=100)
                #plt.ylim(0, 50) 
                #plt.savefig(fname_ylim, bbox_inches='tight',dpi=100)
            elif len(contigs) == 1:
                fig, ax = plt.subplots(nrows = 1, figsize=(10, 2))
                col = "blue"
                ctg = contigs[0]
                ax.plot(sum_coverage_per_pos[ctg], label=ctg, color=col, linewidth=0.5)
                ax.legend()
                plt.savefig(out_png, bbox_inches='tight',dpi=100)
            else: # 0 contigs
                plt.savefig(out_png, bbox_inches='tight',dpi=1, figsize=(10,10))
                #plt.savefig(fname_ylim, bbox_inches='tight',dpi=1)
            plt.close()

            print(str(datetime.now())+" :: Stats :: [reference coverage] finished drawing reference coverage plot.", flush=True, file = log)

        except Exception as e: 
            print("Stats :: reference coverage :: iteration "+str(self.it)+" reference coverage image draw failed", flush=True, file=self.log_file)
            print(e, flush=True, file=self.log_file)
        finally:
            if lock.locked():
                with open(self.lock_log, "a+") as llog:
                    print("thread for it. " +str(self.it)+ " releases lock in [reference_coverage]", flush=True, file=llog)
                lock.release()

    def n50(self, assembly, reference):
        """
        pocita N50 score - min. dlzka contigu, pre ktory vsetky dlhsie/rovnako dlhe dokopy pokryvaju 50% assembly 
        filtered assembly
        reference length (sum)
        NEKONTROLUJE ci su zrovnane k referencii
        """
        if os.path.getsize(assembly) == 0:
            #print("assembly size is 0")
            return -1 #no assembly

        #print("assembly size is non-zero")
        with pysam.FastaFile(reference) as ref:
            ref_lengths= list(ref.lengths)
        reference_length = sum(ref_lengths)

        os.system("samtools faidx "+assembly)
        asm_lengths = []
        with pysam.FastaFile(assembly) as asm:
            asm_lengths= list(asm.lengths)
        sorted_asm_lengths = sorted(asm_lengths, reverse=True)
        n_num = 50
        fraction = n_num/100
        s = 0
        for l in sorted_asm_lengths:
            if s+l > reference_length*fraction:
                return l
            else: # s+l <= reference_length*fraction
                s += l
        return l

#################################################################
    ##### WITH or without ref ####
#################################################################
    
    def start_end_times(self, start,end):
        self.start = start
        self.end = end
        self.duration = end-start

    def dotplot_with_coverage(self,paf_asm_to_asm, coverage_threshold, save_path,lock_log, lock, it, bam_reads_to_asm=None):
        print(str(datetime.now())+" :: Stats :: [dotplot with coverage] drawing dotplot for "+paf_asm_to_asm+", save to: "+save_path, flush=True, file=self.log_file)
        
        
        if bam_reads_to_asm!=None:
            #index bam
            os.system("samtools index "+bam_reads_to_asm)

        #read paf
        contigs = {}  #dictionary to store data for each contig & query
        ref_contigs = []
        assembly_contigs = []
        ref_contig_lengths = []
        assembly_contig_lengths = []
        target_starts = [] #for sorting
        with open(paf_asm_to_asm, 'r') as f:
            for line in f:
                parts = line.split('\t')
                contig_id = parts[5]  # target id (ref. assembly contig)
                contig_length = int(parts[6])
                query_id = parts[0]
                query_length = int(parts[1])
                query_start = int(parts[2])  
                target_start = int(parts[7]) 
                query_end = int(parts[3]) 
                target_end = int(parts[8])
                strand = parts[4]
                
                if contig_id not in ref_contigs:
                    contigs[contig_id] = {}
                    ref_contigs.append(contig_id)
                    ref_contig_lengths.append(contig_length)
                if not query_id in assembly_contigs:
                    assembly_contigs.append(query_id)
                    assembly_contig_lengths.append(query_length)
                    target_starts.append(target_start) #first one for assembled contig
                if not query_id in contigs[contig_id]:
                    contigs[contig_id][query_id]={'x_starts':[], 'y_starts':[], 'x_ends':[], 'y_ends':[], 'strands':[]}
                contigs[contig_id][query_id]['x_starts'].append(target_start)
                contigs[contig_id][query_id]['y_starts'].append(query_start)
                contigs[contig_id][query_id]['x_ends'].append(target_end)
                contigs[contig_id][query_id]['y_ends'].append(query_end)
                contigs[contig_id][query_id]['strands'].append(strand)
        
        num_ref_contigs = len(ref_contigs)
        num_assembly_contigs = len(assembly_contigs)

        if num_assembly_contigs > 0 and num_ref_contigs > 0:
            #sort row contigs by start position in ref
            tuples_to_sort = zip(assembly_contigs, assembly_contig_lengths, target_starts)
            sorted_tuples = sorted(tuples_to_sort, key=lambda x: x[2], reverse=True)
            assembly_contigs, assembly_contig_lengths, target_starts = zip(*sorted_tuples) 
            #sort column contigs by lengths
            tuples_to_sort = zip(ref_contigs, ref_contig_lengths)
            sorted_tuples = sorted(tuples_to_sort, key=lambda x: x[1], reverse=True)
            ref_contigs, ref_contig_lengths = zip(*sorted_tuples)

        try:
            lock.acquire()
            with open(lock_log, "a+") as llog:
                print("thread for it. " +str(it)+ " acquired lock in [refstats dotplot with coverage]", flush=True, file=llog)

            plt.clf()

            if num_assembly_contigs > 0 and num_ref_contigs > 0:
                if bam_reads_to_asm!=None:
                    num_plot_rows = num_assembly_contigs + 2 #add 2 subplots for each contig for coverage
                else:
                    num_plot_rows = num_assembly_contigs
                fig, axes = plt.subplots(num_plot_rows, num_ref_contigs, figsize=(15, 1.5*num_plot_rows), squeeze=False, sharex='col', sharey='row', gridspec_kw={"width_ratios":ref_contig_lengths})
                for row in range(num_assembly_contigs):
                    for col in range(num_ref_contigs):
                        ref_contig = ref_contigs[col]
                        ref_contig_length = ref_contig_lengths[col]
                        asm_contig = assembly_contigs[row]
                        asm_contig_length = assembly_contig_lengths[row]
                        if asm_contig in contigs[ref_contig]:
                            num_lines = len(contigs[ref_contig][asm_contig]['x_starts'])
                            for i in range(num_lines):
                                x1 = contigs[ref_contig][asm_contig]['x_starts'][i]
                                x2 = contigs[ref_contig][asm_contig]['x_ends'][i]
                                y1 = contigs[ref_contig][asm_contig]['y_starts'][i]
                                y2 = contigs[ref_contig][asm_contig]['y_ends'][i]
                                strand = contigs[ref_contig][asm_contig]['strands'][i]
                                axes[row][col].set_xlim([0, ref_contig_length])
                                axes[row][col].set_ylim([0, asm_contig_length])
                                if strand == '+':
                                    axes[row][col].plot([x1,x2],[y1,y2], color='blue', linewidth=2)
                                else: #strand == '-'
                                    axes[row][col].plot([x1,x2],[y2,y1], color='red', linewidth=2)
                        #axes[row][col].set_xlabel(ref_contig)
                        if col == 0:
                            axes[row][col].set_ylabel(asm_contig, rotation=0, labelpad=-100)

                if bam_reads_to_asm!=None:
                    #read bam, caluclate coverage
                    bam = pysam.AlignmentFile(bam_reads_to_asm, 'rb')
                    for col in range(num_ref_contigs):
                        #calculate coverage for the current contig
                        contig_id = ref_contigs[col]
                        coverage = [0] * ref_contig_lengths[col]
                        ax = axes[num_assembly_contigs][col]
                        ax_ylim = axes[num_assembly_contigs+1][col]

                        """
                        for read in bam.fetch(contig_id):
                            if not read.is_unmapped:
                                start_x = read.reference_start
                                end_x = read.reference_end
                                for i in range(start_x,end_x):
                                    coverage[i]+=1
                        """
                        A,C,G,T = bam.count_coverage(contig_id)
                        coverage = [a+c+g+t for a,c,g,t in zip(A,C,G,T)]
                        
                        #fill
                        ax.fill_between([i for i in range(ref_contig_lengths[col])], coverage, [0]*ref_contig_lengths[col], color='green', alpha=0.5)
                        ax_ylim.fill_between([i for i in range(ref_contig_lengths[col])], coverage, [0]*ref_contig_lengths[col], color='green', alpha=0.5)

                        #plot coverage in graph at the bottom
                        ax.plot(coverage, color="green")
                        ax_ylim.plot(coverage, color="green")
                        

                        ax_ylim.set_ylim(0,coverage_threshold)
                        
                        #ax.set_xlabel(contig_id)
                        if col == 0:
                            ax.set_ylabel('coverage', rotation=0, labelpad=-100)
                            ax_ylim.set_ylabel("0-"+str(coverage_threshold)+'\ncoverage', rotation=0, labelpad=-100)

                        ax_ylim.set_xlabel(contig_id)

                plt.subplots_adjust(wspace=0, hspace=0)
                #plt.tight_layout()
                plt.savefig(save_path, bbox_inches='tight',dpi=100)

                plt.close()
        except Exception as e:
            print("error drawing dotplot")
        finally:
            if lock.locked():
                with open(lock_log, "a+") as llog:
                    print("thread for it. " +str(it)+ " releases lock in [refstats dotplot with coverage]", flush=True, file=llog)
                lock.release()

    def print_stats(self):
        json_string = json.dumps(self.dict)
        with open(self.stats_file, "a+") as outfile: 
        #    json.dump(self.dict, outfile)
            print(json_string, file = outfile)
        print(str(datetime.now())+" :: Stats :: [print_stats] wrote stats to ", self.stats_file, flush=True, file=self.log_file)
        print(str(datetime.now())+" :: Stats :: [print_stats] summary: "+str(self.dict["summary"]), flush=True, file=self.log_file)


    def print_html(self):
        if len(self.dict) == 0:
            self.build_dict(self.log_file)
            self.files_processed([])
        html_header = "<!DOCTYPE html>\n<html>\n<head>\n<title>Iteration "+str(self.it)+"</title>\n</head>\n<body>\n"
        processed_files = len(self.dict["processed_files"])
        html = "<H1>Iteration: "+str(self.it)+"</H1>"
        html += '<a href="last_it.html">load latest</a></br>\n'
        html += '<a href="it_'+str(max(self.it-1,0))+'.html">previous</a>\n'
        html += '<a href="it_'+str(self.it+1)+'.html">next</a> </br>\n'
        html += "processed files: "+str(processed_files)+"</br>\n"
        html += "processed reads:"+str(processed_files*self.num_reads_per_file,)+"</br>\n"
        html += "num of files in input directory:"+str(self.dict["all_files"])+"</br>\n"
        html += "catch-up:"+str(self.dict["catch-up"])+"</br>\n"
        html += "start time: "+str(self.start)+"</br>\n"
        html += "end time: "+str(self.end)+"</br>\n"
        html += "duration: "+str(self.end-self.start)+"</br>\n"
        stats_end = datetime.now()
        html += "stats end: "+str(stats_end)+"</br>\n"
        html += "stats duration: "+str(stats_end-self.end)+"</br>\n"
        html += "</br>subsample/assembly stats: </br>\n"
        html += "# contigs: "+str(self.num_contigs)+"</br>\n"
        html += "# subsampled reads: "+str(self.num_reads)+"</br>\n"
        html += "# unmapped reads: "+str(self.num_unmapped_reads)+"</br>\n"
        html += "sum length of contigs: "+str(self.assembly_length)+"</br>\n"
        html += "sum length of reads: "+str(self.total_reads_len)+"</br>\n"
        if self.assembly_length>0:
            html += "avg coverage (crude) [sum of read lengths/assembly length]: "+str(round(self.total_reads_len/self.assembly_length,2))+"</br>\n"
        if self.have_ref:
            html += "N50 (raw assembly): "+str(self.dict["summary"]["N50"])+"</br>"

        html += "</br></br> \n"

        #html += '<a href="../quast/quast_'+str(self.it)+'">Quast</a> \n </br> \n'

        html += "<details> <summary>Contig lengths</summary>"
        if not self.num_contigs == 0:
            html += '<img style="max-width:100%;" src="../without_reference/contig_lengths_graph_per_it/contig_lengths_'+str(self.it)+'.png"> \n </br> \n'
        else:
            html += '<p>contig lengths graph for this iteration is irrelevant, since we have 0 contigs</p> \n'
        html += "</details>"

        """
        html += "<details> <summary>Assemblies alignment (prev. it. to this)</summary>"
        if not self.it == 0:
            image_path = "../without_reference/prev_asm_to_this/prev_asm_to_this_dotplot_"+str(self.it)+".png"
            if os.path.exists(image_path):
                html += '<img style="max-width:100%;" src="'+image_path+'" > \n </br> \n'
            else:
                html += '<p>file not found</p>\n'
        #    for i in range(self.paf_figures+1):
        #        html += '<img style="" src="../without_reference/prev_asm_to_this/map_pafviz_'+str(self.it)+'_'+str(i)+'.png"> \n </br> \n'
        else:
            html += '<p>paf visualization is irrelevant for 0th iteration, since there is no previous assembly</p> \n'
        html += "</details>"
        """

        html += "<details> <summary>avg coverage per-contig (bar graph)</summary>"
        if not self.num_contigs == 0:
            image_path = "../without_reference/coverage_graph_per_it/coverage_bar_"+str(self.it)+".png"
            if os.path.exists(self.stats_dir+"without_reference/coverage_graph_per_it/coverage_bar_"+str(self.it)+".png"):
                html += '<img style="max-width:100%;" src="'+image_path+'"> \n </br> \n'
            else:
                html += '<p>'+image_path+' : file not found</p>\n'
        else:
            html += '<p>avg coverage graph for this iteration is irrelevant, since we have 0 contigs</p> \n'
        html += "</details>"

        html += "<details> <summary>Per-base coverage graph for each contig</summary>"
        if not self.num_contigs == 0:
            image_path = "../without_reference/coverage_graph_per_it/coverage_"+str(self.it)+".png"
            if os.path.exists(self.stats_dir+"without_reference/coverage_graph_per_it/coverage_"+str(self.it)+".png"):
                html += '<img style="max-width:100%;" src="'+image_path+'" > \n </br> \n'
            else:
                html += '<p>'+image_path+' : file not found</p>\n'
        else:
            html += '<p>per-base coverage graph for this iteration is irrelevant, since we have 0 contigs</p> \n'
        html += "</details>"

        if self.have_ref:
            html+="<details> <summary>reference coverage with contigs</summary>" 
            image_path = "../with_reference/reads_to_ref_coverage/ref_cov_"+str(self.it)+".png"
            if os.path.exists(self.stats_dir+"with_reference/reads_to_ref_coverage/ref_cov_"+str(self.it)+".png"):
                html += '<img style="max-width:100%;" src="'+image_path+'" > \n </br> \n'
                html+="</br><p>assembled contigs are shown as colored lines at the bottom of the coverage graphs</p>"
            else:
                html+="<p>"+image_path+" : file not found </p>"
            html+="</details>"
            html+="<details> <summary>assembly to reference dotplot (from minidot)</summary>"
            pdf_path = "../with_reference/assembly_to_ref_dotplots/dotplot_"+str(self.it)+".pdf"
            if os.path.exists(self.stats_dir+"with_reference/assembly_to_ref_dotplots/dotplot_"+str(self.it)+".pdf"):
                html+="<a href='"+pdf_path+"'>pdf dotplot from minidot</a>"
            html+="</details>"

        html_end = "</body>\n</html>"

        html_file = self.stats_dir+"html/it_"+str(self.it)+".html"
        all_it_html_body = self.stats_dir+"html/all_body.html"
        all_it_html = self.stats_dir+"html/all.html"

        with open(html_file, "a+") as out:
            print(html_header+html+html_end, file=out) 

        
        with open(all_it_html_body, "a+") as out:
            print(html, file=out)

        with open(all_it_html, "w+") as out:
            inp = open(all_it_html_body, "r")
            inner_html = inp.read()
            inp.close()

            all_html_start = "<!DOCTYPE html>\n<html>\n<head>\n<title>Iterations 0-"+str(self.it)+"(all)</title>\n</head>\n<body>\n"
            print(all_html_start+inner_html+html_end, file=out)
        

        #latest_link_path = self.out_dir+"last_it.html"
        #latest_link_val = self.stats_dir+"html/it_"+str(self.it)+".html"
        latest_link_path_stats = self.stats_dir+"html/last_it.html"
        latest_link_val_stats = "it_"+str(self.it)+".html"


        #if os.path.exists(latest_link_path):
        #    os.unlink(latest_link_path)
        if os.path.exists(latest_link_path_stats):
            os.unlink(latest_link_path_stats)
        #os.symlink(latest_link_val, latest_link_path)
        os.symlink(latest_link_val_stats,latest_link_path_stats)






    ### cleanup ###

    def clean_files(self):
        os.unlink(self.subsample_fastq_copy)
        os.unlink(self.subsample_to_assembly_bam)
        os.unlink(self.subsample_to_assembly_bam+".bai")
        os.unlink(self.subsample_to_ref_bam)
        os.unlink(self.subsample_to_ref_bam+".bai")
        self.log_file.close()
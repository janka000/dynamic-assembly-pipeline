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
from Stats import Stats
import traceback
import random
from threading import Thread

class FastStats(Stats):
    """
    computing stats for completed assembly, uses paf alignments only, as creating bam per-base alignments is too slow
    """
    def __init__(self, in_dir, out_dir, stats_dir, subsample_fastq, coverage_threshold, minimap, miniasm_dir, reference = None):
        super().__init__(in_dir, out_dir, stats_dir, subsample_fastq, coverage_threshold, minimap, miniasm_dir, reference = reference)
        self.debug_or_have_time = False

    def create_no_ref_necessary_files(self,log):
        print("\n"+str(datetime.now())+" :: Stats :: [create_necessary_files] creating files for stats without reference ", flush=True, file=log)
        def create_subsample_to_asm_paf(assembly,subsample):
            self.minimap.map_reads_to_assembly_paf(subsample,assembly,self.subsample_to_assembly_paf)
        def create_paf_asm_to_asm(prev_asm,this_asm):
            if self.it > 0: #previous does not exist otherwise
                self.minimap.asm_to_asm_paf(this_asm, prev_asm, self.paf_asm_to_asm)
        print(str(datetime.now())+" :: Stats :: [create_no_ref_necessary_files] creating subsample to assembly paf ", flush=True, file=log)
        create_subsample_to_asm_paf(self.assembly_fasta,self.subsample_fastq_copy)
        print(str(datetime.now())+" :: Stats :: [create_no_ref_necessary_files] creating assembly to assembly paf ", flush=True, file=log)
        create_paf_asm_to_asm(self.prev_it_assembly, self.assembly_fasta)

    def create_ref_necessary_files(self,log):
        print("\n"+str(datetime.now())+" :: Stats :: [create_necessary_files] creating files for stats with reference ", flush=True, file=log)
        def create_subsample_to_ref_paf(ref,subsample):
            self.minimap.map_reads_to_assembly_paf(subsample,ref,self.subsample_to_ref_paf)
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
        print(str(datetime.now())+" :: Stats :: [create_ref_necessary_files] creating subsample to ref paf ", flush=True, file=log)
        create_subsample_to_ref_paf(self.reference, self.subsample_fastq_copy)
        print(str(datetime.now())+" :: Stats :: [create_ref_necessary_files] creating assembly to ref paf ", flush=True, file=log)
        create_assembly_to_ref_paf(self.assembly_fasta, self.reference)
        print(str(datetime.now())+" :: Stats :: [create_ref_necessary_files] creating assembly to ref bed (from paf) ", flush=True, file=log)
        create_assembly_to_ref_bed(self.assembly_to_ref_paf, self.assembly_to_ref_bed)

    def set_iteration(self, it):
        super().set_iteration(it)
        self.subsample_to_ref_paf = self.stats_dir+"with_reference/bam_subsample_reads_to_ref/subsample_to_ref_"+str(it)+".paf"
        self.subsample_to_assembly_paf = self.stats_dir+"stats_subsample_to_assembly_"+str(it)+".paf"
        self.subsample_to_ref_bam = None
        self.subsample_to_assembly_bam = None

    def alignment_stats(self):
        self.alignment_stats_from_paf(self.subsample_to_assembly_paf, self.assembly_fasta)
    
    def count_aligned_lengths_per_contig(self,paf):
            aligned_lengths = {}
            alignment_starts = {}
            with open(paf, "r") as pf:
                for line in pf:
                    #print("paf line: ",line)
                    split_line = line.split()
                    #6	string	Target sequence name
                    contig = split_line[5]

                    #8	int	Target start on original strand (0-based)  
                    target_start = int(split_line[7])   
                    #9	int	Target end on original strand (0-based)
                    target_end = int(split_line[8])   

                    alignment_len = abs(target_end - target_start)
                    if not contig in aligned_lengths:
                        aligned_lengths[contig] = []
                        alignment_starts[contig] = []
                    aligned_lengths[contig].append(alignment_len)
                    alignment_starts[contig].append(target_start)
            return aligned_lengths,alignment_starts

    def alignment_stats_from_paf(self,paf, assembly_fa):
        self.num_aligned_reads = {} #number of reads for each contig
        self.avg_min_max_coverages = {} #avg, min and max coverage for a position in contig
        self.coverage_per_base = {}
        self.sum_coverage_per_pos = {}
        aligned_lengths,alignment_starts = self.count_aligned_lengths_per_contig(paf)
        if os.path.getsize(assembly_fa) == 0:
            return
        with pysam.FastaFile(assembly_fa) as ref:
            ref_lengths= list(ref.lengths)
            ref_names = list(ref.references)
            nl = sorted(zip(ref_names,ref_lengths), reverse=True)
        for contig_tup in nl:
            ctg = contig_tup[0]
            clen = contig_tup[1]
            if ctg in aligned_lengths:
                sum_lengths = sum(aligned_lengths[ctg])
                self.num_aligned_reads[ctg]=len(aligned_lengths[ctg])
            else:
                sum_lengths = 0
                self.num_aligned_reads[ctg]=0
            avg = sum_lengths / clen
            #self.sum_coverage_per_pos[ctg]=[avg for i in range(contig_tup[1])] #dummy "coverage per base"
            self.sum_coverage_per_pos[ctg],minc,maxc = self.count_paf_coverage(aligned_lengths[ctg],alignment_starts[ctg],clen)
            self.avg_min_max_coverages[ctg]=[avg,minc,maxc]

    def count_paf_coverage(self, al_lengths, al_starts, clen):
        cov_arr = [0] * clen
        minc = float('inf')
        maxc = 0

        for al_len, al_st in zip(al_lengths, al_starts):
            al_end = al_st + al_len
            for i in range(al_st, al_end):
                cov_arr[i] += 1
                if cov_arr[i] < minc:
                    minc = cov_arr[i]
                if cov_arr[i] > maxc:
                    maxc = cov_arr[i]

        minc = minc if minc != float('inf') else 0

        return cov_arr, minc, maxc

    def ref_subsample_coverage(self, lock,log):
        self.reference_coverage(self.assembly_to_ref_bed, self.subsample_to_ref_coverage_png, self.subsample_to_ref_paf, lock,log)

    def reference_coverage(self, bed, out_png, paffile, lock, log):
        print(str(datetime.now())+" :: Stats :: [reference coverage] generating plotting data.", flush=True, file = log)
        
        with pysam.FastaFile(self.reference) as ref:
            contig_lengths= list(ref.lengths)
            contigs = list(ref.references)

        sum_coverage_per_pos = {}

        aligned_lengths,alignment_starts = self.count_aligned_lengths_per_contig(paffile)

        for contig_id,clen in zip(contigs,contig_lengths):
            if contig_id in aligned_lengths:
                #avg = sum(aligned_lengths[contig_id])/clen
                #sum_coverage_per_pos[contig_id] = [avg for i in range(clen)]
                sum_coverage_per_pos[contig_id],minc,maxc = self.count_paf_coverage(aligned_lengths[contig_id],alignment_starts[contig_id],clen)
            else:
                sum_coverage_per_pos[contig_id] = [0 for i in range(clen)]

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
            traceback.print_exc() 
        finally:
            if lock.locked():
                with open(self.lock_log, "a+") as llog:
                    print("thread for it. " +str(self.it)+ " releases lock in [reference_coverage]", flush=True, file=llog)
                lock.release()

    def clean_files(self):
        os.unlink(self.subsample_fastq_copy)
        os.unlink(self.subsample_to_assembly_paf)
        os.unlink(self.subsample_to_ref_paf)
        self.log_file.close()


"""
        Col	Type	Description
        1	string	Query sequence name                      ...(reads)
        2	int	Query sequence length
        3	int	Query start (0-based; BED-like; closed)
        4	int	Query end (0-based; BED-like; open)
        5	char	Relative strand: "+" or "-"
        6	string	Target sequence name                     ...(assembly)
        7	int	Target sequence length
        8	int	Target start on original strand (0-based)     ..begigning of the alignemnet of the QUERY IN CONTIG (i.e. start of the aligned part of the CONTIG - the part that is covered by the read..)
        9	int	Target end on original strand (0-based)       ..and the end ..open interval?
        10	int	Number of residue matches
        11	int	Alignment block length
        12	int	Mapping quality (0-255; 255 for missing)
"""
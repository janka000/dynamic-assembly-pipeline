import os, sys
import gzip
import pysam
import random

import time
from datetime import datetime
from Helpers import Helpers
from ReadsSubsetBase import ReadsSubsetBase

class ReadsSubset(ReadsSubsetBase):
    def __init__(self, threshold, strategy, assembly_fasta, it):
        super().__init__(threshold,strategy,assembly_fasta)
        
        self.it = it
        self.mapped_to_contig = [[] for c in range(self.num_contigs)] #[read id, start, end] (for overhangs)
        self.mapped_from_prev_subsample_for_contig = [0 for c in range(self.num_contigs)] #total length of mappings for contig from prev. subsample

        
    def parse_prev_subsample_paf(self, paf): #reads from previous subsample aligned to new assembly - so that we know how many new reads to add
        with gzip.open(paf,'rt') as pf:      
            for line in pf:        
                paf_row = line.split()
                contig = self.ctg_name_to_num[paf_row[5]] #number
                read_start = int(paf_row[7])
                read_end = int(paf_row[8])-1
                read_length = read_end-read_start #mapped read part
                self.mapped_from_prev_subsample_for_contig[contig] += read_length

            for contig in range(self.num_contigs):
                print(str(datetime.now())+" :: ReadsSubset :: [parse_prev_subsample_paf] :: "+str(self.mapped_from_prev_subsample_for_contig[contig])+" sum of mappings from previous subsample for contig "+self.contig_names[contig])

    def add_mapped_to_contig(self,mapping,read_id,contig):
        self.mapped_to_contig[contig].append(mapping) #format [read_name, start, end]
        if not read_id in self.mapped_to_contig_read_ids[contig]: #without read duplicates!
            self.mapped_to_contig_read_ids[contig].append(read_id)  

    def all_reads_for_contig(self, contig):
        self.mapped_to_contig_subsampled[contig] = self.mapped_to_contig_read_ids[contig]
        
    def subsamples_for_contigs(self): #assuming that each contig has at least threshold-coverage reads mapped (otherwise should not be created by mianiasm)
        for contig_t in self.contigs_from_longest_to_shortest: # od najdlhsieho !!
            contig = contig_t[0]
            self.choose_reads_for_contig(contig)
        
    def choose_reads_for_contig(self,contig):   
        if self.subsample_strategy == "probabilistic":
            print(str(datetime.now())+" :: ReadSubset :: [choose_reads_for_contig] subsampling reads for contig "+self.contig_names[contig]+" using probabilistic stategy")
            self.subsample_for_contig_probabilistic(contig)
            print(str(datetime.now())+" :: ReadSubset :: [choose_reads_for_contig] num of subsampled reads for contig "+self.contig_names[contig]+" of length "+str(self.contig_lengths[contig])+":"+str(len(self.mapped_to_contig_subsampled[contig])))
        elif self.subsample_strategy == "take_all":
            self.all_reads_for_contig(contig)
        else:
            print(str(datetime.now())+" :: ReadSubset :: [choose_reads_for_contig] ERROR: invalid subsample strategy (please choose one of probabilistic/deterministic/take_all)", file=sys.stderr)
            raise Exception("invalid subsampling strategy")
    
    def subsample_for_contig_probabilistic(self, contig): 
        num_of_mapped_reads = len(self.mapped_to_contig_read_ids[contig])
        print(str(datetime.now())+" :: ReadSubset :: [subsample_for_contig_probabilistic] num of mapped reads for contig "+self.contig_names[contig]+": "+str(num_of_mapped_reads))
        contig_length = self.contig_lengths[contig]
        prob_more_that_one = 0
        
        #consider reads form longer contigs - already selected for longer, and mapped to this
        map_sum_from_longer_contigs = self.cov_from_longer_contigs(contig) #calculate coverage generated from reads subsampled for lerger contigs
        
        #consider reads subsampled in previous iterations
        if self.it > 0:
            map_sum_from_prev_it = self.mapped_from_prev_subsample_for_contig[contig] #coverage generated from mappings for reads that already were in subsample from previous iteration
        else:
            map_sum_from_prev_it = 0

        cov_from_longer_contigs = map_sum_from_longer_contigs/self.contig_lengths[contig]
        cov_from_prev_it = map_sum_from_prev_it/self.contig_lengths[contig]
        threshold = self.coverage_threshold-cov_from_longer_contigs-cov_from_prev_it
        seen_mapped_reads = len(self.seen_reads_for_contig[contig])
        num_unseen_mapped_reads = num_of_mapped_reads - seen_mapped_reads
        unseen_mapped_reads = list(set(self.mapped_to_contig_read_ids[contig]) - set(self.seen_reads_for_contig[contig]))

        if threshold <= 0: #uz mame dostatocne pokrytie contigu, netreba prenho dalsie ready
            return

        #potrebujeme dalsie ready:
        for read in unseen_mapped_reads:
            read_name = read
            mapped_length = self.read_length_dict_for_contig[contig][read_name] # sucet vsetkych zarovnanych usekov pre dany read
            if not read_name in self.seen_reads_for_contig[contig]:
                prob = (contig_length*threshold*mapped_length)/(self.sum_read_lengths_for_contigs_sq[contig]-map_sum_from_longer_contigs) #if there is a contig, num of reads and read length should aways be non-zero.. 
                if random.random() < prob: #choose read with given probability
                    self.mapped_to_contig_subsampled[contig].append(read_name)
                    if prob > 1:
                        #print(str(datetime.now())+" :: ReadSubset :: [subsample_for_contig_probabilistic] !!! warning: choice probability of a read was > 1 "+str(prob)+" ..coverage threshold: "+str(self.coverage_threshold)+", num. of mapped reads: ", str(num_of_mapped_reads))
                        prob_more_that_one+=1
        print(str(datetime.now())+" :: ReadSubset :: [subsample_for_contig_probabilistic] "+str(prob_more_that_one)+" reads chosen with 'prob.' > 1")
        
    def cov_from_longer_contigs(self,contig):
        cov = 0 #pokrytie sposobene readmi zarovnanymi k dalsim contigom
        for ctg_t in self.contigs_from_longest_to_shortest:
            ctg = ctg_t[0]
            if ctg == contig: #od zaciatku pola contigov usporiadanych od najdlhsieho, chceme skoncit na tom, ktory nas zaujima
                return cov
            #(else) ctg je dlhsi ako contig 
            for read_id in self.mapped_to_contig_subsampled[ctg]: #prejdi subsmapolvane ready dlhsieho contigu
                if read_id in self.mapped_to_contig_read_ids[contig]: #read je zarovnany aj k contigu co nas zaujima
                    cov += self.read_length_dict_for_contig[contig][read_id] #pripocitaj sucet dlzok zarovnanych usekov toho readu k cov
                    self.seen_reads_for_contig[contig].add(read_id)
        
    def write_fastq(self, reads_dir, subsample_fastq, unmapped_fastq, log_read_names, prev_log_read_names): #write the selected reads to a new fastq for miniasm
        subsampled_reads = self.get_all_subsampled_reads()
        print(str(datetime.now())+" :: ReadSubset :: [write_fastq] NUM OF SUBSAMPLED READS:"+str(len(subsampled_reads)))
        #print(str(subsampled_reads))
        all_mapped_reads = self.get_all_mapped_reads()
        mapped_count = len(all_mapped_reads)
        subsampled_count = 0
        unmapped_count = 0
        new_unmapped_fastq = unmapped_fastq+".new"
        if prev_log_read_names != log_read_names: #in first it 0 = 0
            os.system("cp "+prev_log_read_names+" "+log_read_names)
        with open(subsample_fastq, "a+") as subsample_fq, open(new_unmapped_fastq,"w+") as unmapped_fq:
            with open(log_read_names, "a+") as read_log:
                for f in os.listdir(reads_dir):
                    fastq = reads_dir+f
                    print(str(datetime.now())+" :: ReadSubset :: [write_fastq] adding reads from "+fastq)
                    with pysam.FastxFile(fastq) as in_fq:
                        for entry in in_fq:
                            #print(entry.name)
                            if entry.name in subsampled_reads: #was subsampled
                                subsample_fq.write(str(entry)+'\n')    
                                read_log.write(entry.name+'\n') #iteration log
                                subsampled_count+=1
                            elif not entry.name in all_mapped_reads: #was not mapped to any contig yet
                                unmapped_fq.write(str(entry)+'\n')  
                                unmapped_count+=1
        reads_in_subsample = Helpers.count_reads_in_fastq(subsample_fastq)
        print(str(datetime.now())+" :: ReadSubset :: [write_fastq] wrote subsample to subsample.fastq, contains "+str(reads_in_subsample)+" reads,"+str(subsampled_count)+ \
        " NEW subsampled reads", flush=True)
        os.system("mv "+new_unmapped_fastq+" "+unmapped_fastq)
        reads_in_unmapped = Helpers.count_reads_in_fastq(unmapped_fastq)
        print(str(datetime.now())+" :: ReadSubset :: [write_fastq] wrote unmapped reads to unmapped.fastq, contains "+str(reads_in_unmapped)+" unmapped reads",flush=True)
    
    def get_all_subsampled_reads(self):
       reads = []
       for contig in self.mapped_to_contig_subsampled:
           for read in contig:
                reads.append(read)
       return reads
       
    def get_all_mapped_reads(self):
        reads = []
        for contig in self.mapped_to_contig_read_ids:
            for read_id in contig:
                reads.append(read_id)
        return reads 

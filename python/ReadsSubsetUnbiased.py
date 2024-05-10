import os, sys
import gzip
import pysam
import random

import time
from datetime import datetime
from Helpers import Helpers
from ReadsSubsetBase import ReadsSubsetBase

class ReadsSubsetBCommon(ReadsSubsetBase):
    def __init__(self, threshold, strategy, assembly_fasta, prob_dict_file):
        super().__init__(threshold,strategy,assembly_fasta)
        
        self.prob_dict = {}
        self.sorted_prob_tuples = []
        self.prob_dict_file = prob_dict_file
        self.contig_read_probs = {}
        self.sorted_prob_tuples_for_contig = []
        
        self.nx_p_for_contig = [0 for c in range(self.num_contigs)]

    def subsample(self,paf,link_dir,subsample_fastq, unmapped_fastq, out_txt):
        self.parse_paf(paf) #paf vyrobeny v self.minimap_to_asm()
        self.load_prob_dict()
        self.generate_p_for_new_reads()
        self.sort_dict()
        self.write_p_for_reads() #print dict to file
        self.probs_for_contig_reads()
        self.find_nx_p_subsample()
        self.write_fastq(link_dir, subsample_fastq, unmapped_fastq, out_txt)
    
    def find_nx_p_subsample(self):
        raise NotImplementedError
    
    def write_p_for_reads(self):
        print(str(datetime.now())+" :: ReadsSubset :: [write_p_for_reads] start")
        with open(self.prob_dict_file, "w+") as reads_prob:
            for read in self.prob_dict:
                print(read+"\t"+str(self.prob_dict[read]), file = reads_prob)
        print(str(datetime.now())+" :: ReadsSubset :: [write_p_for_reads] end")

    def load_prob_dict(self):
        self.prob_dict = {}
        with open(self.prob_dict_file, "r") as reads_prob:
            for line in reads_prob:
                line_split = line.strip().split("\t")
                self.prob_dict[line_split[0]]=float(line_split[1]) #prob_dict[read_id]=p
        print(str(datetime.now())+" :: ReadsSubset :: [load_prob_dict]  - dict size: ",len(self.prob_dict.keys()))


    def add_mapped_to_contig(self,read,read_id,contig):
        if not read_id in self.mapped_to_contig_read_ids[contig]: #without read duplicates!
            self.mapped_to_contig_read_ids[contig].append(read_id)  

    def get_all_subsampled_reads(self):
       reads = set() #lebo read sa moze zarovnat aj k viacerym contigom a nechceme ho v subsample viackrat.. bud spaja 2 contigy (ale nie je ich dost co by ich spajali), alebo mame 2 contigy, ktore sa prekryvaju (t.j. zarovnali by sa k sebe nemlou castou..) 
       num_reads_with_duplicates = 0
       for reads_list_for_contig in self.mapped_to_contig_subsampled:
           for read in reads_list_for_contig:
                num_reads_with_duplicates += 1
                reads.add(read)
       print(str(datetime.now())+" :: ReadsSubset :: [get_all_subsampled_reads] sum of subsampled reads for all contigs:",num_reads_with_duplicates," num of reads without duplicates:",len(reads))
       return reads
       
    def get_all_mapped_reads(self):
        reads = set() #lebo read sa moze zarovnat aj k viacerym contigom
        for contig_reads in self.mapped_to_contig_read_ids:
            for read_name in contig_reads:
                reads.add(read_name)
        return reads
    
    def write_fastq(self, reads_dir, subsample_fastq, unmapped_fastq, log_read_names): #write the selected reads to a new fastq for miniasm
        subsampled_reads = self.get_all_subsampled_reads()
        print(str(datetime.now())+" :: ReadSubset :: [write_fastq] NUM OF SUBSAMPLED READS:"+str(len(subsampled_reads)))
        #print(str(subsampled_reads))
        all_mapped_reads = self.get_all_mapped_reads()
        print(str(datetime.now())+" :: ReadSubset :: [write_fastq] NUM OF MAPPED READS:"+str(len(all_mapped_reads)))
        subsampled_count = 0
        unmapped_count = 0
        new_subsample_fastq = subsample_fastq+".new" #new subsampled reads
        new_unmapped_fastq = unmapped_fastq+".new"
        with open(new_subsample_fastq, "w+") as subsample_fq, open(new_unmapped_fastq,"w+") as unmapped_fq:
            with open(log_read_names, "w+") as read_log:
                for f in os.listdir(reads_dir):
                    fastq = reads_dir+f
                    print(str(datetime.now())+" :: ReadSubset :: [write_fastq] adding reads from "+fastq+" - contains ",Helpers.count_reads_in_fastq(fastq),"reads")
                    with pysam.FastxFile(fastq) as in_fq:
                        for entry in in_fq:
                            if entry.name in subsampled_reads: #was subsampled
                                subsample_fq.write(str(entry)+'\n')    
                                read_log.write(entry.name+'\n') #iteration log
                                subsampled_count+=1
                            elif not entry.name in all_mapped_reads: #was not mapped to any contig yet
                                unmapped_fq.write(str(entry)+'\n')  
                                #read_log.write(entry.name+'\n') #iteration log ..tieto nechceme v subsample.. 
                                unmapped_count+=1      
        os.system("mv "+new_subsample_fastq+" "+subsample_fastq)
        reads_in_subsample = Helpers.count_reads_in_fastq(subsample_fastq)
        print(str(datetime.now())+" :: ReadSubset :: [write_fastq] wrote subsample to subsample.fastq, contains "+str(reads_in_subsample)+" reads", flush=True)
        os.system("mv "+new_unmapped_fastq+" "+unmapped_fastq)
        reads_in_unmapped = Helpers.count_reads_in_fastq(unmapped_fastq)
        print(str(datetime.now())+" :: ReadSubset :: [write_fastq] wrote unmapped reads to unmapped.fastq, contains "+str(reads_in_unmapped)+" unmapped reads")


    def sort_dict(self):
        print(str(datetime.now())+" :: ReadSubset :: [sort_dict] start ")
        self.sorted_prob_tuples = sorted(self.prob_dict.items(), key=lambda x:x[1], reverse=True)
        #print("sorted prob tuples:", self.sorted_prob_tuples)
        self.sorted_prob_tuples_for_contig = [[] for i in range(self.num_contigs)]
        for contig in range(self.num_contigs):
            for tuple in self.sorted_prob_tuples:
                if tuple[0] in self.mapped_to_contig_read_ids[contig]:
                    self.sorted_prob_tuples_for_contig[contig].append(tuple)
            #print("sorted prob tuples for contig:", str(self.sorted_prob_tuples_for_contig[contig]), file=sys.stderr)
        print(str(datetime.now())+" :: ReadSubset :: [sort_dict] end ")
    
    def generate_p_for_new_reads(self):
        print(str(datetime.now())+" :: ReadSubset :: [generate_p_for_new_reads] start ")
        new_reads = 0
        reads_seen = 0
        for contig in range(self.num_contigs):
            for read_id in self.mapped_to_contig_read_ids[contig]:
                if not read_id in self.prob_dict:
                    new_reads += 1
                    self.prob_dict[read_id] = random.uniform(0, 1)
                else:
                    reads_seen += 1
        print(str(datetime.now())+" :: ReadsSubset :: [generate_p_for_new_reads] - num of new reads: "+str(new_reads)+" num of seen reads: "+str(reads_seen))

    def probs_for_contig_reads(self):
        for contig in range(self.num_contigs):
            count = 0
            self.contig_read_probs[contig] = {}
            for read_id in self.mapped_to_contig_read_ids[contig]:
                if not read_id in self.contig_read_probs[contig]:
                    self.contig_read_probs[contig][read_id] = self.prob_dict[read_id]
                    count += 1
            print(str(datetime.now())+" :: ReadsSubset :: [probs_for_contig_reads] num of reads for contig ",self.contig_names[contig]," before subsample:",str(count))


class ReadsSubsetB(ReadsSubsetBCommon):
    def __init__(self, threshold, strategy, assembly_fasta, prob_dict_file):
        super().__init__(threshold,strategy,assembly_fasta,prob_dict_file)  

    def find_nx_p_subsample(self):
        print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p_subsample] - contigs from longest to shortest:", self.contigs_from_longest_to_shortest)
        for contig_t in self.contigs_from_longest_to_shortest: # od najdlhsieho !!
            contig = contig_t[0]
            #use self.sorted_prob_tuples_for_contig[contig]..
            print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p_subsample] calculating coverage from longer contigs")
            map_sum_from_longer_contigs = self.cov_from_longer_contigs(contig) #calculate coverage generated from reads subsampled for lerger contigs
            cov_from_longer_contigs = map_sum_from_longer_contigs/self.contig_lengths[contig]
            target_sum = max(0,self.coverage_threshold * self.contig_lengths[contig] - map_sum_from_longer_contigs) #read lengths sum we want
            print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p_subsample] calculating prefix sums for contig")
            contig_prefix_sums, first_greater = self.prefix_sums_for_contig(contig,target_sum)           
            self.nx_p_for_contig[contig] = self.sorted_prob_tuples_for_contig[contig][first_greater][1]
            print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p] nx_p for contig (contig lenght:",self.contig_lengths[contig],")", self.contig_names[contig],":", self.nx_p_for_contig[contig],"; sum read length for prob (prefix sum [m]):",contig_prefix_sums[first_greater])
            self.read_length_dict_for_contig[contig][self.sorted_prob_tuples_for_contig[contig][0][0]]
            print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p] coverage from reads from longer contigs:",cov_from_longer_contigs,"(mappings length sum:", map_sum_from_longer_contigs,")")
            self.create_subsample_for_contig(contig)

    def prefix_sums_for_contig(self,contig,target_sum):
        contig_prefix_sums = []
        prob_tuples = self.sorted_prob_tuples_for_contig[contig] #(read_id, prob)
        #print("prob_tuples[0]:", str(prob_tuples[0]), file=sys.stderr)
        first_read_id = prob_tuples[0][0]
        contig_prefix_sums.append(self.read_length_dict_for_contig[contig][first_read_id]) #dlzka zarovnanych usekov pre prvy read + zarovnania pre ready vybrate pre ine contigy (nemalo by sa stat, ze neexistuje - ak vznikol contig, aspon jeden read tam musel byt..)
        for i in range(1,len(prob_tuples)):
            read = prob_tuples[i][0]
            read_length = self.read_length_dict_for_contig[contig][read]
            #print("i-1 = "+str(i-1),"len(read_ids)", str(len(read_ids)) ,file=sys.stderr)
            if not read in self.seen_reads_for_contig[contig]:
                contig_prefix_sums.append(contig_prefix_sums[i-1]+read_length)
            else:
                contig_prefix_sums.append(contig_prefix_sums[i-1]) #prefix sum unchanged (already counted)
            if contig_prefix_sums[-1]>=target_sum:    
                return contig_prefix_sums, len(contig_prefix_sums)-1
        return contig_prefix_sums, len(contig_prefix_sums)-1
    
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
    
    def create_subsample_for_contig(self, contig):
        self.mapped_to_contig_subsampled[contig] = []
        count = 0
        sum_lengths = 0
        for tuple in self.sorted_prob_tuples_for_contig[contig]:
            read_name = tuple[0]
            if self.prob_dict[read_name] >= self.nx_p_for_contig[contig]:
                if not read_name in self.seen_reads_for_contig[contig]: #add only new reads
                    self.mapped_to_contig_subsampled[contig].append(read_name)
                    count += 1
                    sum_lengths += self.read_length_dict_for_contig[contig][read_name]
        print(str(datetime.now())+" :: ReadsSubset :: [create_subsample] num of reads subsampled for contig ",self.contig_names[contig],"(length:",self.contig_lengths[contig],") :", str(count), " sum of lengths of selected reads ", str(sum_lengths), "=> coverage:", sum_lengths/self.contig_lengths[contig])


class ReadsSubsetBs(ReadsSubsetBCommon):
    def __init__(self, threshold, strategy, assembly_fasta, prob_dict_file):
        super().__init__(threshold,strategy,assembly_fasta,prob_dict_file)

    def find_nx_p_subsample(self):
        print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p_subsample] - contigs from longest to shortest:", self.contigs_from_longest_to_shortest)
        for contig_t in self.contigs_from_longest_to_shortest: # od najdlhsieho, ale tu by to malo byt jedno
            contig = contig_t[0]
            #use self.sorted_prob_tuples_for_contig[contig]..
            print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p_subsample] finished calculating prefix sums for contig")
            target_sum = self.coverage_threshold * self.contig_lengths[contig] #read lengths sum we want
            print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p_subsample] calculating prefix sums for contig")
            contig_prefix_sums, first_greater = self.prefix_sums_for_contig(contig,target_sum)
            self.nx_p_for_contig[contig] = self.sorted_prob_tuples_for_contig[contig][first_greater][1]
            print(str(datetime.now())+" :: ReadsSubset :: [find_nx_p] nx_p for contig (contig lenght:",self.contig_lengths[contig],")", self.contig_names[contig],":", self.nx_p_for_contig[contig],"; sum read length for prob (prefix sum [m]):",contig_prefix_sums[first_greater])
            self.read_length_dict_for_contig[contig][self.sorted_prob_tuples_for_contig[contig][0][0]]
            self.create_subsample_for_contig(contig)

    def prefix_sums_for_contig(self,contig,target_sum):
        contig_prefix_sums = []
        prob_tuples = self.sorted_prob_tuples_for_contig[contig] #(read_id, prob)
        #print("prob_tuples[0]:", str(prob_tuples[0]), file=sys.stderr)
        first_read_id = prob_tuples[0][0]
        contig_prefix_sums.append(self.read_length_dict_for_contig[contig][first_read_id]) #dlzka zarovnanych usekov pre prvy read + zarovnania pre ready vybrate pre ine contigy (nemalo by sa stat, ze neexistuje - ak vznikol contig, aspon jeden read tam musel byt..)
        for i in range(1,len(prob_tuples)):
            read = prob_tuples[i][0]
            read_length = self.read_length_dict_for_contig[contig][read]
            contig_prefix_sums.append(contig_prefix_sums[i-1]+read_length)
            if contig_prefix_sums[-1]>=target_sum:    
                return contig_prefix_sums, len(contig_prefix_sums)-1
        return contig_prefix_sums, len(contig_prefix_sums)-1
    
    def create_subsample_for_contig(self, contig):
        self.mapped_to_contig_subsampled[contig] = []
        count = 0
        sum_lengths = 0
        for tuple in self.sorted_prob_tuples_for_contig[contig]:
            read_name = tuple[0]
            if self.prob_dict[read_name] >= self.nx_p_for_contig[contig]:
                self.mapped_to_contig_subsampled[contig].append(read_name)
                count += 1
                sum_lengths += self.read_length_dict_for_contig[contig][read_name]
        print(str(datetime.now())+" :: ReadsSubset :: [create_subsample] num of reads subsampled for contig ",self.contig_names[contig],"(length:",self.contig_lengths[contig],") :", str(count), " sum of lengths of selected reads ", str(sum_lengths), "=> coverage:", sum_lengths/self.contig_lengths[contig])


class ReadsSubsetLongReadB(ReadsSubsetB):
    def __init__(self, threshold, strategy, assembly_fasta, prob_dict_file):
        super().__init__(threshold,strategy,assembly_fasta,prob_dict_file)

    def generate_p_for_new_reads(self):
        new_reads = 0
        reads_seen = 0
        for contig in range(self.num_contigs):
            for read_id in self.mapped_to_contig_read_ids[contig]:
                if not read_id in self.prob_dict:
                    new_reads += 1
                    self.prob_dict[read_id] = self.read_length_dict_for_contig[contig][read_id] #READ LENGTH
                else:
                    reads_seen += 1
        print(str(datetime.now())+" :: ReadsSubset :: [generate_p_for_new_reads] - num of new reads: "+str(new_reads)+" num of seen reads: "+str(reads_seen))


class ReadsSubsetLongReadBs(ReadsSubsetBs):
    def __init__(self, threshold, strategy, assembly_fasta, prob_dict_file):
        super().__init__(threshold,strategy,assembly_fasta,prob_dict_file)

    def generate_p_for_new_reads(self):
        new_reads = 0
        reads_seen = 0
        for contig in range(self.num_contigs):
            for read_id in self.mapped_to_contig_read_ids[contig]:
                if not read_id in self.prob_dict:
                    new_reads += 1
                    self.prob_dict[read_id] = self.read_length_dict_for_contig[contig][read_id] #READ LENGTH
                else:
                    reads_seen += 1
        print(str(datetime.now())+" :: ReadsSubset :: [generate_p_for_new_reads] - num of new reads: "+str(new_reads)+" num of seen reads: "+str(reads_seen))

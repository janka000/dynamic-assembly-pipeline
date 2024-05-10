import os, sys
import gzip
import pysam
import random

import time
from datetime import datetime
from Helpers import Helpers
from ReadsSubset import ReadsSubset

class ReadsSubsetMapToAssembly(ReadsSubset):
    def __init__(self, threshold, strategy, assembly_fasta, it):
        super().__init__(threshold, strategy, assembly_fasta, it)

        #containment
        self.overhangs = [[] for c in range(self.num_contigs)]
        
    def coverage_for_contig(self, contig):
        if self.it > 0:
            map_sum_from_prev_it = self.mapped_from_prev_subsample_for_contig[contig] #coverage generated from mappings for reads that already were in subsample from previous iteration
        else:
            map_sum_from_prev_it = 0
        cov_from_prev_it = map_sum_from_prev_it/self.contig_lengths[contig]
        return cov_from_prev_it

    #-------------------------------------
    # map to assembly loop methods
    #------------------------------------    

    def parse_paf_for_map_to_assembly(self, paf):
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
        with gzip.open(paf,'rt') as pf:      
            for line in pf:        
                paf_row = line.split()
                read_name = paf_row[0]
                qs = int(paf_row[2])
                qe = int(paf_row[3])-1
                ql = int(paf_row[1])
                ts = int(paf_row[7])
                te = int(paf_row[8])-1
                tl = int(paf_row[6])
                read = [read_name,qs,qe,ql,ts,te,tl] #[read name, query_start, query_end, query_len, target_start, target_end, target_len]
                read_id = paf_row[0]
                contig = self.ctg_name_to_num[paf_row[5]] #number
                self.add_mapped_to_contig(read,read_id,contig)

                read_start = int(paf_row[7])
                read_end = int(paf_row[8])-1
                read_length = read_end-read_start #mapped read part
                if not read_id in self.read_length_dict_for_contig[contig]:
                    self.read_length_dict_for_contig[contig][read_id] = read_length
                else:
                    self.read_length_dict_for_contig[contig][read_id] += read_length
            for contig in range(self.num_contigs):
                print(str(datetime.now())+" :: ReadsSubset :: [parse_paf] :: "+str(len(self.mapped_to_contig_read_ids[contig]))+" reads mapped to contig "+self.contig_names[contig])


    def get_all_overhangs(self):
        reads = []
        for contig_overhangs_list in self.overhangs:
            reads+=contig_overhangs_list
        return reads

    #---
    #main map to assembly loop methods
    #---
    def find_contig_overhangs(self):
        #all_mapped = self.get_all_mapped_reads()
        for contig in range(len(self.mapped_to_contig)):
            coverage = self.coverage_for_contig(contig)
            thr = int(self.coverage_threshold+0.5*self.coverage_threshold) #50% above threshold
            take_all = coverage < thr
            print(str(datetime.now())+" :: ReadSubset :: [find_contig_overhangs] coverage for contig "+self.contig_names[contig]+" is "+str(coverage)+" - threshold is "+str(self.coverage_threshold) + " - take all reads: "+str(take_all))
            num_mapped_reads = len(set([read[0] for read in self.mapped_to_contig[contig]]))
            print(str(datetime.now())+" :: ReadSubset :: [find_contig_overhangs] num. mappings: "+str(len(self.mapped_to_contig[contig]))+", num mapped reads: "+str(num_mapped_reads))
            contig_overhangs = set()
            if take_all:
                for read in self.mapped_to_contig[contig]:
                    contig_overhangs.add(read[0])
                print(str(datetime.now())+" :: ReadSubset :: [find_contig_overhangs] found "+str(len(contig_overhangs))+" reads for "+self.contig_names[contig]+" (taking all)")
            else: #overhangs only
                """
                for read in self.mapped_to_contig[contig]:
                    if (read[1]==0 or read[2]>=self.contig_lengths[contig]-1):
                        contig_overhangs.add(read[0])
                """
                left_overhangs = 0
                right_overhangs = 0
                right_overhangs_for_contig = []
                left_overhangs_for_contig = []
                for read in self.mapped_to_contig[contig]:
                    read_name, query_start, query_end, query_len, target_start, target_end, target_len = read
                    if query_start > target_start:
                        left_overhangs += 1
                        left_overhangs_for_contig.append(read)
                    elif target_len - target_end < query_len - query_end:
                        right_overhangs += 1
                        right_overhangs_for_contig.append(read)
                sorted_left_overhangs_for_contig = sorted(left_overhangs_for_contig, reverse=True, key=lambda x: abs(x[3]))
                sorted_right_overhangs_for_contig = sorted(right_overhangs_for_contig,  reverse=True, key=lambda x: abs(x[3]))
                print(str(datetime.now())+" :: ReadSubset :: [find_contig_overhangs] num. of left overhangs: "+str(left_overhangs)+", num. of right overhangs: "+str(right_overhangs))
                #mapped but not contained (TODO: add 3,4 form paf to also tell if there really is overlap, i.e. it is not the case when the read starts/ends exactly on the start/end of the contig)
                #if (read[1]==0 or read[2]>=self.contig_lengths[contig]-1):  #no need to check again
                thr = int(self.coverage_threshold+0.5*self.coverage_threshold) #slightly above threshold
                for i in range(thr): #add thr overhangs for each end of the contig (sorted from longest to shortest)
                    if i < len(sorted_left_overhangs_for_contig):
                        contig_overhangs.add(sorted_left_overhangs_for_contig[i][0]) #add read name to set
                    if i < len(sorted_right_overhangs_for_contig):
                        contig_overhangs.add(sorted_right_overhangs_for_contig[i][0]) #add read name to set
                print(str(datetime.now())+" :: ReadSubset :: [find_contig_overhangs] found "+str(len(contig_overhangs))+" reads for "+self.contig_names[contig]+" (taking overhangs only)")
            self.overhangs[contig] = list(contig_overhangs)
    
    def write_contig_overhangs_and_unmapped(self, reads_dir, subsample_fastq):
        overhangs = self.get_all_overhangs()
        all_mapped = self.get_all_mapped_reads()
        #overhangs = all_mapped #test: try to keep all mapped reads
        overhang_count = 0
        unmapped_count = 0
        with open(subsample_fastq, "a+") as out_fq:
            for f in os.listdir(reads_dir):
                fastq = reads_dir+f
                print(str(datetime.now())+" :: ReadSubset :: [write_contig_overhangs_and_unmapped] adding reads from "+fastq)
                with pysam.FastxFile(fastq) as in_fq:
                    for entry in in_fq:
                        if entry.name in overhangs: #overhang     
                            out_fq.write(str(entry)+'\n')
                            overhang_count+=1
                        elif not entry.name in all_mapped: #not mapped
                            out_fq.write(str(entry)+'\n')
                            unmapped_count+=1
        print(str(datetime.now())+" :: ReadSubset :: [write_contig_overhangs_and_unmapped] found "+str(overhang_count)+" overhangs and "+str(unmapped_count)+" new unmapped reads")
    #-------------------------------------
    # end of map to assembly loop methods
    #------------------------------------ 

#common in all readsubset classes
import os, sys
import gzip
import pysam
import random
import time
from datetime import datetime
from Helpers import Helpers

class ReadsSubsetBase:
    def __init__(self, threshold, strategy, assembly_fasta):
        #read from .fasta assembly
        self.num_contigs = 0
        self.contig_lengths = []
        self.contig_names = []
        self.ctg_name_to_num = {}
        self.contigs_from_longest_to_shortest = []
        self.get_contig_stats(assembly_fasta)
        
        #store all reads mapped to contig
        self.mapped_to_contig_read_ids = [[] for c in range(self.num_contigs)] #read names only
        self.read_length_dict_for_contig = [{} for c in range(self.num_contigs)] #for each contig a dictionary of reads aligned to it, with SUM of lengths of aligned regions for that contig   
        self.mapped_to_contig_subsampled = [[] for c in range(self.num_contigs)] #only read names
        self.subsample_position_coverages = [[0 for i in range(self.contig_lengths[c])] for c in range(self.num_contigs)]
        self.seen_reads_for_contig = [set() for c in range(self.num_contigs)] #reads mapped to contig that were already mapped to longer one (avoid duplicates)
        self.sum_read_lengths_for_contigs_sq = [0 for c in range(self.num_contigs)]

        self.subsample_strategy = strategy
        self.coverage_threshold = threshold #min. subsample coverage 
        
    def get_contig_stats(self, assembly_fasta):   
        while not os.path.exists(assembly_fasta): #wait while assembly does not exist
            time.sleep(1)
            print(str(datetime.now())+" :: ReadSubset :: [get_contig_stats] waiting for assembly...", flush=True)
        #if not os.path.exists(assembly_fasta+".fai"): #samtools faidx will not create index on empty file
        if os.stat(assembly_fasta).st_size == 0: #fasta is empty
            self.contig_lengths = []
            self.num_contigs = 0
            self.contig_names = []
            self.ctg_name_to_num = {}
            return
        else: #fasta exists and is not empty
            while not os.path.exists(assembly_fasta+".fai"):  #wait for index - should be created from non-empty file..
                time.sleep(1)
                print(str(datetime.now())+" :: ReadSubset :: [get_contig_stats] waiting for index...")
        with pysam.FastaFile(assembly_fasta) as fa:
            self.contig_lengths = list(fa.lengths)
            self.num_contigs = fa.nreferences
            self.contig_names = list(fa.references)
        i=0
        for contig in self.contig_names:
            self.ctg_name_to_num[contig]=i
            i+=1

        contigs_with_lengths = [(int(i), int(self.contig_lengths[i])) for i in range(self.num_contigs)]
        self.contigs_from_longest_to_shortest = sorted(contigs_with_lengths, key=lambda x: x[1], reverse=True)
        
        
    def parse_paf(self, paf): #paf is created only from new reads aligned to assembly,assembly was created using both new and old subsampled reads
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
                read = [paf_row[0], int(paf_row[7]), int(paf_row[8])-1] #[read name, start, end]
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
                self.sum_read_lengths_for_contigs_sq[contig]+=read_length**2
            for contig in range(self.num_contigs):
                print(str(datetime.now())+" :: ReadsSubset :: [parse_paf] :: "+str(len(self.mapped_to_contig_read_ids[contig]))+" reads mapped to contig "+self.contig_names[contig])

    def add_mapped_to_contig(self,mapping,read_id,contig):
        raise NotImplementedError
    
    def write_fastq(self, reads_dir, subsample_fastq, unmapped_fastq, log_read_names): #write the selected reads to a new fastq for miniasm
        raise NotImplementedError

    def get_all_subsampled_reads(self):
        raise NotImplementedError
       
    def get_all_mapped_reads(self):
        raise NotImplementedError


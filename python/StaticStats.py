import os
import pysam
import time
import json
from datetime import datetime
import random
import string

class StaticStats(object):
    
### with reference ###
    @staticmethod
    def align_assembly_to_ref(assembly_fa, reference, out_paf):
        """
        assembly_fa : assembly form iteration, filtered or not
        reference
        out_paf : where to write minimap output
        """
        print("["+str(datetime.now())+"] align assembly to ref")
        os.system("minimap2 -x map-ont -t 8 --secondary=no "+reference+" "+assembly_fa+" 2> /dev/null > "+out_paf) # -c to print cigar strings
    
    
    @staticmethod
    def percent_aligned_length(paf, reference_length):
        """
        kolko percent z celkovej dlzky assembly je sum(dlzky zarovnani)
        assembly_fa : assembly form iteration, filtered or not
        reference
        paf : alignemnt of assembly_fa to reference
        """
        print("["+str(datetime.now())+"] percent_aligned_length")
        aligned_lengths = StaticStats.count_aligned_lengths(paf)
        percent = sum(aligned_lengths)/reference_length*100
        return percent

    @staticmethod
    def filter_out_unaligned_contigs(assembly_fa, not_filtered_paf, filtered_assembly_fa, filtered_paf):
        """
        odstranit contigy, kotre sa nezarovnali k referencii (resp. zarovnania su kratke (<1000))
        assembly_fa : not filtered
        referecne
        paf : alignemnt of assembly_fa to reference
        filtered_assembly_fa : output file
        """
        print("["+str(datetime.now())+"] filter_out_unaligned_contigs")
        threshold = 1000

        contigs_ok = set() #is ok, if alignment of size > threshold exists

        with open(not_filtered_paf,"r") as nfpaf:
            for line in nfpaf:
                split_line = line.strip().split()
                #1	string	Query sequence name
                contig = split_line[0]
                #8	int	Target start on original strand (0-based)
                #9	int	Target end on original strand (0-based)
                alignment_length = abs(int(split_line[8])-int(split_line[7]))
                if alignment_length > threshold: 
                    contigs_ok.add(contig)

        with open(not_filtered_paf,"r") as nfpaf:
            with open(filtered_paf, "w+") as fpaf:
                for line in nfpaf:
                    #1	string	Query sequence name
                    contig = line.strip().split()[0]
                    if contig in contigs_ok:
                        print(line.strip(), file=fpaf)

        write_next_line = False
        with open(assembly_fa, "r") as fa_in:
            with open(filtered_assembly_fa, "w+") as fa_out:
                for line in fa_in:
                    if line[0]==">":
                        contig = line[1:].strip()
                        if contig in contigs_ok:
                            print(line.strip(), file = fa_out) #write contig name
                            write_next_line = True
                        else:
                            write_next_line = False
                    elif write_next_line:
                        print(line.strip(), file = fa_out) #write sequence
        

    def n50(filtered_assembly, reference_length):
        """
        pocita N50 score - min. dlzka contigu, pre ktory vsetky dlhsie/rovnako dlhe dokopy pokryvaju 50% assembly 
        filtered assembly
        reference length (sum)
        """
        print("["+str(datetime.now())+"] n50")
        if os.path.getsize(filtered_assembly) == 0:
            return 0 #no assembly

        os.system("samtools faidx "+filtered_assembly)
        asm_lengths = []
        with pysam.FastaFile(filtered_assembly) as asm:
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
    
    def count_aligned_lengths(paf):
        print("["+str(datetime.now())+"] count_aligned_lengths")
        aligned_lengths = []
        with open(paf, "r") as pf:
            for line in pf:
                #print("paf line: ",line)
                split_line = line.split()
                #3	int	Query start (0-based; BED-like; closed)
                query_start = int(split_line[2])
                #4	int	Query end (0-based; BED-like; open)
                query_end = int(split_line[3])
                alignment_len = abs(query_end - query_start)
                aligned_lengths.append(alignment_len)
        return aligned_lengths
    
    def corrN50(filtered_paf, reference_length):
        """
        corrected N50 - podobne ako N50, ale namiesto dlzok contigov berieme dlzky zarovnani
        """
        print("["+str(datetime.now())+"] corrN50")
        aligned_lengths = StaticStats.count_aligned_lengths(filtered_paf)
        print("["+str(datetime.now())+"] corrN50 ..")
        n_num = 50
        fraction = n_num/100
        if len(aligned_lengths)==0:
            return 0 #no asssmebly
        sorted_aln_lengths = sorted(aligned_lengths, reverse=True)
        s = 0
        for l in sorted_aln_lengths:
            if s+l > reference_length*fraction:
                return l
            else: # s+l <= reference_length*fraction
                s += l
        return l
    
    def percent_covered_by_zero_or_more_than_one_contig(reference, filtered_paf):
        """
        % of assembly length covered by more than 1 contig and by 0 contigs
        filtered_assembly : assembly form iteration, filtered 
        reference
        paf : alignemnt of assembly_fa to reference
        """
        print("["+str(datetime.now())+"] percent_covered_by_zero_or_more_than_one_contig")
        ref_contigs = []
        contig_lengths = {}
        total_ref_length = 0
        with pysam.FastaFile(reference) as ref:
            ref_contigs = list(ref.references)
            contig_lengths_arr = list(ref.lengths)
            for i in range(len(ref_contigs)):
                contig_lengths[ref_contigs[i]] = contig_lengths_arr[i]
                total_ref_length += contig_lengths_arr[i]

        counts = {}
        for contig in ref_contigs:
            counts[contig] = [0 for i in range(contig_lengths[contig])]

        with open(filtered_paf, "r") as pf:
            for line in pf:
                split_line = line.strip().split()
                #6	string	Target sequence name
                contig = split_line[5]
                #8	int	Target start on original strand (0-based)
                ref_start = int(split_line[7])
                #9	int	Target end on original strand (0-based)
                ref_end =  int(split_line[8])
                
                for i in range(ref_start,ref_end):
                    counts[contig][i] += 1

        pos_more_than_one = 0
        pos_zero = 0
        for contig in ref_contigs:
            for p in counts[contig]:
                if p > 1: 
                    pos_more_than_one +=1
                elif p == 0:
                    pos_zero +=1
        
        percent_zero = pos_zero / total_ref_length *100
        percent_more_than_one = pos_more_than_one / total_ref_length *100

        return percent_zero, percent_more_than_one

### with/out reference ### (applicable if no reference given)
    @staticmethod
    def sum_contig_lengths(assembly_fa):
        """
        assembly fa: not filtered if no reference given, filtered or not if reference is given
        """
        print("["+str(datetime.now())+"] sum_contig_lengths")
        if os.path.getsize(assembly_fa) == 0:
            return 0 #no assembly

        ref_lengths = []
        with pysam.FastaFile(assembly_fa) as ref:
            ref_lengths= list(ref.lengths)
        return sum(ref_lengths)


### reads stats ###

    @staticmethod
    def reads_list_to_fastq(reads_list, fastq_all, fastq_out):
        """
        create fastq form list of read ids
        reads_list: list of ids in txt, one id on one line
        fastq_all : fastq file that contains all the reads in reads list (fastq to filter from)
        fastq_out : where to write output fastq
        """ 
        print("["+str(datetime.now())+"] reads_list_to_fastq")
        os.system("cat "+fastq_all+" | grep -F -f "+reads_list+" -A 3 - | grep -v '^--$' > "+fastq_out)
        
    @staticmethod
    def reads_to_assembly_bam(reads,assembly,bam_out):
        """
        align reads to assembly, outut in bam format
        reads: fastq with reads to align
        assembly : asssembly or reference (if given)
        bam_out: where to write minimap output
        """
        print("["+str(datetime.now())+"] reads_to_assembly_bam")
        minimap_command = "minimap2 -ax map-ont -t 8 --secondary=no "+assembly+" "+reads+" 2> /dev/null | samtools view -S -b -o - | samtools sort - -o "+bam_out
        os.system(minimap_command)
        os.system("samtools index "+bam_out)


    @staticmethod
    def coverge_per_base(assembly, bam):
        """
        count coverage per base array (merged for whole reference)
        assembly : assembly or reference (if given)
        bam from minimap, reads aligned to reference/assembly
        """
        print("["+str(datetime.now())+"] coverge_per_base")
        if os.path.getsize(assembly) == 0:
            return [0] #no assembly


        with pysam.FastaFile(assembly) as ref:
            contig_names = list(ref.references)
        contig_coverages = []
        with pysam.AlignmentFile(bam, "rb", check_sq=False) as aln:
            for ctg in contig_names:
                #contig_coverages[ctg] = []
                #num_aligned_reads[ctg] = al.count(contig=ctg)
                covA, covC, covG, covT = aln.count_coverage(contig=ctg,quality_threshold=0)
                len_cov = len(covA)
                sum_coverage_per_pos_ctg = [covA[i]+covC[i]+covG[i]+covT[i] for i in range(len_cov-1)]
                contig_coverages += sum_coverage_per_pos_ctg
        return contig_coverages
    
    """
    @staticmethod
    def merge_coverage_arrays(per_base_coverage_for_contigs):
        data = []
        for contig in per_base_coverage_for_contigs:
            data+=contig
        return data
    """
    
    @staticmethod
    def assembly_avg_coverage(per_base_coverage):
        print("["+str(datetime.now())+"] assembly_avg_coverage")
        mean = sum(per_base_coverage)/len(per_base_coverage)
        return mean

    @staticmethod
    def assembly_coverage_variance(per_base_coverage, mean):
        """
        per_base_coverage: [c1_{1},...c1_{contig_length}]+[c2_{1},...c2_{contig_length}]+...
        """
        print("["+str(datetime.now())+"] assembly_coverage_variance")
        sum = 0
        for c in per_base_coverage:
            sum += (c - mean)**2
        variance = sum / len(per_base_coverage)
        return variance
    
    @staticmethod
    def count_bases_in_fastq(fastq):
        """
        fastq : subsample or all processed reads
        """
        print("["+str(datetime.now())+"] count_bases_in_fastq")
        line_num = 0
        num_bases = 0
        with open(fastq, 'r') as fp:
            for line in fp:
                if line_num % 4 == 1:
                    num_bases += len(line.strip())
                line_num += 1
        return num_bases

    @staticmethod
    def get_list_processed(in_dir, stats_json):
        print("["+str(datetime.now())+"] get_list_processed")
        def parse_json(json_file):
            content = ""
            with open(json_file, "r") as js:
                content = js.read()
            dict = json.loads(content)
            file_list = dict["processed_files"]
            return file_list
        files = parse_json(stats_json)
        ret_files = [in_dir+f for f in files]
        return ret_files

#### put it all together ###
    @staticmethod
    def compute_numeric_stats_for_it(fastq_seen_reads, subsample_read_ids, assembly, reference=None):
        """
        compute all the numbers for files from ONE ITERATION
        fastq_seen_reads : fastq with all reads seen in before and in that iteration
        subsample_read_ids : read ids of reads in subsample
        assembly : path to file with assembly
        reference : path to file with reference
        """
        print("["+str(datetime.now())+"] ========== compute_numeric_stats_for_it ==========")

        h = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
        timestr = time.strftime("%Y%m%d-%H%M%S")+h

        no_ref_stats = None
        ref_stats = None

        def coverage_stats(fq,asm):
            bamaln = "/tmp/dynamic_assembly__staticstats_"+timestr+"_bamaln.bam"
            StaticStats.reads_to_assembly_bam(fq,asm,bamaln)
            if not os.path.exists(bamaln):
                print("bam alignment does not exist!!")
                return 0,0
            per_base_coverage = StaticStats.coverge_per_base(asm, bamaln)
            avg_coverage = StaticStats.assembly_avg_coverage(per_base_coverage)
            coverage_variance = StaticStats.assembly_coverage_variance(per_base_coverage, avg_coverage)
            os.unlink(bamaln)
            os.unlink(bamaln+".bai")
            return avg_coverage,coverage_variance

        #can be computed without reference:
        subsample_fastq = "/tmp/dynamic_assembly__staticstats_"+timestr+"_subsample.fastq"
        StaticStats.reads_list_to_fastq(subsample_read_ids, fastq_seen_reads, subsample_fastq)
        assembly_length = StaticStats.sum_contig_lengths(assembly)
        sum_subsample_read_lengths = StaticStats.count_bases_in_fastq(subsample_fastq)
        sum_seen_reads_lengths = StaticStats.count_bases_in_fastq(fastq_seen_reads)
        asm_avg_coverage_all = 0 if assembly_length == 0 else sum_seen_reads_lengths/assembly_length
        asm_coverage_variance_all = 0 #coverage_stats(fastq_seen_reads, assembly)
        asm_avg_coverage_subsample = 0 if assembly_length == 0 else sum_subsample_read_lengths/assembly_length
        asm_coverage_variance_subsample = 0 #coverage_stats(subsample_fastq, assembly)

        no_ref_stats = [assembly_length, \
                        sum_subsample_read_lengths, sum_seen_reads_lengths, \
                        asm_avg_coverage_all, asm_coverage_variance_all, \
                            asm_avg_coverage_subsample, asm_coverage_variance_subsample]

        if reference is not None: # with reference
            paf = "/tmp/dynamic_assembly__staticstats_"+timestr+"_asm_to_ref.paf"
            reference_length = StaticStats.sum_contig_lengths(reference)
            StaticStats.align_assembly_to_ref(assembly, reference, paf)
            all_percent_aligned_length = StaticStats.percent_aligned_length(paf, reference_length)
            ref_avg_coverage_all, ref_coverage_variance_all = sum_seen_reads_lengths/reference_length, 0 #coverage_stats(fastq_seen_reads, reference)
            ref_avg_coverage_subsample, ref_coverage_variance_subsample = sum_subsample_read_lengths/reference_length, 0 #coverage_stats(subsample_fastq, reference)
            
            #najdi zarovnane
            filtered_paf = "/tmp/dynamic_assembly__staticstats_"+timestr+"_filtered_asm_to_ref.paf"
            filtered_assembly = "/tmp/dynamic_assembly__staticstats_"+timestr+"_filtered_assemly.fa"
            StaticStats.filter_out_unaligned_contigs(assembly, paf, filtered_assembly, filtered_paf)
            
            #zo zarovnanych
            filtered_assembly_length = StaticStats.sum_contig_lengths(filtered_assembly)
            filtered_percent_aligned_length = StaticStats.percent_aligned_length(filtered_paf, reference_length)
            n50 = StaticStats.n50(filtered_assembly, reference_length)
            corrN50 = StaticStats.corrN50(filtered_paf, reference_length)
            percent_cov_by_zero_ctg, percent_cov_by_more_than_one_ctg = StaticStats.percent_covered_by_zero_or_more_than_one_contig(reference, filtered_paf)

            ref_stats = [reference_length, filtered_assembly_length, \
                         all_percent_aligned_length, filtered_percent_aligned_length, \
                         ref_avg_coverage_all, ref_coverage_variance_all, \
                         ref_avg_coverage_subsample, ref_coverage_variance_subsample, \
                         n50, corrN50, \
                         percent_cov_by_zero_ctg, percent_cov_by_more_than_one_ctg]

        os.unlink(subsample_fastq)
        os.unlink(paf)
        os.unlink(filtered_paf)
        if not os.path.getsize(filtered_assembly) == 0:
            os.unlink(filtered_assembly+".fai")
        os.unlink(filtered_assembly)
        return no_ref_stats, ref_stats            
            
    @staticmethod
    def compute_numeric_stats_for_all(run_in_dir, run_out_dir, it_start, it_end, result_tsv, reference=None): 
        print("computing stats "+result_tsv+" ...")
        with open(result_tsv, "w+") as res:
            header = "iteration\tassembly_length\tsum_subsample_read_lengths\tsum_seen_reads_lengths\tasm_avg_coverage_all\tasm_coverage_variance_all\tasm_avg_coverage_subsample\tasm_coverage_variance_subsample\treference_length\tfiltered_assembly_length\tall_percent_aligned_length\tfiltered_percent_aligned_length\tref_avg_coverage_all\tref_coverage_variance_all\tref_avg_coverage_subsample\tref_coverage_variance_subsample\tn50\tcorrN50\tpercent_cov_by_zero_ctg\tpercent_cov_by_more_than_one_ctg"
            print(header, file = res)
        for i in range(it_start, it_end+1):
            print("computing stats for iteration "+str(i))
            json_stats_file = run_out_dir+"stats/without_reference/json_summary_per_it/stats_"+str(i)+".json"
            if not os.path.exists(json_stats_file):
                print("iteration "+str(i)+" does not exist")
                break
            it_files = StaticStats.get_list_processed(run_in_dir, json_stats_file)
            
            h = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
            timestr = time.strftime("%Y%m%d-%H%M%S")+h
            
            fastq_seen_reads = "/tmp/dynamic_assembly__staticstats_seen_reads_"+str(i)+timestr+".fastq"
            os.system("cat "+(" ").join(it_files)+" > "+fastq_seen_reads)
            subsample_read_ids = run_out_dir+"subsample_filtered/subsampling_filtered_"+str(i)+".txt"
            assembly = run_out_dir+"assemblies/assembly_"+str(i)+".fasta"
            no_ref_stats_for_it, ref_stats_for_it = StaticStats.compute_numeric_stats_for_it(fastq_seen_reads, subsample_read_ids, assembly, reference=reference)
            row = str(i)+"\t"+("\t").join([str(s) for s in no_ref_stats_for_it])+"\t"+("\t").join([str(s) for s in ref_stats_for_it])
            with open(result_tsv,"a+") as res:
                print(row, file=res)
            os.unlink(fastq_seen_reads)

"""
minimap <target> <query>
paf:
1	string	Query sequence name
2	int	Query sequence length
3	int	Query start (0-based; BED-like; closed)
4	int	Query end (0-based; BED-like; open)
5	char	Relative strand: "+" or "-"
6	string	Target sequence name
7	int	Target sequence length
8	int	Target start on original strand (0-based)
9	int	Target end on original strand (0-based)
10	int	Number of residue matches
11	int	Alignment block length
12	int	Mapping quality (0-255; 255 for missing)
"""

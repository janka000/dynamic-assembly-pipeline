import pysam
import os
import time
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class ReferenceStats(object):

    @staticmethod
    def dotplot_with_coverage(paf_asm_to_asm, coverage_threshold, save_path, lock, lock_log, it, log, bam_reads_to_asm=None):
        print(str(datetime.now())+" :: ReferenceStats :: [dotplot with coverage] drawing dotplot for "+paf_asm_to_asm+", save to: "+save_path, flush=True, file=log)
        
        
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

    @staticmethod
    def reference_coverage(ref, it, subsample_fastq, bed, out_png, out_bam, minimapobj, lock, lock_log, log):
        def minimap(final_assembly, fastq, out_bam):
            minimapobj.map_reads_to_assembly(fastq,final_assembly,out_bam)
            return 0
        print(str(datetime.now())+" :: ReferenceStats :: [reference coverage] aligning subsample to ref.", flush=True, file = log)
        if not os.path.exists(out_bam):
            minimap(ref,subsample_fastq,out_bam)
        else:
            print("bam file already exists, skipping",flush=True, file = log)
        print(str(datetime.now())+" :: ReferenceStats :: [reference coverage] generating plotting data.", flush=True, file = log)
        
        bam = pysam.AlignmentFile(out_bam, 'rb')
        contigs = bam.references
        contig_lengths = bam.lengths
        sum_coverage_per_pos = {}
        for contig_id in contigs:
            A,C,G,T = bam.count_coverage(contig_id)
            coverage = [a+c+g+t for a,c,g,t in zip(A,C,G,T)]
            sum_coverage_per_pos[contig_id] = coverage
       
        if lock.acquire():
            with open(lock_log, "a+") as llog:
                print("thread for it. " +str(it)+ " acquired lock in [reference_coverage]", flush=True, file=llog)
        else:
            with open(lock_log, "a+") as llog:
                print("thread for it. " +str(it)+ " lock.acquire() returned false in [reference_coverage]", flush=True, file=llog)

        try:
            plt.clf()
            if len(contigs) > 1:
                fig, ax = plt.subplots(nrows = len(contigs), figsize=(10, len(contigs)*5))
                #print(ax)
                i = 0
                labels = contigs
                lengths = contig_lengths
                labels_lengths = zip(labels, lengths)
                labels_lengths = sorted(labels_lengths, key=lambda x: x[1], reverse = True)
                labels, lengths = zip(*labels_lengths)
                sorted_contigs = labels
                for ctg in sorted_contigs:
                    #col = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
                    col = "blue"
                    #ma = moving_average(self.sum_coverage_per_pos[ctg],3) #too slow
                    ax[i].plot(sum_coverage_per_pos[ctg], label=ctg, color=col, linewidth=0.5)
                    ax[i].legend()
                    i+=1
                #if len(self.contig_lengths)<=15: #draw legend only if there are not too many contigs
                #    plt.legend() 
                plt.savefig(out_png, bbox_inches='tight',dpi=100)
                #plt.ylim(0, 50) 
                #plt.savefig(fname_ylim, bbox_inches='tight',dpi=100)
            elif len(contigs) == 1:
                fig, ax = plt.subplots(nrows = 1, figsize=(10, 5))
                col = "blue"
                ctg = contigs[0]
                ax.plot(sum_coverage_per_pos[ctg], label=ctg, color=col, linewidth=0.5)
                ax.legend()
                plt.savefig(out_png, bbox_inches='tight',dpi=100)
            else: # 0 contigs
                plt.savefig(out_png, bbox_inches='tight',dpi=1, figsize=(10,10))
                #plt.savefig(fname_ylim, bbox_inches='tight',dpi=1)
            plt.close()

        except Exception as e: 
            print("Stats :: draw_coverage :: iteration "+str(it)+" reference coverage image draw failed", flush=True, file=log)
            print(e, flush=True, file=log)
        finally:
            if lock.locked():
                with open(lock_log, "a+") as llog:
                    print("thread for it. " +str(it)+ " releases lock in [reference_coverage]", flush=True, file=llog)
                lock.release()


    @staticmethod
    def fastq_to_bigwig(fastq, final_assembly,out_bigwig,out_bam,genome_sizes, minimapobj, log, it=0):
        print(str(datetime.now())+" :: ReferenceStats :: [fastq_to_bigwig] creating bigwig from fastq: "+fastq, flush=True, file=log)
        #timestr = time.strftime("%Y%m%d-%H%M%S")
        #out_paf = "/tmp/dyn_assembly_"+timestr+"_reads_list_to_bigwig.paf"
        out_bedgraph = "/tmp/dyn_assembly_"+str(it)+"_reads_list_to_bigwig.bedgraph"
        bedgraph_sorted = "/tmp/dyn_assembly_"+str(it)+"_reads_list_to_bigwig.sorted.bedgraph"

        def minimap(final_assembly, fastq, out_bam):
            minimapobj.map_reads_to_assembly(fastq,final_assembly,out_bam)
            #return os.system("minimap2 -ax map-ont "+final_assembly+" "+fastq+" 2> /dev/null | samtools view -S -b -o - | samtools sort - -o "+out_bam)
            return 0

        def bedgraph(out_bam, out_bedgraph, bedgraph_sorted):
            res = os.system("bedtools genomecov -ibam "+out_bam+" -bga -split > "+out_bedgraph)
            res += os.system("cat "+out_bedgraph+" | sort -k1,1 -k2,2n > "+bedgraph_sorted)
            return res
        
        def bedgraph_to_bigwig(bedgraph_sorted,genome_sizes,out_bigwig):
            return os.system("bedGraphToBigWig "+bedgraph_sorted+" "+genome_sizes+" "+out_bigwig)

        def clean(out_bedgraph,bedgraph_sorted):
            # do not delete bam for now.. 
            res = os.system("rm "+out_bedgraph)
            res += os.system("rm "+bedgraph_sorted)
            if res > 0:
                raise Exception("clean after reads list to bigwig failed")
        
        print(str(datetime.now())+" :: ReferenceStats :: [fastq_to_bigwig] running minimap",flush=True, file=log)
        if minimap(final_assembly, fastq, out_bam) > 0:
            print("!!! minimap in fastq_to_bigwig retruned non-zero code", file=log)
            #raise Exception("minimap in reads_list_to_bigwig retruned non-zero code")
        print(str(datetime.now())+" :: ReferenceStats :: [fastq_to_bigwig] running bedgraph",flush=True, file=log)
        if bedgraph(out_bam, out_bedgraph, bedgraph_sorted) > 0:
            print("!!! bedgraph in fastq_to_bigwig retruned non-zero code", file=log)
            #raise Exception("bedgraph in reads_list_to_bigwig retruned non-zero code")
        print(str(datetime.now())+" :: ReferenceStats :: [fastq_to_bigwig] running bedgraph to bigwig",flush=True, file=log)
        if bedgraph_to_bigwig(bedgraph_sorted,genome_sizes,out_bigwig) > 0:
            print("!!! bedgraph_to_bigwig in fastq_to_bigwig retruned non-zero code", flush=True, file=log)
            #raise Exception("bed graph to bigwig in reads_list_to_bigwig retruned non-zero code")
        print(str(datetime.now())+" :: ReferenceStats :: [fastq_to_bigwig] cleaning bedgraph files",flush=True, file=log)
        clean(out_bedgraph,bedgraph_sorted)


    @staticmethod
    def map_assembly_to_assembly(final_assembly, partial_assembly, bed, stats_dir,it, log, miniasm_dir, minimapobj):
        print(str(datetime.now())+" :: ReferenceStats :: [map_assembly_to_assembly] assembly to assembly mapping: final - "+final_assembly+"; partial - "+partial_assembly, flush=True, file=log)
        #timestr = time.strftime("%Y%m%d-%H%M%S")
        paf = stats_dir+"with_reference/assembly_to_ref_paf/align_to_ref"+str(it)+".paf"
        dotplot = stats_dir+"with_reference/assembly_to_ref_dotplots/dotplot_"+str(it)+".ps"
        #bam = stats_dir+"with_reference/assembly_to_ref_bam/align_to_ref_"+str(it)+".bam"
        #sorted_bam = stats_dir+"with_reference/assembly_to_ref_bam/align_to_ref_"+str(it)+".sorted.bam"

        def minimap(final_assembly,partial_assembly,paf):
            #return os.system("minimap2 -cx map-ont "+final_assembly+" "+partial_assembly+" 2> /dev/null > "+paf)
            return minimapobj.asm_to_asm_paf(final_assembly,partial_assembly,paf)
        
        #def minimap_bam(final_assembly,partial_assembly,bam,sorted_bam):
        #    minimap = "minimap2 -ax map-ont "+final_assembly+" "+partial_assembly+" 2> /dev/null > "+bam
        #    sort = "samtools sort "+bam+" >"+sorted_bam
        #    index = "samtools index "+sorted_bam
        #    return os.system(minimap) + os.system(sort) + os.system(index)
        
        def paf2bed(paf, bed):
            awk='{print $6 "\t" $8 "\t" $9 "\t" $1 ":l:" $2 ":eoffset:" ($2-$4) "\t" $10 "\t" $5}'
            return os.system("cat "+paf+" | awk -F '\t' '"+awk+"' > "+bed)
        
        def minidot(paf, miniasm_dir, dotplot):
            minidot_ps = miniasm_dir+"minidot -f3 "+paf+" > "+dotplot
            #minidot = miniasm_dir+"minidot -f12 "+paf+" | ps2pdf -dEPSCrop - "+dotplot+".pdf"
            return os.system(minidot_ps)

        def clean(paf):
            return os.system("rm "+paf)

        if minimap(final_assembly, partial_assembly, paf) > 0:
            print("!!! minimap in map_assembly_to_assembly retruned non-zero code", flush=True, file=log)
            #raise Exception("minimap in map_assembly_to_assembly retruned non-zero code")

        #if minimap_bam(final_assembly,partial_assembly,bam) > 0:
        #    print("!!! minimap_bam in map_assembly_to_assembly returned non-zero code", file=log)

        if paf2bed(paf,bed) > 0:
            print("!!! paf2bed in map_assembly_to_assembly retruned non-zero code", flush=True, file=log)
            #raise Exception("paf2bed in map_assembly_to_assembly retruned non-zero code")

        if minidot(paf, miniasm_dir, dotplot) > 0 :
            print("!!! error creating dotplot", flush=True, file=log)

        #clean(paf)

    @staticmethod
    def genome_sizes(fasta, out):
        os.system("faSize -detailed "+fasta+" > "+out)
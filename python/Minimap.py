import os
import time
from datetime import datetime

class Minimap:
    def __init__(self, out_dir, link_dir, threads, K, ftracker):
        self.ftracker = ftracker
        self.in_dir = self.ftracker.get_dir()
        self.out_dir = out_dir  
        self.link_dir = link_dir
        self.processing_files = []
        self.paf_ava = out_dir+"minimap_ava.paf.gz"
        self.paf_assembly = out_dir+"minimap_assembly.paf.gz"
        self.paf_assembly_subsample = out_dir+"minimap_assembly_prev_it_subsample.paf.gz"
        self.assembly = out_dir+"assembly.fasta"
        self.log_ava = out_dir+"minimap_ava.log"
        self.log_assembly = out_dir+"minimap_assembly.log"
        self.threads = threads
        self.K = K

    def asm_to_asm_paf(self,target, query, paf): #stats
        return os.system("minimap2 -x map-ont -t "+str(self.threads)+" -K "+self.K +"  --secondary=no "+target+" "+query+" 2> /dev/null > "+paf)

    def map_reads_to_assembly_paf(self,reads,assembly,paf): #stats
        m = os.system("minimap2 -x map-ont -t "+str(self.threads)+" -K "+self.K +"  --secondary=no "+assembly+" "+reads+" 2> /dev/null > "+paf)
        return m
    
    def map_reads_to_assembly(self,reads,assembly,bam): #fast stats
        m = os.system("minimap2 -a -x map-ont -t "+str(self.threads)+" -K "+self.K +"  --secondary=no "+assembly+" "+reads+" 2> /dev/null | samtools view -S -b -o - | samtools sort - -o "+bam)
        ind = os.system("samtools index "+bam)
        return (m + ind)
        
    def run_minimap_to_assembly(self, remove_subsample = True):  
        print(str(datetime.now())+" :: Minimap :: [run_minimap_to_assembly] mapping reads to assembly", flush=True)
        while len(os.listdir(self.link_dir)) == 0: 
            print(str(datetime.now())+" :: Minimap :: [run_minimap_to_assembly] link dir empty before Minimap to assembly run.. should not happen ?! ", flush=True)
            time.sleep(1)
        
        if remove_subsample:
        #remove subsample.fastq link - we do not want to subsample the reads that were subsampled in previous iterations
            subsample_link = self.link_dir+"subsample.fastq"
            #subsample_fastq = os.readlink(subsample_link) #save the subsample.fastq path
            if os.path.exists(subsample_link): #does not exist in first iteration
                os.unlink(subsample_link) #remove link 
                print(str(datetime.now())+" :: Minimap :: [run_minimap_to_assembly] removed subsample link", flush=True)
                
        print(str(datetime.now())+" :: Minimap :: [run_minimap_to_assembly] files in link dir: ", os.listdir(self.out_dir+"link/"))
        
        if os.system("cat "+self.link_dir+"* > "+self.out_dir+"tmp__map_to_assembly.fastq") > 0:
            raise Exception("cat before minimap (to assembly) returned non-zero code")
        os.system("echo '\n[new to assembly]\n' >> "+self.log_assembly)
        cmd = "minimap2 -x map-ont --secondary=no -t "+str(self.threads)+" -K "+self.K +" "+self.assembly+" "+self.out_dir+"tmp__map_to_assembly.fastq 2>> "+self.log_assembly+" | gzip -1 > "+self.paf_assembly
        print(str(datetime.now())+" :: Minimap :: [run_minimap_to_assembly] mapping reads to assembly: "+cmd)      
        if os.system(cmd) > 0:
            raise Exception("minimap (to assembly) returned non-zero code")
        os.system("rm "+self.out_dir+"tmp__map_to_assembly.fastq")
        
        #os.symlink(subsample_fastq, subsample_link) #no need to reacreate the link

    def run_minimap_old_subsample_to_new_assembly(self):
        print(str(datetime.now())+" :: Minimap :: [run_minimap_old_subsample_to_new_assembly] mapping reads from prev it. to assembly", flush=True)
        os.system("echo '\n[subsample to assembly]\n' >> "+self.log_assembly)
        cmd = "minimap2 -x map-ont --secondary=no -t "+str(self.threads)+" -K "+self.K +" "+self.assembly+" "+self.out_dir+"subsample.fastq 2>> "+self.log_assembly+" | gzip -1 > "+self.paf_assembly_subsample
        print(str(datetime.now())+" :: Minimap :: [run_minimap_old_subsample_to_new_assembly] mapping reads from prev it. to assembly: "+cmd)      
        if os.system(cmd) > 0:
            raise Exception("minimap (to assembly) returned non-zero code")

    
    def run_minimap_all_vs_all(self):
        print(str(datetime.now())+" :: Minimap :: [run_minimap_all_vs_all] all-vs-all read mapping", flush=True)
        while len(os.listdir(self.link_dir)) == 0: 
            print(str(datetime.now())+" :: Minimap :: [run_minimap_all_vs_all] link dir empty before Minimap all-vs-all run.. should not happen ?! ", flush=True)
            time.sleep(1)
        print(str(datetime.now())+" :: Minimap :: [run_minimap_all_vs_all] files in link dir: ", os.listdir(self.out_dir+"link/"))
        if os.system("cat "+self.link_dir+"* > "+self.out_dir+"tmp__map_ava.fastq") > 0:
            raise Exception("cat before minimap a-v-a returned non-zero code")
        os.system("echo '\n' >> "+self.log_ava)
        cmd = "minimap2 -x ava-ont -t "+str(self.threads)+" -K "+self.K +" "+self.out_dir+"tmp__map_ava.fastq "+self.out_dir+"tmp__map_ava.fastq 2>> "+self.log_ava+" | gzip -1 > "+self.paf_ava
        print(str(datetime.now())+" :: Minimap :: [run_minimap_all_vs_all] running minimap (ava): "+cmd)
        if os.system(cmd) > 0:
            raise Exception("minimap a-v-a returned non-zero code")
        os.system("rm "+self.out_dir+"tmp__map_ava.fastq")

        


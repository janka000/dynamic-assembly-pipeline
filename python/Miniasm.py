import os
import time
from datetime import datetime

class Miniasm:
    def __init__(self, out_dir, link_dir, base_dir, min_unitig_reads, min_coverage):
        self.out_dir = out_dir   
        self.link_dir = link_dir
        self.minimap_paf_ava = out_dir+"minimap_ava.paf.gz"
        self.gfa = out_dir+"miniasm.gfa"
        self.fasta = out_dir+"assembly.fasta"
        self.all_reads = out_dir+"reads.fastq"
        self.log = out_dir+"miniasm.log"
        self.miniasm_path = base_dir+"miniasm/miniasm"
        self.min_unitig_reads = min_unitig_reads
        self.min_coverage = min_coverage
    
    def run_miniasm(self, it):
        while len(os.listdir(self.link_dir)) == 0: 
            print(str(datetime.now())+" :: Miniasm :: [run_miniasm] link dir empty before Miniasm run.. should not happen ?! ", flush=True)
            time.sleep(1)
        cmd = self.miniasm_path+" -f "+self.out_dir+"tmp__miniasm.fastq -e "+str(self.min_unitig_reads)+" -c "+ str(self.min_coverage)+" -p paf_reads_ug -a "+self.out_dir+"filtered.paf -j "+self.out_dir+"miniasm_filtered/miniasm_filtered_"+str(it)+".txt "+self.minimap_paf_ava+" > "+self.gfa+" 2>> "+self.log
        print(str(datetime.now())+" :: Miniasm :: [run_miniasm] running miniasm: "+cmd, flush=True)
        if os.system("cat "+self.link_dir+"* > "+self.out_dir+"tmp__miniasm.fastq") > 0:
             raise Exception("cat in miniasm returned non-zero code")
        os.system("echo '\niteration:"+str(it)+"\n' >> "+self.log)
        if os.system(cmd) > 0:
             raise Exception("miniasm returned non-zero code")
        os.system("rm "+self.out_dir+"tmp__miniasm.fastq")
        
    def gfa_to_fasta(self):
        os.system('perl -lane \'print ">$F[1]\n$F[2]" if $F[0] eq "S"\' '+self.gfa+' > '+self.fasta)
        if (os.path.exists(self.fasta) and os.stat(self.fasta).st_size > 0): #exists and not empty 
                os.system("samtools faidx "+ self.fasta)

    def log_assembly(self,it):
        print(str(datetime.now())+" :: Miniasm :: [log_assembly] copying created assembly to "+self.out_dir+"assemblies/assembly_"+str(it)+".fasta")
        os.system('cp '+self.fasta+' '+self.out_dir+"assemblies/assembly_"+str(it)+".fasta")

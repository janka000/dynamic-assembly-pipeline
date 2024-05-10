import os
from datetime import datetime

class Quast:
    def __init__(self, stats_dir, out_dir, reference=None):
        self.stats_dir = stats_dir
        self.reference = reference
        self.quast_dir = self.stats_dir+"quast/"
        self.out_dir = out_dir
        self.assembly_dir = out_dir+"assemblies/"
        if not os.path.exists(self.quast_dir):
            os.makedirs(self.quast_dir)

    def run(self, it, stats_log):
        self.assembly = self.assembly_dir+"assembly_"+str(it)+".fasta"
        out_dir = self.quast_dir+"quast_"+str(it)+"/"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        if self.reference == None: #no reference
            quast_cmd = "quast.py --min-identity 70 --nanopore --use-all-alignments --no-check --ambiguity-mode permissive "+self.assembly+" -o "+out_dir+" > "+out_dir+"quast.log 2> "+out_dir+"quast.err"
            print(str(datetime.now())+" :: Quast :: [run] running quast (withouot reference) :", quast_cmd, file=stats_log)
            os.system(quast_cmd)
            os.system("cp "+self.assembly+" "+out_dir+"assembly_it_"+str(it)+".fasta")

        else:
            quast_cmd = "quast.py "+self.assembly+" -o "+out_dir+" -r "+self.reference+" > "+out_dir+"quast.log 2> "+out_dir+"quast.err"
            print(str(datetime.now())+" :: Quast :: [run] quast reference: "+self.reference, file=stats_log)
            print(str(datetime.now())+" :: Quast :: [run] running quast (with reference):", quast_cmd, file=stats_log)
            os.system(quast_cmd)
            os.system("cp "+self.assembly+" "+out_dir+"assembly_it_"+str(it)+".fasta")
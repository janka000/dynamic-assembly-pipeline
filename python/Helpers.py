import os
import time
from datetime import datetime
import re

class Helpers(object):

    """
    @staticmethod
    def preprocess_fastq_read(read_array):#[id riadok, read, +, quality]
        #@6bfee659-58a2-4870-ad8f-f37b1f37279f runid=37ed0c95ec68e8222ec9d629d5e1715ef411b250 read=101 ch=394 start_time=2018-10-10T09:41:10Z
        print(read_array[0])
        read_id = read_array[0].strip().split()[0][1:]
        print(read_id)
        write_string = read_array[0]+read_array[1]+read_array[2]+read_array[3]
        print(read_array[0]+read_array[2])
        return read_id, write_string
    """

    @staticmethod
    def alphanumsort(array):
        convert = lambda text: int(text) if text.isdigit() else text 
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
        return sorted(array, key = alphanum_key)

    @staticmethod
    def count_reads_in_fastq(fastq):
        num_lines = 0
        count = -1
        with open(fastq, 'r') as fp:
            for count, line in enumerate(fp):
                pass
        num_lines = count+1
        num_reads = num_lines/4
        return num_reads


    

    
from PipelineRunner import PipelineRunner
from PipelineRunnerA import PipelineRunnerA
from PipelineRunnerB import PipelineRunnerB
from PipelineRunnerC import PipelineRunnerC
import threading
import sys
import traceback
import argparse
import os
import resource

def parse_args():
    #parse args
    parser = argparse.ArgumentParser(description='dynamic assembly pipeline')
    #dirs
    parser.add_argument('-i','--in_dir', help='input dir (where are the reads from the run), should end with /', required=True)
    parser.add_argument('-o','--out_dir', help='where will be the files created by the pipeline (should end with /)', required=True)
    parser.add_argument('-bd','--base_dir', help='base directory (repository root) (should end with /)', required=True)
    #pipeline
    parser.add_argument('-s','--subsample_strategy', help='subsample strategy - probabilistic/take_all', default="probabilistic")
    parser.add_argument('-v','--pipeline_version', default='A', choices=['A', 'B', 'C'], help='pipeline version: A (probabilistic subsampling in each step, append to subample) / B ("unbiased subample", rewrite subsample in each iteration) / C (version with map-to-assembly loop)')
    parser.add_argument('-t','--coverage_threshold', type=int, default = 30, help='subsampling covarage threshold')
    parser.add_argument('-l','--mapping_loops', type=int, default = 0, help='num of map-to-assembly lopps before new assebly is created (ignored in A, B pipeline versions)')
    parser.add_argument('-us','--unbiased_simple', action="store_true", default = False, help='unbiased subsampling (ignered in A,C pipeline versions)')
    parser.add_argument('-plr','--prefer_longer_reads', action="store_true", default = False, help='biased subsampling, strongly prefering longer reads (ignered in A,C pipeline versions)')
    parser.add_argument('-e','--max_loops', type=int, default = 10000, help='maximal number of runs  of main (mapping+assmbly) loop - program will end after max_loops iterations')
    parser.add_argument('-w','--max_waiting_count', type=int, default = 20, help='the program will end after max_waiting_count*wait_time of time if no new files present')
    parser.add_argument('-wt','--wait_time', type=int, default = 20, help='interval in which we check whether new files are present, if no new files seen at the start of the pipeline iteration')
    parser.add_argument('-b','--batch_size', type=int, default = 10, help='max number of files to be processed at once')
    #minimap
    parser.add_argument('-mt','--minimap_threads', type=int, default = 8, help='num of minimap threads')
    parser.add_argument('-mk','--minimap_K', default="500M", help='minimap K parameter (number of bases loaded into memory to process in a mini-batch)')
    #miniasm
    parser.add_argument('-mmur','--miniasm_min_unitig_reads', type=int, default=4, help='min reads needed for unitig in miniasm')
    parser.add_argument('-mmc','--miniasm_min_coverage', type=int, default=3, help='minasm coverage threshold for read filtering')

    #debug/track files
    parser.add_argument('-trfa', '--track_finalassembly', default=None, help='final (correct) assembly')
    parser.add_argument('--have_assembly',action="store_true", default=False, help="we have a reference.. ")
    parser.add_argument('--min_disk_space', action="store_true", default=False, help = "remove all non-essential files at the end of pipeline run")


    args = parser.parse_args()

    return args

args = parse_args()
try:
    pipeline = None
    if len(os.listdir(args.out_dir))>5: #stdout log, stderr log, memprofiler,..
        print("output directory is not empty, exiting")
        exit(1)
    if args.pipeline_version == 'A':
        pipeline = PipelineRunnerA(args)
    elif args.pipeline_version == 'B':
        pipeline = PipelineRunnerB(args)
    elif args.pipeline_version == 'C':
        pipeline = PipelineRunnerC(args)
    else:
        print("unrecognised pipeline version, using A")
        pipeline = PipelineRunnerA(args)
    pipeline.run()
except KeyboardInterrupt:
    print("keyboad interrupt catched, exiting")
    sys.exit(0)
except Exception as e:
    print("exception catched")
    print(e)
    traceback.print_exc() 
    sys.exit(1)
finally:
    threads_alive = threading.enumerate()
    main_thread = threading.current_thread()
    print("waiting for threads...")
    for t in threads_alive:
        if t is main_thread:
            continue
        if t.is_alive():
            print("joining ", t.getName())
            t.join()
    print("all threads joined")
    if pipeline is not None:
        pipeline.clean_at_the_end()
    #TODO: handle stuff so that pipeline run can be restarted.. 
    res_usage = resource.getrusage(resource.RUSAGE_SELF)
    with open(args.out_dir+"resource_usage.txt","w+") as ru:
        print(res_usage, file = ru)
    sys.exit(0)

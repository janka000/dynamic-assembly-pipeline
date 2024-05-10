# source files overview
## asg.c
- assembly graph

## common.c
- constants

## asm.c
- generating unitigs

## dotter.c
- minidot (plotting paf files)

## hit.c
- reads filtering

## main.c
- initialization, 
- processing of the command line argumetnts
- main script that handles running of all the steps in order (using methods from the other files)

## paf.c
- manipulation with the .paf file - open, close, parse, read

## sdict.c
- hash table implementation

## sys.c
- timer
- ..

# algorithm steps and  methods used in them
(from the main.c file)

## Step 0: removing contained reads
- hit.c/ma_hit_no_cont(...)
    - find "contained" reads in the PAF file

## Step 1: reading read mappings
- hit.c/ma_hit_read(...)
    - creates (and returns) an array of MAPPINGS from PAF file, stores the READS from them into the hash table and counts the mappings stored
## Step 2: 1-pass (crude) read selection
if (!no_first)				(no_first == provided argument -1  => skip 1-pass read selection)
- stage>=2
    - hit.c/ma_hit_sub(...)  - filtering based on coverage and identity (creates sub array with longest regions with sufficient coverage)
       - dynamic assebmly note: we should **keep the reads filtered based on coverage..** (they may be and probably will be used in assembly when we have more reads) 
    - hit.c/ma_hit_cut(...)  - filtering using the sub array, fixing positions,..    
- stage >=3
    - hit.c/ma_hit_flt(...)  -  convert hits to arcs, check for containments, keep overlaps only, exclude hits associated with deleted regions

## Step 3: 2-pass (fine) read selection
if(!no_second)				(no_second == provided argument -2, skip 2-pass read selection)
- stage>=4
    - ma_hit_sub,                   
    - ma_hit_cut                  
    - … if !no_first:			(no_first = provided argument -1 = skip 1-pass read selection)
        - ma_sub_merge
- stage>=5
    - ma_hit_contained - remove containments.. 


print_subs, print_hits

**-p STR option** ... -> the filtered subset:
- -p STR    Output information and format [ug]. Possible STR values include 
    - bed: post-filtered read regions in the BED format; 
    - paf: mappings between post-filtered reads;   !!! we want this
    - sg: read overlap graph in the GFA format; 
    - ug: unitig graph in the GFA format.


## Step 4: graph cleaning
- sg = ma_sg_gen
### Step 4.1: transitive reduction
- stage>=6
- asg_arc_del_trans(sg,..)


### Step 4.2: initial tip cutting and bubble popping
- stage>=7
- asg_cut_tip(sg,...)
- asg_pop_bubble(sg, ..)



### Step 4.3: cutting short overlaps (\[opt.n_rounds\] rounds in total)
- stage>=9
``` 
for (i = 0; i <= opt.n_rounds; ++i){     		
    float r = opt.min_ovlp_drop_ratio + (opt.max_ovlp_drop_ratio - opt.min_ovlp_drop_ratio) / opt.n_rounds * i;
    if (asg_arc_del_short(sg, r) != 0) {
        asg_cut_tip(sg, opt.max_ext);
        asg_pop_bubble(sg, opt.bub_dist);
    }
}
```

- asg_arc_del_short(sg,r)
- asg_cut_tip(sg, opt.max_ext)
- asg_pop_bubbles(sg, opt.bub_dist)

<em> 
param.:

- opt.n_rounds 		
    - argument -n , rounds of short overlap removal 
- opt.min_ovlp_drop_ratio, opt.max_ovlp_drop_ration  
    - argument -r …max and min overlap drop ratio [%.2g,%.2g]
- opt.max_ext
    -	-e , small unitig threshold [%d]\n"
- opt.bub_dist
    -	-d, max distance for bubble popping
</em>

### Step 4.4: removing short internal sequences and bi-loops
- stage>=10
- asg_cut_internal(sg, 1)
- asg_cut_biloop(sg, opt.max_ext)
- asg_cut_tip(sg, opt.max_ext)
- asg_pop_bubble(sg, opt.bub_dist)

### Step 5: generating unitigs
- if outformat == “ug”
    - ug = ma_ug_gen(sg)
    - if (fn_reads) ma_ug_seq(ug, d, sub, fn_reads )
    - ma_ug_print(ug, d, sub, stdout)
- else
    - ma_sg_print(sg, d, sub, stdout)

### free allocated memory

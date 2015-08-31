#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''Re-cluster sequencing reads according to the given seed OTUs'''

import os
import sys
import glob
import argparse
import logging
import textwrap
import time,traceback
import pipes
import itertools
import subprocess
from multiprocessing import Pool
from collections import defaultdict
import numpy as np
import extract_reads

write_devnull = open(os.devnull,'w')

def view_sam(sam,header=True):
    '''view sam file
    '''
    if header == True:
        p = subprocess.Popen(['samtools','view','-h','-F1804',sam],stdout=subprocess.PIPE,stderr=write_devnull)    
    else:
        p = subprocess.Popen(['samtools','view','-F1804',sam],stdout=subprocess.PIPE,stderr=write_devnull)
    return p.stdout

def merge_sam(samlist,header=False):
    '''merge a list of sam files
    '''
    header_flag = itertools.chain([header],itertools.repeat(False,len(samlist)-1))
    handles = itertools.imap(lambda f,h:view_sam(f,h),samlist,header_flag)
    return itertools.chain.from_iterable(handles)
    

def convert_sam_to_bam(fa,sam,bam,cores):
    '''convert sam to bam
    '''
    cmd = ['samtools','view','-@',str(cores),'-F1804','-bt',fa+'.fai',sam]
    f = open(bam,'w')
    subprocess.call(cmd,stdout=f)
    f.close()
    
    cmd = ['samtools','sort','-@',str(cores),bam,bam+'.sort']
    subprocess.call(cmd)

    cmd = ['mv',bam+'.sort.bam',bam]
    subprocess.call(cmd)

    cmd = ['samtools','index',bam]
    subprocess.call(cmd)

    
def build_a_seed_otu_dict(fasta,otu_id,out_dir,mapper,verbose):
    '''extract a otu seq, save to out_dir, and build the mapper-associated dict'''
    if verbose:
        logging.info('build the %s dict for otu %s'%(mapper,otu_id))
    otu_seq = os.path.join(out_dir,'%s.fasta'%otu_id)
    # extract otu sequence
    with open(otu_seq,'r') as f:
        subprocess.call(['samtools','faidx',fasta,otu_id],stdout=f,stderr=write_devnull)
    subprocess.call(['samtools','faidx',otu_seq],stderr=write_devnull)

    # build the dict
    if mapper == 'bwa':
        subprocess.call(['bwa','index','-p',otu_seq,otu_seq],stdout=write_devnull,stderr=write_devnull)
    elif mapper == 'bowtie2':
        subprocess.call(['bowtie2-build',otu_seq,otu_seq],stdout=write_devnull,stderr=write_devnull)
    elif mapper == 'bowtie':
        subprocess.call(['bowtie-build',otu_seq,otu_seq],stdout=write_devnull,stderr=write_devnull)

def build_seed_otus_dict(fasta,otus,out_dir,out_file,mapper):
    '''extract otu seqs, save to out_dir, and build the mapper-associated dict
    '''
    out = os.path.join(out_dir,out_file)
    os.system('samtools faidx %s %s > %s 2>/dev/null' % (fasta,otus,out))
    #with open(out,'w') as f:
    #    subprocess.call(['samtools','faidx',fasta,otus],stdout=f,stderr=write_devnull)
    subprocess.call(['samtools','faidx',out],stdout=write_devnull,stderr=write_devnull)

    # build the dict
    if mapper == 'bwa':
        subprocess.call(['bwa','index','-p',out,out],stdout=write_devnull,stderr=write_devnull)
    elif mapper == 'bowtie2':
        subprocess.call(['bowtie2-build',out,out],stdout=write_devnull,stderr=write_devnull)
    elif mapper == 'bowtie':
        subprocess.call(['bowtie-build',out,out],stdout=write_devnull,stderr=write_devnull)
    elif mapper == 'smalt':
        subprocess.call(['smalt','index',out,out],stdout=write_devnull,stderr=write_devnull)


def map_data_to_otu(otu,data,mapper,map_args,out_prefix,verbose):
    '''map sequencing reads to a seed otu reference'''
    if verbose:
        logging.info('map sequencing reads to otu %s' % os.path.splitext(os.path.basename(otu))[0])
    # mapping
    for sample,runs in data.iteritems():
        if mapper == 'bowtie2':
            d = []
            for run in runs:
                if len(run) == 2:
                    d.extend(['-U',run[1]])
                elif len(run) == 3:
                    d.extend(['-1',run[1],'-2',run[2]])
            p = pipes.Template()
            p.append('bowtie2 %s -x %s %s -S %s.%s.sam >&/dev/null'%(map_args,\
                      otu,' '.join(d),\
                      out_prefix,sample),'--')
            f = p.open(None,'w')
            f.close()
        elif mapper == 'bowtie':
            d = ['bowtie',map_args,otu,'-1',run[1],'-2',run[2],'%s.%s.sam'%(out_prefix,sample)]
            subprocess.call(d,stdout=write_devnull,stderr=write_devnull)
        elif mapper == 'bwa':
            for i,run in enumerate(runs):
                p = pipes.Template()
                p.append('bwa %s %s %s 1>%s.%s.%d.sam 2>/dev/null'%(map_args,\
                         otu,' '.join(run),\
                         out_prefix,sample,i),'--')
                f = p.open(None,'w')
                f.close()
            with open('%s.%s.sam'%(out_prefix,sample),'w') as f:
                for row in merge_sam(['%s.%s.%d.sam' % (out_prefix,sample,i) for i in xrange(len(runs))]):
                    f.write(row)

def map_data_to_seed_otus(seed_otus,data,mapper,map_args,out_dir,out_prefix,cores):
    '''map sequencing data to seed otus'''
    samlist = []
    # mapping
    for sample,runs in data.iteritems():
        for run in runs:
            command = [mapper,map_args]
            if mapper == 'bwa':
                command.extend(['-t',str(cores)])
            elif mapper == 'bowtie2':
                command.extend(['-p',str(cores)])
            elif mapper == 'bowtie':
                command.extend(['-p',str(cores)])
            elif mapper == 'smalt':
                command.extend(['-n',str(cores)])
            if mapper == 'bwa':
                command.extend([seed_otus])
            elif mapper == 'bowtie2':
                command.extend(['-x',seed_otus])
            elif mapper == 'bowtie':
                command.extend([seed_otus])
            elif mapper == 'smalt':
                command.extend([seed_otus])
            if len(run) == 2:
                if mapper == 'bwa':
                    command.extend([run(1)])
                elif mapper == 'bowtie2':
                    command.extend(['-U',run[1]])
                elif mapper == 'smalt':
                    command.extend([run(1)])
            elif len(run) == 3:
                if mapper == 'bwa':
                    command.extend([run[1],run[2]])
                elif mapper == 'bowtie2':
                    command.extend(['-1',run[1],'-2',run[2]])
                elif mapper == 'bowtie':
                    command.extend(['-1',run[1],'-2',run[2]])
                elif mapper == 'smalt':
                    command.extend([run[1],run[2]])
            rn = run[0]
            if mapper == 'bwa':
                command.extend(['>',os.path.join(out_dir,'%s.%s.%s.sam'%(out_prefix,sample,rn))])
            elif mapper == 'bowtie2':
                command.extend(['-S',os.path.join(out_dir,'%s.%s.%s.sam'%(out_prefix,sample,rn))])
            elif mapper == 'bowtie':
                command.extend([os.path.join(out_dir,'%s.%s.%s.sam'%(out_prefix,sample,rn))])
            elif mapper == 'smalt':
                command.extend(['>',os.path.join(out_dir,'%s.%s.%s.sam'%(out_prefix,sample,rn))])
            command.extend(['>&/dev/null'])
            # execute
            logging.debug('map run %s of sample %s to seed otus with %d cpus' % (rn,sample,cores))
            t0 = time.time()
            os.system(' '.join(command))
            logging.debug('elapsed time on mapping run %s of sample %s is %f minutes' % (rn,sample,(time.time()-t0)/60.))
            samlist.append(os.path.join(out_dir,'%s.%s.%s.sam'%(out_prefix,sample,rn)))
    return samlist
    
    
def file_size(f):
    if not os.path.exists(f):
        return 0
    statinfo = os.stat(f)
    return statinfo.st_size
    
def extract_marker_gene_reads(bams,out_prefix,cores):
    N = min(len(bams),cores)
    extract_reads.extract_reads(bamlist=bams,cores=N,prefix=out_prefix)
    # merge single reads
    cmd = ['cat']
    for i in xrange(N):
        hs = '%s.%d.single.fastq' % (out_prefix,i)
        if file_size(hs)>0:
            cmd.append(hs)
    if len(cmd)>1:
        f = open(out_prefix+'.single.fastq','w')
        subprocess.call(cmd,stdout=f)
        f.close()
    # merge pe1 reads
    cmd = ['cat']
    for i in xrange(N):
        h1 = '%s.%d.pe1.fastq' % (out_prefix,i)
        if file_size(h1)>0:
            cmd.append(h1)
    if len(cmd)>1:
        f = open(out_prefix+'.pe1.fastq','w')
        subprocess.call(cmd,stdout=f)
        f.close()
    # merge pe2 reads
    cmd = ['cat']
    for i in xrange(N):
        h2 = '%s.%d.pe2.fastq' % (out_prefix,i)
        if file_size(h2)>0:
            cmd.append(h2)
    if len(cmd)>1:
        f = open(out_prefix+'.pe2.fastq','w')
        subprocess.call(cmd,stdout=f)
        f.close()
    # remove temporary file
    cmd = ['rm']
    for i in xrange(N):
        hs = '%s.%d.single.fastq' % (out_prefix,i)
        h1 = '%s.%d.pe1.fastq' % (out_prefix,i)
        h2 = '%s.%d.pe2.fastq' % (out_prefix,i)
        if os.path.exists(hs):
            cmd.append(hs)
        if os.path.exists(h1):
            cmd.append(h1)
        if os.path.exists(h2):
            cmd.append(h2)
    subprocess.call(cmd)
    

def recluster_data(fasta_file,otu_file,bam_file,mapper,map_args,cores,start_from):
    '''main code interface'''
    # make a directory storing seed OTU sequences
    otu_dir = os.path.join(os.getcwd(),"0_otu_dir")
    if not os.path.exists(otu_dir):
        os.mkdir(otu_dir)
    
    # load seed otus
    seed_otus = list()
    with open(otu_file,'r') as f:
        for line in f:
            items = line.rstrip().split()
            seed_otus.append(items[0])

    # build the mapping index dict
    otu_file = 'seed_otus.fasta'
    if start_from == 0:
        logging.info('build the %s index file'%mapper)
        t0 = time.time()
        build_seed_otus_dict(fasta_file,' '.join(seed_otus),otu_dir,otu_file,mapper)
        logging.info('elapsed time on building the %s index file is %f minutes'%(mapper,(time.time()-t0)/60.))

    # load bam files
    bam_files = list()
    with open(bam_file) as f:
        for line in f:
            bam_files.append(line.strip())
#
    # load data file
#    sample_data = defaultdict()
#    sample_record = np.loadtxt(opts.data,dtype=[('run','S1024'),('sample','S1024'),('mode','S1024'),\
#                                                ('single','S1024'),('pe1','S1024'),('pe2','S1024')],\
#                                                comments=None)
#    for d in sample_record:
#        if d['run'].startswith('#'):
#            continue
#        if d['mode'] == 'single':
#            sample_data.setdefault(d['sample'],[]).extend([(d['run'],d['single'])])
#        else:
#            sample_data.setdefault(d['sample'],[]).extend([(d['run'],d['pe1'],d['pe2'])])
#
    # make a directory storing gene-associated reads
    read_dir = os.path.join(os.getcwd(),"1_read_dir")
    if not os.path.exists(read_dir):
        os.mkdir(read_dir)

    # extract marker gene reads
    sample_data = defaultdict()
    if start_from <= 1:
        extract_marker_gene_reads(bam_files,os.path.join(read_dir,'marker_gene_reads'),cores)
    single_reads = os.path.join(read_dir,'marker_gene_reads.single.fastq')
    pe1_reads = os.path.join(read_dir,'marker_gene_reads.pe1.fastq')
    pe2_reads = os.path.join(read_dir,'marker_gene_reads.pe2.fastq')
    if os.path.exists(single_reads):
        sample_data.setdefault('all_samples',[]).extend([('single',single_reads)])
    if os.path.exists(pe1_reads):
        sample_data.setdefault('all_samples',[]).extend([('paired',pe1_reads,pe2_reads)])

    # make a directory storing mapping results
    map_dir = os.path.join(os.getcwd(),"2_map_dir")
    if not os.path.exists(map_dir):
        os.mkdir(map_dir)

    # map sequencing reads to seed otus
    map_prefix = 'to_seed_otus'
    if start_from <= 2:
        logging.info('use %s to map reads onto %d seed otus with %d cores'%(mapper,len(seed_otus),cores))
        t0 = time.time()
        samlist = map_data_to_seed_otus(os.path.join(otu_dir,otu_file),sample_data,mapper,map_args,\
                                        map_dir,map_prefix,cores)
        logging.info('elapsed time on mapping reads onto seed otus is %f minutes'%((time.time()-t0)/60.))
    else:
        samlist = []
        for sample,runs in sample_data.iteritems():
            for run in runs:
                samlist.append('%s.%s.%s.sam' % (map_prefix,sample,run[0]))

    # merge sample mappings
    map_sam_file = map_prefix+".all.sam"
    with open(map_sam_file,'w') as f:
        for row in merge_sam(samlist,header=True):
            f.write(row)

    # covert sam to bam
    map_bam_file = map_prefix+".all.bam"
    convert_sam_to_bam(os.path.join(otu_dir,otu_file),map_sam_file,map_bam_file,cores)


if __name__ == "__main__":
    try:
        # set up command-line parser
        parser = argparse.ArgumentParser(description=globals()['__doc__'],epilog=textwrap.dedent('''\
                                         Examples\n\
                                         --------\n\
                                         1. recluster_data_to_seed_otus \n\
                                         \n'''),formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('fasta',help='gene reference genome file',metavar="GENE_FASTA")
        parser.add_argument('otu',help='a list of seed OTUs, one OTU one line',metavar="SEED_OTUS")
        parser.add_argument('bams',help='a list of bam files, one file one line',metavar="BAM_FILES")
        #parser.add_argument('data',help='a formatted file storing sample datas',metavar="SAMPLE_RUNS")
        parser.add_argument('-c','--core',help='number of computing cores (default: 1)',default=1,dest='cores',type=int)
        parser.add_argument('-m','--mapper',help='mapping program, support bowtie,bowtie2,bwa,smalt (default: bowtie2)',\
                            dest='mapper',default='bowtie2',type=str)
        parser.add_argument('-A','--map_args',help='mapping program arguments (default: --local)',\
                            dest='map_args',default='--local',type=str)
        parser.add_argument('-S','--start_from',help='''start program from stage,\
0: build mapping index,\
1: extract gene-associated reads,\
2: map reads to seed OTUs,\
(default: 0)''',\
                            type=int,default=0,dest='start_from')
        parser.add_argument('-v',help='verbose level, 0:quiet, 1:info, 2:debug (default: 0)',\
                            type=int,default=0,dest='verbose')
        opts = parser.parse_args()

        # set up logging format
        if opts.verbose == 1:
            logging.basicConfig(format="[%(asctime)s] %(levelname)s : %(message)s", level=logging.INFO)
        elif opts.verbose == 2:
            logging.basicConfig(format="[%(asctime)s] %(levelname)s : %(message)s", level=logging.DEBUG)

        if opts.start_from == 0:
            logging.info('start from building the mapping index')
        elif opts.start_from == 1:
            logging.info('start from mapping sequencing reads')

        # TODO main code
        t_start = time.time()
        recluster_data(opts.fasta,opts.otu,opts.bams,\
                       opts.mapper,opts.map_args,\
                       opts.cores,opts.start_from)
        # complete
        logging.info('total elapsed time is %f minutes' % ((time.time()-t_start)/60.))
        sys.exit(0)
    except KeyboardInterrupt,e:
        raise e
    except SystemExit,e:
        raise e
    except Exception,e:
        logging.exception('ERROR, UNEXCEPTED EXCEPTION')
        logging.exception(e)
        traceback.print_exc()
        os._exit(1)

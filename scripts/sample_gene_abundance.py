#!/usr/bin/env python
'''
Mapping-based sample-wise abundance estimation.
'''

import os
import sys
import logging
import argparse
from subprocess import call,Popen,PIPE
from multiprocessing import Pool
import time
import textwrap,traceback
from collections import defaultdict
import numpy as np
from scipy.optimize import minimize
import csv
import extract_reads

write_devnull = open(os.devnull,'w')

def parse_sample_metadata(sample_file):
    sample_data = defaultdict(list)
    with open(sample_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            run,smpl,mode,single,pe1,pe2 = line.strip().split()
            if mode=='single':
                sample_data.setdefault(smpl,[]).extend([(run,single)])
            else:
                sample_data.setdefault(smpl,[]).extend([(run,pe1,pe2)])
    return sample_data

def extract_sample_gene_read(params):
    smpl,bam = params
    extract_reads.extract_reads(bamlist=[bam],prefix=smpl,cores=1)
    return None

def file_size(f):
    if not os.path.exists(f):
        return 0
    statinfo = os.stat(f)
    return statinfo.st_size

def parse_sample_metadata2(sample_file):
    global opts
    sample_list = list()
    sample_bam = defaultdict(list)
    sample_data = defaultdict(list)
    params = list()
    with open(sample_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            smpl,bam = line.strip().split()
            sample_bam.setdefault(smpl,[]).extend([bam])
            params.extend([(smpl,bam)])
            sample_list.extend([smpl])
    # extract reads
    if opts.verbose:
        logging.info('extract 16S rRNA gene associated reads')
    map(extract_sample_gene_read,params)
    # collect reads
    for smpl in sample_list:
        single = smpl+'.0.single.fastq'
        pe10 = smpl+'.0.pe1.fastq'
        pe20 = smpl+'.0.pe2.fastq'
        pe1 = smpl+'.pe1.fastq'
        pe2 = smpl+'.pe2.fastq'
        pes = smpl+'.pes.fastq'
        if file_size(single)>0:
            sample_data.setdefault(smpl,[]).extend([(smpl,single)])
        if file_size(pe10)>0:
            cmd = ['sickle','pe','-t','sanger','-f',pe10,'-r',pe20,'-o',pe1,'-p',pe2,'-s',pes,'-q','10']
            call(cmd)
            sample_data.setdefault(smpl,[]).extend([(smpl,pe1,pe2)])
    return sample_data

def parse_gene_count(gene_file):
    gene_count = defaultdict(int)
    res = Popen(['grep','^>',gene_file],stdout=PIPE).communicate()[0].split('\n')
    for gene_des in res:
        if gene_des.startswith('>'):
            gene,count = gene_des.replace('>','').split()
            gene_count[gene] = int(count)
    return gene_count
    
def make_bowtie2_index(gene_file):
    cmd = ['bowtie2-build',gene_file,os.path.join(os.getcwd(),os.path.basename(gene_file))]
    call(cmd,stderr=write_devnull,stdout=write_devnull)
    return os.path.join(os.getcwd(),os.path.basename(gene_file))

def sam_to_bam(gene_file,prefix):
    cmd = ['samtools','view','-F4','-hbt',gene_file+'.fai',prefix+'.sam']
    with open(prefix+'.bam','w') as f:
        call(cmd,stdout=f,stderr=write_devnull)
    
    cmd = ['samtools','sort',prefix+'.bam',prefix+'.sort']
    call(cmd,stderr=write_devnull)
    
    cmd = ['mv',prefix+'.sort.bam',prefix+'.bam']
    call(cmd)

    cmd = ['samtools','index',prefix+'.bam']
    call(cmd,stderr=write_devnull)

    cmd = ['rm',prefix+'.sam']
    call(cmd)


def sgc(params):
    bam,gene = params
    p1 = Popen(['samtools','view','-F4',bam,gene],stdout=PIPE,stderr=write_devnull)
    p2 = Popen(['wc','-l'],stdin=p1.stdout,stdout=PIPE,stderr=write_devnull)
    p1.stdout.close()
    res = p2.communicate()[0].strip()
    return (gene,int(res))

def sample_gene_count(sample_data,gene_file,gene_bt2idx):
    global opts
    if not os.path.exists(gene_file+'.fai'):
        call(['samtools','faidx',gene_file])

    # load gene names
    gene_set = []
    with open(gene_file+'.fai') as f:
        for line in f:
            gene_set.append(line.split()[0])

    # scan through samples
    sample_count = defaultdict()
    for i,smpl in enumerate(sample_data):
        t0 = time.time()
        if opts.verbose:
            logging.info('processing '+smpl+' '+str(i+1)+'/'+str(len(sample_data)))
        prefix = '%s_to_%s'%(smpl,os.path.basename(gene_file))
        sam = '%s.sam'%prefix
        cmd = ['bowtie2','--local','-p',str(opts.cores),'-x',gene_bt2idx]
        for run in sample_data[smpl]:
            if len(run)==2:
                cmd.extend(['-U',run[1]])
            else:
                cmd.extend(['-1',run[1],'-2',run[2]])
        cmd.extend(['-S',sam])
        # mapping
        if opts.verbose:
            logging.info('mapping reads using %d cores'%opts.cores)
        call(cmd,stderr=write_devnull)
        # sam to bam
        sam_to_bam(gene_file,prefix)
        # scan through gene
        if opts.verbose:
            logging.info('scan through genes')
        smpl_gc = defaultdict(int)
        params = []
        for gene in gene_set:
            params.extend([(prefix+'.bam',gene)])
        pool = Pool(min([opts.cores,len(params)]))
        R = []
        r = pool.map_async(sgc,params,callback=R.extend)
        r.wait()
        pool.close()
        pool.join()
        for gene,gc in R:
            smpl_gc[gene] = gc
        # save
        sample_count[smpl] = smpl_gc
        if opts.verbose:
            logging.info('elapsed time is %f minutes' % ((time.time()-t0)/60.))
    return sample_count
        
def sample_estimation(sample_file,gene_file):
    global opts
    CWD = os.getcwd()
    # make temporary directory
    TmpDir = os.path.join(os.getcwd(),'4_sample_gene_abundance')
    if os.path.exists(TmpDir):
        call(['rm','-rf',TmpDir])
    os.mkdir(TmpDir)
    os.chdir(TmpDir)

    # TODO load sample metadata
    if opts.verbose:
        logging.info('parse sample metadata')
    sample_data = parse_sample_metadata2(sample_file)

    # TODO make bowtie2 index
    if opts.verbose:
        logging.info('make bowtie2 index of '+gene_file)
    gene_bt2idx = make_bowtie2_index(gene_file)
 
    # TODO sample-wise gene read count
    if opts.verbose:
        logging.info('compute sample-wise gene read count')
    sample_count = sample_gene_count(sample_data,gene_file,gene_bt2idx)

    # gene set
    gene_set = []
    with open(gene_file+'.fai') as f:
        for line in f:
            gene_set.append(line.split()[0])
    
    # TODO write to disk file
    with open(os.path.join(CWD,'sample_gene_count.txt'),'w') as f:
        fieldnames = ['sample'] + gene_set
        writer = csv.DictWriter(f, fieldnames=fieldnames,delimiter=' ')
        writer.writeheader()
        for smpl,gc in sample_count.iteritems():
            row = {'sample':smpl}
            for g,c in gc.iteritems():
                row[g] = c
            writer.writerow(row)

    # back to parent directory
    os.chdir(CWD)
    # remove temporary directory
    call(['rm','-rf',TmpDir])

if __name__=="__main__":
    try:
        t_start = time.time()
        parser = argparse.ArgumentParser(description=globals()['__doc__'],epilog=textwrap.dedent('''\
                                         Examples\n\
                                         --------\n\
                                         1. sample_gene_abundance sample_info 16S_contig -v\n\n \
                                         '''),\
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('sample_info',help='a metadata record of samples',metavar='METADATA')
        parser.add_argument('gene_contig',help='the assembled 16S gene contigs',metavar='CONTIG')
        parser.add_argument('-c','--cores',help='computing cores [1]',default=1,type=int,metavar='INT')
        parser.add_argument('-v',dest='verbose',action='store_true',default=False,help='verbose output')
        # parse options
        opts = parser.parse_args()

        # set logging format
        logging.basicConfig(format='[%(asctime)s] %(levelname)s : %(message)s', level=logging.INFO)

        # main routine
        sample_estimation(os.path.abspath(opts.sample_info),\
                          os.path.abspath(opts.gene_contig))

        # complete
        if opts.verbose:
            logging.info('elapsed time is %.5f minutes' % ((time.time()-t_start)/60.))

        sys.exit(0)
    except KeyboardInterrupt,e:
        raise e
    except SystemExit,e:
        raise e
    except Exception,e:
        logging.exception('ERROR, UNEXPECTED EXCEPTION')
        logging.exception(e)
        traceback.print_exc()
        os._exit(1)

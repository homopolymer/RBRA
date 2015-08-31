#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''Extract out reads from SAM/BAM file(s)'''

import os
import sys
import time
import logging
import argparse
import textwrap
import traceback
import subprocess
import itertools
import numpy as np
from mpi4py import MPI
from multiprocessing import Pool


write_devnull = open(os.devnull,'w')

rc = {'A':'T','a':'t',\
      'C':'G','c':'g',\
      'G':'C','g':'c',\
      'T':'A','t':'a',\
      'N':'N','n':'n'}

def read_is_paired(flag):
    return ((flag&int('0x1',0))>0)

def read_is_unmapped(flag):
    return ((flag&int('0x4',0))>0)

def paired_is_unmapped(flag):
    return ((flag&int('0x8',0))>0)

def read_is_reverse_complement(flag):
    return ((flag&int('0x10',0))>0)

def read_is_paired_first(flag):
    return ((flag&int('0x40',0))>0)

def reverse_complement(seq):
    rseq = [rc[b] for b in seq[::-1]]
    return ''.join(rseq)
    
def view_bam(bam, roi=None):
    '''view a bam file with/without roi'''
    logging.debug('view %s' % bam)
    command = ['samtools','view','-F4']
    
    command.append(bam)
    
    if roi is not None and len(roi)>0:
        if isinstance(roi,(list,tuple)):
            command.extend(roi)
        else:
            command.append(roi)

    return subprocess.Popen(command,stdout=subprocess.PIPE,stderr=write_devnull).stdout

def write_read_to_file(handle,name,seq,qual):
    handle.write('@%s\n'%name)
    handle.write('%s\n'%seq)
    handle.write('+\n')
    handle.write('%s\n'%qual)


def extract_reads_of_a_bamlist(bamlist,roi=None,hs=None,h1=None,h2=None):
    '''extract sequence reads from a list of bam files'''
    logging.debug('extract reads from %d bam files' % len(bamlist))
    handles = itertools.imap(lambda b,r:view_bam(b,r),\
                             bamlist,itertools.repeat(roi))
    # output reads
    for read in itertools.chain.from_iterable(handles):
        item = read.rstrip().split()       
        paired  = read_is_paired(int(item[1]))
        mapped  = not read_is_unmapped(int(item[1]))
        paired_mapped = not paired_is_unmapped(int(item[1]))
        if read_is_reverse_complement(int(item[1])):
            seq = reverse_complement(item[9])
            qual = ''.join([x for x in item[10][::-1]])
        else:
            seq = item[9]
            qual = item[10]
        if paired:
            if mapped | paired_mapped:
                if item[6]=="=" or item[6]==item[3] or item[6]=="*" or item[3]=="*":
                    if read_is_paired_first(int(item[1])):
                        write_read_to_file(h1,'%s/1'%item[0],seq,qual)
                    else:
                        write_read_to_file(h2,'%s/2'%item[0],seq,qual)
        elif mapped:
            write_read_to_file(hs,item[0],seq,qual)
            
def sub_extract_reads(bamlist=list(),roilist=None,core=None,prefix=None):
    hs = open('%s.%d.single.fastq'%(prefix,core),'w')
    h1 = open('%s.%d.pe1.fastq'%(prefix,core),'w')
    h2 = open('%s.%d.pe2.fastq'%(prefix,core),'w')
    extract_reads_of_a_bamlist(bamlist,roilist,hs,h1,h2)
    hs.close()
    h1.close()
    h2.close()
    return 0
    
def sub_extract_reads_wrapper(args):
    return sub_extract_reads(*args)

def extract_reads(bamlist=list(),roilist=None,cores=1,prefix=None):
    logging.info('dump reads from %d bam files'%len(bamlist))
    t0 = time.time()
    N = min(len(bamlist),cores)
    R = []
    pool = Pool(N)
    r = pool.map_async(sub_extract_reads_wrapper,\
                       [([bamlist[i]],roilist,i,prefix) for i in xrange(len(bamlist))],callback=R.extend)
    r.wait()
    pool.close()
    pool.join()
    logging.info('elapsed time on dumping reads from %d bam files is %f minutes' % (len(bamlist),(time.time()-t0)/60.))
    

if __name__ == "__main__":
    try:
        # set command-line parser
        parser = argparse.ArgumentParser(description=globals()['__doc__'],epilog=textwrap.dedent('''\
                                         Examples\n\
                                         --------\n\
                                         1. extract_reads\n\
                                         \n'''),formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-b','--bam',help='SAM/BAM file(s)',nargs='*',type=str,default=None)
        parser.add_argument('-f','--fofn',help='file of listing SAM/BAM files, one-file-one-line',type=str,default=None)
        parser.add_argument('-r','--roi',help='region(s) of interesting',nargs='*',type=str,default=None)
        parser.add_argument('-R',help='file listing regions of interesting',default=None,type=str,dest="ROI")
        parser.add_argument('-p',help='number of threading (default: 1)',type=int,default=1,dest='num_thread')
        parser.add_argument('prefix',help='prefix name of output file(s)',type=str)
        parser.add_argument('-v',help='verbose level, 0:quiet, 1:info, 2:debug (default: 0)',\
                            default=1, type=int, dest='verbose')
        opts = parser.parse_args()

        # set logging
        if opts.verbose == 1:
            logging.basicConfig(format="[%(asctime)s] : %(levelname)s : %(message)s", level=logging.INFO)
        elif opts.verbose == 2:
            logging.basicConfig(format="[%(asctime)s] : %(levelname)s : %(message)s", level=logging.DEBUG)

        # set timer
        t0 = time.time()

        # main code
        if True:
            bamlist = list()
            if opts.bam is not None:
                bamlist = opts.bam
            if opts.fofn is not None:
                with open(opts.fofn,'r') as f:
                    for line in f:
                        bamlist.append(line.rstrip())
            roilist = list()
            if opts.roi is not None:
                roilist = opts.roi
            if opts.ROI is not None:
                with open(opts.ROI,'r') as f:
                    for line in f:
                        roilist.append(line.rstrip().split()[0])
            
        extract_reads(bamlist, roilist, opts.num_thread, opts.prefix)

        # complete
        logging.info('total elapsed time is %f minutes' % ((time.time()-t0)/60.))
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

    

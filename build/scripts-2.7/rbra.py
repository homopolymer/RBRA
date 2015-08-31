#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
RBRA: Reference-based Ribosome Assembly
---------------------------------------
It is to assemble full-length 16S rRNA genes of strains
using metagenomics sequencing data.

'''

import os
import sys
import logging
import argparse
import subprocess
from multiprocessing import Pool
import time,random
import textwrap,traceback

class SeqData:
    def __init__(self,_sample="",_bam=""):
        self.SampleData = _sample
        self.BamFile = _bam

    def __repr__(self):
        return 'SampleData=%s, BamFile=%s'%(self.SampleData,self.BamFile)

    def __str__(self):
        return 'SampleData=%s, BamFile=%s'%(self.SampleData,self.BamFile)


class GeneData:
    def __init__(self,_gene="",_index="",_tree="",_tax=""):
        self.GeneSeq = _gene
        self.GeneIndex = _index
        self.GeneTree = _tree
        self.GeneTax = _tax

    def __repr__(self):
        return 'GeneSeq=%s, GeneIndex=%s, GeneTree=%s'%(self.GeneSeq,self.GeneIndex,self.GeneTree)
    
    def __str__(self):
        return 'GeneSeq=%s, GeneIndex=%s, GeneTree=%s'%(self.GeneSeq,self.GeneIndex,self.GeneTree)


def parse_data_info(data_info):
    SampleData = ""
    BamFile = ""
    GeneSeq = ""
    GeneIndex = ""
    GeneTree = ""
    with open(data_info) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('SampleData'):
                SampleData = line.split('=')[1]
                SampleData = SampleData.strip()
            elif line.startswith('BamFiles'):
                BamFile = line.split('=')[1]
                BamFile = BamFile.strip()
            elif line.startswith('GeneSeq'):
                GeneSeq = line.split('=')[1]
                GeneSeq = GeneSeq.strip()
            elif line.startswith('GeneIndex'):
                GeneIndex = line.split('=')[1]
                GeneIndex = GeneIndex.strip()
            elif line.startswith('GeneTree'):
                GeneTree = line.split('=')[1]
                GeneTree = GeneTree.strip()
            elif line.startswith('GeneTax'):
                GeneTax = line.split('=')[1]
                GeneTax = GeneTax.strip()
    return (SeqData(SampleData,BamFile), GeneData(GeneSeq,GeneIndex,GeneTree,GeneTax))


def compute_depth_of_seq_data(GeneData,SeqData,OutFile):
    global opts,ExecDir

    if opts.verbose:
        logging.info('profile sequencing depths')

    N = 0;
    with open(SeqData.BamFile) as f:
        for line in f:
            if len(line.strip())>0:
                N += 1
    N = min([N,opts.cores])
    cmd = ['mpirun','-np',str(N),os.path.join(ExecDir,'coverage_all_samples.py')]

    if opts.verbose:
        cmd.append('-v')

    cmd.extend([SeqData.BamFile,GeneData.GeneIndex])

    f = open(OutFile,'w')
    subprocess.call(cmd,stdout=f)

def compute_gene_abundance(GeneData,DepthFile,OutFile):
    global opts,ExecDir

    if opts.verbose:
        logging.info('compute gene abundance')

    cmd = [os.path.join(ExecDir,'gene_abundance.py')]

    if opts.verbose:
        cmd.append('-v')
    
    cmd.extend([DepthFile,GeneData.GeneIndex])
  
    f = open(OutFile,'w')
    subprocess.call(cmd,stdout=f)


def find_seed_genes(GeneData,GeneAbundance,OutFile):
    global opts,ExecDir

    if opts.verbose:
        logging.info('find seed genes')

    cmd = [os.path.join(ExecDir,'find_seed_otus.py')]

    if opts.verbose:
        cmd.append('-v')

    if len(GeneData.GeneTax)>0:
        cmd.extend(['-T',GeneData.GeneTax])

    cmd.extend(['-s',str(opts.gene_sim),'-d',str(opts.gene_abun)])
    cmd.extend([GeneData.GeneTree,GeneAbundance])

    f = open(OutFile,'w')
    subprocess.call(cmd,stdout=f)

def recluster_data(GeneData,SeqData,SeedGene):
    global opts,ExecDir
    
    if opts.verbose:
        logging.info('map gene reads to seed genes')

    cmd = [os.path.join(ExecDir,'recluster_data_to_seed_otus.py'),'-c',str(opts.cores)]

    if opts.verbose:
        cmd.extend(['-v',str(1)])

    cmd.extend([GeneData.GeneSeq,SeedGene,SeqData.BamFile])
    subprocess.call(cmd)

def strain_call_core(params):
    global opts
    cmd,out = params
    if opts.verbose:
        logging.info(' '.join(cmd))
    with open(out,'w') as f:
        subprocess.call(cmd,stdout=f)
    return None


def strain_call(TmpDir):
    global opts,ExecDir
    seed_gene = []
    with open(os.path.join(TmpDir,'0_otu_dir/seed_otus.fasta.fai')) as f:
        for line in f:
            field = line.split()
            seed_gene.extend(['%s:1-%s'%(field[0],field[1])])
    ScDir = os.path.join(TmpDir,'3_straincall_results')
    if not os.path.exists(ScDir):
        os.mkdir(ScDir)
    params = []
    for roi in seed_gene:
        cmd = [os.path.join(ExecDir,'StrainCall'),\
               '-r',str(roi),'-q',str(opts.map_qual),\
               '-D',str(opts.max_depth),'-I',str(opts.max_ins),\
               '-l',str(opts.read_len),'-t',str(opts.tau),\
               '-d',str(opts.diff_rate),'-w',str(5000),\
               os.path.join(TmpDir,'0_otu_dir/seed_otus.fasta'),\
               os.path.join(TmpDir,'to_seed_otus.all.bam')]
        params.extend([(cmd,os.path.join(ScDir,'%s.fa'%roi))])
 
    R = []
    pool = Pool(opts.cores)
    r = pool.map_async(strain_call_core,params,callback=R.extend)   
    r.wait()
    pool.close()
    pool.join()
    
    cmd = ['cat']
    for roi in seed_gene:
        cmd.extend([os.path.join(ScDir,'%s.fa'%roi)])
    with open(os.path.join(TmpDir,'%s.fa'%opts.prefix),'w') as f:
        subprocess.call(cmd,stdout=f)

#def gene_taxonomy_annotation(GeneData,TmpDir):
#    global opts
    
def rbra_pipe(GeneData,SeqData):
    global opts,ExecDir

    WorkDir = os.getcwd()
    TmpDir = os.path.join(WorkDir,"RBRA_work_dir_%s"%(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')))
    if os.path.exists(TmpDir):
        subprocess.call(['rm','-rf',TmpDir])
    os.mkdir(TmpDir);
    os.chdir(TmpDir);
    
    # TODO  1. compute depth and breadth of sequencing data
    GeneDepth = os.path.join(TmpDir,'gene_depth.txt')
    compute_depth_of_seq_data(GeneData,SeqData,GeneDepth)

    # TODO  2. compute gene abundance
    GeneAbun = os.path.join(TmpDir,'gene_abundance.txt')
    compute_gene_abundance(GeneData,GeneDepth,GeneAbun)

    # TODO  3. define seed gene
    SeedGene = os.path.join(TmpDir,'seed_gene.txt')
    find_seed_genes(GeneData,GeneAbun,SeedGene)
    
    # TODO  4. re-cluster sequencing data
    recluster_data(GeneData,SeqData,SeedGene)

    # TODO  5. strain-level assembly
    strain_call(TmpDir)
    
    # TODO mv result outside
    subprocess.call(['cp',os.path.join(TmpDir,'%s.fa'%opts.prefix),os.path.join(WorkDir,'%s.fa'%opts.prefix)])

    # back to WorkDir and delete TmpDir
    os.chdir(WorkDir)
    if os.path.exists(TmpDir):
        subprocess.call(['rm','-rf',TmpDir])


if __name__=="__main__":
    try:
        t_start = time.time()
        # TODO: set command-line parser
        parser = argparse.ArgumentParser(description=globals()['__doc__'],epilog=textwrap.dedent('''\
                                         Examples\n\
                                         --------\n\
                                         1. rbra [options] data_info -v\n\n \
                                         '''),\
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('data_info',help='data collection',metavar="DATA")
        parser.add_argument('-c','--cores',help='number of computing cores [1]',\
                            default=1,type=int,dest='cores',metavar='INT')
        # command-line options for StrainCall
        parser.add_argument('-D','--max-depth',help='downsample data to the specified depth [800]',\
                            default=800,type=int,dest='max_depth',metavar='INT')
        parser.add_argument('-q','--map-qual',help='only include reads with mapping quality >= INT [3]',\
                            default=3,type=int,dest='map_qual',metavar='INT')
        parser.add_argument('-i','--max-ins',help='only include reads with insert length <= INT [13]',\
                            default=13,type=int,dest='max_ins',metavar='INT')
        parser.add_argument('-l','--read-len',help='only include reads with length >= INT [70]',\
                            default=70,type=int,dest='read_len',metavar='INT')
        parser.add_argument('-t','--tau',help='only include strains with abundance level >= FLT [0.02]',\
                            default=0.02,type=float,dest='tau',metavar='FLT')
        parser.add_argument('-d','--diff-rate',help='only include strains with difference rate >= FLT [0.02]',\
                            default=0.02,type=float,dest='diff_rate',metavar='FLT')
        parser.add_argument('-g','--gene-sim',help='used in finding seed gene, merge genes with similarity >= FLT [0.93]',\
                            default=0.93,type=float,dest='gene_sim',metavar='FLT')
        parser.add_argument('-a','--gene-abun',help='used in finding seed gene, only include genes \n\
with abundance >= INT [13]',default=13,type=int,dest='gene_abun',metavar='INT')
        parser.add_argument('-p','--prefix',help='output filename prefix [16s_gene_assembly]',\
                            default='16s_gene_assembly',dest='prefix',metavar='STR')
        # other options
        parser.add_argument('-v',dest='verbose',action='store_true',default=False,help='verbose output')
        # parse options
        opts = parser.parse_args()

        # TODO: set logging format
        logging.basicConfig(format='[%(asctime)s] %(levelname)s : %(message)s', level=logging.INFO)

        # TODO enter to the pipeline
        SeqData,GeneData = parse_data_info(opts.data_info)       

        ExecDir = os.path.dirname(sys.argv[0])
        ExecDir = os.path.abspath(ExecDir)

        rbra_pipe(GeneData,SeqData)

        # TODO: complete
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


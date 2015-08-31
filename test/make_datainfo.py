#!/usr/bin/env python

import subprocess
import os,sys

# data directory
base_dir = os.path.dirname(os.path.abspath(__file__))

# bam file
bam_file = os.path.join(base_dir,'SRR606249_to_gg_99_otus.bam')

# gene sequence
gene_seq = os.path.join(base_dir,'99_otus.fasta')
# gene taxonomy
gene_tax = os.path.join(base_dir,'99_otu_taxonomy.txt')
# gene tree
gene_tree = os.path.join(base_dir,'99_otus_unannotated.tree')

# decompress gzip file
if os.path.exists('%s.gz'%gene_seq):
    sys.stderr.write('%s\n'%(' '.join(['gzip','-d','%s.gz'%gene_seq])))
    subprocess.call(['gzip','-d','%s.gz'%gene_seq])

if os.path.exists('%s.tar.gz'%gene_seq):
    sys.stderr.write('%s\n'%(' '.join(['gzip','-d','%s.fai.gz'%gene_seq])))
    subprocess.call(['gzip','-d','%s.fai.gz'%gene_seq])

if os.path.exists('%s.gz'%gene_tax):
    sys.stderr.write('%s\n'%(' '.join(['gzip','-d','%s.gz'%gene_tax])))
    subprocess.call(['gzip','-d','%s.gz'%gene_tax])

if os.path.exists('%s.gz'%gene_tree):
    sys.stderr.write('%s\n'%(' '.join(['gzip','-d','%s.gz'%gene_tree])))
    subprocess.call(['gzip','-d','%s.gz'%gene_tree])

# make bams.fofn
bams_fofn = os.path.join(base_dir,'bams.fofn')
with open(bams_fofn,'w') as f:
    f.write('%s\n'%bam_file)

# write to data_info.txt
with open(os.path.join(base_dir,'data_info.txt'),'w') as f:
    f.write('BamFiles=%s\n'%bams_fofn)
    f.write('SampleData=%s\n'%os.path.join(base_dir,'sample.txt'))
    f.write('GeneSeq=%s\n'%gene_seq)
    f.write('GeneTax=%s\n'%gene_tax)
    f.write('GeneTree=%s\n'%gene_tree)
    f.write('GeneIndex=%s.fai\n'%gene_seq)

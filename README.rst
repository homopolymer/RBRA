****
RBRA
****

Created by Feng Zeng Summer 2015.

Copyright 2015 Feng Zeng. All rights reserved.

========
Overview
========

RBRA is a **pipeline** for assembling the **strain-level** full-length **16S rRNA** genes.  It is designed for the large-scale metagenomic studies which contain tens to hundreds of samples.  RBRA uses the partial order graph (POG) and streaming Dirichlet process (DP) mixture clustering techniques to decode the composition of highly similar microorganisms. It is *lightweight*, *ultrafast*, *accurate* and *easy-to-use*.

==========
Dependency
==========

* Python >= 2.7.10
* mpi4py >= 1.3.1
* ETE2 Toolkit >= 2.3.1
* Scipy >= 0.15.1
* Numpy >= 1.9.2
* GCC >= 4.9.1
* Bowtie2 >= 2.2.4
* Samtools >= 1.2
* Bedtools >= 2.24.0

============
Installation
============

* Clone the repository::

    $ git clone https://github.com/homopolymer/RBRA.git

* Build and install to $PWD/bin::

    $ python setup.py install --prefix=$PWD

============
16S Database 
============

RBRA depends on the sequences, taxonomic annotations and phylogenetic tree of 16S rRNA genes.  GreenGenes (99_otus) contains 203,452 16S rRNA genes, and provides the separate files for sequence, taxonomy and tree.  GreenGenes is available at ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/.  While GreenGenes is constantly used during our development, RBRA works for other annotation databases, e.g. SILVA, if all of sequence, taxonomy and tree files are provided.

=====
Input
=====

The input of RBRA is a metadata file.  An example of metadata is given in the <test> directory.

A metadata file, e.g. data_info.txt, contains the following contents:

* Sequencing data
    1) BamFiles=<FILE>, a file listing the mapping files of data against 16S rRNA reference, one BAM file per line. ::

        $ ls -A1 *.bam | xargs realpath > bams.fofn

    2) SampleData=<FILE>, a tab-delimited file recording sample data.  An example is::

         #RUN         SAMPLE          MODE      SINGLE    PAIRED1              PAIRED2
         SRR606249    OAK_MOCK_WGS    paired    na        SRR606249_1.fastq    SRR606249_2.fastq

* Gene data
    1) GeneSeq=<FILE>, a FASTA file storing the reference sequences of 16S rRNA genes.
    2) GeneIndex=<FILE>, the index file of FASTA file.
    3) GeneTree=<FILE>, a Newick file recording the phylgenetic tree of 16S rRNA genes.
    4) GeneTax=<FILE>, a taxonomic annotation of 16S rRNA genes.

=====
Usage
=====

::

    $ python rbra.py [-c cores] [-v] [-p output_prefix] data_info.txt

====
Test
====

The directory <test> contains data and script to test whether RBRA installs and works correctly. ::

    $ cd test
    $ python make_datainfo.py
    $ ../bin/rbra.py -c 20 -v data_info.txt

================
Development Team
================

* Feng Zeng, Xiamen University
* Zicheng Wang, Tsinghua University
* Ting Chen, Tsinghua University


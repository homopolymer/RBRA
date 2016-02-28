****
RBRA
****

Created by Feng Zeng Summer 2015.

Copyright 2015 Feng Zeng. All rights reserved.

========
Overview
========

RBRA is a **pipeline** for assembling the **strain-level** full-length **16S rRNA** genes.  It is designed for the large-scale metagenomic studies which contain tens to hundreds of samples.  RBRA uses the partial order graph (POG) and streaming Dirichlet process (DP) mixture clustering techniques to decode the composition of highly similar microorganisms. It is *lightweight*, *ultrafast*, *accurate* and *easy-to-use*. RBRA works on both Illumina and Roche/454 sequencing data.

==========
Dependency
==========

* Python >= 2.7.10
* ETE2 Toolkit >= 2.3.1
* Scipy >= 0.15.1
* Numpy >= 1.9.2
* GCC >= 4.9.1
* Bowtie2 >= 2.2.4
* Samtools >= 1.2
* Bedtools >= 2.24.0
* Sickle >= 1.33
* Scikit-bio == 0.4.0

============
Installation
============

* Set up environment variable
  
    export CXX=<PATH TO G++4.9>

* Clone the repository::

    $ git clone https://github.com/homopolymer/RBRA.git

* Build and install to $PWD/bin::

    $ python setup.py install --prefix=$PWD

============
16S Database 
============

RBRA depends on the sequences, taxonomic annotations and phylogenetic tree of 16S rRNA genes.  GreenGenes (99_otus) contains 203,452 16S rRNA genes, and provides the separate files for sequence, taxonomy and tree.  GreenGenes is available at ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/.  While GreenGenes is constantly used during our development, RBRA works for other annotation databases, e.g. SILVA, if all of sequence, taxonomy and tree files are provided.

To download GreenGenes::
    
    $ bash download_greengenes.sh

=====
Input
=====

The input of RBRA is a metadata file.  An example of metadata is given in the <test> directory.

A metadata file, e.g. data_info.txt, contains the following contents:

* Sequencing data
    1) BamFiles=<FILE>, a file listing the mapping files of data against 16S rRNA reference, one BAM file per line. ::

        $ ls -A1 *_to_gg_99_otus.bam | xargs realpath > bams.fofn

* Gene data
    1) GeneSeq=<FILE>, a FASTA file storing the reference sequences of 16S rRNA genes.
    2) GeneIndex=<FILE>, the FASTA index file.
    3) GeneTree=<FILE>, a Newick file recording the phylgenetic tree of 16S rRNA genes.
    4) GeneTax=<FILE>, a taxonomic annotation of 16S rRNA genes.
    5) GeneAlign=<FILE>, a FASTA file of 16S sequence alignments.

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
    $ ../bin/rbra.py -c 20 -v data_info.txt 2>&1 | tee log.txt
    $ less log.txt

==========
Change Log
==========
* September 7, 2015, RBRA v0.3.1 released

    1) Add a filtration to StrainCall to filter out noisy assemblies
    
* September 7, 2015, RBRA v0.3.0 released

    1) Add coverage fraction into seed gene reference finding
    2) Add a script to count 16S assembly abundance per sample
    
* September 2, 2015, RBRA v0.2.0 released


    1) Replace mpi4py by multiprocessing
    2) Change mawk to awk
    3) Clean code

* August 31, 2015, RBRA v0.1.0 released

================
Development Team
================

* Feng Zeng, Xiamen University
* Zicheng Wang, Tsinghua University
* Ting Chen, Tsinghua University


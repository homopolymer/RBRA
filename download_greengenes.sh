#!/usr/bin/env bash

mkdir -p GreenGenes
cd GreenGenes

curl -C - -O ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/rep_set/99_otus.fasta
curl -C - -o 99_otus_aligned.fasta ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/rep_set_aligned/99_otus.fasta
curl -C - -O ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt
curl -C - -O ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/trees/99_otus_unannotated.tree 

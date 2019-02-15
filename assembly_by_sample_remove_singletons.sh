#!/bin/bash

## MOTUR analysis pipeline for actinobacterial diversity in Leontopodium tissues
## Created by Thomas Rattei and modified by Jacky Hess

mkdir -p /scratch/edelweiss
rm -f work
ln -s /scratch/edelweiss work

for ffastq in ../../data/*.1.fq ; do
  SAMPLE=$(basename $ffastq | sed 's/.1.fq//')
  ln -s -f /proj/metagenomes/edelweiss/data/$SAMPLE.1.fq work/$SAMPLE.1.fq
  ln -s -f /proj/metagenomes/edelweiss/data/$SAMPLE.2.fq work/$SAMPLE.2.fq
  egrep "primer|$SAMPLE" ../../data/16S_V3V4.oligos >work/$SAMPLE.oligos
  echo "make.contigs(ffastq=work/$SAMPLE.1.fq, rfastq=work/$SAMPLE.2.fq, oligos=work/$SAMPLE.oligos, bdiffs=1, pdiffs=4, processors=16)" | mothur
done
grep -c ">" work/*fasta >work/16S_V3V4.assembly.stats
cat work/*.contigs.groups >work/16S_V3V4.contigs.groups
cat work/*.trim.contigs.fasta >work/16S_V3V4.contigs.fasta
rm -f work/16S_V3V4_*
for file in data/*gz; do
  pigz -dc $file >work/$(basename $file | sed 's/.gz$//')
done
echo """summary.seqs(fasta=work/16S_V3V4.contigs.fasta, group=work/16S_V3V4.contigs.groups)
screen.seqs(fasta=work/16S_V3V4.contigs.fasta, group=work/16S_V3V4.contigs.groups, maxambig=0, minlength=380, maxlength=430)
unique.seqs(fasta=work/16S_V3V4.contigs.good.fasta)
count.seqs(name=work/16S_V3V4.contigs.good.names, group=work/16S_V3V4.contigs.good.groups)
align.seqs(fasta=work/16S_V3V4.contigs.good.unique.fasta, reference=work/silva.nr_v123.align, processors=16)
summary.seqs(fasta=work/16S_V3V4.contigs.good.unique.align, count=work/16S_V3V4.contigs.good.count_table)
screen.seqs(fasta=work/16S_V3V4.contigs.good.unique.align, count=work/16S_V3V4.contigs.good.count_table, summary=work/16S_V3V4.contigs.good.unique.summary, start=6428, end=23444, maxhomop=8)
summary.seqs(fasta=current, count=current)
filter.seqs(fasta=work/16S_V3V4.contigs.good.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=work/16S_V3V4.contigs.good.unique.good.filter.fasta, count=work/16S_V3V4.contigs.good.good.count_table)
pre.cluster(fasta=work/16S_V3V4.contigs.good.unique.good.filter.unique.fasta, count=work/16S_V3V4.contigs.good.unique.good.filter.count_table, diffs=2)
chimera.vsearch(fasta=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.fasta, count=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
remove.seqs(fasta=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
summary.seqs(fasta=current, count=current)
classify.seqs(fasta=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=work/trainset9_032012.pds.fasta, taxonomy=work/trainset9_032012.pds.tax, cutoff=80)
remove.lineage(fasta=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=current, count=current)
summary.seqs(fasta=current, count=current)
dist.seqs(fasta=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.20)
cluster(column=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"""| mothur

## filter out singleton clusters at this point

echo "remove.rare(list=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03, nseqs=1)
make.shared(list=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.pick.list, count=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, label=0.03)
classify.otu(list=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.pick.list, count=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, taxonomy=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)
dist.seqs(fasta=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, output=lt)
clearcut(phylip=work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist)""" | mothur
egrep '^tax|^3' work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tax.summary >read_summary.txt
egrep '^tax|^2' work/16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.tax.summary >otu_summary.txt
#rm -f work/silva* work/trainset* work/*dist work/*align work/*map
#tar -I /archive/bin/pbzip2 -cf data.tar.bz2 work/*
#rm -rf /scratch/edelweiss work

## Subset data to only contain Actinobacteria

echo "get.lineage(constaxonomy=16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.pick.0.03.cons.taxonomy, list=16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.pick.list, taxon=Actinobacteria, label=0.03)" | mothur

## Get representative FASTA sequences for each OTU to run alignments
echo "get.oturep(phylip=16S_V3V4.contigs.good.unique.good.filter.unluster.pick.pick.phylip.dist, list=16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.pick.0.03.pick.list, fasta=16S_V3V4.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)" | mothur
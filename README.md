#scaffold3C

scaffold3C is a bioinformatics analysis tool which attempts to determine the order (and potentially orientation) of genome assembly contigs using Hi-C linkage information. To that end, a genome or metagenome sequencing project will require both conventional shotgun (for assembly) and Hi-C sequencing data..

Finding the order and orientation of contigs within a chromosome is accomplished by transforming the problem of scaffolding into a Travelling Salesman Problem (TSP) and solving this using the local search algorithm [LKH](http://www.akira.ruc.dk/~keld/research/LKH/).

Currently scaffold3C takes as input the results of the genome-clustering tool [bin3C](https://github.com/cerebis/bin3C).

##Installation

###Dependencies

####proxigenomics_toolkit

 


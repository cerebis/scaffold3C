# scaffold3C

scaffold3C is a bioinformatics analysis tool which attempts to determine the order (and potentially orientation) of genome assembly contigs using Hi-C linkage information. To that end, a genome or metagenome sequencing project will require both conventional shotgun (for assembly) and Hi-C sequencing data..

Finding the order and orientation of contigs within a chromosome is accomplished by transforming the problem of scaffolding into a Travelling Salesman Problem (TSP) and solving this using the local search algorithm [LKH](http://www.akira.ruc.dk/~keld/research/LKH/).

Currently scaffold3C takes as input the results of the genome-clustering tool [bin3C](https://github.com/cerebis/bin3C).

## Installation

### Dependencies

- [proxigenomics_toolkit](https://github.com/cerebis/proxigenomics_toolkit)

### Installation steps

Numpy and Cython must be installed first. Note: NumPy's version is restricted due a compatibility issue within a dependency.

It is also highly recommended that a virtual environment is used as there are version requirements on some packages and the final executable of scaffold3C will be created in ```bin/```. 

1. Create a clean Python 2.7 environment

```$bash
mkdir scaffold3C
cd scaffold3C
virtualenv -p python2.7 .
```

2. Install NumPy < 1.15 and Cython.
   
```$bash
bin/pip install "numpy<1.15" cython
```

3. Install scaffold3C from github

```$bash
pip install git+https://github.com/cerebis/scaffold3C
```
 
Once pip completes installation, you will find an entry-point for scaffold3C at ```bin/scaffold3C```.

## Command-line interface

The complete description of the command-line interface ```bin/scaffold3C -h```

```$bash
usage: scaffold3C [-h] [-V] [-s SEED] [-v] [--clobber] [--log LOG]
                  [--max-image MAX_IMAGE] [--min-reflen MIN_REFLEN]
                  [--min-signal MIN_SIGNAL] [--min-size MIN_SIZE]
                  [--min-extent MIN_EXTENT] [--min-ordlen MIN_ORDLEN]
                  [--dist-method {inverse,neglog}] [--only-large]
                  [--skip-plotting] [--fasta FASTA]
                  MAP CLUSTERING OUTDIR

Create a 3C fragment map from a BAM file

positional arguments:
  MAP                   Contact map
  CLUSTERING            Clustering solution
  OUTDIR                Output directory

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         Show the application version
  -s SEED, --seed SEED  Random integer seed
  -v, --verbose         Verbose output
  --clobber             Clobber existing files
  --log LOG             Log file path [OUTDIR/scaffold3C.log]
  --max-image MAX_IMAGE
                        Maximum image size for plots [4000]
  --min-reflen MIN_REFLEN
                        Minimum acceptable reference length [1000]
  --min-signal MIN_SIGNAL
                        Minimum acceptable trans signal [5]
  --min-size MIN_SIZE   Minimum cluster size for ordering [5]
  --min-extent MIN_EXTENT
                        Minimum cluster extent (kb) for ordering [50000]
  --min-ordlen MIN_ORDLEN
                        Minimum length of sequence to use in ordering [2000]
  --dist-method {inverse,neglog}
                        Distance method for ordering [inverse]
  --only-large          Only write FASTA for clusters longer than min_extent
  --skip-plotting       Skip plotting the contact map
  --fasta FASTA         Alternative location of source FASTA from that
                        supplied during clustering
```

## Using scaffold3C

To scaffold a metagenomic assembly, you will first need to analyse the assembly with the Hi-C based genome binning tool [bin3C](https://github.com/cerebis/bin3C).

### Ordering only

After bin3C has analyzed a metagenome, you will have both a contact map and clustering result. These two files form the input to scaffold3C. A standard run will attempt to find the order of each MAG over a minimum total extent and cluster size (default: 50kbp and 5 contigs).

Eg. Ordering MAGs for a Hi-C dataset generated from two enzymatic digestions (Sau3AI and MluCI). 
```$bash
bin/bin3C mkmap -v --seed 1234 -e Sau3AI -e MluCI --min-reflen 5000 [fasta] [bam_file] [out_dir]

bin/bin3C clustering -v --seed 1234 [contact_map] [out_dir]

bin/scaffold3C -v --seed 1234 [contact_map] [clustering] [out_dir]
``` 

### Order and orientation

Inferring the orientation of contigs requires that bin3C create what I have termed a tip-based contact map. This procedure tracks only reads which map within a limited region at the two ends of each contig (the tips). For sanity, tip sizes are constrained to be no larger than the minimum acceptable contig length (tip-size <= min-reflen).

At present, for better performance it is recommended that tip-size be no smaller than 5000bp and minimum contig length be 10kbp.

scaffold3C will automatically detect the presence of a tip-based contact map and solve for both order and orientation.

Eg. Finding both order and orientation for the same dataset above. 
```$bash
bin/bin3C mkmap -v --seed 1234 -e Sau3AI -e MluCI --min-reflen 10000 --tip-size 5000 [fasta] [bam_file] [out_dir]

bin/bin3C clustering -v --seed 1234 [contact_map] [out_dir]

bin/scaffold3C -v --seed 1234 [contact_map] [clustering] [out_dir]
```

## Algorithm details

<under construction>

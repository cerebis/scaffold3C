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
 
Once pip completes installation, you will find an entry-point for scaffold3C at ```./bin/scaffold3C```.

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

To order the contigs of each MAG identified by bin3C, you will need to let bin3C at the time of map creation to produce an additional extent map. This is accomplished by specifying a `bin_size` during the `mkmap` stage. 

An extent map is a contact map where the bins are not simply whole sequences, but rather are measured in base-pair extent. In this way, a bin size of 5kbp will produce 20 bins across a 100kbp contig. 

### Choosing a bin size

Keep in mind that the size of the contact map grows with the square of the number of bins. Although scaffold3C and bin3C make use of sparse matrices, memory constraints can still play a factor. This is particularly the case for complex and deeply sequenced metagenomic assemblies, where there can be >100,000 contigs of significant length.

Further, making bins too small will weaken the Hi-C signal by spreading it over too wide an area. Standard bin3C clustering runs perform fine on contigs down to 1000bp, while the more demanding task of scaffolding requires more signal. A bin_size of 5-10 kbp is sufficiently smal and will not exclude smaller contigs which are still larger than the minimum contig length.

An extent map will be computing by specifying a bin_size during bin3C map creation.

```$bash
bin3C mkmap --bin_size 10000 -e MluCI -e Sau3AI [fasta] [bam_file] [out_dir]
```

Basic usage of scaffold3C after a bin3C analysis is completed.

```$bash
bin/scaffold3C [contact_map] [clustering] [out_dir]
```

#### Orientation

To infer orientation requires a further option `tip_size` to be supplied at the time of map creation in bin3C. This procedure roughly doubles the storage requirements of the extent map and splits the Hi-C signal in two parts, tracking the ends of each contig independently.


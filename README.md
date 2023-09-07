# PicklSeq
A python command line tool to turn a fastq and fasta to a pickle file for ML or other stuff.
Currently, only mac is fully supported, unix based os should also work as well but no tests have been done on it.

# Installation
1. Clone this repo.
2. install minimap2: https://github.com/lh3/minimap2
3. install samtools: http://www.htslib.org/
4. install chopper: https://github.com/wdecoster/chopper

tip: Use homebrew for Mac

Version used during development:
Python 3.11.3
minimap2 2.26-r1175
samtools 1.17 Using htslib 1.17
chopper 0.5.0

# Run
Run picklseq.py to trim sequences to a pickl file
```
usage: picklseq.py [-h] [-f FILE] [-t TYPE] [-m MINLENGTH] [-M MAXLENGTH] [-q QUALITY] [-T THREADS] [-c MATCHCOUNT]
                   [-p] [-k]

A script to trim DNA sequences

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  The file name of the fastq
  -t TYPE, --type TYPE  The type of the DNA sequence
  -m MINLENGTH, --minlength MINLENGTH
                        min length of acceptable DNA sequence
  -M MAXLENGTH, --maxlength MAXLENGTH
                        max length of acceptable DNA sequence
  -q QUALITY, --quality QUALITY
                        quality of Chopper reads
  -T THREADS, --threads THREADS
                        number of CPU threads
  -c MATCHCOUNT, --matchcount MATCHCOUNT
                        number of the number of CIGAR matches
  -p, --pad             pad the start of the DNA sequence when pos is not 1
  -k, --keep            keep subprocess file output
```

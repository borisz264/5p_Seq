# 5p-seq

for mapping of 5p-seq reads to a yeast genome using STAR

## Getting Started
*   Download all of the dependencies and make sure they are in your PATH or PYTHONPATH respectively
*   Assemble your genome FASTA files into a single folder.
*   Assemble FASTA files of noncoding RNA into a different folder. it is recommended that rRNA have 'rRNA' int he sequence names, and tRNAs have 'tRNA' or 'MT' in the sequence names.
*   Make sure the FASTA files are uncompressed and end in .fa (sorry STAR compiler doesn't seem to do compression)
*   gzip your fastq files and put them in one folder
*   make a settings file for your experiment. An example file is included, and an annotated version is included, as comments are not allowed in json files.

*   If you change a parameter and restart the pipeline, it is best to delete all of the output, or at least all of the output past the affected point. Otherwise the pipeline will use your old intermediates.

## running pipeline
`python fpseq_main.py /path/to/settings.json.txt`


This will run the basic read trimming and mapping

#### Optional Parameters
* **--threads**   Specify number of threads to use. Default is 8.


## Prerequisites
*   I highly recommend 32 GB of RAM for a human or mouse genome. Mapping can be done with less (with a speed tradeoff). This requirement scales with genome size, not dataset size.

#### Dependencies
*   skewer read adaptor trimmer downloaded and added to path (https://github.com/relipmoc/skewer) (0.2.2 tested)
*   STAR RNA-seq aligner downlaoded and added to path (tested versionSTAR 020201)
*   samtools (http://www.htslib.org/ version 1.5 tested) compiled and added to path (make install may do this automatically)
*   FASTX-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/) tested with 0.0.14
*   UMI-tools (https://github.com/CGATOxford/UMI-tools/blob/master/doc/QUICK_START.md)

#### Python dependencies (for python 2.7) installed with pip install (on mac, a local install with --user may be recommended)
*   simplejson (3.11.1)
*   numpy (1.13.1)
*   scipy (0.19.1)


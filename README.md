# sexassign
The `sexassign.py` script uses counts of reads mapping to the X chromosome
and autosomes to determine the genetic sex of an individual.  Under the
assumption that the number of reads obtained should reflect the chromosome
copy numbers and chromosome lengths, two binomial models are constructed,
one for males and one for females.  A likelihood ratio test is then used
to distinguish between the two models.  This is effective with as few as
5000 mapped shotgun reads.

# Prerequisites
The `sexassign.py` script requires python 2.7, numpy, scipy, and matplotlib.

# Usage
Read counts must first be extracted from mapped BAM files using
`samtools idxstats`, and `sexassign.py` can then be used output a
tab-separated table of sex assignments, and to plot read dosages.

```
$ for i in *.bam; do samtools idxstats ${i} > ${i}.idxstats; done
$ sexassign.py -w -o sexes.pdf *.idxstats > sexes.txt
```

Further information can be obtained from the help output.

```
$ sexassign.py  -h
usage: sexassign.py [-h] [-w] [-o OPDF] [-n MIN_READS] [-l MIN_LENGTH]
                    [-X CHRX_LIST] [-e EXCLUDE_LIST] [-i INCLUDE_LIST] [-p]
                    infiles [infiles ...]

sex determination from read dosages

positional arguments:
  infiles               input file(s)

optional arguments:
  -h, --help            show this help message and exit
  -w, --wide            plot widescreen ratio (16x9) [False]
  -o OPDF, --opdf OPDF  output filename [out.pdf]
  -n MIN_READS, --min-reads MIN_READS
                        exclude samples with fewer reads than this [5000]
  -l MIN_LENGTH, --min-length MIN_LENGTH
                        exclude contigs shorter than this [1000000]
  -X CHRX_LIST, --chrX-list CHRX_LIST
                        file containing list of X linked contigs, if not
                        specified, defaults to X, chrX, ChrX
  -e EXCLUDE_LIST, --exclude-list EXCLUDE_LIST
                        file containing list of contigs to exclude, if not
                        specified, defaults to Y, chrY, ChrY, M, MT, Mt
  -i INCLUDE_LIST, --include-list INCLUDE_LIST
                        file containing list of contigs to include, if not
                        specified, defaults to all contigs
  -p, --pca             Plot PCA of read dosages
```

### Usage
```
cat NCLscan.result | ./ncl_exon_number_extractor.py Annotation.gtf.gz [--show-all] [--expand] [--detail] > OUTPUT
```

### Input
The first 6 columns of the input file should be the followings:
```
(1) Chromosome name of the donor side (5'ss) 
(2) Junction coordinate of the donor side
(3) Strand of the donor side
(4) Chromosome name of the acceptor side (3'ss) 
(5) Junction coordinate of the acceptor side
(6) Strand of the acceptor side
```

### Output
The following 9 columns are **appended**:
```
(1) ID of the longest transcript (donor)
(2) ID of the longest transcript (accepter)
(3) Exon number (donor)
(4) Exon number (accepter)
(5) Length of flanking intron (donor)
(6) Length of flanking intron (accepter)
(7) Total length of intermediate exons
(8) Total number of exons of the transcript (donor)
(9) Total number of exons of the transcript (accepter)
```

In `--detail` mode, there are 4 additional columns:
```
(10) Is protein-coding transcript (donor)
(11) Is protein-coding transcript (acceptor)
(12) Total length of the transcript (donor)
(13) Total length of the transcript (acceptor)
```

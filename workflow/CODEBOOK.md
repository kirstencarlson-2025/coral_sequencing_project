# Code notebook
Kirsten Carlson
<br>
<br>Coding notebook for recording code used, troubleshooting steps, brainstorming, etc.
<br>
```bash
Date: 
Goal:     
Code Used:
```

```bash
Date: 3/15/2026
<br>Goal: Using MarineOmics.github workflow, use SeqArray in R to evaluate missingness.
<br>Code Used: (code documentation in /scratch/kcarls36/projects/data/alignments/sint_align/discosnp/k25_D5/SeqArray_K25_D5.Rmd)
<br>Before converting vcf to gds with Seq Array, in command line:
<br>- sort the clustered vcf
<br>- rename headers (discnosnp names each sample G1-Gx in the order of the fof.txt - create a reheader.txt from fof and reheader for correct sample IDs)
<br>- simplify vcf for faster processing (-x INFO,FORMAT,CONTIG)
<br>- index final vcf
<br>
<br>Use SeqArray to convert final vcf to gds format (I requested 32G and 4:00:00 in R session on hpc. Process took 1:15:00)
```

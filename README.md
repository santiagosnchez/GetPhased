# GetPhasedGenome
RANDOMLY sorts out IUPAC ambiguities into different haplotypes

## IMPORTANT NOTE:
Do not use these sequences for phylogenetics!!
They should only be used for site-based summary statistics.

The program gets phased sequences that result from running read-based phasers, such as [HapCompass](http://www.brown.edu/Research/Istrail_Lab/hapcompass.php). FASTQ files should be first parsed by OrderFromSamtools.pl (found in this repository). You also need the original alignments (not phased), which are used to retrieve non-polymorphic regions. In addition, the script can randomly resolve polymorphisms, which can be used for instance, to estimate frequency-based diversity statistics. 

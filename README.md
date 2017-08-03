# GetPhased
Sorts out IUPAC ambiguities into different haplotypes

The program gets phased sequences that result from running read-based phasers, such as [HapCompass](http://www.brown.edu/Research/Istrail_Lab/hapcompass.php). FASTQ files should be first parsed by [OrderFromSamtools.pl] (https://github.com/santiagosnchez/OrderFromSamtools). You also need the original alignments (not phased), which are used to retrieve non-polymorphic regions. In addition, the script can randomly resolve polymorphisms, which can be used for instance, to estimate frequency-based diversity statistics.

## IMPORTANT NOTE:
Do not use RANDOMLY resolved sequences for phylogenetics!!
They should only be used for site-based summary statistics.

## Installation

    git clone https://github.com/santiagosnchez/GetPhased
    cd GetPhased
    chmod +x GetPhased.pl
    sudo cp GetPhased.pl /usr/local/bin

## 

# GetPhased
Sorts out IUPAC ambiguities into different haplotypes

The program gets phased sequences that result from running read-based phasers, such as [HapCompass](http://www.brown.edu/Research/Istrail_Lab/hapcompass.php). FASTQ files should be first parsed by [OrderFromSamtools.pl](https://github.com/santiagosnchez/OrderFromSamtools). You also need the original alignments (not phased), which are used to retrieve non-polymorphic regions. In addition, the script can randomly resolve polymorphisms, which can be used for instance, to estimate frequency-based diversity statistics.

## IMPORTANT NOTE:
Do not use RANDOMLY resolved sequences for phylogenetics!!
They should only be used for site-based summary statistics.

## Installation

    git clone https://github.com/santiagosnchez/GetPhased
    cd GetPhased
    chmod +x getPhased.pl
    sudo cp getPhased.pl /usr/local/bin

## Running the code

Use the `-h` flag for more details:

    perl getPhased.pl -h
    
    Try:
    perl getPhased.pl -hap hapfile.fasta    [ FASTA alignemnt with haploid sequences; result from OrderFromSamtools.pl ]
                      -aln alignment.fasta  [ FASTA alignment with original sequences; i.e. with ambiguities ]
                      -out outfile.fasta    [ Your result will be stored here ]
                      -res                  [ Optional argument for randomly resolving ambiguities, in case no haploid sequence were found ]

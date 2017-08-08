#!/usr/bin/perl
# ©Santiago Sanchez-Ramirez, University of Toronto

no warnings;
my $hapfile;
my $alnfile;
my $outfile;
my $res=0;

if (grep { /^-he{0,1}l{0,1}p{0,1}$/ } @ARGV){
	die "
Try:
perl getPhased.pl -hap hapfile.fasta    [ FASTA alignemnt with haploid sequences; result from OrderFromSamtools.pl ]
                  -aln alignment.fasta  [ FASTA alignment with original sequences; i.e. with ambiguities ]
                  -out outfile.fasta    [ Your result will be stored here ]
                  -res                  [ Optional argument for randomly resolving ambiguities, in case no haploid sequences were found ]\n\n";
}

if (grep { /^-res$/ } @ARGV){
	$res=1;
}

if (my ($indHap) = grep { $ARGV[$_] =~ /^-hap$/ } 0 .. $#ARGV){
	$hapfile = $ARGV[$indHap+1];
} else {
	die "-hap flag not found.\n";
}

if (my ($indAln) = grep { $ARGV[$_] =~ /^-aln$/ } 0 .. $#ARGV){
	$alnfile = $ARGV[$indAln+1];
} else {
	die "-aln flag not found.\n";
}

if (my ($indOut) = grep { $ARGV[$_] =~ /^-out$/ } 0 .. $#ARGV){
	$outfile = $ARGV[$indOut+1];
} else {
	die "-out flag not found.\n";
}

my %alnhash=();
my %haphash=();
my @head=();
my @seqs=();

open(ALN, "<", $alnfile);
while(<ALN>){
	next if (/^$/);
	if (/^>/){
		chomp($_);
		push @head, $_;
	} else {
		chomp($_);
		push @seqs, $_;
	}
}

if (scalar(@head) == scalar(@seqs)){
	for my $i (0 .. $#head){
		$seqs[$i] =~ tr/a-z/A-Z/;
		$alnhash{$head[$i]} = $seqs[$i];
	}
} else {
	die "Unequal number of headers and sequences in $alnfile. FASTAs must be in sequential format.\n";
}

my @head=();
my @seqs=();

open(HAP, "<", $hapfile);
while(<HAP>){
	next if (/^$/);
	if (/^>/){
		chomp($_);
		push @head, $_;
	} else {
		chomp($_);
		push @seqs, $_;
	}
}

if (scalar(@head) == scalar(@seqs)){
	for my $i (0 .. $#head){
		$seqs[$i] =~ tr/a-z/A-Z/;
		$seqs[$i] =~ s/\?//g;
		$haphash{$head[$i]} = $seqs[$i];
	}
} else {
	die "Unequal number of headers and sequences in $hapfile. FASTAs must be in sequential format.\n";
}

my @seqnumb=();

open(OUT, ">", $outfile);
foreach(sort {$a cmp $b} keys %alnhash){
	if (exists $haphash{$_}){
		@hapseq = split '', $haphash{$_};
		@alnseq = split '', $alnhash{$_};
		my @phased1=();
		my @phased2=();
		for my $i (0 .. $#hapseq){
			if (($hapseq[$i] eq 'N') and ($alnseq[$i] ne 'N')){
				if ($res == 0){
					push @phased1, $alnseq[$i];
					push @phased2, $alnseq[$i];
				} else {
					if ($alnseq[$i] eq 'Y'){
						push @phased1, 'C';
						push @phased2, 'T';
					}
					elsif ($alnseq[$i] eq 'R'){
						push @phased1, 'A';
						push @phased2, 'G';
					}
					elsif ($alnseq[$i] eq 'W'){
						push @phased1, 'T';
						push @phased2, 'A';
					}
					elsif ($alnseq[$i] eq 'K'){
						push @phased1, 'T';
						push @phased2, 'G';
					}
					elsif ($alnseq[$i] eq 'M'){
						push @phased1, 'C';
						push @phased2, 'A';
					}
					elsif ($alnseq[$i] eq 'S'){
						push @phased1, 'C';
						push @phased2, 'G';
					} else {
						push @phased1, $alnseq[$i];
						push @phased2, $alnseq[$i];
					}
				}		
			}
			elsif (($hapseq[$i] eq $alnseq[$i]) and ($res == 0)){
				push @phased1, $alnseq[$i];
				push @phased2, $alnseq[$i];
			}
			elsif (($hapseq[$i] eq $alnseq[$i]) and ($res == 1)){
				if (($alnseq[$i] eq 'Y') and ($hapseq[$i] eq 'Y')){
					push @phased1, 'C';
					push @phased2, 'T';
				}
				elsif (($alnseq[$i] eq 'R') and ($hapseq[$i] eq 'R')){
					push @phased1, 'A';
					push @phased2, 'G';
				}
				elsif (($alnseq[$i] eq 'M') and ($hapseq[$i] eq 'M')){
					push @phased1, 'A';
					push @phased2, 'C';
				}
				elsif (($alnseq[$i] eq 'W') and ($hapseq[$i] eq 'W')){
					push @phased1, 'A';
					push @phased2, 'T';
				}
				elsif (($alnseq[$i] eq 'K') and ($hapseq[$i] eq 'K')){
					push @phased1, 'G';
					push @phased2, 'T';
				}
				elsif (($alnseq[$i] eq 'S') and ($hapseq[$i] eq 'S')){
					push @phased1, 'C';
					push @phased2, 'G';
				} else {
					push @phased1, $alnseq[$i];
					push @phased2, $alnseq[$i];
				}
			}
			elsif (($alnseq[$i] eq 'Y') and ($hapseq[$i] eq 'C')){
				push @phased1, 'C';
				push @phased2, 'T';
			}
			elsif (($alnseq[$i] eq 'Y') and ($hapseq[$i] eq 'T')){
				push @phased1, 'T';
				push @phased2, 'C';
			}
			elsif (($alnseq[$i] eq 'R') and ($hapseq[$i] eq 'A')){
				push @phased1, 'A';
				push @phased2, 'G';
			}
			elsif (($alnseq[$i] eq 'R') and ($hapseq[$i] eq 'G')){
				push @phased1, 'G';
				push @phased2, 'A';
			}
			elsif (($alnseq[$i] eq 'S') and ($hapseq[$i] eq 'G')){
				push @phased1, 'G';
				push @phased2, 'C';
			}
			elsif (($alnseq[$i] eq 'S') and ($hapseq[$i] eq 'C')){
				push @phased1, 'C';
				push @phased2, 'G';
			}
			elsif (($alnseq[$i] eq 'K') and ($hapseq[$i] eq 'G')){
				push @phased1, 'G';
				push @phased2, 'T';
			}
			elsif (($alnseq[$i] eq 'K') and ($hapseq[$i] eq 'T')){
				push @phased1, 'T';
				push @phased2, 'G';
			}
			elsif (($alnseq[$i] eq 'W') and ($hapseq[$i] eq 'T')){
				push @phased1, 'T';
				push @phased2, 'A';
			}
			elsif (($alnseq[$i] eq 'W') and ($hapseq[$i] eq 'A')){
				push @phased1, 'A';
				push @phased2, 'T';
			}
			elsif (($alnseq[$i] eq 'M') and ($hapseq[$i] eq 'A')){
				push @phased1, 'A';
				push @phased2, 'C';
			}
			elsif (($alnseq[$i] eq 'M') and ($hapseq[$i] eq 'C')){
				push @phased1, 'C';
				push @phased2, 'A';
			}
			else {
				push @phased1, $alnseq[$i];
				push @phased2, $alnseq[$i];
			}
		}
		if (scalar(@phased1) < scalar(@alnseq)){
			if ($res == 0){
				push @phased1, @alnseq[$#phased1+1 .. $#alnseq];
				push @phased2, @alnseq[$#phased2+1 .. $#alnseq];
			} else {
				push @phased1, @alnseq[$#phased1+1 .. $#alnseq];
				push @phased2, @alnseq[$#phased2+1 .. $#alnseq];
				map { s/Y/C/g; s/R/A/g; s/W/A/g; s/K/G/g; s/S/C/g; s/M/C/g; } @phased1;
				map { s/Y/T/g; s/R/G/g; s/W/T/g; s/K/T/g; s/S/G/g; s/M/A/g; } @phased2;
			}
		}
		if ((scalar(@phased1) == scalar(@alnseq)) && (scalar(@phased2) == scalar(@alnseq))){
			if ($res == 0){
				if (@phased1 ~~ @phased2){
					print OUT $_ . "\n" . join('', @phased1) . "\n";
					push @seqnumb, 1;
				} else {
					print OUT $_ . "_a\n" . join('', @phased1) . "\n";
					print OUT $_ . "_b\n" . join('', @phased2) . "\n";
					push @seqnumb, 2;
				}
			} else {
				print OUT $_ . "_a\n" . join('', @phased1) . "\n";
				print OUT $_ . "_b\n" . join('', @phased2) . "\n";
				push @seqnumb, 2;
			}
		} else {
			die "Phased and original sequences differ in length\n";
			`rm $outfile`;
		}
	} else {
		if ($res == 0){
			print OUT $_ . "\n" . $alnhash{$_} . "\n";
			push @seqnumb, 1;
		} else {
			$seqA = $alnhash{$_};
			$seqB = $alnhash{$_};
			$seqA =~ s/Y/C/g;
			$seqB =~ s/Y/T/g;
			$seqA =~ s/R/A/g;
			$seqB =~ s/R/G/g;
			$seqA =~ s/W/A/g;
			$seqB =~ s/W/T/g;
			$seqA =~ s/K/G/g;
			$seqB =~ s/K/T/g;
			$seqA =~ s/S/C/g;
			$seqB=~ s/S/G/g;
			$seqA =~ s/M/C/g;
			$seqB =~ s/M/A/g;
			print OUT $_ . "_a" . "\n" . $seqA . "\n";
			print OUT $_ . "_b" . "\n" . $seqB . "\n";
			push @seqnumb, 2;
		}
	}
}

my $x = scalar(keys %alnhash);
my $y = eval join '+', @seqnumb;

print "Initial number of sequences $x\nFinal number of sequences $y\n";



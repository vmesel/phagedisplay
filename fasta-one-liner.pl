#!/usr/bin/perl
# Programmer: Elton Vasconcelos (02/Aug/2013)
# Script that place the entire sequence in one single line in each sequence inside the FASTA file
# Usage: perl fasta-single-line.pl [infile.fasta] >outfile.fasta
#################################################################################################

open (FILE, "$ARGV[0]") or die ("Can't open file $ARGV[0]\n");
my $line = <FILE>;
chomp($line);
my $seq;
my $headline;

while ($line ne "") {
    if ($line =~ m/^>.*/) {
        $headline = $line;
        $line = <FILE>;
        chomp($line);
        until ($line =~ m/^>/ || $line eq "") {
            $seq .= $line;            #putting the whole chromosome sequence in one single line
            $line = <FILE>;
            chomp($line);
        }
        print("$headline\n$seq\n");
        $seq = "";
    }
    else {        #Just in the case the FASTA file is not starting with a headline (^>.+)
        $line = <FILE>;
        chomp($line);
    }
}

#!/usr/bin/perl
#Subset samtools mpileup for single nucleotide polymorphisms. This script will ignore indels and non polymorphic positions.
#1) This script will first identify positions that are covered by sequence reads in both of the parental genotypes.
#2) Next it will identify positions that are polymorphic between the two parental genotypes.
#3) Finally it will scan through the example DH RIL and only output positions that are segregating between the two original parents.
#4) The output data will also contain a new column for genotype calls. Parent 1 = 0, Parent 2 = 1.

#This script requires 3 mpileup files - 1 for each parent, and 1 for the focal DH RIL. It also requires the user to input the output file name.

use strict; use warnings;
die "usage: 00.extract_SNPs_from_mpileup.pl <Parent1 mpileup> <Parent2 mpileup> <DH mpileup> <output file>\n" unless @ARGV==4;

open (P1, "<",$ARGV[0]) or die "Cannot open PARENT1 MPILEUP file\n"; #Open mpileup for Parent 1
open (P2, "<",$ARGV[1]) or die "Cannot open PARENT2 MPILEUP file\n"; #Open mpileup for Parent 2
open (DH, "<",$ARGV[2]) or die "Cannot open DH MPILEUP file\n"; #Open mpileup for DH RIL
open (OUT, ">",$ARGV[3]) or die "Cannot open output file\n"; #Open output file


##########
#1) Identify positions that have genotype information in both parents
my %p1; #Store the ID of positions covered in Parent 1
my %p2; #Store the ID of positions covered in Parent 2

#Extract positions covered in Parent 1
while (<P1>) {
  chomp;
  my ($chr,$pos,$ref,$cov,$alt,$qual) = split "\t",$_; #Store information for each line of the mpileup file. Chromosome, position, Reference base, Read coverage, Alternative base, Base quality

  if ($alt =~ m/^[ATGC.,]$/) { #Check if the alternative base matches A,T,G,C or . or ,. The comma and the period denote positions where the sequenced read matches the reference base. This line also ensures that the base in question is only 1 nucleotide long, removing indels.
    #If the alternative base column is a period or comma, store the reference base in %p1.
    if ($alt =~ m/^[.,]$/) {$p1{$chr."_".$pos} = $ref} #Unique IDs are stored as Chr_pos.
    #If the alternative base column is an A,T,G, or C, store the alternative base in %p1.
    if ($alt =~ m/^[ATGC]$/) {$p1{$chr."_".$pos} = $alt}
  }
}

#Extract positions covered in Parent 1
while (<P2>) {
  chomp;
  my ($chr,$pos,$ref,$cov,$alt,$qual) = split "\t",$_; #Store information for each line of the mpileup file. Chromosome, position, Reference base, Read coverage, Alternative base, Base quality

  if ($alt =~ m/^[ATGC.,]$/) { #Check if the alternative base matches A,T,G,C or . or ,. The comma and the period denote positions where the sequenced read matches the reference base. This line also ensures that the base in question is only 1 nucleotide long, removing indels.
    #If the alternative base column is a period or comma, store the reference base in %p2.
    if ($alt =~ m/^[.,]$/) {$p2{$chr."_".$pos} = $ref}
    #If the alternative base column is an A,T,G, or C, store the alternative base in %p2.
    if ($alt =~ m/^[ATGC]$/) {$p2{$chr."_".$pos} = $alt}
  }
}

#Count the number of positions covered in each parent and output that information
my $p1pos = scalar keys %p1;
my $p2pos = scalar keys %p2;

print $p1pos." positions covered in Parent 1\n";
print $p2pos." positions covered in Parent 2\n";



##########
#2) Identify postions that are polymorphic between the 2 parents by comparing %p1 and %p2
my $overlap; #Keep track of the number of overlapping positions
my $snpcount; #Keep track of the number of overlapping positions that are polymorphic
my %snp; #Store information for polymorphic positions

foreach my $p (keys %p1) { #For each unique id in %p1
  if (exists $p2{$p}) { #Ask if that unique id is also present in %p2
    $overlap++; #If so, count the position as overlapping

    #For overlapping positions, identify positions that do not share the same sequence
    if ($p1{$p} !~ m/$p2{$p}/) {
      $snpcount++; #Keep track of the number of polymorphic positions
      $snp{$p} = $p1{$p}.$p2{$p}; #Store the nucleotide sequence of P1 and P2 for each unqiue SNP id.
    }
  }
}

print $overlap." positions overlap between parents\n";
print $snpcount. " positions are polymorphic between parents\n";


#########
#3) and 4) Output positions covered in the DH RIL that were determined to be polymorphic in the original parents. Also, output genotype calls for each position. Parent 1 = 0, Parent 2 = 1.
my $dhsnpcount; #Keep track of the number of polymorphic positions detected in the DH RIL

while (<DH>) {
  chomp;
  my ($chr,$pos,$ref,$cov,$alt,$qual) = split "\t",$_; #Store information for each line of the mpileup file. Chromosome, position, Reference base, Read coverage, Alternative base, Base quality

  #Only keep positions that are polymorphic in the parents
  if (exists $snp{$chr."_".$pos}) { #Check if the current position was determined to be polymorphic between the 2 parents (%snp)

    if ($alt =~ m/^[ATGC.,]$/) { #Make sure the the current position is only 1 nucleotide long (remove indels in DH mpileup)
      $dhsnpcount++; #Count SNPs detected in the DH
      print OUT $chr."\t".$pos."\t".$ref."\t".$cov."\t".$alt."\t"; #Print DH RIL mpileup information for polymorphic positions to output file

      #4) For each position, identify if the DH RIL is carrying the Parent 1 allele (0) or the Parent 2 allele (1)
      my $geno; #Create variable to store genotype
      if ($alt =~ m/^[.,]$/) {$geno = $ref} #If the alternative bases is a period or a comma, the genotype is the reference allele
      if ($alt =~ m/^[ATGC]$/) {$geno = $alt} #If the alternative bases is A,T,G,or C, the genotype is the alternative base

      if ($geno =~ m/$p1{$chr."_".$pos}/) {print OUT "0"} #If the genotype matches the Parent 1 allele, output 0.
      if ($geno =~ m/$p2{$chr."_".$pos}/) {print OUT "1"} #If the genotype matches the Parent 2 allele, output 1.
      print OUT "\n"; #Output a newline
    }
  }
}

print "Number of SNPs detected in DH: ".$dhsnpcount."\n";


close (P1); #Close Parent 1 mpileup
close (P2); #Close Parent 2 mpileup
close (DH); #Close DH RIL mpileup
close (OUT); #Close output file

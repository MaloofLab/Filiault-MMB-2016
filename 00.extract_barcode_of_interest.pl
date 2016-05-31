#!/usr/bin/perl
#Extract sequences with a given 5' barcode from a raw FASTQ file
#Provide the input FASTQ sequence, the desired barcode sequence (capital letter), and an output file name
use strict; use warnings;
die "usage: 00.extract_barcode_of_interest.pl <FASTQ file> <barcode sequence> <OUTPUT file>\n" unless @ARGV==3;

open (FASTQ, "<",$ARGV[0]) or die "Cannot open FASTQ file\n"; #Open FASTQ file
open (OUT, ">",$ARGV[2]) or die "Cannot open output file\n"; #Open empty output file

my $barcode = $ARGV[1]; #Get barcode sequence

##########
#These variables will be used to keep track of the line number.
#FASTQ file format is dependent on 4 lines of information and these lines iterate throughout the file.
my $linecount;
my $id1; #Line 1 - Sequence ID. Begins with the @ character
my $seq; #Line 2 - Actual Sequence read
my $id2; #Line 3 - Sequence ID repeated. Begins with the + character
my $qual; #Line 4 - Quality score of each sequenced based in line 2.


##########
#Loop through the raw FASTQ file.
#Store each iteration of the file (every for lines) and then test if the sequence matches the desired barcode. If so, print out the current FASTQ iteration.
while (<FASTQ>) {
  chomp;
  $linecount++; #Keep track of the line number
  if ($linecount==1) {$id1 = $_} #Store Line 1
  if ($linecount==2) {$seq = $_} #Store Line 2
  if ($linecount==3) {$id2 = $_} #Store Line 3
  if ($linecount==4) {
    $qual = $_; #Store Line 3

    #Check if current sequence (line 2) matches the 5' barcode of interest
    if ($seq =~ m/^$barcode/){
      #If the sequence contains the desired barcode, remove the barcode and KpnI restriction enzyme site CATG (5+4 nt)
      my $shortseq = substr($seq,9); #Remove first 9 bases from current sequence read
      my $shortqual = substr($qual,9); #Remove first 9 bases from current quality read
      print OUT $id1."\n".$shortseq."\n".$id2."\n".$shortqual."\n"; #Print FASTQ iteration with shorted sequence and quality reads to output file.
    }
    #Reset line counter to zero and empty variables with stored data before the next iteration of the FASTQ file
    $linecount = 0;
    $id1 = ();
    $seq = ();
    $id2 = ();
    $qual = ();
  }
}

close (FASTQ); #Close input FASTQ file
close (OUT); #Close output file

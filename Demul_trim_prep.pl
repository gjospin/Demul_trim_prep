#!/usr/bin/env perl
use strict;
use warnings;
######################################
#
# Authors : Aaron Darling and Guillaume Jospin
#
#
#
#
# demultiplex_dualBC.pl <illumina_directory> <barcode_list>
#
# Illumina_directory : Directory where the illumina files are located for the dual barcode demultiplex
#                      For a MiSeq run, files should look like :
#                      XTDM1_NoIndex_L001_R1_001.fastq.gz
#                      XTDM1_NoIndex_L001_R2_001.fastq.gz
#                      XTDM1_NoIndex_L001_R3_001.fastq.gz
#                      XTDM1_NoIndex_L001_R4_001.fastq.gz
# Barcode_list : list of the barcode name and barcode.  <barcode_label><tab><barcode>
#
# Printing to summary.txt the read counts for each barcode
# Mismatched_barcode : one of the two mates had a barcode that was not recognized from the list
# BC1.BC2 : Whenever a barcode pair was not listed in the name mapping file
#
#
######################################
## Use FLASH <read1> <read2> -m 30 -M 70 -x 0.25 -p 33 -o <merged file> -r 250 -f 250 -s 25
## Change -p 33 to whatever type calls for in the script.
use IO::Zlib;
use Getopt::Long;
my ( $full_print, $reverse, $fragment_length, $skip_merge );
my $read_length       = 250;
my $quality_threshold = 20;
my $fragment_std      = 20;
my $min_overlap       = 10;
my $max_overlap       = 70;
my $mismatch_ratio    = 0.25;
my $phred_value       = 33;
GetOptions(
			"full-print"       => \$full_print,
			"reverse"          => \$reverse,
			"q=i"              => \$quality_threshold,
			"read-len=i"       => \$read_length,
			"frag-len=i"       => \$fragment_length,
			"frag-std=i"       => \$fragment_std,
			"min-overlap=i"    => \$min_overlap,
			"max-overlap=i"    => \$max_overlap,
			"mismatch-ratio=i" => \$mismatch_ratio,
			"skip-merge"       => \$skip_merge,
			"phred=i"          => $phred_value,
);
my $usage = "Wrong number of arguments\nUsage:\ndemultiplex_dualBC.pl <options> <illumina_directory> <mapping_file> <output_directory> <filename_core>\n";
die("$usage") if @ARGV != 4;

#reading barcodes
#my %barcode_forward = ();
my %barcode_rev   = ();
my $output_dir    = $ARGV[2];
my $out_file_core = $ARGV[3];
##if the $output_dir does not exists create it.
print STDERR "Creating $output_dir\n" unless -e $output_dir;
`mkdir -p $output_dir`                unless -e $output_dir;
##if the $output_dir does not exists create it.
print STDERR "Creating $output_dir/interleaved_fastq\n" unless -e "$output_dir/interleaved_fastq";
`mkdir -p $output_dir/interleaved_fastq`                unless -e "$output_dir/interleaved_fastq";
##if the $output_dir does not exists create it.
print STDERR "Creating $output_dir/raw_fastq\n" unless -e $output_dir."/raw_fastq";
`mkdir -p $output_dir/raw_fastq`                unless -e $output_dir."/raw_fastq";
##if the $output_dir does not exists create it.
print STDERR "Creating $output_dir/qiime_ready\n" unless -e $output_dir."/qiime_ready";
`mkdir -p $output_dir/qiime_ready`                unless -e $output_dir."/qiime_ready";

print STDERR "Reading barcodes and sample mapping\n";
open( INBC, $ARGV[1] );
my $header = <INBC>;
$header =~ s/^#//;
my @cols                     = split( /\s+/, $header );
my %barcode_forward          = ();
my %barcode_reverse          = ();
my %mapping                  = ();
my %output_filehandles_1     = ();
my %output_filehandles_2     = ();
my %output_filehandles_full  = ();
my %output_filehandles_qiime = ();
my %sample_read_count        = ();
my %mapping_file             = ();
my %barcode_names_forward = ();
my %barcode_names_reverse = ();
while (<INBC>) {
	chomp($_);
	my @line = split( /\t/, $_ );
	for ( my $i = 0; $i < scalar(@line); $i++ ) {
		my $key   = $cols[$i];
		my $value = $line[$i];
		$mapping_file{ $line[0] }{$key} = $value;
	}
	$barcode_names_forward{$mapping_file{ $line[0] }{"BarcodeSequence"}} = $mapping_file{$line[0]}{"BarcodeName"} unless exists $barcode_names_forward{$mapping_file{ $line[0] }{"BarcodeSequence"}} ;
	$barcode_names_reverse{$mapping_file{ $line[0] }{"ReverseBarcode"}} = $mapping_file{$line[0]}{"ReverseName"} unless exists $barcode_names_reverse{$mapping_file{ $line[0] }{"ReverseBarcode"}};
	# insert all single-error barcodes
	for ( my $i = 0; $i < length( $mapping_file{ $line[0] }{"BarcodeSequence"} ); $i++ ) {
		my @chars = ( "A", "C", "G", "T", "N" );
		my $s = $mapping_file{ $line[0] }{"BarcodeSequence"};
		foreach my $ck (@chars) {
			substr( $s, $_, 1 ) =~ s/[ACGT]/$ck/ for $i;
			print STDERR "Barcode collision! $s => ".$mapping{ $line[0] }{"BarcodeSequence"}." was already defined as $barcode_forward{$s}!!\n"
			  if defined $barcode_forward{$s}
			  && $s                   ne $mapping_file{ $line[0] }{"BarcodeSequence"}
			  && $barcode_forward{$s} ne $mapping_file{ $line[0] }{"BarcodeSequence"};
			$barcode_forward{$s} = $mapping_file{ $line[0] }{"BarcodeSequence"};
		}
	}
	for ( my $i = 0; $i < length( $mapping_file{ $line[0] }{"ReverseBarcode"} ); $i++ ) {
		my @chars = ( "A", "C", "G", "T", "N" );
		my $s = $mapping_file{ $line[0] }{"ReverseBarcode"};
		foreach my $ck (@chars) {
			substr( $s, $_, 1 ) =~ s/[ACGT]/$ck/ for $i;
			print STDERR "Barcode collision! $s => ".$mapping_file{ $line[0] }{"ReverseBarcode"}." was already defined as $barcode_reverse{$s}!!\n"
			  if defined $barcode_reverse{$s} && $s ne $mapping_file{ $line[0] }{"ReverseBarcode"};
			$barcode_reverse{$s} = $mapping_file{ $line[0] }{"ReverseBarcode"};
		}
	}
	$mapping{ $mapping_file{ $line[0] }{"BarcodeSequence"} }{ $mapping_file{ $line[0] }{"ReverseBarcode"} } = $line[0];
	$output_filehandles_full{ $line[0] } = new IO::Zlib;
	$output_filehandles_full{ $line[0] }->open( "$output_dir/interleaved_fastq/$out_file_core"."_$line[0].fastq.gz", "wb9" );
	$output_filehandles_qiime{ $line[0] } = new IO::Zlib;
	$output_filehandles_qiime{ $line[0] }->open( "$output_dir/qiime_ready/$out_file_core"."_$line[0].fastq.gz", "wb9" );
	open( $output_filehandles_1{ $line[0] }, ">$output_dir/raw_fastq/$out_file_core"."_$line[0]"."_1".".fastq" ) unless $skip_merge;
	open( $output_filehandles_2{ $line[0] }, ">$output_dir/raw_fastq/$out_file_core"."_$line[0]"."_2".".fastq" ) unless $skip_merge;

}
close(INBC);

#foreach my $fw(keys %mapping){
#	foreach my $rv (keys %{$mapping{$fw}}){
#		print "$fw\t$rv\t$mapping{$fw}{$rv}\n";
#	}
#}
## FastQ files are assumed to be gzipped.  .gz filenames
my %bc_count = ();
my @files    = <$ARGV[0]/*_R1_*.fastq.gz>;

foreach my $file (@files) {
	print STDERR "Processing $file\n";
	$file =~ m/^(\S+)_\S\S_(\d+).fastq.gz/;
	my $core  = $1;
	my $index = $2;
	open( my $TYPETEST, "zcat $file |" );
	my $type = get_sequence_input_type($TYPETEST);

	#print STDERR "TYPE :".$type->{qtype}."\n";
	close($TYPETEST);
	my $read_type = $type->{qtype};
	open( READ1,  "zcat $file |" );
	open( READ2,  "zcat $core"."_R4_$index.fastq.gz |" );
	open( INDEX1, "zcat $core"."_R2_$index.fastq.gz |" );
	open( INDEX2, "zcat $core"."_R3_$index.fastq.gz |" );
	my @read1  = ();
	my @read2  = ();
	my @index1 = ();
	my @index2 = ();
	## read the reads and indices from the input files.
	while (1) {
		for ( my $i = 0; $i < 4; $i++ ) {
			if ($reverse) {
				$read1[$i]  = <READ2>;
				$read2[$i]  = <READ1>;
				$index1[$i] = <INDEX2>;
				$index2[$i] = <INDEX1>;
			} else {
				$read1[$i]  = <READ1>;
				$read2[$i]  = <READ2>;
				$index1[$i] = <INDEX1>;
				$index2[$i] = <INDEX2>;
			}
		}
		## stop reading if we reached the end of the files
		last if !defined( $read1[0] );
		## figure out what barcodes we are dealing with
		my $BC1 = "";
		my $BC2 = "";
		my $i1  = $index1[1];
		my $i2  = $index2[1];
		chomp($i1);
		chomp($i2);
		if ( exists $barcode_forward{$i1} ) {
			$BC1 = $barcode_forward{$i1};
		}
		if ( exists $barcode_reverse{$i2} ) {
			$BC2 = $barcode_reverse{$i2};
		}
		my $code;
		##print "BC1 : $BC1\t BC2\t $BC2\n";
		if ( $BC1 eq "" || $BC2 eq "" ) {
			$bc_count{mismatch} = 1 if !exists $bc_count{mismatch};
			$bc_count{mismatch}++ if exists $bc_count{mismatch};
			$code = "mismatched_barcode";
		} else {
			#next unless exists $mapping{$BC1}{$BC2};    ##barcode combo is not part of our samples, skip.
			if(!exists $mapping{$BC1}{$BC2}){
				$mapping{$BC1}{$BC2} = $barcode_names_forward{$BC1}.".".$barcode_names_reverse{$BC2};
			}
			$code = $mapping{$BC1}{$BC2} ;
			
			$bc_count{$code} = 1 if !exists $bc_count{$code};
			$bc_count{$code}++ if exists $bc_count{$code};
			
		}
		if ( !exists $output_filehandles_1{$code} ) {
			print "new CODE : $code\n";
			$output_filehandles_full{$code} = new IO::Zlib;
			$output_filehandles_full{$code}->open( "$output_dir/$out_file_core"."_$code.fastq.gz", "wb9" );
			open( $output_filehandles_1{$code}, ">$output_dir/$out_file_core"."_$code"."_1".".fastq" ) unless $skip_merge;
			open( $output_filehandles_2{$code}, ">$output_dir/$out_file_core"."_$code"."_2".".fastq" ) unless $skip_merge;
			$output_filehandles_qiime{ $code } = new IO::Zlib;
			$output_filehandles_qiime{ $code }->open( "$output_dir/qiime_ready/$out_file_core"."_$code.fastq.gz", "wb9" );
		}
		next unless exists $mapping{$BC1}{$BC2};
		my $READ_HANDLE_1     = $output_filehandles_1{ $mapping{$BC1}{$BC2} };
		my $READ_HANDLE_2     = $output_filehandles_2{ $mapping{$BC1}{$BC2} };
		my $READ_HANDLE_FULL  = $output_filehandles_full{ $mapping{$BC1}{$BC2} };
		my $READ_HANDLE_QIIME = $output_filehandles_qiime{ $mapping{$BC1}{$BC2} };
		## quality trim the reads
		qtrim_read( read => \@read1, quality => $quality_threshold, readtype => $type );
		qtrim_read( read => \@read2, quality => $quality_threshold, readtype => $type );
		## print the trimmed reads to an intermediate file if specified in $full_print

		## merge the reads if possible
		## my @merged_read = align_and_merge_reads(read1=> \@read1, read2=> \@read2 );

		## print the merged reads to an intermediate file if specified in $full_print
		## discard the reads if it wasn't possible to merge
		$read1[0] = clean_line( line => $read1[0], num => 1 );
		$read2[0] = clean_line( line => $read2[0], num => 2 );

		## Change the read header to accomodate for barcoding.

		## print the reads to their respective files
		print $READ_HANDLE_FULL @read1 if exists $output_filehandles_full{ $mapping{$BC1}{$BC2} };
		print $READ_HANDLE_FULL @read2 if exists $output_filehandles_full{ $mapping{$BC1}{$BC2} };
		print $READ_HANDLE_1 @read1    if exists $output_filehandles_1{ $mapping{$BC1}{$BC2} };
		print $READ_HANDLE_2 @read2    if exists $output_filehandles_2{ $mapping{$BC1}{$BC2} };
		##my @qiime_read1 = convert_to_qiime_read(read_array => \@read1, sample => $mapping{$BC1}{$BC2});
		##my @qiime_read2 = convert_to_qiime_read(read_array => \@read1, sample => $mapping{$BC1}{$BC2});
		##print $READ_HANDLE_QIIME @read1 if exists $output_filehandles_qiime{ $mapping{$BC1}{$BC2} };
		##print $READ_HANDLE_QIIME @read2 if exists $output_filehandles_qiime{ $mapping{$BC1}{$BC2} };
	}
}
##close files and flush IO buffers
foreach my $handle ( keys(%output_filehandles_1) ) {
	$output_filehandles_1{$handle}->close();
	next unless  -e "$output_dir/raw_fastq/$out_file_core"."_$handle"."_1.fastq" && -e "$output_dir/raw_fastq/$out_file_core"."_$handle"."_2.fastq";
	#print "HANDLE : $handle\n";
	#print "PHRED VALUE : $phred_value\n";
	my $options   = "-m $min_overlap -M $max_overlap -p $phred_value -s $fragment_std -r $read_length -x $mismatch_ratio";
	my $flash_cmd =
	   "flash_250 $output_dir/raw_fastq/$out_file_core"."_$handle"
	  ."_1.fastq $output_dir/raw_fastq/$out_file_core"."_$handle"
	  ."_2.fastq $options -d $output_dir/qiime_ready -o $out_file_core".".$handle";

	print "RUNNING : $flash_cmd\n";
	system($flash_cmd)
	  unless -z "$output_dir/raw_fastq/$out_file_core"."_$handle"."_1.fastq" && -z "$output_dir/raw_fastq/$out_file_core"."_$handle"."_2.fastq";
	open( OUT_QIIME, ">$output_dir/qiime_ready/$out_file_core.fasta" );
	print STDERR "Reading $output_dir/qiime_ready/$out_file_core.extendedFrags.fastq\n";
	open( INMERGED,  "$output_dir/qiime_ready/$out_file_core.".$handle.".extendedFrags.fastq" );
	my @read = ();
	while (<INMERGED>) {
		$read[0] = $_;
		$read[1] = <INMERGED>;
		$read[2] = <INMERGED>;
		$read[3] = <INMERGED>;
		my @qiime_read1 = convert_to_qiime_read( read_array => \@read, sample => $handle );
		print OUT_QIIME @qiime_read1;
		@read = ();
	}

	close(INMERGED);
	close(OUT_QIIME);
}

## Write a summary of the demultiplex
open( OUT, ">summary.txt" );
foreach my $code ( sort { $bc_count{$a} cmp $bc_count{$b} } keys %bc_count ) {
	print OUT $code."\t".$bc_count{$code}."\n";
}
close(OUT);

sub convert_to_qiime_read {
	my %args   = @_;
	my @read   = @{ $args{read_array} };
	my $sample = $args{sample};
	my @return_array;
	my $read_seq  = $read[1];
	chomp($read_seq);
	$read_seq = reverse($read_seq);
	$read_seq =~ tr/ACGTacgt/TGCAtgca/;
	$sample_read_count{$sample} = 0 unless exists $sample_read_count{$sample};
	$sample_read_count{$sample}++;
	$sample =~ s/_/./g; #qiime does not like _ in the sample names.
	$return_array[0] = ">$sample"."_$sample_read_count{$sample}\n";
	$return_array[1] = $read_seq."\n";
	
	#print "@read\n@return_array";

	return @return_array;
}

sub clean_line {
	my %args = @_;
	$args{line} =~ s/ /:/g;
	chomp $args{line};
	$args{line} .= "/".$args{num}."\n";
	my @line = split( /:/, $args{line} );
	splice( @line, 7, 1 );
	$args{line} = join( ':', @line );
	return $args{line};
}

=head2 qtrim_read

trims a fastq read to a particular quality score using Heng Li's algorithm from bwa.
code based on SGA's implementation.

=cut

sub qtrim_read {
	my %args     = @_;
	my $read     = $args{read};
	my $q        = $args{quality};
	my $readtype = $args{readtype};

	$q += 33 if $readtype->{qtype} eq "phred33";
	$q += 64 if $readtype->{qtype} eq "phred64";

	# Perform a soft-clipping of the sequence by removing low quality bases from the
	# 3' end using Heng Li's algorithm from bwa

	my $seq = @$read[1];
	chomp $seq;
	my $qq = @$read[3];
	chomp $qq;
	my @qual          = split( //, $qq );
	my $endpoint      = 0;                  # not inclusive
	my $max           = 0;
	my $i             = length($seq) - 1;
	my $terminalScore = ord( $qual[$i] );

	# Only perform soft-clipping if the last base has qual less than $q
	return if ( $terminalScore >= $q );

	my $subSum = 0;
	while ( $i >= 0 ) {
		my $ps    = ord( $qual[$i] );
		my $score = $q - $ps;
		$subSum += $score;
		if ( $subSum > $max ) {
			$max      = $subSum;
			$endpoint = $i;
		}
		$i--;
	}

	# Clip the read
	@$read[1] = substr( $seq, 0, $endpoint )."\n";
	@$read[3] = substr( @$read[3], 0, $endpoint )."\n";
}

=head2 get_sequence_input_type

Checks whether input is FastA, FastQ, which quality type (33 or 64), and DNA or AA
Returns a hash reference with the values 'seqtype', 'format', and 'qtype' populated.

=cut

sub get_sequence_input_type {
	my $FILE = shift;
	my %type;
	my $counter    = 0;
	my $maxfound   = 0;
	my $dnacount   = 0;
	my $line_count = 0;
	$type{seqtype} = "dna";
	$type{format}  = "unknown";
	$type{qtype}   = "none";
	my $allcount = 0;
	my $sequence = 1;
	my $minq     = 255;    # minimum fastq quality score (for detecting phred33/phred64)
	my @lines;

	while ( my $line = <$FILE> ) {

		#print "$line\n";
		if ( $line =~ /^>/ ) {
			$maxfound = $counter > $maxfound ? $counter : $maxfound;
			$counter = 0;
			$type{format} = "fasta" if $type{format} eq "unknown";
		} elsif ( $line =~ /^@/ || $line =~ /^\+/ ) {
			$counter  = 0;
			$sequence = 1;
			$type{format} = "fastq" if $type{format} eq "unknown";

			#print STDERR "File format detected is fastq\n";
		} elsif ( $line =~ /^\+/ ) {
			$sequence = 0;
			$type{format} = "fastq" if $type{format} eq "unknown";
		} elsif ( $type{format} eq "fastq" && !$sequence ) {

			#print STDERR "Checking quality scores\n";

			# check whether qualities are phred33 or phred64
			for my $q ( split( //, $line ) ) {
				$minq = ord($q) if ord($q) < $minq;
			}
		} elsif ($sequence) {

			#$sequence =~ s/[-\.]//g;    #removing gaps from the sequences
			$line =~ s/[-\.]//g;
			$counter  += length($line) - 1;
			$dnacount += $line =~ tr/[ACGTUNacgtun]//;
			$allcount += length($line) - 1;
		}
		$line_count++;
		last if ( $line_count > 10 );
		last if ( $counter > 100000 );
	}

	$maxfound = $counter > $maxfound ? $counter : $maxfound;
	$type{seqtype} = "protein" if ( $dnacount < $allcount * 0.75 );
	$type{seqtype} = "dna"
	  if ( $type{format} eq "fastq" );    # nobody using protein fastq (yet)
	$type{qtype} = "phred64" if $minq < 255;
	$type{qtype} = "phred33" if $minq < 64;
	$type{paired} = 0;                    # TODO: detect interleaved read pairing
	return \%type;
}

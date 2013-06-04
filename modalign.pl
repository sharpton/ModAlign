#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::AlignIO;
use IPC::System::Simple qw(capture $EXITVAL);


my ( $in_seqs, $ffdb );
my $format = 'fasta';    #input sequence format
my $clust_thresh = 95;  
my $cover_thresh = 0.9;
my $is_P         = "T";
my $aln_format   = 'stockholm'; #output alignment format

GetOptions(
    "i=s"  => \$in_seqs,
    "d=s"  => \$ffdb,
    "f:s"  => \$format,
    "c:i"  => \$clust_thresh,
    "l:s"  => \$cover_thresh,
    "p:s"  => \$is_P,         #either T or F
    "of:s" => \$aln_format,   #fasta or stockholm (default)
    );

unless( $in_seqs ){
    die "You must specify an input sequence file using the -i option!\n";
}
unless( $ffdb && -d $ffdb ){
    die "You must specify an existing directory location in which I can dump output files!\n";
}

#MAIN
#Start the work here
print "Starting a new modalign run:\n";
print "modalign.pl -i $in_seqs -d $ffdb -f $format -c $clust_thresh -l $cover_thresh -p $is_P -of $aln_format\n";
my $cluster_file  = $ffdb . "clusters.tab";
run_blastclust( $in_seqs, $cluster_file, $is_P, $clust_thresh, $cover_thresh );
my $reps_file = $ffdb . "reps.fa";
select_reps( $in_seqs, $cluster_file, $reps_file, $format );
my $align_out = $ffdb . "reps.mfa";
run_muscle( $reps_file, $align_out );
my $hmm_file  = $ffdb . "rep_model.hmm";
run_hmmbuild( $align_out, $hmm_file, $is_P );
my $all_align = $ffdb . "aligned_seqs.stk";
run_hmmalign( $in_seqs, $hmm_file, $all_align );

if( $aln_format eq 'fasta' ){
    my $fasta_all_align = $ffdb . "aligned_seqs.mfa";
    stockholm_to_fasta( $all_align, $fasta_all_align );
}


#SUBROUTINES
sub run_blastclust{
    my( $in_seqs, $cluster_file, $is_P, $clust_thresh, $cover_thresh ) = @_;
    my @args = ( "-i $in_seqs", "-o $cluster_file", "-p $is_P", "-L $cover_thresh", "-S $clust_thresh", "-b T" );
    print( "Running blastclust:\n" );
    print( "blastclust @args\n" );
    my $results = capture( "blastclust @args" );
    if( $EXITVAL != 0 ){
        warn("Error running blastclust on $in_seqs: $results, $EXITVAL\n");
	exit(0);
    }
    return $results;    
}

sub select_reps{
    my( $in_seqs, $cluster_file, $reps_file, $format ) = @_;
    print( "Selecting representatives\n");
    my $seqs = Bio::SeqIO->new( -file => $in_seqs, -format => $format );
    my $out  = Bio::SeqIO->new( -file => ">$reps_file", -format => $format );
    open( CLUSTS, $cluster_file ) || die "Can't open $cluster_file for read: $!\n";
    CLUST: while(<CLUSTS>){
	chomp $_;
	my @ids = split( ' ', $_ );
	my $id = $ids[0];
	while( my $seq = $seqs->next_seq() ){
	    my $seqid = $seq->display_id();
	    if( $seqid eq $id ){
		$out->write_seq( $seq );
		next CLUST;
	    }
	}
    }
    close CLUSTS;
    return 1;
}

sub run_muscle{
    my( $inseqs, $output ) = @_;
    my @args = ( "-in $inseqs", "-out $output" );
    print( "Running muscle:\n" );
    print( "muscle @args\n" );
    my $results = capture( "muscle @args" );
    if( $EXITVAL != 0 ){
        warn("Error running muscle on $inseqs: $results\n");
	exit(0);
    }
    return $results;    
}

sub run_hmmbuild{
    my( $in_aln, $out_hmm, $is_P ) = @_;
    my @args = ();
    if( $is_P eq "T" ){
	@args = ( "$out_hmm", "$in_aln" );
    }
    else{
	@args = ( "--informat afa", "--dna", "$out_hmm", "$in_aln" );
    }
    print( "Running hmmbuild:\n" );
    print( "hmmbuild @args\n" );
    my $results  = capture( "hmmbuild @args" );
    if( $EXITVAL != 0 ){
        warn("Error running hmmbuild on $in_aln: $results\n");
	exit(0);
    }
    return $results;
}

sub run_hmmalign{
    my ( $inseqs, $hmm, $output ) = @_;
    my @args = ( "-o $output", "$hmm", "$inseqs" );
    print( "Running hmmalign:\n" );
    print( "hmmalign @args\n" );
    my $results = capture( "hmmalign @args" );
    if( $EXITVAL != 0 ){
	warn( "Error running hmmalign on $inseqs: $results\n" );
	exit(0);
    }
    return $results;
}

sub stockholm_to_fasta{
    my ( $all_align, $fasta_all_align ) = @_;
    my $in  = Bio::AlignIO->new( -file => $all_align,          -format => 'stockholm' );
    my $out = Bio::AlignIO->new( -file => ">$fasta_all_align", -format => 'fasta' );
    while( my $aln = $in->next_aln() ){
	$out->write_aln( $aln );
    }
}

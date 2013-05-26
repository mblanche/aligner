use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path remove_tree);

use Getopt::Long;
use Pod::Usage;

our ($man,$help,$debug,@prog_args,$Force);
our (@dirs,$ebwt,$tx_idx);
our ($cpu);

## Careful... This is not a great idea... 
our (@excluded); 
our $aligner = 'tophat';

MAIN:{
  init();
  testIndex();
  my %res = getSample2Files();
  qsub_align(%res);
  exit 0;
}

## delete simlink if exists
# unlink "$bam_dir/${exp}.bam" if (-l "$bam_dir/${exp}.bam");
# symlink("$ENV{PWD}/$aln_dir/$exp/accepted_hits.bam","$bam_dir/${exp}.bam");
# print STDERR "Indexing the bam file $bam_dir/${exp}.bam\n";
# system("samtools index $bam_dir/${exp}.bam");

sub testIndex{
  my $ori_tx_idx = $tx_idx;
  my $ori_ebwt = $ebwt;

  ## Expand tilde to full path, 
  $tx_idx =~ s/~/$ENV{HOME}/;
  $ebwt =~ s/~/$ENV{HOME}/;
  
  unless (-e "$tx_idx.gff"){
    if (exists $ENV{TOPHAT_TX_INDEXES}){
      $ENV{TOPHAT_TX_INDEXES} =~ s/~/$ENV{HOME}/;
      $tx_idx = -e "$ENV{TOPHAT_TX_INDEXES}/$tx_idx.gff" ? 
	"$ENV{TOPHAT_TX_INDEXES}/$tx_idx" : undef;
    } else {
      $tx_idx = undef;
    }
  }
  
  unless (-e "$ebwt.1.bt2"){
    if (exists $ENV{BOWTIE2_INDEXES}){
      $ENV{BOWTIE2_INDEXES} =~ s/~/$ENV{HOME}/;
      $ebwt = -e "$ENV{BOWTIE2_INDEXES}/$ebwt.1.bt2" ? 
	"$ENV{BOWTIE2_INDEXES}/$ebwt" : undef;
    } else {
      $ebwt = undef;
    }
  }
  
  die "Can't find $ori_tx_idx\n" unless $tx_idx;
  die "Can't find $ori_ebwt\n" unless $ebwt;
  
}

sub getSample2Files {
  my %res2;
  
  for my $dir (@dirs){
    my ($csv,@rest) = glob "$dir/*.csv";
    die " there is more than one CSV file in $dir->[0]" if @rest;
    
    open FH, $csv;
    my($filename, $dir, $suffix) = fileparse($csv);
    my %res;
    while (<FH>){
      next if $. == 1;
    chomp;
      my @line = split /,/;
      $res{$line[0]} = $line[2];
    }
    my @fastq = glob "${dir}*[ACGT].fastq.gz";
    
    my %excluded = map{$_,1} @excluded if @excluded;
    
    @fastq = grep{!exists $excluded{$_}} @fastq;
    
    for my $path (@fastq){
      my($file) = fileparse($path);
      
      if (!exists $res{$file}){
        print STDERR "$file is not associated with a sample name\n";
      } else {
        push @{$res2{$res{$file}}},$path
      }
    }
  }
  return %res2
}

sub qsub_align{
  my %res = @_;

  my $aln_dir = "TopHat_aln";
  my $sge_out = "SGE_out";
  make_path($aln_dir,$sge_out);
  
  for my $exp (keys %res) {
    print STDERR "Aligning $exp\n";
    
    my $align = "tophat";
    $align   .= " -p$cpu";
    $align   .= " --transcriptome-index=$tx_idx";
    $align   .= " -o $aln_dir/$exp";
    $align   .= " $ebwt";
    my $fastq = join(",",@{$res{$exp}});
        
    my $com = "qsub -j y -o $sge_out -q test-mar22 -N _$exp -V -cwd -b y '$align $fastq'";
    print "$com\n\n";
    
    next if $debug;
    my $res = system $com;
    print $res,"\n";
  }
}

sub init {
  my $filename = basename($0);
  @prog_args = ($filename,@ARGV);

  GetOptions('debug'         => \$debug,
	     'help|?'        => \$help,
	     'man'           => \$man,
	     
	     'directories=s{,}'      => \@dirs,
	     'bowtie_index=s'        => \$ebwt,
	     'transcriptome-index=s' => \$tx_idx,
	     
	     'excluded=s{,}'            => \@excluded,
	     'cpu=i'     => \$cpu,

	    ) or pod2usage(1);

  if ($help || !@dirs){
    pod2usage(-exitstatus => 0, -verbose => 2);
  }elsif ($man){
    pod2usage(-exitstatus => 0, -verbose => 0);
  }

  $cpu ||= 8;
  $tx_idx ||=  "$ENV{HOME}/Genomes/bowtie_indexes/tophat_tx/Drosophila_melanogaster.BDGP5.71.min.tc";
  $ebwt ||= "$ENV{HOME}/Genomes/bowtie_indexes/bwt2/Drosophila_melanogaster.BDGP5.71.min";

  print STDERR "##########RUNNING IN DEBUGING MODE##########\n" if $debug;
}


__END__

example:

perl alignSeq_SGS.pl --dir test --bowtie Drosophila_melanogaster.BDGP5.71.min  --trans Drosophila_melanogaster.BDGP5.71.min.tc

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path remove_tree);

use threads;
use threads::shared;
use Thread::Queue;

use Getopt::Long;
use Pod::Usage;

our ($man,$help,$debug,@prog_args,$Force);
our (@dirs,$ebwt,$tx_idx);
our ($cpu);

## Careful... This is not a great idea... 
our (@excluded); 
our $aligner = 'tophat';


init();
print STDERR "##########RUNNING IN DEBUGING MODE##########\n" if $debug;

#syncGenomes();
my %res = getSample2Files();
qsub_aling(%res);

my $aln_dir = "TopHat_aln";
my $bam_dir = "bam";
make_path($aln_dir,$bam_dir);






 
## delete simlink if exists
# unlink "$bam_dir/${exp}.bam" if (-l "$bam_dir/${exp}.bam");
# symlink("$ENV{PWD}/$aln_dir/$exp/accepted_hits.bam","$bam_dir/${exp}.bam");
# print STDERR "Indexing the bam file $bam_dir/${exp}.bam\n";
# system("samtools index $bam_dir/${exp}.bam");

sub syncGenomes {

  my $src = '/Users/lab_project/Genomes/bowtie_indexes';
  my $dest = 'lepus:Genomes';

  my $res = system "rsync -va $src $dest";
  die "Couldn't sync the indexes\n" unless $res == 0;
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
  my $aligner_param;
  $aligner_param .= ' -p8';
  #$aligner_param .= ' -g1';
  $aligner_param .= ' --transcriptome-index='.$tx_idx;
  
  
  for my $exp (keys %res) {
    
    print STDERR "Aligning $exp\n";
    my $com = "$aligner $aligner_param -o $exp $ebwt ".join(",",@{$res{$exp}});
    print STDERR "$com\n\n";
    next if $debug;
    #system $com;
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

  #$db   ||= 'd_melanogaster_fb5_22';
  #$tdf_db ||= 'dmel5_22';
  $cpu ||= 8;
  
}


__END__

example:

perl alignSeq_SGS.pl --dir test --bowtie Drosophila_melanogaster.BDGP5.71.min --trans Drosophila_melanogaster.BDGP5.71.min.tc

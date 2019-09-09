#!/usr/bin/perl -w
#
# take file and blastx it against a protein db and form an alignment
# Assumption: file contains ortholog sequences of a gene

use File::Basename;
use File::Temp;
$tempdir = File::Temp::tempdir( CLEANUP => 1);

use XML::LibXML;
my $parser = XML::LibXML->new(); 

use Getopt::Long;
&Getopt::Long::Configure qw(pass_through);

my $db = "";
my $query = "";
my $name = "";
my $special = "";
my $skeleton = "";
my $seq = {};
my $rhit = {};
my @alignment = ();
my @gap_pos = ();
my $cutoff = 50;
my $frame = 1;
my $identity = 65;
my $doframe = 1;
my $selano = 0;

if( $ARGV[0] eq "" || $ARGV[0] eq "-h" ){
  print "Usage: $0 -query <fasta file of orthologs> -db <reference db> -special <special orgname> -frame <frame no> -identity <%identity in blast> -cutoff <% gap in cleaning>\n";
  exit;
}

GetOptions (
  'db=s' => \$db,
  'query=s' => \$query,
  'identity=f' => \$identity,
  'special=s' => \$special,
  'cutoff=f' => \$cutoff,
  'doframe!' => \$doframe,
  'frame=i' => \$frame,
  'selano!' => \$selano,
);

#------------------------------------------------------------------------------------------

# convert fasta to hash
open (Q, $query) || die "Cannot open $query: $!\n";
$seq = fasta2hash(<Q>);
close Q;

# replace any known gaps with N (as blast ignores gaps)
open T, ">$tempdir/fasta";
my $seqN = "";
my $length = 0;
foreach (keys %$seq) {
  $seqN = $seq->{$_}; 
  $seqN =~ s/-/N/g ;
  $seq->{$_} = $seqN;
  print T ">$_\n";
  print T "$seqN\n";
}
close T;

# format the db
system "cp $db $tempdir";
my ($filename) = fileparse($db);
system "formatdb -i $tempdir/$filename";

# run blastx with the filter off -b parameter to restrict no. of hits to 1 doesnt apply to xml format
system "blastall -p blastx -i $tempdir/fasta -d $tempdir/$filename -e 0.00001 -m 7 -F F  > $tempdir/blastx_xml";
my $doc = $parser->parse_file("$tempdir/blastx_xml");

# findnodes returns an array, findvalue returns a string (even when there are multiple values)
# Predicates: Hit[1] get the first hit for the node, Hsp[1] gets the first Hsp for the Hit
my @nodes = $doc->findnodes('/BlastOutput/BlastOutput_iterations/Iteration');
foreach $node (@nodes){
  $name = $node->findvalue('./Iteration_query-def');
  @child =  $node->nonBlankChildNodes;
  foreach $child (@child){
    next if( $child->nodeName ne "Iteration_hits" );
    $rhit->{$name}->{subject} = $child->findvalue('./Hit[1]/Hit_def');
    @hsps = $child->findnodes('./Hit[1]//Hsp');
    foreach my $hsp (@hsps){
      if( $doframe ){
        if( $hsp->findvalue('./Hsp_query-frame') == $frame ){
          if( $hsp->findvalue('./Hsp_identity')/$hsp->findvalue('./Hsp_align-len') >= $identity/100 ){
            push @{$rhit->{$name}->{qstart}}, $hsp->findvalue('./Hsp_query-from');
            push @{$rhit->{$name}->{qend}}, $hsp->findvalue('./Hsp_query-to');
            push @{$rhit->{$name}->{qframe}}, $hsp->findvalue('./Hsp_query-frame');
            push @{$rhit->{$name}->{qseq}}, $hsp->findvalue('./Hsp_qseq');
            push @{$rhit->{$name}->{hseq}}, $hsp->findvalue('./Hsp_hseq');
            push @{$rhit->{$name}->{hstart}}, $hsp->findvalue('./Hsp_hit-from');
            push @{$rhit->{$name}->{hend}}, $hsp->findvalue('./Hsp_hit-to');
            push @hit_end, $hsp->findvalue('./Hsp_hit-to');
            push @{$rhit->{$name}->{align_length}}, $hsp->findvalue('./Hsp_align-len');
            push @{$rhit->{$name}->{identity}}, $hsp->findvalue('./Hsp_identity');
            push @{$rhit->{$name}->{positive}}, $hsp->findvalue('./Hsp_positive');
          }
        }
      } else {
        if( $hsp->findvalue('./Hsp_identity')/$hsp->findvalue('./Hsp_align-len') >= $identity/100 ){
          push @{$rhit->{$name}->{qstart}}, $hsp->findvalue('./Hsp_query-from');
          push @{$rhit->{$name}->{qend}}, $hsp->findvalue('./Hsp_query-to');
          push @{$rhit->{$name}->{qframe}}, $hsp->findvalue('./Hsp_query-frame');
          push @{$rhit->{$name}->{qseq}}, $hsp->findvalue('./Hsp_qseq');
          push @{$rhit->{$name}->{hseq}}, $hsp->findvalue('./Hsp_hseq');
          push @{$rhit->{$name}->{hstart}}, $hsp->findvalue('./Hsp_hit-from');
          push @{$rhit->{$name}->{hend}}, $hsp->findvalue('./Hsp_hit-to');
          push @hit_end, $hsp->findvalue('./Hsp_hit-to');
          push @{$rhit->{$name}->{align_length}}, $hsp->findvalue('./Hsp_align-len');
          push @{$rhit->{$name}->{identity}}, $hsp->findvalue('./Hsp_identity');
          push @{$rhit->{$name}->{positive}}, $hsp->findvalue('./Hsp_positive');
        } 
      }
    }
  }
}

if( scalar(@hit_end) ){
  # get the earliest start position of hit protein among all hits
  #push @hit_start, $_->textContent foreach($doc->findnodes('//Hit[1]//Hsp[1]/Hsp_hit-from'));
  #@hit_start = sort {$a <=> $b} @hit_start;
  #$hit_start = $hit_start[0]; # the earliest start position of the hit protein
  
  # get the greatest end position of hit protein among all hits
  @hit_end = sort {$b <=> $a} @hit_end;
  $hit_end = $hit_end[0]; # the last position of the hit protein
  
  # make the skeleton sequence to contruct the alignment on
  $i = 0;
  $skeleton = "";
  while( $i< $hit_end*3 ){
    $skeleton .= "-";
    $i++; 
  }
} else {
  print "Alignment is empty!\n";
  exit;
}

$process = 0;
if( $special ne ""){
  foreach (keys %$rhit){
    if( $_ =~ /$special/i ){
      if( $rhit->{$_}->{positive}->[0]/$rhit->{$_}->{align_length}->[0] == 1 || 
          $rhit->{$_}->{identity}->[0]/$rhit->{$_}->{align_length}->[0] == 1 ){
        $process = 1;
      } else {
        $process = 0;
      }
    }
  }
} else {
  $process = 1;
} 

if( $process ){
SEQ:foreach $seq_name (keys %$rhit){
    $rhit->{$seq_name}->{dna} = $skeleton;
  
    foreach my $num (0..$#{$rhit->{$seq_name}->{qstart}}){
      my $dna = substr($seq->{$seq_name}, $rhit->{$seq_name}->{qstart}->[$num] - 1, $rhit->{$seq_name}->{qend}->[$num]-$rhit->{$seq_name}->{qstart}->[$num]+1);

      # insert gaps into DNA substring, for - and X
      my @gap_index = ();
      my @aa = split( //, $rhit->{$seq_name}->{qseq}->[$num] );
      foreach my $gap (0..$#aa){
        if( $aa[$gap] eq "-" ){
          substr($dna, $gap*3, 0, "---"); 
        } elsif ( $aa[$gap] eq "X" ){
          substr($dna, $gap*3, 3, "---"); 
        } elsif( $aa[$gap] eq "*" ){
          if( $selano ){
          } else {
            next SEQ;
          }
          #substr($dna, $gap*3, 3, "---");
        }
      }
      @aa = split( //, $rhit->{$seq_name}->{hseq}->[$num] );
      foreach my $gap (reverse 0..$#aa){
        if( $aa[$gap] eq "-" ){
          substr($dna, $gap*3, 3, ""); 
        }
      }
      substr($rhit->{$seq_name}->{dna}, ($rhit->{$seq_name}->{hstart}->[$num]- 1)*3, ($rhit->{$seq_name}->{hend}->[$num] - $rhit->{$seq_name}->{hstart}->[$num]+1)*3 , $dna);
    }
  
    if( $special ne "") { 
     if( length($rhit->{$seq_name}->{dna}) == length($skeleton) ){
       push @alignment, ">$seq_name";
       push @alignment, $rhit->{$seq_name}->{dna};
     }
    } else {
      push @alignment, ">$seq_name";
      push @alignment, $rhit->{$seq_name}->{dna};
    }
  }
} else {
  print "$special sequence does not match 100% to database $special sequence\n";
  system "rm -rf $tempdir";
  exit;
}

# trim alignment - if more than 90% of the alignment is a gap at position x
# then delete the codon of position x from whole alignment
# Also, if more than 90% of a sequence is gaps, remove the sequence from the alignment.
my %count;
my $count;
my $removed = 0;
if( scalar(@alignment) ){

  # remove sequences with more than 90% gaps in the entire alignment
#  for(my $i=$#alignment; $i > 0; $i -= 2){ # reverse index so it doesnt affect array index when removing elements
#    $count = 0;
#    @chars = split //, $alignment[$i];
#    foreach $j (0..$#chars){
#      if( $chars[$j] eq "-" ){
#        $count++; # no. of gaps
#        if( $count{$j} ){ # freq of gaps at position j
#          $count{$j}++;
#        } else {
#          $count{$j} = 1;
#        } 
#      }
#    }
#    if( $count >= ($cutoff*0.01*(length($alignment[$i]))) ){
#      splice( @alignment, $i-1, 2 );
#      $removed++;
#    }
#  }
#  # remove positions with more than 90% gaps in the entire alignment
#  foreach $seq (@alignment){
#    if( $seq =~ /^>/ ){} else{ 
#      @pos = sort{ $b <=> $a } keys %count;
#      for(my $i = 0; $i < scalar(@pos); $i+=3){ # reverse index - so it doesn't affect the seq index
##print $count{$pos[$i]}, " $pos[$i]\n", $count{$pos[$i+1]}," $pos[$i+1]\n", $count{$pos[$i+2]}," $pos[$i+2]\n";
#        if( $count{$pos[$i]}-$removed >= ($cutoff*0.01*(scalar(@alignment)/2)) ||
#            $count{$pos[$i+1]}-$removed >= ($cutoff*0.01*(scalar(@alignment)/2)) ||
#            $count{$pos[$i+2]}-$removed >= ($cutoff*0.01*(scalar(@alignment)/2)) ){ 
#          substr($seq, $pos[$i+2], 3) = ""; 
#        }
#      }
#    }
#  }
#
  my $alignment = join ("\n", @alignment);
  $ali = fasta2hash($alignment);
  foreach (sort keys %$ali){
    print ">$_\n";
    print $ali->{$_}, "\n";
  }
} else {
  print "Alignment is empty!\n"
}

system "rm -rf $tempdir"

sub fasta2hash {
  # try to chomp newlines so we can take any format
  # try to be smart about strings vs. arrays
  my @a = @_;
  my @b = ();
  my $seq;
  my $name;
  my $a;
  my $i;
  my $test;
  foreach $a (@a) {
    push @b, split (/\n/, $a);
  }
  foreach $a (@b) {
    next if $a =~ /^#/;
    chomp $a;
    if ($a =~ s/^>//) {
      $i = 1;
      $name = $a;
      $test = $name;
      while (defined $seq->{$test}) {
        $test = $name . $i;
        $i++;
      }
      $name = $test;
      $seq->{$name} = "";
    } elsif (defined $name) {
      $seq->{$name} .= $a;
#      undef $name;
    }
  }
  return $seq;
}


#!/usr/bin/env perl

my $version="1.0.0";
my $author="Sebastian Mackowiak";
my $last_change="15.01.2018";

use strict;
use warnings;

use Getopt::Std;
use Term::ANSIColor;

my %options;
getopts("a:A:m:uf:F:bo:B:G:c:sO:Sdp:CeMqQI:t:",\%options);

my $R="\x1b[31m";
my $G="\x1b[32m";
my $Y="\x1b[33m";
my $B="\x1b[34m";
my $M="\x1b[35m";
my $C="\x1b[36m";
my $E="\x1b[0m";

sub col{
    return "$_[0]$_[1]$E";
}

my $usage="
Usage:
${B}exonize_template.pl$C -m$E mapping_file$C -a$E annotation_file [options]

[options]
".col($G,'-f')." [char] -> feature to use if gff file given, eg, exon, 5utr,3utr, intron ...
".col($G,'-F')." [char] -> feature to ignore if gff file given as option -A, eg, exon, 5utr,3utr, intron ...
".col($G,'-u')."        -> use only unique mapping reads 
".col($G,'-C')."        -> if your mapping file does not have the NH tag then supply this to make a file with id counts (attention, you must have output all alignments by the aligner)
".col($G,'-b')."        -> use only mapping positions with exactly one gene annotation  
".col($G,'-o')." [char] -> output file name, otherwise default file will be used 
".col($G,'-B')." [char] -> output mapped reads of id=options{'B'} to file obam and obam2
".col($G,'-G')."        -> file that maps refseq ids to gene_symbols 
".col($G,'-c')." [char] -> chromosome to process, omit if all chromosomes shall processed
".col($G,'-s')."        -> do not sort input file cause it is already sorted
".col($G,'-S')."        -> do strand specific assignments
".col($G,'-O')." [char] -> file with feature precedence, eg exon > lincRNA > intron ... ## order can vary  
".col($G,'-d')."        -> read multiplicity is ignored, so if id has _x100 it will be counted as 1 only  
".col($G,'-p')." [char] -> print reads mapped to category other to file given as option 'p'
".col($G,'-e')."        -> output file contains also per sample counts given by first tag at fasta ID separated by '_' 
".col($G,'-q')."        -> count number of read mapping positions in gene instead of read abundance purely
".col($G,'-Q')."        -> output startpos and read abundance per gene instead of mapping positions
".col($G,'-I')." [chars]-> list of samples to be ignored, comma separated 
".col($G,'-t')." [ncT]  -> use only pool files which have the particular indicator- n:nonclipped,c:clip2,0 for clip1

";

if($options{'p'}){
    open OTHER,">$options{'p'}" or die "Could not open file $options{'p'}\n";
}

my %order=();
my %features=();
if($options{'O'}){
    open IN,$options{'O'} or die "File with order not found\n";
    my $count=1000;
    while(<IN>){
        chomp;
        if(/^(\d+)\s+(\S+)/){
            $order{$1}=$2;
            $features{$2}=1;
        }else{
            $count++;
            $order{$count}=$_;
            $features{$_}=1;
        }
    }
    close IN;
}

my %tags=();

my $gchr;

$gchr=$options{'c'} if($options{'c'});

my %a2b=();
if($options{'G'}){
    open IN,$options{'G'} or die "No file with remappings given\n$usage";
    while(<IN>){
        chomp;
        my @l=split();
        if($a2b{$l[0]}){
            if($a2b{$l[0]} =~ /Mir/ and $l[1] !~ /Mir/){
                #print STDERR "$l[0]\t$a2b{$l[0]}\t$l[1]     next \n";
                next;
            }elsif($a2b{$l[0]} !~ /Mir/ and $l[1] =~ /Mir/){
                #print STDERR "$l[0]\t$a2b{$l[0]}\t$l[1]\n";
            }
        }
        $a2b{$l[0]}=$l[1];
    }
    close IN;
}



if(not $options{'a'}){ die "No annotation file given by option $B a $E \n$usage";}
my $anno=$options{'a'};

my $anno_ex;


if(not $options{'m'}){ die "No mapping file given by option $B m $E \n$usage";}
my $map=$options{'m'};

my %h=();
my %counts=();
my %countsq=();
my %ignore=();
my $uniq=0;
my $ambig=1;

if($options{'I'}){
    %ignore=map{$_ => 1} split(",",$options{'I'});
}

my $feature='exon';
$feature=$options{'f'} if($options{'f'});

my $feature2='exon';
$feature2=$options{'F'} if($options{'F'});

if(defined $options{'u'}){
    $uniq=1;
}

if(defined $options{'b'}){
    $ambig=0;
}

## window must be the size of the longest region we look at. eg longest transcript, intron, exon etc depending on the feature we check
my $window=100;

if(defined $options{'w'}){
    $window=$options{'w'};
}

######## options output

my @opts=qw(a A m u f F b o B G M);

print STDERR "\n\noptions used\n\ta: Annotation file                         = $B $options{'a'} $E\n";
my $ig='';
if($options{'A'}){ $ig=$options{'A'};}
print STDERR "\tA: Annotation file with features to ignore = $B $ig $E\n" if($options{'A'});
print STDERR "\tm: read mapping file                       = $B $options{'m'} $E\n";
print STDERR "\tu: only uniq mappings                      = $B $uniq $E\n";
print STDERR "\tf: feature to use in GTF file              = $B $feature $E\n";
print STDERR "\tF: feature2 to use in GTF file             = $B $feature2 $E\n" if($options{'A'});
print STDERR "\tb: allow multiple genes at the same locus  = $B $ambig $E\n";
print STDERR "\ts: input file sorted already               = $B 1 $E\n" if($options{'s'});
print STDERR "\tS: reads are strand specific               = $B 1 $E\n" if($options{'S'});
print STDERR "\tO: preference file for genomic features    = $B $options{'O'} $E\n" if($options{'O'});
print STDERR "\tM: ignore genes with Gm[INT] id            = $B $options{'M'} $E\n" if($options{'M'});
print STDERR "\te: use sample tags                         = $B $options{'e'} $E\n" if($options{'e'});


my $of='';
if($options{'o'}){$of=$options{'o'};}

print STDERR "\to: output file                             = $B $of $E\n";

my $or=0;
if(exists $options{'B'}){ $or=1;$or.=" $options{'B'}" if(defined $options{'B'});}

print STDERR "\tB: output mapped reads                     = $B $or$E\n";

my $rfg=0;
if($options{'G'}){$rfg=1;}


print STDERR "\tG: refseq 2 geneSymbol file                = $B $rfg $E\n\n\n";

## we can read in multiple annotation files in GTF or BED with information of skipping features ...
my %h2=();
if($options{'A'}){
    my @l=split(",",$options{'A'});
    foreach my $f(@l){
        $anno_ex=$f;
        read_anno($anno_ex,\%h2,$feature2,1);
    }
}

my @l=split(",",$anno);
foreach my $f(@l){
    read_anno($f,\%h,$feature,0);
    print STDERR "Reading annotation file $f done\n";
}

sub read_anno{
    my ($anno,$h,$feature,$skip)=@_;


    my @l=();
    my $id;

    my $bed=0;
    if($anno =~ /.bed$/){
        $bed=1;
    }elsif($anno =~ /.g[tf]f[123]*$/){
    }else{
        die "Annotation file ending is not gtf/gff or bed\nMake sure you use one of these files\nOther files are not supported"
    }

    open IN,$anno or die "No annotation file found\n";
    my $cnt=0;
    my $beg=0;
    my $end=0;
    my $lid;
    my $lend;
    my @dmp=split(",",$feature);
    my %fs;

    while(<IN>){
        @l=split();
        next if($gchr and $gchr ne $l[0]);
        if(not $bed){
            if($feature eq 'all'){
                $id = $l[2];	
                if($id =~ /(\S+)_\d+/){$id=$1;}
                ## if we gave a feature order list then ignore all features that are not in the list
                next if($options{'O'} and !$features{$id});

                $counts{$id}{'global'}{'t'}=0 if(not $skip);
                $cnt++;

                # print STDERR $cnt,"\n" if($cnt % 10000 == 0);
                $beg=$l[3];
                $end=$l[4];

                push(@{$$h{$l[0]}{$beg}{$end}{'f'}},$id);	
                #experimental, hash it
                #$$h{$l[0]}{$beg}{$end}{'f'}{$id}++;
                push(@{$$h{$l[0]}{$beg}{$end}{'s'}},$l[6]);       ## this is the strand information we want to use	

                #die "@{$$h{$l[0]}{$beg}{$end}{'f'}}\t$$h{$l[0]}{$beg}{$end}{'s'}\t$_\n";

                if($end-$beg > $window){
                    $window=$end-$beg+1;
                    $lid=$id;
                    $lend=$end;
                }


            }else{
                #next if(!$fs{$l[2]});
                next if($l[1] ne $feature);
                if(/transcript_id "([()A-Z\.a-z-0-9_]+)"/){
                    $id=$1;
                    if($id =~ /(\S+)_dup\d+/){ $id = $1;}
                    if($options{'G'} and $a2b{$id}){ $id=$a2b{$id}; }
                    next if($id =~ /^Gm\d+$/i and $options{'M'});
                    $counts{$id}{'global'}{'t'}=0 if(not $skip);
                    $cnt++;
                    #			print STDERR $cnt,"\n" if($cnt % 10000 == 0);
                    $beg=$l[3];
                    $end=$l[4];
                    push(@{$$h{$l[0]}{$beg}{$end}{'f'}},$id);	
                    $$h{$l[0]}{$beg}{$end}{'s'}=$l[6];       ## this is the strand information we want to use	
                    if($end-$beg > $window){
                        $window=$end-$beg+1;
                        $lid=$id;
                        $lend=$end;
                    }
                }else{
                    die "no id found $_\n";
                }
            }
        }else{
            $id=$l[3];
            if($l[3] =~ /(\S+)_dup\d+/){ $id = $1;}
            if($options{'G'} and $a2b{$id}){ $id=$a2b{$id}; }
            next if($id =~ /^Gm\d+$/i and $options{'M'});
            $counts{$id}{'global'}{'t'}=0 if(not $skip);
            $cnt++;
            #print STDERR $cnt,"\n" if($cnt % 10000 == 0);
            $beg=$l[1]+1; ## since bed is 0 based we add plus one cause the sam file we are reading is one based!!!
            $end=$l[2];
            push(@{$$h{$l[0]}{$beg}{$end}{'f'}},$id);	
            $$h{$l[0]}{$beg}{$end}{'s'}=$l[5];       ## this is the strand information we want to use	
            if($end-$beg > $window){
                $window=$end-$beg+1;
                $lid=$id;
                $lend=$end;
            }
        }
    }
    print STDERR $cnt," entries read\n";
    if(not $cnt){ print STDERR "no features read\nplease define a feature from your file to be taken into account\nthe default is set to exon\nto change it either say -f all or -f intron or -f miRNA etc ...\n";exit;}
    close IN;
}


my %nh;
my $mfiles=0;

my $map_fh='';

if($map =~ /.sam/){
#open (MAP, "cat $map |sort -k3,3 -k4,4n |") or die "can not open $map\n";
    $map_fh=$map;
    make_NH_tags($map,\%nh) if($options{'C'});		
}elsif($map =~ /.bam$/){
#open (MAP, "samtools view $map|");#|sort -k3,3 -k4,4n |") or die "can not open $map file or no samtools installed\n";
    $map_fh=$map;
    make_NH_tags($map,\%nh) if($options{'C'});		
}else{
    print STDERR "Processing mapping file $map now\n";
    open T,"$map" or die "$map file could not be opened\n";
    $mfiles=1;
}

##### if our config file contains a bunch of files then we process them file wise
if($mfiles){
    ## this routine is different from the original since we fuse all input files and output only 1 file!!!
    my $ofile=$options{'o'};
    if($options{'t'}){
        $ofile.=$options{'t'};
    }

    open OUT,">$ofile" or die "Could not create output file $ofile\n";
    while(<T>){
        chomp;
        my $mapfile=$_;
        if($options{'t'}){
            if($options{'t'} eq 'T'){
                next if($mapfile =~ /pool\d[nc]_vs/);
            }else{
                next if($mapfile !~ /pool\d$options{'t'}_vs/);
            }
        }

        make_NH_tags($mapfile,\%nh) if($options{'C'});		
        print STDERR "Processing mapping file $mapfile now X\n";
        my $outh= \*OUT;
        counting($mapfile,$outh);
        
    }
    printcounts(\%counts,\%countsq,\%tags,\*OUT);
    close OUT;
}else{
    my $outh= \*STDOUT;
    if($options{'o'}){
        my $ofile=$options{'o'};
        open OUT,">$ofile" or die "Could not create outfile $ofile\n";
        $outh=\*OUT;
    }
    counting($map_fh,$outh);
    printcounts(\%counts,\%countsq,\%tags,$outh);
    #my ($counts,$countsq,$tags,$out_fh)
    close $outh if($options{'o'});
}

exit;


sub counting{
    my ($map_fh,$out_fh)=@_;
    foreach my $k(keys %counts){
        if(not $counts{$k}{'global'}{'t'}){  
            $counts{$k}{'global'}{'t'}= 0;
        } 
    }
    my $oldchr='';
    my @a;
    my @a_ex;
    my $prev_min=0;
    my $prev_min_ex=0;
    my $cc=0;

    if($map_fh =~ /sam/){	
        if(not $options{'s'} and $map_fh !~ /sorted$/){
            open (MAP, "cat $map_fh |sort -k3,3 -k4,4n |") or die "can not open $map_fh\n";
        }else{
            open (MAP, "cat $map_fh |") or die "can not open $map_fh\n";
        }	
    }elsif($map_fh =~ /.bam$/){
        open (MAP, "samtools view $map_fh|sort -k3,3 -k4,4n |") or die "can not open $map_fh file or no samtools installed\n";
    }else{
        die "\n---File $map_fh cannot be read by this tool\n";
    }
    if($options{'B'}){
        my $oi='';
        if($map_fh =~ /(\w{2,3}_[123])_including_ERCC92/){
            $oi=$1;
        }elsif($map_fh =~ /(\w{1,2}_[123])/){
            $oi=$1;
        }else{
            die "could not match index\n";
        }

        open OBAM,">$oi.obam" or die "Cannot create obam outfile\n";
        open OBAM2,">$oi.obam2" or die "Cannot create obam2 outfile\n";
    }
    my $processed=0;
    while(my $line = <MAP>){
        next if($line =~ /^@/);
        $cc++;
        if($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S.+)$/){
            my $query=$1;
            my $tag ='';
            my $dmp;
            ## get sample tag here only if we have option e specified
            if($options{'e'}){  
                ($tag,$dmp)=split('_',$query);
                if($ignore{$tag}){
                    next;
                }    
                $tags{$tag}++;
            }


            my $mult=1;
            my $flag=$2;
            next if($flag == 4);

            my $mstrand='+';
            ## that should work but we need to check still
            if($flag & 16){ $mstrand = '-';}

            ## next if unmapped
            my $chr=$3;

            next if($gchr and $gchr ne $chr);
            $processed++;	
            print STDERR "$processed\r";
            my $map_pos=$4; ## one base end end coord is included
            my $cigar=$6;
            my $rest=$7;
            

            die "$cc\n$line" if($flag !~ /^\d+$/);
            next if($flag == 4);

            my $mappings=1;
            if($rest=~/NH:i:(\S+)/){
                $mappings=$1;
            }elsif($rest =~ /NH/){
                print "problem with mappings $map_fh with line $cc\n$line\n"; 
                return;

            }else{
                if($options{'C'}){
                    $mappings=$nh{$query};
                }
            }

            ## get multiplicity of read here
            ## lets soften the multiplicity here
            ##if(not $options{'d'} and $query =~ /_x(\d+)$/)
            if(not $options{'d'} and $query =~ /_x(\d+)/){
                $mult=$1;
            }


## these guys are not in the gtf file but the mappings file 
            if(not $h{$chr}){
                print OTHER $line if ($options{'p'});
                $counts{'other'}{'global'}{'t'}+=($mult/$mappings);
                $counts{'other'}{$tag}{'t'}+=($mult/$mappings) if($tag);
                if($mappings == 1){
                    $counts{'other'}{'global'}{'u'}++;
                    $counts{'other'}{$tag}{'u'}++ if($tag);
                }else{
                    $counts{'other'}{'global'}{'m'}++;
                    $counts{'other'}{$tag}{'m'}++ if($tag);
                }
                next;
            }

## ignore multimappers
            next if($uniq and $mappings > 1);

## check read here now
            my $best_middle=();
            my $best_lng=0;

            my $startpos=$map_pos;
            ## in case we have a read mapping to the antisense strand then the actual start position is start+read-len -1 => This means also that mismatches occur very likely in the beginning
            ## of the read in the sam file since quality is worst here ... ## however, this is heuristically solved right now
            if($flag & 16){
                my @tmp=split(/\s+/,$rest);
                $startpos+=length($tmp[3])-1;
            }



            while($cigar=~/(\d+)(\w)/g){
                my $lng=$1;
                my $type=$2;
                if($type eq "M"){
                    if($lng>$best_lng){
                        $best_middle=int($map_pos+(0.5*$lng));
                    }
                }
                $map_pos=$map_pos+$lng;
            }	

## shouldnt be necessary since it is ordered in input already
            if($oldchr eq ''){
                $oldchr = $chr;
                $prev_min=0;
                @a=sort {$a <=> $b} keys %{$h{$chr}};	
                @a_ex=sort {$a <=> $b} keys %{$h2{$chr}} if($anno_ex);	
            }

            if($oldchr ne $chr){
                $oldchr = $chr;
                $prev_min=0;
                @a=sort {$a <=> $b} keys %{$h{$chr}};
                @a_ex=sort {$a <=> $b} keys %{$h2{$chr}} if ($anno_ex);	
            }

## optimize here, we have sorted input so we dont need to look 
            if($anno_ex){
                my $mins_ex=bins(\@a_ex,$best_middle,$prev_min_ex);
#$prev_min_ex=$mins_ex;
                my $skip=0;
                if($mins_ex == 0 and $best_middle < $a_ex[$mins_ex]){
                }else{
## now we check if 
                    my $i=$mins_ex;
                    my @out=();

## check all indices in question lower than the closest one until we reach certain size limit,,,eg if window is 300, query is 3000 then only ids will be returned were query-start < 300
#running time O(n)
                    if(not $a_ex[$i]){
                        print STDERR "exiting here $i $#a_ex   \n";
                        exit 0;
                    }
                    while( $best_middle-$a_ex[$i] < $window){
## die $best_middle,"\t",keys %{$h{$chr}{$a[$i]}};
                        foreach my $k(keys %{$h2{$chr}{$a_ex[$i]}}){

                            if($best_middle <= $k){ ## check if end position is bigger than query
                                push(@out,$i);
                                last;
                            }
                        }
                        $i--;
                        last if($i < 0);
                    }
## check out indices closer now and get us a weight
# running time O(n)
                    $skip=0;
                    foreach my $i(@out){
###               chr   startpos
                        for my $k(keys %{$h2{$chr}{$a_ex[$i]}}){
## if end position is less then query we add it to out array
                            if($best_middle <=$k){
                                $skip=1;
                            }
                        }
                    }
                }
                next if($skip);
            }

            my $mins=bins(\@a,$best_middle,$prev_min);
            $prev_min=$mins; ## previous min index
            ##$prev_min=0;
            ## doesnt really matter somehow ... difference in mouse genome parsing is in milli seconds :-)
            if($mins == 0 and $best_middle < $a[$mins]){
                print OTHER $line if ($options{'p'});
                $counts{'other'}{'global'}{'t'}+=($mult/$mappings);
                $counts{'other'}{$tag}{'t'}+=($mult/$mappings) if($tag);
                if($mappings ==1){
                    $counts{'other'}{'global'}{'u'}++;
                    $counts{'other'}{$tag}{'u'}++ if($tag);
                }else{
                    $counts{'other'}{'global'}{'m'}++;
                    $counts{'other'}{$tag}{'m'}++ if($tag);
                }
                next;
            }

## now we check if 
            my $i=$mins;
            my @out=();

## check all indices in question lower than the closest one until we reach certain size limit,,,eg if window is 300, query is 3000 then only ids will be returned were query-start < 300
#running time O(n)
            while( $best_middle-$a[$i] < $window){
## die $best_middle,"\t",keys %{$h{$chr}{$a[$i]}};
                foreach my $k(keys %{$h{$chr}{$a[$i]}}){
                    if($best_middle <= $k){ ## check if end position is bigger than query
                        push(@out,$i);
                        last;
                    }
                }
                $i--;
                last if($i < 0);
            }
            if(scalar @out ==0){
                print OTHER $line if ($options{'p'});
                $counts{'other'}{'global'}{'t'}+=($mult/$mappings);
                $counts{'other'}{$tag}{'t'}+=($mult/$mappings) if($tag);
                if($mappings == 1){
                    $counts{'other'}{'global'}{'u'}++;
                    $counts{'other'}{$tag}{'u'}++ if($tag);
                }else{
                    $counts{'other'}{'global'}{'m'}++;
                    $counts{'other'}{$tag}{'m'}++ if($tag);
                }
                next;
            }
## check out indices closer now and get us a weight
# running time O(n)

            my $weight=0; ## number of featured annotations at this position (number of isoforms)
            my $weight2=0; ## number of gene-models annotated at this position (number of genes)
            my @outl=(); ## ids which are in range and hit go in here
            my %whash=();
            my %collect=();
            foreach my $i(@out){
                ###         chr   startpos
                for my $k(keys %{$h{$chr}{$a[$i]}}){
                    ## if end position is more then query we add it to out array
                    if($best_middle <=$k){          ### yhiah
                        if($options{'O'}){
                            ## inspect all features at this end pos
                            for(my $ff=0;$ff < scalar @{$h{$chr}{$a[$i]}{$k}{'f'}};$ff++){
                                if($options{'S'}){
                                    if($h{$chr}{$a[$i]}{$k}{'s'}[$ff] eq $mstrand){ ## check for strand and if correct count it
                                        $collect{$h{$chr}{$a[$i]}{$k}{'f'}[$ff]}++;  ## count feature if strand is correct
                                    }
                                }else{
                                    $collect{$h{$chr}{$a[$i]}{$k}{'f'}[$ff]}++;  ## count feature if with disregard of the strand
                                }
                            }
                        }else{
                            foreach my $v(@{$h{$chr}{$a[$i]}{$k}{'f'}}){
#->				#### check here for feature order
                                $whash{$v}++;
                                $weight+=scalar @{$h{$chr}{$a[$i]}{$k}{'f'}};
                                push(@outl,@{$h{$chr}{$a[$i]}{$k}{'f'}});
                            }
                        }
                        ## testing procedure still

                        ## adapt here, this weight works only properly if we just have exons annotated 
                        ## these two were move in the loop above if we do hashing now
                        #$weight+=scalar @{$h{$chr}{$a[$i]}{$k}{'f'}};
                        #push(@outl,@{$h{$chr}{$a[$i]}{$k}{'f'}});
                    }
                }
            }

            if($options{'O'}){
                my $weight2=0;
                for my $no(sort {$a <=> $b} keys %order){
                    if($collect{$order{$no}}){
                        $weight=$collect{$order{$no}};
                        $weight2=$collect{$order{$no}};
                        for(my $i=0;$i< $weight;$i++){
                            push(@outl,$order{$no});
                        }
                    }
                    last if($weight2);
                }
                ## prepare output here
            }else{
                $weight2=scalar keys %whash;
            }

            ## now we add weighed count to each id to get results
            ## foreach my $kkk(keys %whash){print STDERR "$kkk\t$whash{$kkk} $query\n";}
            foreach my $e (@outl){	        
                next if($ambig == 0 and $weight2 > 1);
                $counts{$e}{'global'}{'t'}+=($mult/($weight*$mappings));
                $counts{$e}{$tag}{'t'}+=($mult/($weight*$mappings)) if($tag);
                if($mappings ==1){
                    $counts{$e}{'global'}{'u'}++;
                    $counts{$e}{$tag}{'u'}++ if($tag);
                }else{
                    $counts{$e}{'global'}{'m'}++;
                    $counts{$e}{$tag}{'m'}++ if($tag);
                }
                if($options{'B'}){
                    if($options{'B'} =~ /\w/){

                        if($e eq $options{'B'}){
                            print OBAM $line;
                            print OBAM2 "$e\t",$mult/($weight*$mappings),"\n";
                        }
                    }else{
                        print OBAM $line;
                        print OBAM2 "$e\t",$mult/($weight*$mappings),"\n";

                    }
                }
                if($options{'q'}){
                    if(not exists $countsq{$e}{'global'}{'t'}{$startpos}){
                        $countsq{$e}{'global'}{'mp'}++;
                        #print $line if($e eq 'Mir294');
                    }
                    if($e eq 'Tgif2'){
                        print "$mult\t$weight\t$mappings\t$line\t$map_fh\n";
                    }


                    $countsq{$e}{'global'}{'t'}{$startpos}+=($mult/($weight*$mappings));
                    $countsq{$e}{$tag}{'mp'}++ if(not exists $countsq{$e}{$tag}{'t'}{$startpos});
                    $countsq{$e}{$tag}{'t'}{$startpos}+=($mult/($weight*$mappings)) if($tag);
                    if($mappings ==1){
                        $countsq{$e}{'global'}{'u'}{$startpos}++;
                        $countsq{$e}{$tag}{'u'}{$startpos}++ if($tag);
                    }else{
                        $countsq{$e}{'global'}{'m'}{$startpos}++;
                        $countsq{$e}{$tag}{'m'}{$startpos}++ if($tag);
                    }
                }
            }
        }
    }
    close MAP;
    if($options{'B'}){
        close OBAM;
        close OBAM2;
    }	
}

sub printcounts{
    my ($counts,$countsq,$tags,$out_fh) = (@_);
    if($options{'q'} or $options{'Q'}){
        my @kt=sort keys %$tags;
        print $out_fh "#id\tglobal";

        ## if option given then we output the start pos and abundance for all included cells
        if($options{'Q'}){
        ;
        }else{
            #print all samples in file  
            foreach my $k(@kt){
                print $out_fh "\t$k";
            }
        }
        print $out_fh "\n";

        for my $k(sort {$$countsq{$b}{'global'}{'mp'} <=> $$countsq{$a}{'global'}{'mp'}} keys %$countsq){
            print $out_fh "$k\t$$countsq{$k}{'global'}{'mp'}";
            if($options{'Q'}){
                for my $pos(sort {$a <=> $b} keys %{$$countsq{$k}{'global'}{'t'}}){
                    print $out_fh "\t${pos}:$$countsq{$k}{'global'}{'t'}{$pos}";
                }
                print $out_fh "\n";
            }else{
                foreach my $tag(@kt){
                        if($$countsq{$k}{$tag}){
                            print $out_fh "\t$$countsq{$k}{$tag}{'mp'}";
                        }else{
                             print $out_fh "\t0";
                         }
                }
                print $out_fh "\n";
            }
        }
        exit;
    }


    if($options{'e'}){
        my @kt=sort keys %$tags;
        print $out_fh "#id\tglobal";

        foreach my $k(@kt){
            print $out_fh "\t$k";
        }
        print $out_fh "\n";

        for my $k(sort {$$counts{$b}{'global'}{'t'} <=> $$counts{$a}{'global'}{'t'}} keys %$counts){
            print $out_fh "$k\t$$counts{$k}{'global'}{'t'}";
            foreach my $tag(@kt){
                    if($$counts{$k}{$tag}){
                        print $out_fh "\t$$counts{$k}{$tag}{'t'}";
                    }else{
                         print $out_fh "\t0";
                     }
            }
            print $out_fh "\n";
        }

    }else{
        for my $k(sort {$$counts{$b}{'global'}{'t'} <=> $$counts{$a}{'global'}{'t'}} keys %$counts){
            print $out_fh "$k\t$$counts{$k}{'global'}{'t'}";
            if(not $uniq and $$counts{$k}{'global'}{'m'}){
                if(not $$counts{$k}{'global'}{'u'}){$$counts{$k}{'global'}{'u'}=0;}
                print $out_fh "\t$$counts{$k}{'global'}{'u'}\t$$counts{$k}{'global'}{'m'}";
            }
            print $out_fh "\n";
        }
    }
}

## done reading annotation file

## binary search that gives the array index which holds the value closest to the query but is smaller than the query!
sub bins{
    my ($a,$query,$prev)=(@_);
    my ($imin,$imax,$imid)=($prev, $#$a,0);

    while ($imin <= $imax){
##calculate the midpoint for roughly equal partition
        my $imid = int(($imin+$imax)/2);
        if ($$a[$imid] == $query){
## key found at index imid
            return $imid;
## determine which subarray to search
        }elsif ($$a[$imid] < $query){
##change min index to search upper subarray
            $imin = $imid + 1;
        }else{        
##change max index to search lower subarray
            $imax = $imid - 1;
        }
    }
## key was not found
    return $imin if(not $imin);
    return $imin-1;
}

sub make_NH_tags{
    my ($f,$nh)=@_;
    if(not -f "${f}_mult"){
        if($f !~ /.bam$/){
            system("cut -f1 $f|sort|uniq -c > ${f}_mult");
        }else{
            system("samtools view $f|cut -f1 $f|sort|uniq -c > ${f}_mult");
        }
    }else{
        print STDERR "File ${f}_mult exists\n";
    }

    open IN,"${f}_mult" or die "File ${f}_mult not found\n";

    my @l;

    while(<IN>){
        chomp;
        @l=split();
        next if($l[1] =~ /^@/);
        $$nh{$l[1]}=$l[0];
    }
    close IN;
}


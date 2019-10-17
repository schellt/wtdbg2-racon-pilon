#!/usr/bin/perl

use strict;
use warnings;
use Path::Tiny;
use Cwd 'abs_path';
use IPC::Cmd qw[can_run run];

my $version = "0.4";

sub print_help{
	print STDOUT "\n";
	print STDOUT "wtdbg2-racon-pilon.pl v$version\n";
	print STDOUT "\n";
	print STDOUT "Description:\n";
	print STDOUT "\tAutomatic execution of (wtdbg2,) racon and pilon.\n";
	print STDOUT "\t- The tools wtdbg2 and wtpoa-cns need to be in your \$PATH if longreads should be assembled.\n";
	print STDOUT "\t- For execution of long read polishing your \$PATH should contain minimap2 and racon.\n";
	print STDOUT "\t- If short read polishing with pilon should be executed java, bwa and samtools need to be in\n";
	print STDOUT "\t  your \$PATH, as well -pilon-path and at least one Illumina read file needs to be specified.\n";
	print STDOUT "\t  Paired and unpaired short reads are mapped to the assembly by bwa mem, resulting bam\n";
	print STDOUT "\t  files are sorted by samtools and pilon is executed on the bam files and the assembly.\n";
	print STDOUT "\n";
	print STDOUT "Usage:\n";
	print STDOUT "\twtdbg2-racon-pilon.pl [-l <longreads.fq> -x <longread tech>] [-a <assembly.fa>]\n";
	print STDOUT "\t                      [-p <paired_1.fq>,<paired_2.fq> -u <unpaired.fq>]\n";
	print STDOUT "\n";
	print STDOUT "Mandatory (if -a is not set):\n";
	print STDOUT "\t-l STR\t\t\tFile containing long reads in fastq format\n";
	print STDOUT "\t-x STR\t\t\tLong read sequencing technology for wtdbg2 and minimap2 presets\n";
	print STDOUT "\t\t\t\tValid arguments are: rsII, rs, sequel, sq, nanopore, ont,\n";
	print STDOUT "\t\t\t\tcorrected, ccs\n";
	print STDOUT "Mandatory (if -l and -x is not set):\n";
	print STDOUT "\t-a STR\t\t\tFile containing an assembly in fasta format that should be polished\n";
	print STDOUT "\n";
	print STDOUT "Input/pipeline options: [default]\n";
	print STDOUT "\t-racon-rounds INT\tNumber of racon iterations [3]\n";
	print STDOUT "\t-pilon-rounds INT\tNumber of pilon iterations [3]\n";
	print STDOUT "\t-pilon-path STR\t\tComplete path to the pilon jar file\n";
	print STDOUT "\t-pr STR\t\t\tFile in bed format to restrict pilon polishing to certain regions\n\t\t\t\t[polish all positions]\n";
	print STDOUT "\t-p STR\t\t\tTwo files with paired short reads in fastq format comma sperated\n";
	print STDOUT "\t\t\t\tCan be specified multiple times\n";
	print STDOUT "\t-u STR\t\t\tOne file with unpaired short reads in fastq format\n";
	print STDOUT "\t\t\t\tCan be specified multiple times\n";
	print STDOUT "\t-xmx STR/INT\t\tMaximum java heap size for pilon [current available]\n";
	print STDOUT "\t-t INT\t\t\tNumber of parallel executed processes [1]\n";
	print STDOUT "\t\t\t\tAffects wtdbg2, wtpoa-cns, minimap2, racon, bwa mem, samtools sort \n";
	print STDOUT "\t\t\t\tand pilon.\n";
	print STDOUT "\tPass specific options to tools with:\n";
	print STDOUT "\t-wtdbg-opts [], -wtpoa-opts [], -minimap-opts [], -racon-opts [-u], -bwa-opts [-a -c 10000],\n\t-pilon-opts [--diploid]\n";
	print STDOUT "\tFor example like -wtdbg-opts \'-g 100m\'\n";
	print STDOUT "\n";
	print STDOUT "Output options: [default]\n";
	print STDOUT "\t-o STR\t\t\tOutput directory [.]\n";
	print STDOUT "\t\t\t\tWill be created if not existing\n";
	print STDOUT "\t-pre STR\t\tPrefix of output files [{w|<-a>.}r<-racon-rounds>p<-pilon-rounds>]\n";
	print STDOUT "\t-v\t\t\tPrint executed commands to STDERR [off]\n";
	print STDOUT "\t-dry-run\t\tOnly print commands to STDERR instead of executing [off]\n";
#	print STDOUT "\t-kt\t\tKeep temporary files [off]\n";
	print STDOUT "\n";
	print STDOUT "\t-h or -help\t\tPrint this help and exit\n";
	print STDOUT "\t-version\t\tPrint version number and exit\n";
	exit;
}

sub exe_cmd{
	my ($cmd,$verbose,$dry) = @_;
	if($verbose == 1){
		print STDERR "CMD\t$cmd\n";
	}
	if($dry == 0){
		system("$cmd") == 0 or die "ERROR\tsystem $cmd failed: $?";
	}
}

my $out_dir = abs_path("./");
my $assembly = "";
my $longreads = "";
my $longread_tech = "";
my @paired = ();
my @unpaired = ();
my $threads = 1;
my $prefix = "";
my $verbose = 0;
my $wtdbg2_opts = "";
my $wtpoa_opts = "";
my $minimap_opts = "";
my $racon_opts = "-u ";
my $bwa_opts = "-a -c 10000 ";
my $pilon_opts = "--diploid ";
my $racon_rounds = 3;
my $pilon_rounds = 3;
my $pilon_path = "";
my $keep_tmp = 0;
my $xmx = "";
my $dry = 0;
my $cmd;
my $final_assembly = "";
my $home = `echo \$HOME`;
chomp $home;
my $polish_regions = "";

my $input_error = 0;

if(-f "/proc/meminfo"){
	$xmx = `grep "MemAvailable:" /proc/meminfo | awk '{print \$2"k"}'`;
	chomp $xmx;
}

for (my $i = 0; $i < scalar(@ARGV);$i++){
	if ($ARGV[$i] eq "-o"){
		$out_dir = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-a"){
		if(-l $ARGV[$i+1]){
			my $path =  path($ARGV[$i+1]);
			$assembly = $path->absolute;
		}
		else{
			$assembly = abs_path($ARGV[$i+1]);
		}
	}
	if ($ARGV[$i] eq "-l"){
		$longreads = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-x"){
		$longread_tech = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-p"){
		push(@paired,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-u"){
		push(@unpaired,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-t"){
		$threads = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-xmx"){
		$xmx = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-pre"){
		$prefix = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-v"){
		$verbose = 1;
	}
	if ($ARGV[$i] eq "-bwa-opts"){
		$bwa_opts = $ARGV[$i+1] . " ";	#nonsense flags are skipped from bwa
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-wtdbg-opts"){
		$wtdbg2_opts = $ARGV[$i+1] . " ";
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-wtpoa-opts"){
		$wtpoa_opts = $ARGV[$i+1] . " ";
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-minimap-opts"){
		$minimap_opts = $ARGV[$i+1] . " ";
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-racon-opts"){
		$racon_opts = $ARGV[$i+1] . " ";
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-pilon-opts"){
		$pilon_opts = $ARGV[$i+1] . " ";
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-racon-rounds"){
		$racon_rounds = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-pilon-rounds"){
		$pilon_rounds = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-pilon-path"){
		$pilon_path = abs_path($ARGV[$i+1]);
	}
	if($ARGV[$i] eq "-kt"){
		$keep_tmp = 1;
	}
	if ($ARGV[$i] eq "-dry-run"){
		$dry = 1;
		$verbose = 1;
	}
	if ($ARGV[$i] eq "-pr"){
		$polish_regions = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-h" or $ARGV[$i] eq "-help"){
		print_help;
	}
	if ($ARGV[$i] eq "-version"){
		print STDERR $version . "\n";
		exit;
	}
}

if(scalar(@ARGV) == 0){
	print_help;
}

print STDERR "CMD\t" . $0 . " " . join(" ",@ARGV) . "\n";

if($longreads eq ""){
	if($racon_rounds != 0){
		print STDERR "INFO\tNo long read file specified! Setting racon rounds to 0\n";
		$racon_rounds = 0;
	}
}

if($longreads ne ""){
	if($longread_tech eq ""){
		print STDERR "ERROR\tLong read technology needs to be specified!\n";
		$input_error = 1;
	}
	else{
		if($longread_tech ne "rsII" and $longread_tech ne "rs" and $longread_tech ne "sequel" and $longread_tech ne "sq" and $longread_tech ne "nanopore" and $longread_tech ne "ont" and $longread_tech ne "corrected"and $longread_tech ne "ccs"){
			print STDERR "ERROR\tNo valid long read technology specified!\n";
			$input_error = 1;
		}
	}
}

if($assembly eq "" and $longreads eq ""){
	print STDERR "ERROR\tSpecify either longreads or an assembly!\n";
	$input_error = 1;
}

if($assembly eq ""){
	if(not defined(can_run("wtdbg2"))){
		print STDERR "ERROR\twtdbg2 is not in your \$PATH\n";
		$input_error = 1;
	}
	if(not defined(can_run("wtpoa-cns"))){
		print STDERR "ERROR\twtpoa-cns is not in your \$PATH\n";
		$input_error = 1;
	}
}

if($racon_rounds > 0){
	if(not defined(can_run("minimap2"))){
		print STDERR "ERROR\tminimap2 is not in your \$PATH and racon rounds > 0!\n";
		$input_error = 1;
	}
	if(not defined(can_run("racon"))){
		print STDERR "ERROR\tracon is not in your \$PATH and racon rounds > 0!\n";
		$input_error = 1;
	}
}

if($pilon_rounds > 0){
	if(not defined(can_run("bwa"))){
		print STDERR "ERROR\tbwa is not in your \$PATH and pilon rounds > 0!\n";
		$input_error = 1;
	}
	if(not defined(can_run("samtools"))){
		print STDERR "ERROR\tsamtools is not in your \$PATH and pilon rounds > 0!\n";
		$input_error = 1;
	}
	if(not defined(can_run("java"))){
		print STDERR "ERROR\tjava is not in your \$PATH and pilon rounds > 0!\n";
		$input_error = 1;
	}
	if($pilon_path eq ""){
		print STDERR "ERROR\t-pilon-path is not specified and pilon rounds > 0!\n";
		$input_error = 1;
	}
	if($polish_regions ne "" and not -f "$polish_regions"){
		print STDERR "ERROR\t-pr is not a file!\n";
		$input_error = 1;
	}
}


if(-f "$out_dir"){
	print STDERR "ERROR\tOutput directory $out_dir is already a file!\n";
	$input_error = 1;
}

if($threads !~ m/^\d+$/ or $threads < 1){
	print STDERR "ERROR\tThreads is no integer >= 1!\n";
	$input_error = 1;
}
if($racon_rounds !~ m/^\d+$/ or $racon_rounds < 0){
	print STDERR "ERROR\t-racon-rounds is no integer >= 0!\n";
	$input_error = 1;
}
if($pilon_rounds !~ m/^\d+$/ or $pilon_rounds < 0){
	print STDERR "ERROR\t-pilon-rounds is no integer >= 0!\n";
	$input_error = 1;
}

if($xmx eq "" and $pilon_rounds > 0){
	print STDERR "ERROR\tJava Xmx not specified and /proc/meminfo not present!\n";
	$input_error = 1;
}
if($longreads ne ""){
	if(not -f "$longreads"){
		print STDERR "ERROR\tNo existing long read file specified!\n";
		$input_error = 1;
	}
}

if($assembly ne ""){
	if(not -f "$assembly"){
		print STDERR "ERROR\tNo existing assembly file specified!\n";
		$input_error = 1;
	}
}

if ($input_error == 1){
	print STDERR "ERROR\tInput error detected!\n";
	exit 1;
}

if(not -d "$out_dir"){
	print STDERR "INFO\tCreating output directory $out_dir\n";
	$cmd="mkdir -p $out_dir";
	exe_cmd($cmd,$verbose,$dry);
}

if($prefix eq ""){
	if($assembly eq ""){
		$prefix = "wr" . $racon_rounds . "p" . $pilon_rounds;
	}
	else{
		$a=(split /\//,$assembly)[-1];
		$prefix = $a . ".r" . $racon_rounds . "p" . $pilon_rounds;
	}
	print STDERR "INFO\tSetting Outpufile prefix to $prefix\n";
}

if($longread_tech eq "rsII" or $longread_tech eq "rs" or $longread_tech eq "sequel" or $longread_tech eq "sq"){
	if($minimap_opts eq ""){
		$minimap_opts = "-x map-pb -H ";
	}
	else{
		$minimap_opts = "-x map-pb -H " . $minimap_opts;
	}
}

if($longread_tech eq "nanopore" or $longread_tech eq "ont"){
	if($minimap_opts eq ""){
		$minimap_opts = "-x map-ont ";
	}
	else{
		$minimap_opts = "-x map-ont " . $minimap_opts;
	}
}

my %paired_filter;

foreach(@paired){
	my @pair = split(/,/,$_);
	foreach(@pair){
		if($_ =~ m/^~/){
			$_ =~ s/^~/$home/;	#~ is translated by bash into $HOME. This does not work if there is no space infront. That means if the second file starts with "~" it will not be recognized even though it exists
		}
	}
	if(scalar(@pair) != 2){
		print STDERR "INFO\tNot a pair: $_ - skipping these file(s)\n";
	}
	else{
		my $file_error = 0;
		if(not -f "$pair[0]"){
			print STDERR "INFO\tNo file $pair[0] - skipping pair $_\n";
			$file_error = 1;
		}
		if(not -f "$pair[1]"){
			print STDERR "INFO\tNo file $pair[1] - skipping pair $_\n";
			$file_error = 1;
		}
		if($file_error == 0){
			if(exists($paired_filter{abs_path($pair[0]) . "," . abs_path($pair[1])})){
				print STDERR "INFO\tpair " . abs_path($pair[0]) . "," . abs_path($pair[1]) . " already specified\n";
			}
			else{
				$paired_filter{abs_path($pair[0]) . "," . abs_path($pair[1])} = 1;
			}
		}
	}
}

my %unpaired_filter;

foreach(@unpaired){
	if(not -f "$_"){
		print STDERR "INFO\tNo file $_ - skipping this file\n";
	}
	else{
		if(exists($unpaired_filter{abs_path($_)})){
			print STDERR "INFO\tfile " . abs_path($_) . " already specified\n";
		}
		else{
			$unpaired_filter{abs_path($_)} = 1;
		}
	}
}

if(scalar(keys(%paired_filter)) == 0 and scalar(keys(%unpaired_filter)) == 0){
	print STDERR "INFO\tNo existing short read files specified! Setting pilon rounds to 0\n";
	$pilon_rounds = 0;
}

if($racon_rounds == 0 and $pilon_rounds == 0 and $assembly ne ""){
	print STDERR "ERROR\tAssembly as input (-a) and racon as well as pilon rounds are 0\n";
	print STDERR "ERROR\tNothing to do.\n";
	exit 1;
}

my $wtdbg2_version;
my $wtpoa_cns_version;

if($assembly eq ""){
	$wtdbg2_version = `wtdbg2 -V | awk '{print \$2}'`;
	chomp $wtdbg2_version;

	$wtpoa_cns_version = `wtpoa-cns -V | awk '{print \$2}'`;
	chomp $wtpoa_cns_version;
}

my $minimap_version;
my $racon_version;
if($racon_rounds > 0){
	$minimap_version = `minimap2 --version`;
	chomp $minimap_version;
	
	$racon_version = `racon --version | sed 's/^v//'`;
	chomp $racon_version;
}

my $bwa_version;
my $samtools_version;
my $java_version;
my $pilon_version;
if($pilon_rounds > 0){
	$bwa_version = `bwa 2>&1 | head -3 | tail -1 | sed 's/^Version: //'`;
	chomp $bwa_version;
	
	$samtools_version = `samtools --version | head -1 | sed 's/^samtools //'`;
	chomp $samtools_version;
	
	$java_version = `java -version 2>&1 | grep "version"`;
	chomp $java_version;
	
	$pilon_version = `java -jar $pilon_path --version | sed 's/^Pilon version //;s/ .*\$//'`;
	chomp $pilon_version;
}

my $verbose_word = "";
if($verbose == 0){
	$verbose_word = "No";
}
else{
	$verbose_word = "Yes";
}

my $keep_word = "";
if($keep_tmp == 0){
	$keep_word = "No";
}
else{
	$keep_word = "Yes";
}

print "\n";
print "wtdbg2-racon-pilon.pl v$version\n";
print "\n";
print "Detected tools\n";
print "==============\n";
if($assembly eq ""){
	print "wtdbg:                " . $wtdbg2_version . "\n";
	print "wtpoa-cns:            " . $wtpoa_cns_version . "\n";
}
if($racon_rounds > 0){
	print "minimap:              " . $minimap_version . "\n";
	print "racon:                " . $racon_version . "\n";
}
if($pilon_rounds > 0){
	print "bwa:                  " . $bwa_version . "\n";
	print "samtools:             " . $samtools_version . "\n";
	print "java:                 " . $java_version . "\n";
	print "pilon:                " . $pilon_version . "\n";
}
print "\n";
print "User defined input\n";
print "==================\n";
print "output directory:     " . $out_dir . "\n";
if($assembly ne ""){
	print "assembly:             " . $assembly . "\n";
}
if($longreads ne ""){
	print "long reads:           " . $longreads . "\n";
}
if(scalar(keys(%paired_filter)) > 0){
	print "paired reads:         ";
	print join("\n                      ",keys(%paired_filter)) . "\n";
}
if(scalar(keys(%unpaired_filter)) > 0){
	print "unpaired reads:       ";
	print join("\n                      ",keys(%unpaired_filter)) . "\n";
}
print "number of threads:    " . $threads . "\n";
print "outpufile prefix:     " . $prefix . "\n";
print "verbose:              " . $verbose_word . "\n";
#print "Keep temporary files: " . $keep_word . "\n";
if($assembly eq ""){
	if($wtdbg2_opts ne ""){
		print "wtdbg options:        " . $wtdbg2_opts . "\n";
	}
	if($wtpoa_opts ne ""){
		print "wtpoa-cns options:    " . $wtpoa_opts . "\n";
	}
}
if($racon_rounds > 0){
	if($minimap_opts ne ""){
		print "minimap options:      " . $minimap_opts . "\n";
	}
	if($racon_opts ne ""){
		print "racon options:        " . $racon_opts . "\n";
	}
}
if($pilon_rounds > 0){
	if($bwa_opts ne ""){
		print "bwa options:          " . $bwa_opts . "\n";
	}
	if($pilon_opts ne ""){
		print "pilon options:        " . $pilon_opts . "\n";
	}
	print "java Xmx:             " . $xmx . "\n";
	if($polish_regions ne ""){
		print "polish regions:       " . $polish_regions . "\n";
	}
}
print "racon rounds:         " . $racon_rounds . "\n";
print "pilon rounds:         " . $pilon_rounds . "\n";

if($assembly eq ""){
	$cmd = "wtdbg2 -x $longread_tech $wtdbg2_opts-i $longreads -t $threads -o $out_dir/$prefix > $out_dir/$prefix\_wtdbg2.log 2> $out_dir/$prefix\_wtdbg2.err";
	exe_cmd($cmd,$verbose,$dry);
	
	$cmd = "wtpoa-cns $wtpoa_opts-t $threads -i $out_dir/$prefix.ctg.lay.gz -o $out_dir/$prefix.ctg.fa > $out_dir/$prefix\_wtpoa.log 2> $out_dir/$prefix\_wtpoa.err";
	exe_cmd($cmd,$verbose,$dry);
	
	$final_assembly = "$out_dir/$prefix.ctg.fa";
	
	$assembly = "$out_dir/$prefix.ctg.fa";
}

if($longreads ne ""){
	for(my $i=1; $i < $racon_rounds+1; $i++){
		if($i > 1){
			my $before = $i-1;
			$assembly = "$out_dir/$prefix\_racon_$before.fasta";
		}
		$cmd = "minimap2 $minimap_opts-t $threads $assembly $longreads > $out_dir/$prefix\_minimap2_$i.paf 2> $out_dir/$prefix\_minimap2_$i.err";
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "racon $racon_opts-t $threads $longreads $out_dir/$prefix\_minimap2_$i.paf $assembly > $out_dir/$prefix\_racon_$i.fasta 2> $out_dir/$prefix\_racon_$i.err";
		$assembly = "$out_dir/$prefix\_racon_$i.fasta";
		$final_assembly = "$out_dir/$prefix\_racon_$i.fasta";
		exe_cmd($cmd,$verbose,$dry);
	}
}

my $samtools_threads = $threads - 1;

for(my $j = 1; $j < $pilon_rounds+1; $j++){
	my @paired_bams = ();
	my @unpaired_bams = ();
	if($j > 1){
		my $before = $j-1;
		$assembly = "$out_dir/$prefix\_pilon_$before/$prefix\_pilon_$before.fasta";
	}

	if(not -f $assembly . ".amb" or not -f $assembly . ".ann" or not -f $assembly . ".bwt" or not -f $assembly . ".pac" or not -f $assembly . ".sa"){
		$cmd = "bwa index $assembly > $out_dir/$prefix\_index_$j.log 2> $out_dir/$prefix\_index_$j.err";
		exe_cmd($cmd,$verbose,$dry);
	}
	else{
		print STDERR "INFO\tIndex files for $assembly already existing\n";
	}
	
	my $p = 0;
	foreach(keys(%paired_filter)){
		$p++;
		my ($for,$rev) = split(/,/,$_);
		if($polish_regions ne ""){
			$cmd = "bwa mem $bwa_opts-t $threads $assembly $for $rev 2> $out_dir/$prefix\_bwa_mem_paired_$j.$p.err | samtools view -L $polish_regions -1 -b - > $out_dir/$prefix\_paired_$j.$p.bam";
		}
		else{
			$cmd = "bwa mem $bwa_opts-t $threads $assembly $for $rev 2> $out_dir/$prefix\_bwa_mem_paired_$j.$p.err | samtools view -1 -b - > $out_dir/$prefix\_paired_$j.$p.bam";
		}
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "samtools sort -l 9 -@ $samtools_threads -T $out_dir/$prefix\_paired_$j.$p -o $out_dir/$prefix\_paired_$j.$p.sort.bam $out_dir/$prefix\_paired_$j.$p.bam";
		push(@paired_bams,"$out_dir/$prefix\_paired_$j.$p.sort.bam");
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "rm $out_dir/$prefix\_paired_$j.$p.bam";
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "samtools index $out_dir/$prefix\_paired_$j.$p.sort.bam";
		exe_cmd($cmd,$verbose,$dry);
	}
	my $u = 0;
	foreach(keys(%unpaired_filter)){
		$u++;
		if($polish_regions ne ""){
			$cmd = "bwa mem $bwa_opts-t $threads $assembly $_ 2> $out_dir/$prefix\_bwa_mem_unpaired_$j.$u.err | samtools view -L $polish_regions -1 -b - > $out_dir/$prefix\_unpaired_$j.$u.bam";
		}
		else{
			$cmd = "bwa mem $bwa_opts-t $threads $assembly $_ 2> $out_dir/$prefix\_bwa_mem_unpaired_$j.$u.err | samtools view -1 -b - > $out_dir/$prefix\_unpaired_$j.$u.bam";
		}
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "samtools sort -l 9 -@ $samtools_threads -T $out_dir/$prefix\_unpaired_$j.$u -o $out_dir/$prefix\_unpaired_$j.$u.sort.bam $out_dir/$prefix\_unpaired_$j.$u.bam";
		push(@unpaired_bams,"$out_dir/$prefix\_unpaired_$j.$u.sort.bam");
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "rm $out_dir/$prefix\_unpaired_$j.$u.bam";
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "samtools index $out_dir/$prefix\_unpaired_$j.$u.sort.bam";
		exe_cmd($cmd,$verbose,$dry);
	}
	
	my $paired_bam_files = "";
	foreach(@paired_bams){
		$paired_bam_files = $paired_bam_files . " --frags " . $_;
	}
	$paired_bam_files =~ s/^ //;
	my $unpaired_bam_files = "";
	foreach(@unpaired_bams){
		$unpaired_bam_files = $unpaired_bam_files . " --unpaired " . $_;
	}
	$unpaired_bam_files =~ s/^ //;
	$cmd = "java -Xmx$xmx -jar $pilon_path $pilon_opts--genome $assembly $paired_bam_files $unpaired_bam_files --output $prefix\_pilon_$j --outdir $out_dir/$prefix\_pilon_$j --threads $threads > $out_dir/$prefix\_pilon_$j.log 2> $out_dir/$prefix\_pilon_$j.err";
	exe_cmd($cmd,$verbose,$dry);
	
	if($dry == 0){
		open (IN, '<', "$out_dir/$prefix\_pilon_$j/$prefix\_pilon_$j.fasta") or die "Could not open Inputfile $out_dir/$prefix\_pilon_$j/$prefix\_pilon_$j.fasta\n";
		open (OUT, '>', "$out_dir/$prefix\_pilon_$j/$prefix\_pilon_$j.fasta.tmp") or die "Could not open Outputfile $out_dir/$prefix\_pilon_$j/$prefix\_pilon_$j.fasta.tmp\n";
		while (my $line = <IN>){
			chomp $line;
			if($line =~ m/^>/){
				$line =~ s/.pilon$//;
				print OUT $line . "\n";
			}
			else{
				print OUT $line . "\n";
			}
		}
		close IN;
		close OUT;
	}
	else{
		print STDERR "I would rename headers of $out_dir/$prefix\_pilon_$j/$prefix\_pilon_$j.fasta\n";
	}
	$cmd = "mv $out_dir/$prefix\_pilon_$j/$prefix\_pilon_$j.fasta.tmp $out_dir/$prefix\_pilon_$j/$prefix\_pilon_$j.fasta";
	exe_cmd($cmd,$verbose,$dry);
	
	$final_assembly = "$out_dir/$prefix\_pilon_$j/$prefix\_pilon_$j.fasta";
}

print "Final assembly:       " . $final_assembly . "\n";

exit;

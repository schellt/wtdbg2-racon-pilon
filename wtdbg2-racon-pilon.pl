#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use IPC::Cmd qw[can_run run];

my $version = "0.2";

sub print_help{
	print STDOUT "\n";
	print STDOUT "wtdbg2-racon-pilon.pl v$version\n";
	print STDOUT "\n";
	print STDOUT "Description:\n";
	print STDOUT "\tAutomatic execution of wtdbg2, racon and pilon.\n";
	print STDOUT "\t- The tools wtdbg2 and wtpoa-cns need to be in your \$PATH.\n";
	print STDOUT "\t- For execution of racon your \$PATH should contain racon and minimap2.\n";
	print STDOUT "\t  Afterwards long reads are mapped to the assembly and correction(s) with racon are executed.\n";
	print STDOUT "\t- If pilon should be executed java, bwa and samtools need to be in your \$PATH,\n";
	print STDOUT "\t  as well -pilon-path and at least one Illumina read file need to be specified.\n";
	print STDOUT "\t  Paired and unpaired short reads are mapped to the assembly by bwa mem, resulting bam\n";
	print STDOUT "\t  files are sorted by samtools and pilon is executed on the bam files and the assembly.\n";
	print STDOUT "\n";
	print STDOUT "Usage:\n";
	print STDOUT "\twtdbg2-racon-pilon.pl -l <longreads.fq> -x <longread tech>\n";
	print STDOUT "\t[-p <paired_1.fq>,<paired_2.fq> -u <unpaired.fq>]\n";
	print STDOUT "\n";
	print STDOUT "Mandatory:\n";
	print STDOUT "\t-l STR\t\t\tFile containing long reads in fastq format\n";
	print STDOUT "\t-x STR\t\t\tLong read sequencing technology for wtdbg2 presets\n";
	print STDOUT "\t\t\t\tValid arguments are: rsII, rs, sequel, sq, nanopore, ont,\n";
	print STDOUT "\t\t\t\tcorrected, ccs\n";
	print STDOUT "\n";
	print STDOUT "Input/pipeline options: [default]\n";
	print STDOUT "\t-racon-rounds INT\tNumber of racon iterations [3]\n";
	print STDOUT "\t-pilon-rounds INT\tNumber of pilon iterations [3]\n";
	print STDOUT "\t-pilon-path STR\t\tComplete path to the pilon jar file\n";
	print STDOUT "\t-sr-mapper STR\t\tShort read mapper for correction with pilon [bwa]\n";
	print STDOUT "\t\t\t\tValid arguments are: bwa, ngm\n";
	print STDOUT "\t\t\t\tTo execute different mappers in different pilon rounds supply a\n";
	print STDOUT "\t\t\t\tcomma separated list of mappers. Overrides -pilon-rounds.\n";
	print STDOUT "\t-p STR\t\t\tTwo files with paired short reads in fastq format comma sperated\n";
	print STDOUT "\t\t\t\tCan be specified multiple times\n";
	print STDOUT "\t-u STR\t\t\tOne file with unpaired short reads in fastq format\n";
	print STDOUT "\t\t\t\tCan be specified multiple times\n";
	print STDOUT "\t-xmx STR/INT\t\tMaximum java heap size for pilon [current available]\n";
	print STDOUT "\t-t INT\t\t\tNumber of parallel executed processes [1]\n";
	print STDOUT "\t\t\t\tAffects wtdbg2, wtpoa-cns, minimap2, racon, bwa mem, ngm,\n";
	print STDOUT "\t\t\t\tsamtools sort and pilon.\n";
	print STDOUT "\tPass specific options to tools with:\n";
	print STDOUT "\t-wtdbg-opts, -wtpoa-opts, -minimap-opts, -racon-opts, -bwa-opts, -ngm-opts, -pilon-opts\n";
	print STDOUT "\tFor example like -wtdbg-opts \'-g 100m\'\n";
	print STDOUT "\n";
	print STDOUT "Output options: [default]\n";
	print STDOUT "\t-o STR\t\t\tOutput directory [.]\n";
	print STDOUT "\t\t\t\tWill be created if not existing\n";
	print STDOUT "\t-pre STR\t\tPrefix of output files [wr<racon rounds>p<pilon rounds>]\n";
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
#		print STDERR "I would exe sth\n";
		system("$cmd") == 0 or die "ERROR\tsystem $cmd failed: $?";
	}
}

my $out_dir = abs_path("./");
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
my $racon_opts = "";
my $bwa_opts = "-a -c 10000 ";
my $ngm_opts = "";
my $pilon_opts = "";
my $racon_rounds = 3;
my $pilon_rounds = 3;
my $pilon_path = "";
my $sr_mapper = "bwa";
my @srmapper = ();
my $keep_tmp = 0;
my $xmx = "";
my $dry = 0;
my $cmd;
my $final_assembly = "";
my $home = `echo \$HOME`;
chomp $home;

my $input_error = 0;

if(-f "/proc/meminfo"){
	$xmx = `grep "MemAvailable:" /proc/meminfo | awk '{print \$2"k"}'`;
	chomp $xmx;
}

for (my $i = 0; $i < scalar(@ARGV);$i++){
	if ($ARGV[$i] eq "-o"){
		$out_dir = abs_path($ARGV[$i+1]);
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
	if ($ARGV[$i] eq "-ngm-opts"){
		$ngm_opts = $ARGV[$i+1] . " ";
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
	if ($ARGV[$i] eq "-sr-mapper"){
		$sr_mapper = $ARGV[$i+1];
	}
	if($ARGV[$i] eq "-kt"){
		$keep_tmp = 1;
	}
	if ($ARGV[$i] eq "-dry-run"){
		$dry = 1;
		$verbose = 1;
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

if(not defined(can_run("wtdbg2"))){
	print STDERR "ERROR\twtdbg2 is not in your \$PATH\n";
	$input_error = 1;
}
if(not defined(can_run("wtpoa-cns"))){
	print STDERR "ERROR\twtpoa-cns is not in your \$PATH\n";
	$input_error = 1;
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

if($sr_mapper =~ m/,/){
	@srmapper = split(/,/,$sr_mapper);
	my $sr_mapper_error = 0;
	foreach(@srmapper){
		if($_ ne "bwa" and $_ ne "ngm"){
			print STDERR "ERROR\t$_ is not a valid argument for -sr-mapper!\n";
			$input_error = 1;
			$sr_mapper_error = 1;
		}
	}
	if($sr_mapper_error == 1){
		print STDERR "ERROR\tNo valid short read mapper specified!\n";
		$input_error = 1;
	}
	else{
		$pilon_rounds = scalar(@srmapper);
	}
}
else{
	if($sr_mapper ne "bwa" and $sr_mapper ne "ngm" and $pilon_rounds > 0){
		print STDERR "ERROR\tNo valid short read mapper specified!\n";
		$input_error = 1;
	}
	else{
		@srmapper = (($sr_mapper) x $pilon_rounds);
	}
}

if($pilon_rounds > 0){
	if(not defined(can_run("bwa")) and $sr_mapper =~ /bwa/){
		print STDERR "ERROR\tbwa is not in your \$PATH and pilon rounds > 0!\n";
		$input_error = 1;
	}
	if(not defined(can_run("ngm")) and $sr_mapper =~ /ngm/){
		print STDERR "ERROR\tngm is not in your \$PATH and pilon rounds > 0!\n";
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
}


if(-f "$out_dir"){
	print STDERR "ERROR\tOutput directory $out_dir is already a file!\n";
	$input_error = 1;
}
else{
	if(not -d "$out_dir"){
		print STDERR "INFO\tCreating output directory $out_dir\n";
		$cmd="mkdir -p $out_dir";
		exe_cmd($cmd,$verbose,$dry);
	}
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
if($xmx eq "" and $pilon_rounds > 0){
	print STDERR "ERROR\tJava Xmx not specified and /proc/meminfo not present!\n";
	$input_error = 1;
}
if(not -f "$longreads"){
	print STDERR "ERROR\tNo existing long read file specified!\n";
	$input_error = 1;
}

if ($input_error == 1){
	print STDERR "ERROR\tInput error detected!\n";
	exit 1;
}

if($prefix eq ""){
	$prefix = "wr" . $racon_rounds . "p" . $pilon_rounds;
	print STDERR "INFO\tSetting Outpufile prefix to $prefix\n";
}
if($longread_tech eq "rsII" or $longread_tech eq "rs" or $longread_tech eq "sequel" or $longread_tech eq "sq"){
	if($minimap_opts eq ""){
		$minimap_opts = "-H ";
	}
	else{
		$minimap_opts = "-H " . $minimap_opts
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

if(scalar(keys(%paired_filter)) == 0 and scalar(keys(%unpaired_filter)) == 0 and $pilon_rounds > 0){
	print STDERR "ERROR\tNo existing short read files specified and number of pilon rounds > 0!\n";
	exit 1;
}

my $wtdbg2_version = `wtdbg2 -V | awk '{print \$2}'`;
chomp $wtdbg2_version;

my $wtpoa_cns_version = `wtpoa-cns -V | awk '{print \$2}'`;
chomp $wtpoa_cns_version;

my $minimap_version;
my $racon_version;
if($racon_rounds > 0){
	$minimap_version = `minimap2 --version`;
	chomp $minimap_version;
	
	$racon_version = `racon --version | sed 's/^v//'`;
	chomp $racon_version;
}

my $bwa_version;
my $ngm_version;
my $samtools_version;
my $java_version;
my $pilon_version;
if($pilon_rounds > 0){
	if($sr_mapper =~ /bwa/){
		$bwa_version = `bwa 2>&1 | head -3 | tail -1 | sed 's/^Version: //'`;
		chomp $bwa_version;
	}
	
	if($sr_mapper =~ /ngm/){
		$ngm_version = `ngm 2>&1 | grep "\\[MAIN\\] NextGenMap " | awk '{print \$NF}'`;
		chomp $ngm_version;
	}
	
	$samtools_version = `samtools --version | head -1 | sed 's/^samtools //'`;
	chomp $samtools_version;
	
	$java_version = `java -version 2>&1 | head -1 | awk '{gsub("\\"","",\$NF);print \$NF}'`;
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
print "wtdbg:                " . $wtdbg2_version . "\n";
print "wtpoa-cns:            " . $wtpoa_cns_version . "\n";
if($racon_rounds > 0){
	print "minimap:              " . $minimap_version . "\n";
	print "racon:                " . $racon_version . "\n";
}
if($pilon_rounds > 0){
	if($sr_mapper =~ m/bwa/){
		print "bwa:                  " . $bwa_version . "\n";
	}
	if($sr_mapper =~ m/ngm/){
		print "ngm:                  " . $ngm_version . "\n";
	}
	print "samtools:             " . $samtools_version . "\n";
	print "java:                 " . $java_version . "\n";
	print "pilon:                " . $pilon_version . "\n";
}
print "\n";
print "User defined input\n";
print "==================\n";
print "output directory:     " . $out_dir . "\n";
print "long reads:           " . $longreads . "\n";
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
if($wtdbg2_opts ne ""){
	print "wtdbg options:        " . $wtdbg2_opts . "\n";
}
if($wtpoa_opts ne ""){
	print "wtpoa-cns options:    " . $wtpoa_opts . "\n";
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
	print "short read mapper:    " . $sr_mapper . "\n";
	if($bwa_opts ne "" and $sr_mapper =~ m/bwa/){
		print "bwa options:          " . $bwa_opts . "\n";
	}
	if($ngm_opts ne "" and $sr_mapper =~ m/ngm/){
		print "ngm options:          " . $ngm_opts . "\n";
	}
	if($pilon_opts ne ""){
		print "pilon options:        " . $pilon_opts . "\n";
	}
	print "java Xmx:             " . $xmx . "\n";
}
print "racon rounds:         " . $racon_rounds . "\n";
print "pilon rounds:         " . $pilon_rounds . "\n";

$cmd = "wtdbg2 -x $longread_tech $wtdbg2_opts-i $longreads -t $threads -o $out_dir/$prefix > $out_dir/$prefix\_wtdbg2.log 2> $out_dir/$prefix\_wtdbg2.err";
exe_cmd($cmd,$verbose,$dry);

$cmd = "wtpoa-cns $wtpoa_opts-t $threads -i $out_dir/$prefix.ctg.lay.gz -o $out_dir/$prefix.ctg.fa > $out_dir/$prefix\_wtpoa.log 2> $out_dir/$prefix\_wtpoa.err";
exe_cmd($cmd,$verbose,$dry);

$final_assembly = "$out_dir/$prefix.ctg.fa";

my $assembly = "$out_dir/$prefix.ctg.fa";
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

for(my $j = 0; $j < scalar(@srmapper); $j++){
	my $i = $j+1;
	$sr_mapper = $srmapper[$j];
	my @paired_bams = ();
	my @unpaired_bams = ();
	if($i > 1){
		my $before = $i-1;
		$assembly = "$out_dir/$prefix\_pilon_$before/$prefix\_pilon_$before.fasta";
	}

	if($sr_mapper eq "bwa"){	
		$cmd = "bwa index $assembly > $out_dir/$prefix\_index_$i.log 2> $out_dir/$prefix\_index_$i.err";
		exe_cmd($cmd,$verbose,$dry);
	}
	
	my $p = 0;
	foreach(keys(%paired_filter)){
		$p++;
		my ($for,$rev) = split(/,/,$_);
		if($sr_mapper eq "bwa"){
			$cmd = "bwa mem $bwa_opts-t $threads $assembly $for $rev 2> $out_dir/$prefix\_bwa_mem_paired_$i.$p.err | samtools view -1 -b - > $out_dir/$prefix\_paired_$i.$p.bam";
		}
		if($sr_mapper eq "ngm"){
			my $paired_ngm_opts = $ngm_opts;
			my @p_ngm_opts = split(/ /,$ngm_opts);
			my $change_ngm_opts = 0;
			for(my $i = 0; $i < scalar(@p_ngm_opts); $i++){
				if($p_ngm_opts[$i] eq "-n" or $p_ngm_opts[$i] eq "--topn"){
					$change_ngm_opts = 1;
					splice (@p_ngm_opts,$i,2);
				}
			}
			if($change_ngm_opts == 1){
				print STDERR "INFO\tPaired end mode with -n/--topn > 1 is not supported in ngm. Removing the option.\n";
				$paired_ngm_opts = join(" ",@p_ngm_opts) . " ";
			}
			$cmd = "ngm $paired_ngm_opts-t $threads -r $assembly -1 $for -2 $rev -o $out_dir/$prefix\_paired_$i.$p.sam > $out_dir/$prefix\_ngm_paired_$i.$p.log 2> $out_dir/$prefix\_ngm_paired_$i.$p.err && samtools view -@ $threads -b $out_dir/$prefix\_paired_$i.$p.sam > $out_dir/$prefix\_paired_$i.$p.bam 2> $out_dir/$prefix\_view_paired_$i.$p.err && rm $out_dir/$prefix\_paired_$i.$p.sam";
		}
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "samtools sort -l 9 -@ $threads -T $out_dir/$prefix\_paired_$i.$p -o $out_dir/$prefix\_paired_$i.$p.sort.bam $out_dir/$prefix\_paired_$i.$p.bam";
		push(@paired_bams,"$out_dir/$prefix\_paired_$i.$p.sort.bam");
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "rm $out_dir/$prefix\_paired_$i.$p.bam";
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "samtools index $out_dir/$prefix\_paired_$i.$p.sort.bam";
		exe_cmd($cmd,$verbose,$dry);
	}
	my $u = 0;
	foreach(keys(%unpaired_filter)){
		$u++;
		if($sr_mapper eq "bwa"){
			$cmd = "bwa mem $bwa_opts-t $threads $assembly $_ 2> $out_dir/$prefix\_bwa_mem_unpaired_$i.$u.err | samtools view -1 -b - > $out_dir/$prefix\_unpaired_$i.$u.bam";
		}
		if($sr_mapper eq "ngm"){
			$cmd = "ngm $ngm_opts-t $threads -r $assembly -q $_ -o $out_dir/$prefix\_unpaired_$i.$u.sam > $out_dir/$prefix\_ngm_unpaired_$i.$u.log 2> $out_dir/$prefix\_ngm_unpaired_$i.$u.err && samtools view -@ $threads -b $out_dir/$prefix\_unpaired_$i.$u.sam > $out_dir/$prefix\_unpaired_$i.$u.bam 2> $out_dir/$prefix\_view_unpaired_$i.$u.err && rm $out_dir/$prefix\_unpaired_$i.$u.sam";
		}
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "samtools sort -l 9 -@ $threads -T $out_dir/$prefix\_unpaired_$i.$u -o $out_dir/$prefix\_unpaired_$i.$u.sort.bam $out_dir/$prefix\_unpaired_$i.$u.bam";
		push(@unpaired_bams,"$out_dir/$prefix\_unpaired_$i.$u.sort.bam");
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "rm $out_dir/$prefix\_unpaired_$i.$u.bam";
		exe_cmd($cmd,$verbose,$dry);
		$cmd = "samtools index $out_dir/$prefix\_unpaired_$i.$u.sort.bam";
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
	$cmd = "java -Xmx$xmx -jar $pilon_path $pilon_opts--genome $assembly $paired_bam_files $unpaired_bam_files --output $prefix\_pilon_$i --outdir $out_dir/$prefix\_pilon_$i --threads $threads > $out_dir/$prefix\_pilon_$i.log 2> $out_dir/$prefix\_pilon_$i.err";
	exe_cmd($cmd,$verbose,$dry);
	$final_assembly = "$out_dir/$prefix\_pilon_$i/$prefix\_pilon_$i.fasta";
}

print "Final assembly:       " . $final_assembly . "\n";

exit;

# wtdbg2-racon-pilon.pl v0.4

## Description
__Automatic execution of (wtdbg2,) Racon and Pilon.__

Longread assembly will be executed with wtdbg2 if no assembly file is specified.  
Iterative correction with long reads by mapping them to the assembly with minimap2 and polishing with Racon.  
Iterative correction with short reads by mapping them with bwa mem, (merging and) sorting of bam files with samtools and polishing with Pilon.

## Dependencies

`wtdbg2-racon-pilon.pl` needs the following perl modules and will search for executables in your `$PATH`:

__Perl modules__  
- [Path::Tiny](https://metacpan.org/pod/Path::Tiny)

__Assembly__  
- [wtdbg2](https://github.com/ruanjue/wtdbg2): `wtdbg2`, `wtpoa-cns`

__Long read polishing__  
- [Minimap2](https://github.com/lh3/minimap2): `minimap2`
- [Racon](https://github.com/isovic/racon): `racon`

__Short read polishing__  
- Java: `java`
- [bwa (mem)](https://github.com/lh3/bwa): `bwa`
- [samtools](https://github.com/samtools/samtools): `samtools`
- [Pilon](https://github.com/broadinstitute/pilon): Specify the path to Pilon `jar` file with `-pilon-path`

## Usage

```
wtdbg2-racon-pilon.pl [-l <longreads.fq> -x <longread tech>] [-a <assembly.fa>]
                      [-p <paired_1.fq>,<paired_2.fq> -u <unpaired.fq>]

Mandatory (if -a is not set):
	-l STR			File containing long reads in fastq format
	-x STR			Long read sequencing technology for wtdbg2 presets
				Valid arguments are: rsII, rs, sequel, sq, nanopore, ont,
				corrected, ccs
Mandatory (if -l and -x is not set):
	-a STR			File containing an assembly that should be polished in fasta format

Input/pipeline options: [default]
	-racon-rounds INT	Number of racon iterations [3]
	-pilon-rounds INT	Number of pilon iterations [3]
	-pilon-path STR		Complete path to the pilon jar file
	-pr STR			File in bed format to restrict pilon polishing to certain regions
				[polish all positions]
	-p STR			Two files with paired short reads in fastq format comma sperated
				Can be specified multiple times
	-u STR			One file with unpaired short reads in fastq format
				Can be specified multiple times
	-xmx STR/INT		Maximum java heap size for pilon [current available]
	-t INT			Number of parallel executed processes [1]
				Affects wtdbg2, wtpoa-cns, minimap2, racon, bwa mem, samtools sort
				and pilon.
	Pass specific options to tools with:
	-wtdbg-opts [], -wtpoa-opts [], -minimap-opts [], -racon-opts [-u], -bwa-opts [-a -c 10000],
	-pilon-opts [--diploid]
	For example like -wtdbg-opts '-g 100m'

Output options: [default]
	-o STR			Output directory [.]
				Will be created if not existing
	-pre STR		Prefix of output files [{w|<-a>.}r<-racon-rounds>p<-pilon-rounds>]
	-v			Print executed commands to STDERR [off]
	-dry-run		Only print commands to STDERR instead of executing [off]

	-h or -help		Print this help and exit
	-version		Print version number and exit
```

## Citation
__If you use this tool please cite the dependencies as well:__

- wtdbg2:  
Ruan J, Li H (2019). Fast and accurate long-read assembly with wtdbg2. _bioRxiv_, <https://www.biorxiv.org/content/10.1101/530972v1>
- Minimap2:  
Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. _Bioinformatics_, 34:3094–3100, <https://doi.org/10.1093/bioinformatics/bty191>
- Racon:  
Vaser R, Sović I, Nagarajan N, Šikić M (2017). Fast and accurate de novo genome assembly from long uncorrected reads. _Genome research_, 27(5):737–746, <https://dx.doi.org/10.1101%2Fgr.214270.116>
- bwa mem:  
Li H (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. _arXiv preprint arXiv:1303.3997_, <https://arxiv.org/abs/1303.3997>
- samtools:  
Li H, Handsaker B, Wysoker A, Fennell T, Ruan J et al. (2009). The Sequence Alignment/Map format and SAMtools. _Bioinformatics_, 25(16):2078–2079, <https://doi.org/10.1093/bioinformatics/btp352>
- Pilon:  
Walker BJ, Abeel T, Shea T, Priest M, Abouelliel A et al. (2014). Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. _PLOS ONE_, 9(11):e112963, <https://doi.org/10.1371/journal.pone.0112963>

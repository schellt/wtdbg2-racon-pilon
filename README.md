# wtdbg2-racon-pilon.pl v0.3

## Description
__Automatic execution of (wtdbg2,) racon and pilon.__

Longread assembly will be executed with wtdbg2.
Iterative correction with long reads by mapping to the assembly with minimap2 and polishing with Racon.
Iterative correction with short reads by mapping with either bwa mem or NextGenMap, (merging and) sorting of bam files with samtools and polishing with Pilon.

## Dependencies

`wtdbg2-racon-pilon.pl` will search for the following executables in your `$PATH`:

__Assembly__
- [wtdbg2](https://github.com/ruanjue/wtdbg2): `wtdbg2`, `wtpoa-cns`

__Long read polishing__
- [Racon](https://github.com/isovic/racon): `racon`
- [Minimap2](https://github.com/lh3/minimap2): `minimap2`

__Short read polishing__
- Java: `java`
- [bwa (mem)](https://github.com/lh3/bwa): `bwa`
- [NextGenMap](http://cibiv.github.io/NextGenMap/): `ngm`
- [samtools](https://github.com/samtools/samtools): `samtools`
- [Pilon](https://github.com/broadinstitute/pilon): Specify the path to Pilon `jar` file with `-pilon-path`

## Usage

```
wtdbg2-racon-pilon.pl {-l <longreads.fq> -x <longread tech> | -a <assembly.fa>}
                      [-p <paired_1.fq>,<paired_2.fq> -u <unpaired.fq>]

Mandatory:
	-l STR			File containing long reads in fastq format
	-x STR			Long read sequencing technology for wtdbg2 presets
				Valid arguments are: rsII, rs, sequel, sq, nanopore, ont,
				corrected, ccs
	OR
	-a STR			File containing an assembly that should be polished in fasta format

Input/pipeline options: [default]
	-racon-rounds INT	Number of racon iterations [3]
	-pilon-rounds INT	Number of pilon iterations [3]
	-pilon-path STR		Complete path to the pilon jar file
	-sr-mapper STR		Short read mapper for correction with pilon [bwa]
				Valid arguments are: bwa, ngm
				To execute different mappers in different pilon rounds supply a
				comma separated list of mappers. Overrides -pilon-rounds.
	-p STR			Two files with paired short reads in fastq format comma sperated
				Can be specified multiple times
	-u STR			One file with unpaired short reads in fastq format
				Can be specified multiple times
	-xmx STR/INT		Maximum java heap size for pilon [current available]
	-t INT			Number of parallel executed processes [1]
				Affects wtdbg2, wtpoa-cns, minimap2, racon, bwa mem, ngm,
				samtools sort and pilon.
	Pass specific options to tools with:
	-wtdbg-opts, -wtpoa-opts, -minimap-opts, -racon-opts, -bwa-opts, -ngm-opts, -pilon-opts
	For example like -wtdbg-opts '-g 100m'

Output options: [default]
	-o STR			Output directory [.]
				Will be created if not existing
	-pre STR		Prefix of output files [{w|<-a>.}r<racon rounds>p<pilon rounds>]
	-v			Print executed commands to STDERR [off]
	-dry-run		Only print commands to STDERR instead of executing [off]

	-h or -help		Print this help and exit
	-version		Print version number and exit
```

## Citation
If you use this tool please cite the dependencies.

## Known issues
In some cases Pilon can not polish the assembly with bam files created by NextGenMap.

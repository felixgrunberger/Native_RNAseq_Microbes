Nanopore-based native RNA sequencing provides insights into prokaryotic
transcription, operon structures, rRNA maturation and modifications
================
Felix Grünberger<sup>1</sup>, Robert Knüppel<sup>2</sup>, Michael
Jüttner<sup>2</sup>, Martin Fenk<sup>1</sup>, Andreas
Borst<sup>3</sup>, Robert Reichelt<sup>1</sup>, Winfried
Hausner<sup>1</sup>, Jörg Soppa<sup>3</sup>,
<a href="https://orcid.org/0000-0002-0522-843X">Sébastien
Ferreira-Cerca<sup>2°</sup></a>, and
<a href="https://orcid.org/0000-0002-0570-2517">Dina
Grohmann<sup>1°</sup></a>  

<sup>1</sup> Department of Biochemistry, Genetics and Microbiology,
Institute of Microbiology, Single-Molecule Biochemistry Lab &
Biochemistry Centre Regensburg, University of Regensburg,
Universitätsstraße 31, 93053 Regensburg, Germany

<sup>2</sup> Biochemistry III – Institute for Biochemistry, Genetics and
Microbiology, University of Regensburg, Universitätsstraße 31, 93053
Regensburg, Germany.

<sup>3</sup> Goethe University, Institute for Molecular Biosciences,
Max-von-Laue-Str. 9, D-60438, Frankfurt, Germany

<sup>°</sup> Corresponding authors

  - [Information about this
    repository](#information-about-this-repository)
  - [Data generation](#data-generation)
  - [Data analysis](#data-analysis)
      - [Demultiplexing using
        <span>`poreplex`</span>](#demultiplexing-using-poreplex)
      - [Basecalling using `guppy`](#basecalling-using-guppy)
      - [Mapping using <span>`minimap2`</span>](#mapping-using-minimap2)
  - [Data availability](#data-availability)
      - [Raw sequencing files](#raw-sequencing-files)
      - [Additional data](#additional-data)
  - [License](#license)
  - [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

-----

## Information about this repository

This is the repository for the manuscript “Nanopore-based native RNA
sequencing provides insights into prokaryotic transcription, operon
structures, rRNA maturation and modifications”. It contains a
description of the bioinformatical tools used to process native RNA
sequencing data and the downstream analysis mostly based on custom
[Rscripts](Rscripts).

The repository is currently actively developed.

[![Active
Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)

## Data generation

Libraries for Nanopore sequencing were prepared from poly(A)-tailed RNAs
according to the SQK-RNA001 Kit protocol (Oxford Nanopore, Version:
DRS\_9026\_v1\_revP\_15Dec2016) with minor modifications for barcoded
libraries. In this case, Agencourt AMPure XP magnetic beads (Beckman
Coulter) in combination with 1 µl of RiboGuard RNase Inhibitor (Lucigen)
were used instead of the recommended Agencourt RNAclean XP beads to
purify samples after enzymatic reactions. For the barcoded libraries,
the RTA adapter was replaced by custom adapters described in
<https://github.com/hyeshik/poreplex> and reverse transcription (RT) was
performed in individual tubes for each library. After RT reactions, cDNA
was quantified using the Qubit DNA HS assay kit (Thermo Fisher
Scientific) and equimolar amounts of DNA for the multiplexed samples
were used in the next step for ligation of the RNA Adapter (RMX) in a
single tube. Subsequent reactions were performed according to the
protocols recommended by ONT. The libraries were sequenced on a MinION
using R9.4 flow cells and subsequently, FAST5 files were generated using
the recommended script in MinKNOW.

## Data analysis

Before starting a sequencing run, different running options can be
selected in MinKNOW. In case reads are stored in multi-FAST5-containing
files (e.g. 4000 reads per file), files can be converted to single-read
FAST5 files using the
[ont\_fast5\_api](https://github.com/nanoporetech/ont_fast5_api), as
some workflows
(e.g. [`nanopolish`](https://nanopolish.readthedocs.io/en/latest/) and
[`tombo`](https://nanoporetech.github.io/tombo/)) rely on single-FAST5
files for further analysis.  
After a run, reads are stored in two folders (*fast5\_failed*,
*fast5\_passed*). To prevent actual good reads from beeing discarded we
**included all reads from both folders** in the following steps of the
analysis.  
First, we converted multi-FAST5-files with the `multi_to_single_fast5`
command from the
[ont\_fast5\_api](https://github.com/nanoporetech/ont_fast5_api):

``` bash
# set input_path, save_path, search in all folders for fast5 files, set number of threads
multi_to_single_fast5 \
    --input_path <(path) folder containing multi_read_fast5 files> \
    --save_path <(path) to folder where single_read fast5 files will be output> \
    --recursive <recursively search sub-directories> \
    --threads <number of CPU threads to use>
```

The output will be single-read FAST5 files in the *save\_path* folder
with one subfolder per multi-read input file.

### Demultiplexing using [`poreplex`](https://github.com/hyeshik/poreplex)

Multiplexed libraries (how to is described here:
<https://github.com/hyeshik/poreplex>) can be demultiplexed using
[`poreplex`](https://github.com/hyeshik/poreplex). Following this
approach four direct RNA sequencing libraries can be barcoded, pooled
and sequenced together.  
`Poreplex` can demultiplex the libraries into separate folders with:

``` bash
# trim adapters, basecall using albacore, de-multiplex, create symbolic fast5 link, increase number of working processes, sort reads to folders according to barcodes
poreplex \
    -i <path/to/fast5> \
    -o <path/to/output> \
    --trim-adapter <trim 3′ adapter sequences from FASTQ outputs> \
    --barcoding <sort barcoded reads into separate outputs> \
    --basecall <call the ONT albacore for basecalling on-the-fly> \
    --symlink-fast5 <create symbolic links to FAST5 files in output directories even when hard linking is possible> \
    --parallel <number of worker processes>
```

Reads can be basecalled automatically during demultiplexing using
`albacore`, the outdated basecaller of ONT. As we observed dramatic
differences in the detected read qualities and number of reads that can
be mapped, we chose to only sort the reads based on `poreplex` and used
`guppy` (Version 3.0.3) for basecalling. In this case the `poreplex` can
be shortened to: `poreplex -i <input> -o <output> --barcoding
--parallel`.

### Basecalling using `guppy`

Demultiplexed raw FAST5 files and other raw FAST5 files that have not
been barcoded (raw MinKNOW output) can be basecalled (*translated* in
FASTQ data) using `guppy`, the ONT-developed basecaller (available in
the ONT Community). We used version 3.0.3 for basecalling of all of our
reads:

``` bash
ont-guppy-cpu/bin/guppy_basecaller \
    --flowcell FLO-MIN106 \
    --kit SQK-RNA001 \
    --input $input \
    --save_path $output \
    --recursive \
    --reverse_sequence yes \
    --hp_correct 1 \
    --disable_pings 1 \
    --enable_trimming 0 \
    --cpu_threads_per_caller 4 \
    --calib_detect
```

### Mapping using [`minimap2`](https://github.com/lh3/minimap2)

(Li [2018](#ref-Li2018))

## Data availability

### Raw sequencing files

The raw sequencing data in fast5 format will be submitted to the NCBI
sequence read archive
(<a href="https://www.ncbi.nlm.nih.gov/sra">SRA</a>) under BioProject
accession number XXX.

### Additional data

-----

## License

This project is under the general MIT License - see the
[LICENSE](LICENSE) file for details

## References

<div id="refs" class="references hanging-indent">

<div id="ref-Li2018">

Li, Heng. 2018. “Minimap2: Pairwise alignment for nucleotide sequences.”
*Bioinformatics*. <https://doi.org/10.1093/bioinformatics/bty191>.

</div>

</div>

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

<sup>3</sup>

<sup>°</sup> Corresponding authors

  - [Information about this
    repository](#information-about-this-repository)
  - [Data generation](#data-generation)
  - [Data analysis](#data-analysis)
      - [Demultiplexing using
        `poreplex`](#demultiplexing-using-poreplex)
      - [Basecalling using `guppy`](#basecalling-using-guppy)
  - [Data availability](#data-availability)
      - [Raw sequencing files](#raw-sequencing-files)
      - [Additional data](#additional-data)
  - [License](#license)

<!-- README.md is generated from README.Rmd. Please edit that file -->

-----

## Information about this repository

This is the repository for the manuscript “Nanopore-based native RNA
sequencing provides insights into prokaryotic transcription, operon
structures, rRNA maturation and modifications”

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

### Demultiplexing using `poreplex`

### Basecalling using `guppy`

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

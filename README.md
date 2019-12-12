-   [Nanopore-based native RNA sequencing provides insights into
    prokaryotic transcription, operon structures, rRNA maturation and
    modifications](#nanopore-based-native-rna-sequencing-provides-insights-into-prokaryotic-transcription-operon-structures-rrna-maturation-and-modifications)
    -   [**Authors:** Felix Grünberger<sup>1</sup>, Robert
        Knüppel<sup>2</sup>, Michael Jüttner<sup>2</sup>, Martin
        Fenk<sup>1</sup>, Andreas Borst<sup>3</sup>, Robert
        Reichelt<sup>1</sup>, Winfried Hausner<sup>1</sup>, Jörg
        Soppa<sup>3</sup>, Sébastien Ferreira-Cerca<sup>2°</sup>, and
        Dina
        Grohmann<sup>1°</sup>](#authors-felix-grünberger1-robert-knüppel2-michael-jüttner2-martin-fenk1-andreas-borst3-robert-reichelt1-winfried-hausner1-jörg-soppa3-sébastien-ferreira-cerca2-and-dina-grohmann1)
-   [exampleRPackage](#examplerpackage)
    -   [1](#section)
    -   [2](#section-1)
    -   [3](#section-2)
-   [B](#b)
    -   [B1](#b1)
-   [Motivation](#motivation)
    -   [Why adopt a common standard?](#why-adopt-a-common-standard)
    -   [Which standard to adopt?](#which-standard-to-adopt)
-   [How to Create an R Package](#how-to-create-an-r-package)
    -   [Create a new R package with R
        Studio](#create-a-new-r-package-with-r-studio)
-   [Further Reading](#further-reading)
    -   [Online Resources](#online-resources)
    -   [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->
Nanopore-based native RNA sequencing provides insights into prokaryotic transcription, operon structures, rRNA maturation and modifications
===========================================================================================================================================

#### **Authors:** Felix Grünberger<sup>1</sup>, Robert Knüppel<sup>2</sup>, Michael Jüttner<sup>2</sup>, Martin Fenk<sup>1</sup>, Andreas Borst<sup>3</sup>, Robert Reichelt<sup>1</sup>, Winfried Hausner<sup>1</sup>, Jörg Soppa<sup>3</sup>, Sébastien Ferreira-Cerca<sup>2°</sup>, and Dina Grohmann<sup>1°</sup>

<sup>1</sup> Department of Biochemistry, Genetics and Microbiology,
Institute of Microbiology, Single-Molecule Biochemistry Lab &
Biochemistry Centre Regensburg, University of Regensburg,
Universitätsstraße 31, 93053 Regensburg, Germany  
<sup>2</sup> Biochemistry III – Institute for Biochemistry, Genetics and
Microbiology, University of Regensburg, Universitätsstraße 31, 93053
Regensburg, Germany. <sup>3</sup>

<sup>°</sup> Corresponding authors

exampleRPackage
===============

1
-

2
-

3
-

B
=

B1
--

exampleRPackage is an example R package [available on
GitHub](https://github.com/mvuorre/exampleRPackage).

This is the Git(Hub) repository of an example R package. In our
manuscript (not yet available), we describe why and how researchers
might choose to share their research products[1] as R packages. This
repository is the example described in the manuscript, and can be viewed
online for details of the implementation (that is, the R package’s
source code). The exampleRPackage can also be installed from github
(although as an example package it does not contain anything useful):

    # install.packages("devtools")
    #devtools::install_github("mvuorre/exampleRPackage")

It is also permanently stored on OSF. To install from OSF:

    #temporary_file <- tempfile(fileext = ".tar.gz")
    #download.file("https://osf.io/mqd6f/download", destfile = temporary_file)
    #install.packages(temporary_file, repos = NULL)

The file you are reading now is the package’s README, which describes
how to create R packages with functions, data, and appropriate
documentation. In writing this online tutorial, we relied heavily on
Hadley Wickham’s “R Packages”, which is an excellent source of
information on creating R packages \[@wickham\_r\_2015\].

------------------------------------------------------------------------

Motivation
==========

Lack of reproducibility has been identified as a key limiting factor to
the reliability of scientific research, and researchers are now urged to
focus on factors in their workflow that could improve the
reproducibility of their results
\[@Munafomanifestoreproduciblescience2017\].

> “A research project is computationally *reproducible* if a second
> investigator (including you in the future) can recreate the final
> reported results of the project, including key quantitative findings,
> tables, and figures, given only a set of files and written
> instructions.” \[@KitzesPracticeReproducibleResearch2017\]

The role of standards and common practices has always been important in
ensuring the continuity of one’s work, especially in areas where results
depend on computational work. For example, the BIDS (Brain Imaging Data
Structure) has been introduced as a standardized organization for brain
imaging to facilitate collaborative work on very large data sets and
computationally demanding projects \[@gorgolewski\_brain\_2016\].
However, similar standards have not been established for less
computationally intense areas of behavioral science (although there have
been some attempts, for example the Tier protocol\[^tier\]). \[^tier\]:
<a href="http://www.projecttier.org/tier-protocol/specifications/" class="uri">http://www.projecttier.org/tier-protocol/specifications/</a>

Why adopt a common standard?
----------------------------

Quoting @gorgolewski\_brain\_2016:

-   Minimized curation: Common standards make it possible for
    researchers who were not directly involved in data collection to
    understand and work with the data. This is particularly important to
    ensure that data remain accessible and usable by different
    researchers over time in the following instances:
    -   within a laboratory over time
    -   between labs facilitating collaboration and making combining
        data in multi-center studies easier and less ambiguous
    -   between public databases (i.e., OpenfMRI) allowing for the quick
        ingestion of big data organized according to a common scheme.
-   Error reduction: Errors attributed to the misunderstanding of the
    meaning of a given datum (e.g., when variable names are not
    explicitly stated in the data file and standardized across files).
-   Optimized usage of data analysis software is made possible when the
    metadata necessary for analysis (i.e., details of the task or
    imaging protocol) are easily accessible in a standardized and
    machine- readable way. This enables the application of completely
    automated analysis workflows, which greatly enhances reproducibility
    and efficiency.
-   Development of automated tools for verifying the consistency and
    completeness of datasets is realized. Such tools make it easier to
    spot missing metadata that limit how the data could be analyzed in
    the future.

Which standard to adopt?
------------------------

Instead of suggesting yet another arbitrary standard, we propose that
behavioral scientists could adopt a well-established standard from
statistical software development for their project organizing
principles. Specifically, in this tutorial we describe how to organize
and share one’s data sets, functions, and related materials including
analyses, as packages for the R programming language
\[@RCoreTeamLanguageEnvironmentStatistical2017\].

Why R? R is an easily accessible programming language for statistical
computing and graphics, and is rapidly increasing in popularity in the
behavioral sciences. For example, the APS Observer recently ran a series
of articles promoting R within the psychological science community
\[@yee\_why\_2017\]. Importantly, users can create R packages with no or
minimal coding, because many of the procedures have been implemented in
the RStudio \[@RStudioTeamRStudioIntegratedDevelopment2016\] graphical
interface.

How to Create an R Package
==========================

The outline of the tutorial is as follows:

1.  [Create a new R package with R
    Studio](#create-a-new-r-package-with-r-studio)
    -   With a few button clicks, this automatically sets up the
        underlying software infrastructure
2.  [Describe the package](#describe-the-package)
    -   DESCRIPTION and README files
3.  [Add data to package](#add-data)
    -   Raw data, preprocessing scripts, R data object
4.  [Create and add functions](#create-functions)
    -   \[todo\]
5.  [Document the package](#document-the-package)
    -   Describe the package, its functions, and data, in a machine- and
        human-readable format

After these simple steps, you will have a functional R package on your
computer. We will also go through advanced (optional) steps.

-   [Sharing the R package](#sharing-the-R-package)
    -   Upload it to GitHub so it is easily available to anyone (R user
        or otherwise)
    -   Mint a DOI for citeability and longevity (todo)
    -   Connect to Open Science Framework (todo)
-   [Document data analysis as package
    vignette](#documenting-analysis-as-package-vignette)
    -   Creates a readable file showing how to use the package. For
        example the vignette can describe how the data was (or could be)
        analyzed
-   [Create a website for the
    package](#creating-a-website-for-the-data-package)
    -   Showcase your R package online with a website

You will need one R package (R developer tools) to follow these
instructions:

    #install.packages("devtools")

The **devtools** package \[@wickham\_devtools:\_2017\] contains helpful
functions for creating R packages.

Create a new R package with R Studio
------------------------------------

First, use R Studio to create a new R Project. While creating the
project, make sure to create the project as an R Package:

Creating an R (Package) Project with R Studio sets up the necessary
infrastructure leaving little work for the user. After creating the
package, the project’s files and folders look like this
(`exampleRPackage` is the project’s root folder):

    exampleRPackage/
    ├── man/
    |   └── hello.Rd
    ├── R/
    |   └── hello.R
    ├── DESCRIPTION
    ├── NAMESPACE
    ├── exampleRPackage.Rproj
    ├── .gitignore
    └── .Rbuildignore

`man/` is the “manuals” folder which will have files documenting the
package. `R/` is a folder for R functions. `DESCRIPTION` is a file
describing the package, and `NAMESPACE` its functions.
`exampleRPackage.Rproj` identifies the folder as an R package project.
`.gitignore` and `.Rbuildignore` are hidden files, and specify which
files should be ignored by Git \[@vuorre\_curating\_2017\], and R
package building operations, respectively. These last three files can be
ignored for now.

At this point, you can delete `man/hello.Rd` and `R/hello.R`. These two
files are examples of R function files and R documentation files.

This is already a fully functional R package (but it contains nothing so
it’s useless.) We now need to introduce content, and change some of the
included files, to turn it into an R package.

Further Reading
===============

Online Resources
----------------

-   <a href="http://r-pkgs.had.co.nz/" class="uri">http://r-pkgs.had.co.nz/</a>:
    Website of Hadley Wickham’s R Packages book
-   [Writing an R package from
    scratch](https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/):
    A short and good blog post on how to create minimal R packages
-   [Writing R
    Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html):
    The official R documentation on writing R packages. This is the
    complete and definitive set of instructions on how to write R
    packages. It is almost unreadable in it’s comprehensiveness, and
    unnecessary for small R packages such as the data package described
    here.

References
----------

[1] By “product”, we mean any combination of text (manuscripts), code,
data, stimuli, and other research materials.

*** Resubmission
This is a resubmission. In this corrected version I have:

* Replaced directed quotes in the Description file and the vignettes by undirected quotes.


#####

> devtools::check_rhub()
── R CMD build ────────────────────────────────────
✔  checking for file 'C:\Users\CWERNER\OneDrive - CIMMYT\Desktop\EiB\projects\FieldSimR\fieldsimr/DESCRIPTION' (1.2s)
─  preparing 'FieldSimR': (39.3s)
✔  checking DESCRIPTION meta-information ... 
─  installing the package to build vignettes
✔  creating vignettes (45.7s)
─  checking for LF line-endings in source and make files and shell scripts (25.4s)
─  checking for empty or unneeded directories
─  building 'FieldSimR_1.1.0.tar.gz'
   
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/FieldSimR_1.1.0.tar.gz-6fcc0f4920ba4c19891dccaa447e7dd4
   https://builder.r-hub.io/status/FieldSimR_1.1.0.tar.gz-a220a374f51242eab1aa9086f31dac61
   https://builder.r-hub.io/status/FieldSimR_1.1.0.tar.gz-73386ff1ae0e44e0a7202b41d30be67f
─  Build started
─  Creating new user
─  Downloading and unpacking package file
─  Querying package dependencies
─  Installing package dependencies
─  Running R CMD check
   setting _R_CHECK_FORCE_SUGGESTS_ to false
   setting R_COMPILE_AND_INSTALL_PACKAGES to never
   setting R_REMOTES_STANDALONE to true
   setting R_REMOTES_NO_ERRORS_FROM_WARNINGS to true
   setting _R_CHECK_FORCE_SUGGESTS_ to true
   setting _R_CHECK_CRAN_INCOMING_USE_ASPELL_ to true
   'getOption("repos")' replaces Bioconductor standard repositories, see
   'help("repositories", package = "BiocManager")' for details.
   Replacement repositories:
       CRAN: https://cloud.r-project.org
─  using log directory 'C:/Users/USERLsXIVjzHfc/FieldSimR.Rcheck'
─  using R Under development (unstable) (2023-01-14 r83615 ucrt)
─  using platform: x86_64-w64-mingw32 (64-bit) (739ms)
─  R was compiled by (742ms)
       gcc.exe (GCC) 12.2.0
       GNU Fortran (GCC) 12.2.0
─  running under: Windows Server x64 (build 20348) (1.4s)
─  using session charset: UTF-8
─  using option '--as-cran' (737ms)
✔  checking for file 'FieldSimR/DESCRIPTION'
─  this is package 'FieldSimR' version '1.1.0' (1.5s)
─  package encoding: UTF-8
─  checking CRAN incoming feasibility ... [16s] NOTE
   Maintainer: 'Christian Werner <werner.christian@proton.me>'
   
   Possibly misspelled words in DESCRIPTION:
     Phenotypes (2:38)
✔  checking package namespace information
✔  checking package dependencies (2.2s)
✔  checking if this is a source package (3s)
✔  checking if there is a namespace (2.2s)
✔  checking for executable files (2.9s)
✔  checking for portable file names (3.7s)
✔  checking serialization versions (2.2s)
✔  checking for hidden files and directories (2.2s)
✔  checking whether package 'FieldSimR' can be installed
✔  checking installed package size (1.5s)
✔  checking package directory
✔  checking for future file timestamps
✔  checking 'build' directory
✔  checking DESCRIPTION meta-information (1.5s)
✔  checking top-level files
✔  checking for left-over files (735ms)
✔  checking index information
✔  checking package subdirectories (1.4s)
✔  checking R files for non-ASCII characters (736ms)
✔  checking R files for syntax errors (734ms)
✔  checking whether the package can be loaded
✔  checking whether the package can be loaded with stated dependencies (1.4s)
✔  checking whether the package can be unloaded cleanly (735ms)
✔  checking whether the namespace can be loaded with stated dependencies
✔  checking whether the namespace can be unloaded cleanly (731ms)
✔  checking loading without being on the library search path (740ms)
✔  checking use of S3 registration (789ms)
✔  checking dependencies in R code
✔  checking S3 generic/method consistency
✔  checking replacement functions (1.5s)
✔  checking foreign function calls (788ms)
✔  checking R code for possible problems
✔  checking Rd files (798ms)
✔  checking Rd metadata (1.5s)
✔  checking Rd line widths
✔  checking Rd cross-references
✔  checking for missing documentation entries
✔  checking for code/documentation mismatches (1.5s)
✔  checking Rd \usage sections
✔  checking Rd contents
✔  checking for unstated dependencies in examples (739ms)
✔  checking contents of 'data' directory (1.5s)
✔  checking data for non-ASCII characters
✔  checking LazyData (755ms)
✔  checking data for ASCII and uncompressed saves
✔  checking installed files from 'inst/doc' (749ms)
✔  checking files in 'vignettes'
✔  checking examples
✔  checking for unstated dependencies in vignettes
✔  checking package vignettes in 'inst/doc' (1.5s)
─  checking re-building of vignette outputs ... [99s] OK
─  checking PDF version of manual ... [18s] OK
✔  checking HTML version of manual (1.5s)
✔  checking for non-standard things in the check directory
N  checking for detritus in the temp directory
   Found the following files/directories:
     'lastMiKTeXException'
   
─  Done with R CMD check
─  Cleaning up files and user
    

── FieldSimR 1.1.0: NOTE

  Build ID:   FieldSimR_1.1.0.tar.gz-6fcc0f4920ba4c19891dccaa447e7dd4
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  10m 54.8s ago
  Build time: 7m 55.5s

❯ checking CRAN incoming feasibility ... [16s] NOTE
  Maintainer: 'Christian Werner <werner.christian@proton.me>'
  
  Possibly misspelled words in DESCRIPTION:
    Phenotypes (2:38)

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── FieldSimR 1.1.0: CREATED

  Build ID:   FieldSimR_1.1.0.tar.gz-a220a374f51242eab1aa9086f31dac61
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  10m 54.9s ago


── FieldSimR 1.1.0: CREATED

  Build ID:   FieldSimR_1.1.0.tar.gz-73386ff1ae0e44e0a7202b41d30be67f
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  10m 54.9s ago


*** COMMENTS

2 notes are noted. Phenotypes is spelled correctly.

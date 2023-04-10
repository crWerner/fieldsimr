> devtools::check_rhub()
── R CMD build ────────────────────────────────────────────────────────────────────────
✔  checking for file 'C:\Users\CWERNER\OneDrive - CIMMYT\Desktop\EiB\projects\FieldSimR\fieldsimr/DESCRIPTION' ...
─  preparing 'FieldSimR': (10.5s)
✔  checking DESCRIPTION meta-information ... 
─  installing the package to build vignettes
✔  creating vignettes (44.2s)
─  checking for LF line-endings in source and make files and shell scripts (21.1s)
─  checking for empty or unneeded directories
─  building 'FieldSimR_1.1.0.tar.gz'
   
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/FieldSimR_1.1.0.tar.gz-90d5e28af40b438dbc4871ab0792bd9f
   https://builder.r-hub.io/status/FieldSimR_1.1.0.tar.gz-95316becb5d24589b5bd54b998db7e91
   https://builder.r-hub.io/status/FieldSimR_1.1.0.tar.gz-696f0583f14d41308a0178a7d4f4cc3d
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
─  using log directory 'C:/Users/USERcSMqIYnlhN/FieldSimR.Rcheck'
─  using R Under development (unstable) (2023-01-14 r83615 ucrt)
─  using platform: x86_64-w64-mingw32 (64-bit) (4.1s)
─  R was compiled by
       gcc.exe (GCC) 12.2.0
       GNU Fortran (GCC) 12.2.0
─  running under: Windows Server x64 (build 20348)
─  using session charset: UTF-8
─  using option '--as-cran'
✔  checking for file 'FieldSimR/DESCRIPTION'
─  this is package 'FieldSimR' version '1.1.0'
─  package encoding: UTF-8
   Warning in strsplit(lines[ind], ": ", fixed = TRUE) :
     input string 3 is invalid in this locale
   Warning in strsplit(lines[ind], ": ", fixed = TRUE) :
     input string 4 is invalid in this locale
   Warning in strsplit(lines[ind], ": ", fixed = TRUE) :
     input string 5 is invalid in this locale
─  checking CRAN incoming feasibility ... [16s] NOTE (4.4s)
   Maintainer: 'Christian Werner <werner.christian@proton.me>'
   
   Possibly misspelled words in DESCRIPTION:
     Phenotypes (2:38)
✔  checking package namespace information
✔  checking package dependencies
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files
✔  checking for hidden files and directories
✔  checking for portable file names
✔  checking serialization versions
✔  checking whether package 'FieldSimR' can be installed
✔  checking installed package size
✔  checking package directory
✔  checking for future file timestamps (4.7s)
✔  checking 'build' directory
✔  checking DESCRIPTION meta-information
✔  checking top-level files
✔  checking for left-over files
✔  checking index information
✔  checking package subdirectories
✔  checking R files for non-ASCII characters
✔  checking R files for syntax errors
✔  checking whether the package can be loaded
✔  checking whether the package can be loaded with stated dependencies
✔  checking whether the package can be unloaded cleanly
✔  checking whether the namespace can be loaded with stated dependencies
✔  checking whether the namespace can be unloaded cleanly
✔  checking loading without being on the library search path
✔  checking use of S3 registration
✔  checking dependencies in R code
✔  checking S3 generic/method consistency
✔  checking replacement functions
✔  checking foreign function calls
✔  checking R code for possible problems (8.5s)
✔  checking Rd files
✔  checking Rd metadata
✔  checking Rd line widths
✔  checking Rd cross-references
✔  checking for missing documentation entries
✔  checking for code/documentation mismatches
✔  checking Rd \usage sections
✔  checking Rd contents (6.1s)
✔  checking for unstated dependencies in examples
✔  checking contents of 'data' directory
✔  checking data for non-ASCII characters
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
✔  checking examples (4.4s)
✔  checking for unstated dependencies in vignettes
✔  checking package vignettes in 'inst/doc'
─  checking re-building of vignette outputs ... [89s] OK (1m 26.7s)
─  checking PDF version of manual ... [14s] OK (14.3s)
✔  checking HTML version of manual (4.4s)
✔  checking for non-standard things in the check directory
N  checking for detritus in the temp directory
   Found the following files/directories:
     'lastMiKTeXException'
─  Done with R CMD check
─  Cleaning up files and user
    

── FieldSimR 1.1.0: NOTE

  Build ID:   FieldSimR_1.1.0.tar.gz-90d5e28af40b438dbc4871ab0792bd9f
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  7m 5.8s ago
  Build time: 6m 52.9s

❯ checking CRAN incoming feasibility ... [16s] NOTE
  Maintainer: 'Christian Werner <werner.christian@proton.me>'
  
  Possibly misspelled words in DESCRIPTION:
    Phenotypes (2:38)

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── FieldSimR 1.1.0: CREATED

  Build ID:   FieldSimR_1.1.0.tar.gz-95316becb5d24589b5bd54b998db7e91
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  7m 5.9s ago


── FieldSimR 1.1.0: CREATED

  Build ID:   FieldSimR_1.1.0.tar.gz-696f0583f14d41308a0178a7d4f4cc3d
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  7m 5.9s ago


*** COMMENTS

2 notes are noted. Phenotypes is spelled correctly.

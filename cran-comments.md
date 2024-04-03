
- Tests ran without error and 3 notes.
--> Note that there are no misspelled words in DESCRIPTION.

> devtools::check_rhub() # Used for cran-comments.md
── R CMD build ───────────────────────────────────────
✔  checking for file 'C:\Users\CWERNER\OneDrive - CIMMYT\Desktop\ABI\projects\FieldSimR\fieldsimr/DESCRIPTION' ... 
─  preparing 'FieldSimR': (3.3s)
✔  checking DESCRIPTION meta-information ...
─  installing the package to build vignettes
✔  creating vignettes (33s)
─  checking for LF line-endings in source and make files and shell scripts (1.1s)
─  checking for empty or unneeded directories
─  building 'FieldSimR_1.3.0.tar.gz'
   
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/FieldSimR_1.3.0.tar.gz-2bcafa5ed9204e91a173242d9c7aaf2d
   https://builder.r-hub.io/status/FieldSimR_1.3.0.tar.gz-bcf6c227617c4a3fac781be42361185d
   https://builder.r-hub.io/status/FieldSimR_1.3.0.tar.gz-7b18ca3e8a0246eea66f309ba73cc1f5
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
─  using log directory 'C:/Users/USERDItLhgOHMd/FieldSimR.Rcheck'
─  using R Under development (unstable) (2024-03-01 r86033 ucrt)
─  using platform: x86_64-w64-mingw32
─  R was compiled by
       gcc.exe (GCC) 12.3.0
       GNU Fortran (GCC) 12.3.0
─  running under: Windows Server 2022 x64 (build 20348)
─  using session charset: UTF-8
─  using option '--as-cran'
✔  checking for file 'FieldSimR/DESCRIPTION'
─  this is package 'FieldSimR' version '1.3.0' (2.3s)
─  package encoding: UTF-8
─  checking CRAN incoming feasibility ... [13s] NOTE (2.4s)
   Maintainer: 'Christian Werner <werner.christian@proton.me>'
   
   Possibly misspelled words in DESCRIPTION:
     AlphaSimR (15:88)
     GxE (15:28)
✔  checking package namespace information
✔  checking package dependencies
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files
✔  checking for hidden files and directories
✔  checking for portable file names
✔  checking whether package 'FieldSimR' can be installed (7.2s)
✔  checking installed package size
✔  checking package directory
✔  checking for future file timestamps
✔  checking 'build' directory
✔  checking DESCRIPTION meta-information
✔  checking top-level files
✔  checking for left-over files
✔  checking index information
✔  checking package subdirectories (2.5s)
✔  checking R files for non-ASCII characters
✔  checking R files for syntax errors
✔  checking whether the package can be loaded
✔  checking whether the package can be loaded with stated dependencies
✔  checking whether the package can be unloaded cleanly
✔  checking whether the namespace can be loaded with stated dependencies
✔  checking whether the namespace can be unloaded cleanly
✔  checking loading without being on the library search path (2.3s)
✔  checking use of S3 registration
✔  checking dependencies in R code (2.3s)
✔  checking S3 generic/method consistency
✔  checking replacement functions
✔  checking foreign function calls
✔  checking R code for possible problems (7.1s)
✔  checking Rd files (2.6s)
✔  checking Rd metadata
✔  checking Rd line widths
✔  checking Rd cross-references
✔  checking for missing documentation entries
✔  checking for code/documentation mismatches
✔  checking Rd \usage sections
✔  checking Rd contents (2.4s)
✔  checking for unstated dependencies in examples
✔  checking contents of 'data' directory
✔  checking data for non-ASCII characters
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
✔  checking examples (4.8s)
✔  checking for unstated dependencies in vignettes
✔  checking package vignettes
─  checking re-building of vignette outputs ... [48s] OK (48.5s)
─  checking PDF version of manual ... [12s] OK (12.3s)
✔  checking HTML version of manual (2.7s)
N  checking for non-standard things in the check directory
   Found the following files/directories:
     ''NULL''
N  checking for detritus in the temp directory
   Found the following files/directories:
     'lastMiKTeXException'
   
─  Done with R CMD check
─  Cleaning up files and user
    

── FieldSimR 1.3.0: NOTE

  Build ID:   FieldSimR_1.3.0.tar.gz-2bcafa5ed9204e91a173242d9c7aaf2d
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  6m 31s ago
  Build time: 6m 21.3s

❯ checking CRAN incoming feasibility ... [13s] NOTE
  Maintainer: 'Christian Werner <werner.christian@proton.me>'
  
  Possibly misspelled words in DESCRIPTION:
    AlphaSimR (15:88)
    GxE (15:28)

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 3 notes ✖
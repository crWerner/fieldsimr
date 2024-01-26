## rhub::check_for_cran results

> devtools::check_rhub()
── R CMD build ─────────────────────────────────────────────────────────
✔  checking for file 'C:\Users\CWERNER\OneDrive - CIMMYT\Desktop\ABI\projects\FieldSimR\fieldsimr/DESCRIPTION' ...
─  preparing 'FieldSimR': (10.9s)
✔  checking DESCRIPTION meta-information ... 
─  installing the package to build vignettes
✔  creating vignettes (1m 1s)
─  checking for LF line-endings in source and make files and shell scripts (22.1s)
─  checking for empty or unneeded directories
─  building 'FieldSimR_1.2.0.tar.gz'
   
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/FieldSimR_1.2.0.tar.gz-efef2c49c1404c21954a669bfcce2f6f
   https://builder.r-hub.io/status/FieldSimR_1.2.0.tar.gz-e2750b94eb5643e98cbf8c0004addb65
   https://builder.r-hub.io/status/FieldSimR_1.2.0.tar.gz-4b45906d980148e59748f8f4892a554c
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
─  using log directory 'C:/Users/USERcJswXQFVSH/FieldSimR.Rcheck'
─  using R Under development (unstable) (2023-10-14 r85331 ucrt)
─  using platform: x86_64-w64-mingw32
─  R was compiled by
       gcc.exe (GCC) 12.2.0
       GNU Fortran (GCC) 12.2.0
─  running under: Windows Server 2022 x64 (build 20348)
─  using session charset: UTF-8
─  using option '--as-cran' (43.6s)
✔  checking for file 'FieldSimR/DESCRIPTION'
─  this is package 'FieldSimR' version '1.2.0'
─  package encoding: UTF-8
─  checking CRAN incoming feasibility ... [13s] Note_to_CRAN_maintainers
   Maintainer: 'Christian Werner <werner.christian@proton.me>'
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
✔  checking for future file timestamps
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
✔  checking R code for possible problems
✔  checking Rd files
✔  checking Rd metadata
✔  checking Rd line widths
✔  checking Rd cross-references
✔  checking for missing documentation entries
✔  checking for code/documentation mismatches
✔  checking Rd \usage sections
✔  checking Rd contents
✔  checking for unstated dependencies in examples
✔  checking contents of 'data' directory
✔  checking data for non-ASCII characters
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
✔  checking examples
✔  checking for unstated dependencies in vignettes
✔  checking package vignettes in 'inst/doc'
─  checking re-building of vignette outputs ... [90s] OK (1m 18.2s)
─  checking PDF version of manual ... [12s] OK (10.3s)
✔  checking HTML version of manual (10.5s)
N  checking for non-standard things in the check directory
   Found the following files/directories:
     ''NULL''
N  checking for detritus in the temp directory
   Found the following files/directories:
     'lastMiKTeXException'
   
─  Done with R CMD check
─  Cleaning up files and user
    

── FieldSimR 1.2.0: NOTE

  Build ID:   FieldSimR_1.2.0.tar.gz-efef2c49c1404c21954a669bfcce2f6f
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  8m 1.5s ago
  Build time: 7m 3.5s

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖
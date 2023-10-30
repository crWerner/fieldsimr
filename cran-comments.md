## rhub::check_for_cran results

> devtools::check_rhub() 
── R CMD build ─────────────────────────────────────────────────────────────────────────────────
✔  checking for file 'C:\Users\CWERNER\OneDrive - CIMMYT\Desktop\ABI\projects\FieldSimR\fieldsimr/DESCRIPTION' (401ms)
─  preparing 'FieldSimR': (11.5s)
✔  checking DESCRIPTION meta-information ... 
─  installing the package to build vignettes
✔  creating vignettes (50.9s)
─  checking for LF line-endings in source and make files and shell scripts (21.9s)
─  checking for empty or unneeded directories
─  building 'FieldSimR_1.2.0.tar.gz'
   
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/FieldSimR_1.2.0.tar.gz-d32c1bd9de6f4adfb2a0ba01020f6530
   https://builder.r-hub.io/status/FieldSimR_1.2.0.tar.gz-5320021750fe466b800de88744cdbdca
   https://builder.r-hub.io/status/FieldSimR_1.2.0.tar.gz-9b51c455a9674bc59b62e90c663536d2
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
─  using log directory 'C:/Users/USERTBFaDvUrIi/FieldSimR.Rcheck'
─  using R Under development (unstable) (2023-10-14 r85331 ucrt)
─  using platform: x86_64-w64-mingw32
─  R was compiled by (774ms)
       gcc.exe (GCC) 12.2.0
       GNU Fortran (GCC) 12.2.0
─  running under: Windows Server 2022 x64 (build 20348)
─  using session charset: UTF-8 (1.5s)
─  using option '--as-cran'
✔  checking for file 'FieldSimR/DESCRIPTION'
─  this is package 'FieldSimR' version '1.2.0'
─  package encoding: UTF-8 (740ms)
─  checking CRAN incoming feasibility ... [12s] Note_to_CRAN_maintainers
   Maintainer: 'Christian Werner <werner.christian@proton.me>'
✔  checking package namespace information
✔  checking package dependencies (737ms)
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files
✔  checking for hidden files and directories (988ms)
✔  checking for portable file names
✔  checking serialization versions
✔  checking whether package 'FieldSimR' can be installed (4.4s)
✔  checking installed package size (770ms)
✔  checking package directory
✔  checking for future file timestamps
✔  checking 'build' directory
✔  checking DESCRIPTION meta-information
✔  checking top-level files (753ms)
✔  checking for left-over files (770ms)
✔  checking index information
✔  checking package subdirectories
✔  checking R files for non-ASCII characters
✔  checking R files for syntax errors
✔  checking whether the package can be loaded
✔  checking whether the package can be loaded with stated dependencies (1.4s)
✔  checking whether the package can be unloaded cleanly
✔  checking whether the namespace can be loaded with stated dependencies
✔  checking whether the namespace can be unloaded cleanly
✔  checking loading without being on the library search path (754ms)
✔  checking use of S3 registration (747ms)
✔  checking dependencies in R code (1.5s)
✔  checking S3 generic/method consistency (2.9s)
✔  checking replacement functions
✔  checking foreign function calls
✔  checking R code for possible problems (6s)
✔  checking Rd files
✔  checking Rd metadata
✔  checking Rd line widths (788ms)
✔  checking Rd cross-references
✔  checking for code/documentation mismatches
✔  checking for missing documentation entries
✔  checking Rd \usage sections (776ms)
✔  checking Rd contents
✔  checking for unstated dependencies in examples (770ms)
✔  checking contents of 'data' directory
✔  checking data for non-ASCII characters
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
✔  checking examples (4.4s)
✔  checking for unstated dependencies in vignettes (763ms)
✔  checking package vignettes in 'inst/doc'
─  checking re-building of vignette outputs ... [90s] OK (1m 32.2s)
─  checking PDF version of manual ... [12s] OK (11.1s)
✔  checking HTML version of manual (1.5s)
N  checking for non-standard things in the check directory
   Found the following files/directories:
     ''NULL''
N  checking for detritus in the temp directory
   Found the following files/directories:
     'lastMiKTeXException'
   
─  Done with R CMD check
─  Cleaning up files and user
    

── FieldSimR 1.2.0: NOTE

  Build ID:   FieldSimR_1.2.0.tar.gz-d32c1bd9de6f4adfb2a0ba01020f6530
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  7m 17s ago
  Build time: 7m 10.6s

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── FieldSimR 1.2.0: IN-PROGRESS

  Build ID:   FieldSimR_1.2.0.tar.gz-5320021750fe466b800de88744cdbdca
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  7m 17s ago


── FieldSimR 1.2.0: IN-PROGRESS

  Build ID:   FieldSimR_1.2.0.tar.gz-9b51c455a9674bc59b62e90c663536d2
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  7m 17.1s ago

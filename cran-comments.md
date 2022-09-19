ℹ Check <werner.christian@proton.me> for the results in
  15-30 mins (~04:38 PM).
> devtools::check_rhub()
✔  checking for file 'C:\Users\CWERNER\OneDrive - CIMMYT\Desktop\fsr\fieldsimr/DESCRIPTION' (488ms)
─  preparing 'FieldSimR': (6.3s)
✔  checking DESCRIPTION meta-information ... 
─  checking for LF line-endings in source and make files and shell scripts (15.2s)
─  checking for empty or unneeded directories
   Omitted 'LazyData' from DESCRIPTION
─  building 'FieldSimR_1.0.0.tar.gz'
   
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/FieldSimR_1.0.0.tar.gz-4d47c544430d4dd59a2dde54cc50948d
   https://builder.r-hub.io/status/FieldSimR_1.0.0.tar.gz-2cbfd6ae76c6480f9f8055910a785562
   https://builder.r-hub.io/status/FieldSimR_1.0.0.tar.gz-14f269db7ec54f09a7890b6520e0bc77
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
   Error : Bioconductor does not yet build and check packages for R version 4.3; see
     https://bioconductor.org/install
─  using log directory 'C:/Users/USERYbMpzmFPgc/FieldSimR.Rcheck' (813ms)
─  using R Under development (unstable) (2022-08-15 r82721 ucrt) (976ms)
─  using platform: x86_64-w64-mingw32 (64-bit)
─  using session charset: UTF-8 (850ms)
─  using option '--as-cran'
✔  checking for file 'FieldSimR/DESCRIPTION' (1.7s)
─  this is package 'FieldSimR' version '1.0.0'
─  package encoding: UTF-8 (838ms)
─  checking CRAN incoming feasibility ... [14s] NOTE
   Maintainer: 'Christian Werner <werner.christian@proton.me>'
   
   New submission
   
   Possibly misspelled words in DESCRIPTION:
     phenotypes (21:39)
✔  checking package namespace information
✔  checking package dependencies (832ms)
✔  checking if this is a source package (2.5s)
✔  checking if there is a namespace (1.7s)
✔  checking for executable files (836ms)
✔  checking for hidden files and directories
✔  checking for portable file names (1.8s)
✔  checking serialization versions
✔  checking whether package 'FieldSimR' can be installed (852ms)
✔  checking installed package size
✔  checking package directory (860ms)
✔  checking for future file timestamps
✔  checking DESCRIPTION meta-information (1.1s)
✔  checking top-level files
✔  checking for left-over files (2.5s)
✔  checking index information (836ms)
✔  checking package subdirectories (1.2s)
✔  checking R files for non-ASCII characters
✔  checking whether the package can be loaded (850ms)
✔  checking R files for syntax errors
✔  checking whether the package can be loaded with stated dependencies
✔  checking whether the package can be unloaded cleanly (846ms)
✔  checking whether the namespace can be loaded with stated dependencies (900ms)
✔  checking whether the namespace can be unloaded cleanly
✔  checking loading without being on the library search path
✔  checking use of S3 registration
✔  checking dependencies in R code (1.7s)
✔  checking S3 generic/method consistency
✔  checking replacement functions
✔  checking foreign function calls
✔  checking R code for possible problems (844ms)
✔  checking Rd files (829ms)
✔  checking Rd metadata
✔  checking Rd line widths
✔  checking Rd cross-references (839ms)
✔  checking for missing documentation entries (869ms)
✔  checking for code/documentation mismatches
✔  checking Rd \usage sections
✔  checking Rd contents
✔  checking for unstated dependencies in examples
✔  checking examples
─  checking PDF version of manual ... [26s] OK
✔  checking HTML version of manual (930ms)
✔  checking for non-standard things in the check directory
N  checking for detritus in the temp directory
   Found the following files/directories:
     'lastMiKTeXException'
   
─  Done with R CMD check
─  Cleaning up files and user
    

── FieldSimR 1.0.0: NOTE

  Build ID:   FieldSimR_1.0.0.tar.gz-4d47c544430d4dd59a2dde54cc50948d
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  4m 57.1s ago
  Build time: 4m 52.6s

❯ checking CRAN incoming feasibility ... [14s] NOTE
  Maintainer: 'Christian Werner <werner.christian@proton.me>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    phenotypes (21:39)

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── FieldSimR 1.0.0: IN-PROGRESS

  Build ID:   FieldSimR_1.0.0.tar.gz-2cbfd6ae76c6480f9f8055910a785562
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  4m 57.2s ago


── FieldSimR 1.0.0: IN-PROGRESS

  Build ID:   FieldSimR_1.0.0.tar.gz-14f269db7ec54f09a7890b6520e0bc77
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  4m 57.2s ago



### Resubmission
This is a resubmission. As requested, I have:

* Put package and software names in single quotes in title and description.

* Corrected one typo as kindly noted by the reviewer.

* Uncommented the example in rand_cor_mat.Rd which was accidently commented out.






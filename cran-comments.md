
> check()
══ Documenting ═══════════════════════════════════════════════════════════════
ℹ Updating FieldSimR documentation
ℹ Loading FieldSimR

══ Building ══════════════════════════════════════════════════════════════════
Setting env vars:
• CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
• CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
• CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX14FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX17FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX20FLAGS: -Wall -pedantic -fdiagnostics-color=always
── R CMD build ───────────────────────────────────────────────────────────────
✔  checking for file 'C:\... (415ms)
─  preparing 'FieldSimR': (10.9s)
✔  checking DESCRIPTION meta-information ...
─  installing the package to build vignettes
✔  creating vignettes (39.5s)
─  checking for LF line-endings in source and make files and shell scripts (2.7s)
─  checking for empty or unneeded directories
─  building 'FieldSimR_1.4.0.tar.gz'
   
══ Checking ══════════════════════════════════════════════════════════════════
Setting env vars:
• _R_CHECK_CRAN_INCOMING_REMOTE_               : FALSE
• _R_CHECK_CRAN_INCOMING_                      : FALSE
• _R_CHECK_FORCE_SUGGESTS_                     : FALSE
• _R_CHECK_PACKAGES_USED_IGNORE_UNUSED_IMPORTS_: FALSE
• NOT_CRAN                                     : true
── R CMD check ───────────────────────────────────────────────────────────────
─  using log directory 'C:/Users/CWERNER/AppData/Local/Temp/RtmpkdkNUl/file2ef4257b4a23/FieldSimR.Rcheck' (465ms)
─  using R version 4.2.3 (2023-03-15 ucrt)
─  using platform: x86_64-w64-mingw32 (64-bit)
─  using session charset: UTF-8
─  using options '--no-manual --as-cran'
✔  checking for file 'FieldSimR/DESCRIPTION'
─  this is package 'FieldSimR' version '1.4.0'
─  package encoding: UTF-8
✔  checking package namespace information
✔  checking package dependencies (2.1s)
✔  checking if this is a source package ...
✔  checking if there is a namespace
✔  checking for executable files (1.1s)
✔  checking for hidden files and directories ...
✔  checking for portable file names
✔  checking whether package 'FieldSimR' can be installed (5.8s)
✔  checking installed package size ... 
✔  checking package directory
✔  checking for future file timestamps ... 
✔  checking 'build' directory ...
✔  checking DESCRIPTION meta-information (446ms)
✔  checking top-level files
✔  checking for left-over files ...
✔  checking index information ... 
✔  checking package subdirectories (371ms)
✔  checking R files for non-ASCII characters ... 
✔  checking R files for syntax errors ... 
✔  checking whether the package can be loaded ... 
✔  checking whether the package can be loaded with stated dependencies ... 
✔  checking whether the package can be unloaded cleanly ... 
✔  checking whether the namespace can be loaded with stated dependencies ... 
✔  checking whether the namespace can be unloaded cleanly (433ms)
✔  checking loading without being on the library search path (483ms)
✔  checking dependencies in R code (1.2s)
✔  checking S3 generic/method consistency (667ms)
✔  checking replacement functions ... 
✔  checking foreign function calls (433ms)
✔  checking R code for possible problems (8.6s)
✔  checking Rd files (650ms)
✔  checking Rd metadata ... 
✔  checking Rd line widths ... 
✔  checking Rd cross-references ... 
✔  checking for missing documentation entries ... 
✔  checking for code/documentation mismatches (920ms)
✔  checking Rd \usage sections (962ms)
✔  checking Rd contents ... 
✔  checking for unstated dependencies in examples (534ms)
✔  checking contents of 'data' directory ...
✔  checking data for non-ASCII characters ... 
✔  checking LazyData ...
✔  checking data for ASCII and uncompressed saves ... 
✔  checking installed files from 'inst/doc' ... 
✔  checking files in 'vignettes' ...
✔  checking examples (4.1s)
✔  checking for unstated dependencies in vignettes (369ms)
✔  checking package vignettes in 'inst/doc' ...
─  checking re-building of vignette outputs ... [46s] OK (45.8s)
✔  checking for non-standard things in the check directory
✔  checking for detritus in the temp directory
   
   
── R CMD check results ────────────────────────────────── FieldSimR 1.4.0 ────
Duration: 1m 21.7s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
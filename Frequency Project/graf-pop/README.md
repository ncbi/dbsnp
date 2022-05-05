# GrafPop Source Code

GrafPop source code includes C++ programs and Perl scripts. See xxx (to be added) for GrafPop software documentation.

### Make C++ binary `grafpop`
Under the root directory, execute:
```sh
$ make
```

To regenerate the C++ binary after editing the code, execute:
```sh
$ make clean
$ make
```

### Run medium tests

Test scripts and test cases are placed under medium_testing directory. Test cases are saved in `test_manifest.txt`. Perl script test_grafpop.pl is used for manually running these test cases.
1. Make sure environment variable `PATH` includes current directory `.`.
1. If necessary, set environment variable `GARFPATH` to include the directory where GrafPop binary and Perl scripts are located.
1. Under `medium_testing` directory, execute: `test_grafpop.pl test_manifest.txt`.
1. If source code is updated, update `test_manifest.txt` to add new test cases or modify existing cases, and execute `test_grafpop.pl test_manifest.txt 1` to update the baseline.
1. Check the baseline files and make sure they are all correct, then execute `test_grafpop.pl test_manifest.txt` (without the second parameter) again.
1. Make sure all test cases pass. 
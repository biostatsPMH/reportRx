# reportRx 1.03

* Added a `NEWS.md` file to track changes to the package. (May 20, 2021)

## covsum()

* Specify the number of digits for continuous covariates with `digits` (March 10, 2021)
* Display percentages by row with the argument `percentage`  (April 11, 2021)
* Specify if NAs from the main covariate are included with `include_missing` (May 20, 2021)
* Specify if the full sample should be included in the table with`full` (May 20, 2021)
* The following were added/changed on 25 May 2021
* Added functionality to test for small counts in contingency table and perform Fisher.exact if req'd
* If testcont=T will perform unequal variance t-test for two groups
* Added the option to display the statistical tests with `show.tests`
* Added the option of showing both the IQR and the range on separate lines with `all.tests`
* Added the option to exclude levels from association tests `excludeLevels`
* ---stop may 25 changes

## outTable()

* Can now auto-detect output and print accordingly (25 May, 2021)
* Can now specify digits for rounding (28 June, 2021)
* Fixed bug with special characters in caption breaking Latex

## rm_covsum(), rm_uvsum(), rm_mvsum()

- removed documentation and changed to optional arguments to better match p* functions (25 May, 2021)
- fixed bug in rm_uvsum where tableOnly=T returned an error (22 July, 2021)

## forestplot2

* now works correctly with ordinal models from MASS::polr (20 June 2021)
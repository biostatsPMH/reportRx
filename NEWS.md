# reportRx 1.03

* Added a `NEWS.md` file to track changes to the package. (May 20, 2021)

## covsum()

* bug fix to allow all values of cov in a level of maincov to be missing (5 Aug 2021)
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

* Added an argument align to control column alignment (29 Jul 2021)
* Fixed bug with special characters in caption breaking Latex
* Can now specify digits for rounding (28 June, 2021)
* Can now auto-detect output and print accordingly (25 May, 2021)

## nestTable() 

- new function to combine two columns into a single column with one column serving as header rows

## rm_covsum() 
- removed documentation and changed to optional arguments to better match p* functions (25 May, 2021)

## rm_uvsum()

- fixed bug in rm_uvsum where tableOnly=T returned an error (22 July, 2021)

## mvsum()

- bug fix to handle empty factor levels (3 Aug 2021)
- functionality to allow for interaction terms

## forestplot2

* now works correctly with ordinal models from MASS::polr (20 June 2021)

## rmdBibfile

* bug fix with EOF within quoted string and improved messages and warnings (5 Aug 2021)
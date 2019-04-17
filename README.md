# vpowerin. 
Stata module to calculate various voting power indices.
***
<br>

## Description
`vpowerin` implements dynamic programming algorithms (Kurz 2016) to calculate various voting power indices. In addition, `vpowerin` alternatively implements generating functions and methods of direct enumeration (see [options](#options) below).

You may either calculate the Shapley-Shubik power index, the absolute and the standardized Banzhaf index, or all three indices by specifying the respective options.

`vpowerin` also calculates the effective number of parties if required. Finally, `vpowerin` also estimates all possible minimal winning coalitions.

## Installation
You can install the latest version of `vpowerin` by executing the following code:
```Stata
net install vpowerin, from("https://raw.githubusercontent.com/eckerale/vpowerin/master")
```

## Dependencies
`vpowerin` requires the user-written package `moremata` (http://fmwww.bc.edu/RePEc/bocode/m) by Ben Jann (2005) to be installed.

## Syntax
```Stata
[by varlist:] vpowerin id weight [if] [in] [weight] [, ssi banzhaf mwc(newvar) effective gfunction enumeration generate(newvar) quota(integer) noprint]
```
where *id* is a variable which uniquely identifies each player and *weight* is a variable which indicates each player's weight.


## Options
`ssi` calculates the Shapley-Shubik power index.<br>

`banzhaf` calculates the absolute and the standardized Banzhaf index.<br>

`mwc(newvar)` creates a series of indicator variables with *stub* *newvar* for all minimal winning coalitions indicating whether each player is part of the respective minimal winning coalition.<br>

`effective` calculates the effective number of parties.<br>

`gfunction` estimate voting power indices via generating functions (see, e.g., Chessa 2014 for additional information).<br>

`enumeration` estimate voting power indices via method of direct enumeration.<br>

`generate(newvar)` creates a new variable *newvar*_ssi which returns the according value of the Shapley-Shubik power index. If any of the additional options `banzhaf` or `effective` is specified, an additional variable (*newvar_bi_abs*, *newvar_bi_std* or *newvar_eff*, respectively) is generated which returns the corresponding value for each observation.<br>

`quota(integer)` specifies the decision rule. The default option is 50 percent of the weights + 1.<br>

`noprint` suppresses the output.<br>

## Remarks
`vpowerin` requires the data to be in long format. Use the `reshape` command if your data are in wide format. Also note that `vpowerin` uses absolute weights (e.g., number of seats in the legislature).

## Examples
    . bysort cabinet_id: vpowerin party_id seats

    . vpowerin party_id seats, ssi generate(power)

    . vpowerin party_id seats, effective

    . vpowerin party_id seats, banzhaf


## Saved results
`vpowerin` saves the following in r():

**r(results)**  matrix of results

## References
Chessa, M. 2014. A generating functions approach for computing the Public Good index efficiently. *TOP* 22(2): 658-73. Available from https://doi.org/10.1007/s11750-013-0286-8.

Jann, B. 2005. moremata: Stata module (Mata) to provide various functions. Available from http://ideas.repec.org/c/boc/bocode/s455001.html.

Kurz, S. 2016. Computing the power distribution in the IMF. CoRR. Available from http://arxiv.org/abs/1603.01443.

## Author
A. Ecker<br>
Mannheim Centre for European Social Research, University of Mannheim.<br>
Please email to alejandro.ecker@mzes.uni-mannheim.de if you observe any problems.

## How to cite
Thanks for citing this Stata module as follows:<br>
Ecker, Alejandro. 2019. vpowerin: Stata module to calculate various voting power indices. Available from "https://github.com/eckerale/vpowerin".

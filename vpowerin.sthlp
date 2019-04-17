{smcl}
{* *! version 0.9.1 01 April 2019}{...}
{cmd: help vpowerin}
{hline}


{title:Title}

{pstd}
{bf:vpowerin} {hline 2} Calculate various voting power indices.


{title:Syntax}

{pstd}
[by {varlist}:] {cmd:vpowerin} {it:id} {it:weight} {ifin} {weight} [{cmd:,} {cmdab:ssi:} {cmdab:banz:haf} {cmdab:mwc:(}{it:newvar}{cmd:)} {cmdab:eff:ective} {cmdab:gf:unction} {cmdab:enum:eration} {cmdab:gen:erate(}{it:newvar}{cmd:)} 
{cmdab:quo:ta(}{it:integer}{cmd:)} {cmdab:noprint}]

{p 4 4 2}
where {it:id} is a variable which uniquely identifies each player and {it:weight} is a variable which indicates each player's weight.


{title:Installation}

{p 4 4 2}
You can install the latest version of {bf:vpowerin} by executing the following code:

{p 4 4 2}
{cmd:. net install vpowerin, from("https://raw.githubusercontent.com/eckerale/vpowerin/master")}


{title:Dependencies}

{p 4 4 2} 
{cmd:vpowerin} requires the user-written package moremata ({net describe moremata:http://fmwww.bc.edu/RePEc/bocode/m/moremata.html}) by Ben Jann (2005) to be installed. 


{title:Description}

{p 4 4 2} 
{cmd:vpowerin} implements dynamic programming algorithms (Kurz 2016) to calculate various voting power indices. In addition, {cmd:vpowerin} alternatively implements generating functions and methods of direct enumeration (see options below).
You may either calculate the Shapley-Shubik power index, the absolute and the standardized Banzhaf index, or all three indices by specifying the respective options. {cmd:vpowerin} 
also calculates the effective number of parties if required. Finally, {cmd:vpowerin} also estimates all possible minimal winning coalitions. 


{title:Options}

{p 4 4 2} 
{cmdab:ssi:} calculates the Shapley-Shubik power index.

{p 4 4 2} 
{cmdab:banz:haf} calculates the absolute and the standardized Banzhaf index.

{p 4 4 2} 
{cmdab:mwc:(}{it:newvar}{cmd:)} creates a series of indicator variables with {it:stub} {it:newvar} for all minimal winning coalitions indicating whether each player is part of the respective minimal winning coalition.

{p 4 4 2} 
{cmdab:eff:ective} calculates the effective number of parties.

{p 4 4 2} 
{cmdab:gf:unction} estimate voting power indices via generating functions (see, e.g., Chessa 2014 for additional information).

{p 4 4 2} 
{cmdab:enum:eration} estimate voting power indices via method of direct enumeration.

{p 4 4 2} 
{cmd:generate(}{it:newvar}{cmd:)} creates a new variable {it:newvar}_ssi which returns the according value of the Shapley-Shubik power index. If any of the additional options ({cmdab:banz:haf} or {cmdab:eff:ective})
is specified, an additional variable ({it:newvar}_bi_abs, {it:newvar}_bi_std or {it:newvar}_eff, respectively) is generated which returns the corresponding value for each observation.

{p 4 4 2} 
{cmdab:quo:ta(}{it:integer}{cmd:)} specifies the decision rule. The default option is 50 percent of the weights + 1.

{p 4 4 2} 
{cmd:noprint} suppresses the output.


{title:Remarks}

{p 4 4 2}
{cmd:vpowerin} requires the data to be in long format. Use the {help reshape} command if your data are in wide format. Also note that {cmd:vpowerin} uses absolute weights (e.g., number of seats in the legislature).


{title:Examples}

{p 4 4 2}
{cmd:. bysort cabinet_id: vpowerin party_id seats}

{p 4 4 2}
{cmd:. vpowerin party_id seats, ssi generate(power)}

{p 4 4 2}
{cmd:. vpowerin party_id seats, effective}

{p 4 4 2}
{cmd:. vpowerin party_id seats, banzhaf}


{title:Saved results}

{p 4 4 2}
{cmd:vpowerin} saves the following in {cmd:r()}:

{synoptset 12 tabbed}{...}
{synopt:{cmd:r(results)}}matrix of results{p_end}


{title:References}

{p 4 4 2} Chessa, M. 2014. A generating functions approach for computing the Public Good index efficiently. {it:TOP} 22(2): 658-73. Available from {browse "https://doi.org/10.1007/s11750-013-0286-8"}.

{p 4 4 2} Jann, B. 2005. moremata: Stata module (Mata) to provide various functions. Available from {browse "http://ideas.repec.org/c/boc/bocode/s455001.html"}.

{p 4 4 2} Kurz, S. 2016. Computing the power distribution in the IMF. {it:CoRR}. Available from {browse "http://arxiv.org/abs/1603.01443"}.


{title:Author}

{p 4 4 2}
A. Ecker, Mannheim Centre for European Social Research, University of Mannheim. Please email to {browse "mailto:alejandro.ecker@mzes.uni-mannheim.de":alejandro.ecker@mzes.uni-mannheim.de} if you observe any problems.


{title:How to cite}

{p 4 4 2}
Thanks for citing this Stata module as follows: Ecker, Alejandro. 2019. vpowerin: Stata module to calculate various voting power indices. Available from {browse "https://github.com/eckerale/vpowerin"}.


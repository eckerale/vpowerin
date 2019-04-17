* Define name of program *

capture program drop vpowerin
*! vpowerin v0.9.1 AEcker 01april2019 corrected bug in generating functions
* powinabs v0.9.0 AEcker 01february2019 added dynamic programming algorithm
* powinabs v0.8.1 AEcker 30september2015 fixed bug for majority situations and generating functions
* powinabs v0.8.0 AEcker 06june2015 implemented generating function for MWCs
* powinabs v0.7.2 AEcker 28may2015 improved handling of majority situations and parties with weight = 0
* powinabs v0.7.1 AEcker 28may2015 check for integer weights
* powinabs v0.7.0 AEcker 25may2015 implemented generating functions
* powinabs v0.6.0 AEcker 18may2015 added indicator variables for MWCs
* powinabs v0.5.0 AEcker 18august2014 based on absolute number of seats
* powinabs v0.4.1 AEcker 25july2013 deleted error-message
* powinabs v0.4.0 AEcker 25july2013
* powinabs v0.3.0 AEcker 24april2013
* powinabs v0.2.0 AEcker 18october2012 restricted vote-share to [0,1]
* powinabs v0.1.3 AEcker 12september2012 adding effective number of parties
* powinabs v0.1.2 AEcker 27june2012 making programm byable
* powinabs v0.1.1 AEcker 27june2012

program vpowerin, byable(recall) rclass
	version 10.0
	
	* define syntax *
	syntax varlist(min=2 max=2) [if] [in] [fweight] [,SSI BANZhaf MWC(name) EFFective GFunction ENUMeration GENerate(name) QUOta(integer -99) noPRINT]

	* define marksample *
	marksample touse
	
	* tokenize varlist *
	tokenize `varlist'	
	
	* check for dependencies - moremata *
	qui findfile moremata.hlp
	if "`r(fn)'" == "" {
         di as txt "user-written package moremata needs to be installed first;"
         di as txt "use -ssc install moremata- to do that"
         exit 498
	}
	
	* generate/check output variables *
	if "`generate'" != "" & _byindex() == 1 {
		if "`ssi'" != "" {
			qui generate `generate'_ssi = .
		}
		if "`banzhaf'" != "" {
			qui generate `generate'_bi_abs = .
			qui generate `generate'_bi_std = .
		}
		if "`effective'" != "" {
			qui generate `generate'_eff = .
		}
	}
	
	if "`mwc'" != "" & _byindex() == 1 {
		capture unab varmwc: `mwc'?
		if "`varmwc'" != "" {
			display in red "variable(s) `mwc' already exist"
			exit 110
		}
	}
	
	* specify options *
	if "`generate'" != "" {
		if "`ssi'" != "" {	
			local generatevars `generate'_ssi
		}
		if "`banzhaf'" != "" {
			local generatevars `generatevars' `generate'_bi_abs `generate'_bi_std
		}
		if "`effective'" != "" {
			local generatevars `generatevars' `generate'_eff
		}
	}
	
	if "`generate'" == "" {
		local generate nogenerate
	}
	
	if "`ssi'" == "" {
		local ssi nossi
	}

	if "`banzhaf'" == "" {
		local banzhaf nobanzhaf
	}
	
	if "`mwc'" == "" {
		local mwc nomwc
	}	
	
	if "`effective'" == "" {
		local effective noeffective
	}	

	if "`gfunction'" != "" & "`enumeration'" != "" {
		display in red "cannot specify more than one estimation method"
			exit 184
	}	
	
	if "`gfunction'" == "" & "`enumeration'" == "" {
		local method dynprog
	}	

	if "`gfunction'" != "" {
		local method gfunction
	}	
	
	if "`enumeration'" != "" {
		local method enumeration
	}
		
	* call mata routine *
	mata vpowerinmata("`1'", "`2'", "`generate'", "`ssi'", "`banzhaf'", "`mwc'", "`effective'", "`method'", "`touse'", "`generatevars'", `quota', "`print'")
	
end

version 10.0
mata
	void vpowerinmata(string scalar players, string scalar weights, string scalar generate, string scalar ssi, string scalar banzhaf, string scalar mwc, string scalar effective, string scalar method, string scalar touse, string scalar generatevars, real scalar quota, string scalar print) {
		
		// initial definitions and certifications
			// define scalar/vector/matrix of players and weights
			vec_players = st_data(.,players,touse)
			sca_players = rows(vec_players)
			mat_data = ((1::rows(st_data(.,weights,touse))),st_data(.,weights,touse))
			mat_weights = mat_data[.,2]
			
			// define vector with specified options
			vec_options = (ssi != "nossi"),(banzhaf != "nobanzhaf"),(effective != "noeffective")
			
			// identify players with weight = 0
			if (anyof(mat_data[.,2],0)) {
				printf(" {txt}One or more players have a weight of 0. You may want to exclude them to reduce computation time. \n")
				// if `ENUMeration' option NOT specified ignore players with weight = 0
				if (method != "enumeration") {
					mat_data = select(mat_data, mat_data[.,2]:!=0)
				}
			}
		
			// identify players with non-integer weights
			if (allof(mod(mat_data[.,2],1),0)==0) {
				printf(" {txt}One or more players have non-integer weights. Weights rounded to the closest integer. \n")
				mat_data[.,2] = round(mat_data[.,2])
			}
			
			// identify instances with quota <= max(weights)
			if (any(quota:<=mat_data[.,2]) & method == "gfunction" & quota > 0) {
				printf(" {err}Quota smaller than maximum weight. Method via generating functions not applicable. \n")
				exit(error(498))
			}
			
			// identify instances with quota > sum(weights)
			if (quota > sum(mat_data[.,2]) & quota > 0) {
				printf(" {err}Quota larger than sum of weights. \n")
				exit(error(498))
			}			
			
			// identify missing values for weights
			if (colmissing(st_data(.,weights)) > 0) {
				printf(" {err}Fewer weights than there are players. \n")
				exit(error(498))
			}		
			
			// assure at least one power index/MWCs specified
			if (ssi == "nossi" & banzhaf == "nobanzhaf" & effective == "noeffective" & mwc == "nomwc") {
				printf(" {err}Please specify at least one power index/MWCs to be calculated. \n ")
				exit(error(498))
			}
			
			// define seat of median legislator if `QUOta' option NOT specified
			sca_odd = mod(colsum(mat_data[.,2]),2)
			// if even number of seats: 50% + 1
			if(sca_odd==0) sca_medianseat = (colsum(mat_data[.,2])/2)+1
			// if odd number of seats: 50% + 0.5
			if(sca_odd==1) sca_medianseat = (colsum(mat_data[.,2])/2)+0.5
			if (quota==-99) {
				// identify instances with quota < max(weights)
				majority_party_vec = J(1,2,0)
				if (any(sca_medianseat:<=mat_data[.,2])) {
					majority_party_vec = select(mat_data,mat_data[.,2]:>=sca_medianseat)
					// printf(" {err}Quota smaller than maximum weight. Method via generating functions not applicable. \n")
					// exit(error(498))
				}
			}
			
			// define seat of median legislator if `QUOta' option specified and quota larger than/equal to 50% + 1
			if (quota>0 & quota>=sca_medianseat) {
				sca_medianseat = quota
				// identify instances with quota < max(weights)
				majority_party_vec = J(1,2,0)
				if (any(sca_medianseat:<=mat_data[.,2])) {
					majority_party_vec = select(mat_data,mat_data[.,2]:>=sca_medianseat)
					// printf(" {err}Quota smaller than maximum weight. Method via generating functions not applicable. \n")
					// exit(error(498))
				}
			}
			
			// define seat of median legislator if `QUOta' option specified and quota smaller than 50% + 1
			if (quota>0 & quota<sca_medianseat) {
				sca_medianseat = quota
				// by definition no instances with quota < max(weights)
				majority_party_vec = J(1,2,0)
			}
			
			// define output matrix
			vec_out = (ssi!="nossi"),(banzhaf!="nobanzhaf"):*2,(effective!="noeffective")
			stats = J(sca_players,colsum(vec_out'),.)
			statstring = ""
			display = ""
			counter_stats = 1		
		// initial definitions and certifications end
		
		// calculation of power indexes via dynamic programming algorithm
			if (method == "dynprog") {

				// Shapley-Shubik Index
				if (ssi != "nossi") {
					mat_results_ssi = J(sca_players,1,0)
					if (majority_party_vec[1,1] != 0) {
						mat_results_ssi[majority_party_vec[1,1],1] = 1 
					}
					
					if (majority_party_vec[1,1] == 0) {
						// (empty) matrix with sum of weights of all potential winning coalitions including varying numbers of parties
						mat_out = J(sum(mat_data[,2])+1,rows(mat_data)+1,0)
						mat_out[rows(mat_out),cols(mat_out)] = 1						
						// define lower bound for backward counting of coalitions per sum of weights
						vec_lb = J(1,rows(mat_data),.)
						colsum = mm_colrunsum(mat_data[,2])'
						for (i=1;i<=rows(mat_data);++i) {  
							lb = sca_medianseat+1,(sum(mat_data[,2]):-colsum)[i]
							vec_lb[i] = max(lb) + mat_data[i,2]'
						}
						
						// backward counting of coalitions per sum of weights with different cardinality
						for (i=1;i<=rows(mat_data);++i) {  
							lb = vec_lb[i] - mat_data[i,2]'
							ub = sum(mat_data[,2]) + 1 - mat_data[i,2]'
							if (lb<=ub) {
								mat_add = mat_out[vec_lb[i]..sum(mat_data[,2]) + 1, 2..rows(mat_data) + 1]
								mat_out[lb..ub, 1..rows(mat_data)] = mat_out[lb..ub, 1..rows(mat_data)]:+mat_add 
							}
						}			
						// backward counting of coalitions per sum of weights with different cardinality
						
						// coalitions with weight sum x and cardinality s that contain player i
						for (i=1;i<=rows(mat_data);++i) {  
							lb = sum(mat_data[,2]) - mat_data[i,2]' + 1
							ub = sum(mat_data[,2]) + 1
							if (lb<=ub) {
								mat_p_out = J(sum(mat_data[,2])+1,rows(mat_data)+1,0)
								mat_p_out[lb..ub,] = mat_out[lb..ub,]
								for (x=sum(mat_data[,2])-mat_data[i,2]'+1;x>=sca_medianseat+1;x--) {
									for (j=rows(mat_data);j>=1;--j) {
										mat_p_out[x,j] = mat_out[x,j] - mat_p_out[x+mat_data[i,2]',j+1]
									}
								}
							}
							mat_p_out = mat_p_out[2..sum(mat_data[,2])+1,2..rows(mat_data)+1]
							for (s=0;s<=rows(mat_data)-1;++s) { 
								mat_results_ssi[mat_data[i,1]] = mat_results_ssi[mat_data[i,1]] + (factorial(s) * factorial(rows(mat_data)-s-1) / factorial(rows(mat_data))) * sum(mat_p_out[sca_medianseat..(sca_medianseat + mat_data[i,2]' - 1),s+1])
							} 
							
						}
					}
					// coalitions with weight sum x and cardinality s that contain player i
					
					statstring = statstring + "Shapley-Shubik "
					stats[.,counter_stats] = mat_results_ssi
					display = "%9.3g"
					counter_stats = counter_stats+1				
				}
				// Shapley-Shubik Index
				
				// Normalized and non-normalized Banzhaf Index
				if (banzhaf != "nobanzhaf") {
					mat_results_bi = J(sca_players,2,.)
					if (majority_party_vec[1,1] != 0) {
						mat_results_bi = J(sca_players,2,0)
						mat_results_bi[majority_party_vec[1,1],.] = 1,1
					}				
					
					if (majority_party_vec[1,1] == 0) {
						// (empty) matrix with sum of weights of all potential winning coalitions
						mat_out = J(sum(mat_data[.,2]),1,0)
						mat_out[sum(mat_data[.,2])] = 1
						
						// backward counting of coalitions per sum of weights
						for (i=1;i<=sca_players;++i) { 
							select_vec = J(1,sca_players,1) 
							select_vec[i] = 0
							sca_max = max((sca_medianseat + mat_weights'[i], (sum(mat_data[.,2]) - sum(select(mat_weights', select_vec)))))
							for (x=sca_max;x<=sum(mat_data[.,2]);++x) {
								mat_out[x - mat_weights'[i]] = mat_out[x] + mat_out[x - mat_weights'[i]] 
							}
						}
						// backward counting of coalitions per sum of weights
						
						// coalitions with weight sum x that contain player i
						for (i=1;i<=sca_players;++i) { 
							lb = sum(mat_data[.,2]) - mat_weights'[i] + 1
							ub = sum(mat_data[.,2])
							if (lb<=ub) {
								mat_p_out = J(sum(mat_data[.,2]),1,0)
								mat_p_out[lb..ub] = mat_out[lb..ub]
								for (x=sum(mat_data[.,2])-mat_weights'[i];x>=sca_medianseat;x--) {
									mat_p_out[x] = mat_out[x] - mat_p_out[x+mat_weights'[i]]
								}
								mat_results_bi[i,1] = 1/2^(sca_players-1)*sum(mat_p_out[sca_medianseat..(sca_medianseat + mat_weights'[i] - 1)])
							}
						}
					}
					// coalitions with weight sum x that contain player i
					
					mat_results_bi[,2] = mat_results_bi[,1]:/colsum(mat_results_bi[,1])
					_editmissing(mat_results_bi, 0)
					
					statstring = statstring + " Banzhaf (abs)  Banzhaf (std) "
					stats[.,counter_stats..counter_stats+1] = mat_results_bi
					display = display+"       %9.3g       %9.3g"			
				}
				// Normalized and non-normalized Banzhaf Index
			}
		// calculation of power indexes via dynamic programming algorithm end
		
		// calculation of power indexes via method of direct enumeration
			if (method == "enumeration") {
				
				// Shapley-Shubik Index
				if (ssi != "nossi") {
					counter_ssi = 1																					// counter running from 1 to factorial(number of parties)																	
					info  = cvpermutesetup(mat_data[.,1])															// setup cvpermute-command
					while ((perm=cvpermute(info)) != J(0,1,.)) {
						mat_permute = perm,mat_data[perm,2],runningsum(mat_data[perm,2])							// create matrix with permutation, corresponding seat-share und running sum of seat-shares
						if (counter_ssi==1) {
							vec_ssi_res = select(mat_permute,mat_permute[.,3]:>=sca_medianseat)[1,1]				// vector to identify 'critical party' for each permutation
						}
						else {
							vec_ssi_res = vec_ssi_res\select(mat_permute,mat_permute[.,3]:>=sca_medianseat)[1,1]
						}
						++counter_ssi
					}
			
					statstring = statstring + "Shapley-Shubik "
					stats[.,counter_stats] = mm_freq(vec_ssi_res,1,1::rows(vec_players)):/factorial(rows(mat_data))
					display = "%9.3g"
					counter_stats = counter_stats+1
				}
				// Shapley-Shubik Index
			
				// Normalized and non-normalized Banzhaf Index
				if (banzhaf != "nobanzhaf") {
					counter_bi = 1																				// counter running from 1 to number of winning coalitions
					for (i=1;i<=rows(mat_data);++i) {
						for (k=1;k<=cols(mm_subsets(rows(mat_data),i));++k) {
							if (colsum(mat_data[mm_subsets(rows(mat_data),i)[.,k],2])>=sca_medianseat) {		// identify whether combination is a winning coalition
								if (counter_bi==1) {
									vec_bi_res = select(mm_subsets(rows(mat_data),i)[.,k],((colsum(mat_data[mm_subsets(rows(mat_data),i)[.,k],2])):-mat_data[mm_subsets(rows(mat_data),i)[.,k],2])[.,1]:<sca_medianseat)			// identify critical players
								}
								else {
									vec_bi_res = vec_bi_res\select(mm_subsets(rows(mat_data),i)[.,k],((colsum(mat_data[mm_subsets(rows(mat_data),i)[.,k],2])):-mat_data[mm_subsets(rows(mat_data),i)[.,k],2])[.,1]:<sca_medianseat)	// identify critical players
								}
								++counter_bi
							}
						}	
					}
					
					statstring = statstring + " Banzhaf (abs)  Banzhaf (std) "
					stats[.,counter_stats..counter_stats+1] = mm_freq(vec_bi_res,1,1::rows(vec_players)):/(2^(rows(vec_players)-1)), mm_freq(vec_bi_res,1,1::rows(vec_players)):/(rows(vec_bi_res))
					display = display+"       %9.3g       %9.3g"
				}
				// Normalized and non-normalized Banzhaf Index		
				
			}
		// calculation of power indexes via method of direct enumeration end
		
		// calculation of power indexes via generating functions
			if (method == "gfunction") {
				
				// Shapley-Shubik Index
				if (ssi != "nossi") {
					mat_results_ssi = J(sca_players,1,.)
					if (majority_party_vec[1,1] != 0) {
						mat_results_ssi = J(sca_players,1,0)
						mat_results_ssi[majority_party_vec[1,1],1] = 1 
					}
					
					if (majority_party_vec[1,1] == 0) {
						for (i=1;i<=rows(mat_data);++i) {
							if (i==1) {
								mat_excl_play = J(1,max(mat_data[.,2])+1,0)
							}
							if (i>1) {
								mat_excl_play = mat_excl_play\J(1,max(mat_data[.,2])+1,0)
							}
							mat_excl_play[i,1] = 1
							mat_excl_play[i,mat_data[i,2]+1] = 1
						}
						
						// estimating number of sets in which j players other than i have a sum of weights equal to k
						for (j=1;j<=rows(mat_data);++j) {
							// select players other than i
							vec_deselect = J(1,rows(mat_data),1)
							vec_deselect[1,j] = 0
							mat_select = select(mat_excl_play,vec_deselect')
							mat_select = mat_select[.,1],(mat_select[.,2..cols(mat_select)]:*10)
							
							// number of sets in which j players (rows) have a sum of weights equal to k (columns)
							mat_poly_multi = (1\0),J(2,colsum(mat_data[.,2])-mat_data[j,2], 0)
							
							// consecutively (i.e. by player) multiply multivariate polynomials 
							for(k=1;k<=rows(mat_select)-1;++k) {
								if (k==1) {
									vec_poly = polymult(mat_select[k,.],mat_select[k+1,.])
									if (cols(vec_poly)<cols(mat_poly_multi)) vec_poly = vec_poly,J(1,cols(mat_poly_multi)-cols(vec_poly),0)
									else vec_poly = polymult(mat_select[k,.],mat_select[k+1,.])[1..cols(mat_poly_multi)]
									// update vector for all sets with 1 player
									mat_poly_multi[1,] = (mod(vec_poly, 10^2) :!= 0) :* vec_poly
									// update vector for all sets with 2 players
									mat_poly_multi[2,] = (mod(vec_poly, 10^2) :== 0) :* vec_poly :/ 10
								}
								if (k>1) {
									for (x=k; x>=1; x--) {
										if (x == k) mat_poly_multi = mat_poly_multi\((polymult(mat_poly_multi[x,],mat_select[k+1,.]) :== 10^2) :* 10)[1..cols(mat_poly_multi)]
										if (x <  k) {
											tmp = polymult(mat_poly_multi[x,],mat_select[k+1,.])[1..cols(mat_poly_multi)] :- mat_poly_multi[x,]
											mat_poly_multi[x+1,] = mat_poly_multi[x+1,] :+ (mod(((tmp :== 0) :* mat_poly_multi[x,] :+ tmp), 10^2) :== 0) :* (tmp :/ 10)
										}
										if (x == 1) {
											mat_poly_multi[x,1..cols(mat_select[k+1,.])] = mat_poly_multi[x,1..cols(mat_select[k+1,.])] :+ mat_select[k+1,.]
											mat_poly_multi[1,1] = 1
										}
									}
								}
							}
							mat_poly_multi = mat_poly_multi:/10
							mat_poly_multi[1,1] = 1
							// sum(mat_poly_multi)
							mat_first_prod = ((factorial(1::rows(mat_poly_multi))):*factorial((J(rows(mat_poly_multi),1,rows(mat_data)):+((1::rows(mat_poly_multi)):*-1):+J(rows(mat_poly_multi),1,-1)))):/factorial(rows(mat_data))
							mat_results_ssi[mat_data[j,1],1] = sum(mat_first_prod:*rowsum(mat_poly_multi[.,(sca_medianseat-mat_data[j,2]+1)..(sca_medianseat-1+1)]))
						}
					}
					
					statstring = statstring + "Shapley-Shubik "
					_editmissing(mat_results_ssi, 0)
					stats[.,counter_stats] = mat_results_ssi
					display = "%9.3g"
					counter_stats = counter_stats+1				
				}
				// Shapley-Shubik Index
				
				// Normalized and non-normalized Banzhaf Index	
				if (banzhaf != "nobanzhaf") {
					mat_results_bi = J(sca_players,2,.)
					if (majority_party_vec[1,1] != 0) {
						mat_results_bi = J(sca_players,2,0)
						mat_results_bi[majority_party_vec[1,1],.] = 1,1
					}
					
					if (majority_party_vec[1,1] == 0) {				
						for (i=1;i<=rows(mat_data);++i) {
							if (i==1) {
								mat_excl_play = J(1,max(mat_data[.,2])+1,0)
							}
							if (i>1) {
								mat_excl_play = mat_excl_play\J(1,max(mat_data[.,2])+1,0)
							}
							mat_excl_play[i,1] = 1
							mat_excl_play[i,mat_data[i,2]+1] = 1
						}					
						
						for (j=1;j<=rows(mat_data);++j) {
							vec_deselect = J(1,rows(mat_data),1)
							vec_deselect[1,j] = 0
							mat_select = select(mat_excl_play,vec_deselect')
		
							for(k=1;k<=rows(mat_select)-1;++k) {
								if (k==1) {
									mat_poly = polymult(mat_select[k,.],mat_select[k+1,.])
								}
								if (k>1) {
									mat_poly = polymult(mat_poly,mat_select[k+1,.])
								}
							}
							mat_results_bi[mat_data[j,1],2] = sum(mat_poly[.,(sca_medianseat-mat_data[j,2]+1)..(sca_medianseat-1+1)])
						}
						
						sca_absbi = sum(mat_results_bi[.,2])/(2^(rows(mat_data)-1))
						mat_results_bi[.,2] = mat_results_bi[.,2]:/sum(mat_results_bi[.,2])				
						mat_results_bi[.,1] = mat_results_bi[.,2]:*sca_absbi
						_editmissing(mat_results_bi, 0)
					}
					
					statstring = statstring + " Banzhaf (abs)  Banzhaf (std) "
					stats[.,counter_stats..counter_stats+1] = mat_results_bi
					display = display+"       %9.3g       %9.3g"
				}
				// Normalized and non-normalized Banzhaf Index	
		
			}
		// calculation of power indexes via generating functions end
		
		// Effective number of parties
			if (effective != "noeffective") {
				mat_eff_num_par = J(rows(stats), 1, 1/colsum((mat_data[.,2]:/colsum(mat_data[.,2])):^2))
				
				statstring = statstring + " Eff.# Parties "
				stats[.,cols(stats)] = mat_eff_num_par
				display = display+"       %9.3g"
			}
		// Effective number of parties end	
		
		// Generate Variable(s) 
			if (generate != "nogenerate") {
				st_view(res,.,(tokens(generatevars)),touse)
				res[.,.] = stats
			}
		// Generate Variable(s) end
		
		// Generate MWC indicator-variables
			if (mwc != "nomwc") {
				// sort descending by weights
				mat_data_mwcs = (st_data(.,weights,touse))
				vec_order_desc = order(mat_data_mwcs,-1)				
				mat_weights_mwcs = mat_weights[vec_order_desc,.]
				vec_order_orig = invorder(vec_order_desc)
				
				// select subset of players with weight > 0
				sca_weight_0 = anyof(mat_weights_mwcs,0)
				if (sca_weight_0==1) {
					vec_players_weight_1 = selectindex(mat_weights_mwcs':!=0)
					mat_weights_subset = select(mat_weights_mwcs,mat_weights_mwcs:>0)
					mat_weights_mwcs = mat_weights_subset
				}
								
				// generate vector with prime numbers
				vec_prime = 2,J(1,cols(mat_weights_mwcs')-1,.)
				for (n=2;n<=cols(vec_prime);++n) {
					j = vec_prime[1,n-1] + 1
					sca_max_divisor = ceil(sqrt(j))
					sca_isprime = 0
	
					while (sca_isprime == 0) {
						sca_divisible = 0
						for (k=1;k<=n-1;++k) {
							if (vec_prime[1,k] > sca_max_divisor) {
								continue
								break
							}
							if (mod(j,vec_prime[1,k])==0) {
								sca_divisible = 1
								continue
								break
							}
						}
						if (sca_divisible == 1) {
							j = j + 1
						}
						else {
							vec_prime[1,n] = j
							sca_isprime = 1
						}
					}
				}			
				
				// identify MWCs
				// loop over players: k
				for (k=1;k<=cols(mat_weights_mwcs');++k) {
					// create vector vec_deselect_k to deselect player k
					vec_deselect_k = J(1,cols(mat_weights_mwcs'),1)
					vec_deselect_k[1,k] = 0
					// vector vector_weights_k with weights excluding k
					vector_weights_k = select(mat_weights_mwcs',vec_deselect_k)
					// vector vector_prime_k with primes excluding k
					vector_prime_k = select(vec_prime,vec_deselect_k)
					// scalar q: quota - weight of player k
					sca_q = sca_medianseat - mat_weights_mwcs'[1,k]
					// scalar max_poly: maximum polynomial, i.e. sum of weights of reduced game
					sca_max_poly = rowsum(vector_weights_k)
					// obtain Holler recursive function by looping over players of reduced game
					for (l=1;l<=cols(vector_weights_k);++l) {
						// for first player of reduced game
						if (l==1) {
							// create polynomial vector
							vec_poly = J(1,sca_max_poly+1,0)
							vec_poly[1,1] = 1
							vec_poly[1,vector_weights_k[1,l]+1] = 1
							// create polynomial vector with prime
							vec_poly_prime = vec_poly[.,1],(vec_poly[.,2..cols(vec_poly)]:*vector_prime_k[1,l])
							// create polynomial matrix with prime
							mat_poly_prime = vec_poly_prime
							// create polynomial vector and polynomial vector with prime for l-1
							vec_poly_l_1 = vec_poly
							vec_poly_prime_l_1 = vec_poly_prime
						}
						// for all remaining players of reduced game
						if (l>1) {
							// derive quotient polynomial
							vec_poly_quo = vec_poly[1,1..sca_q]
							vec_poly_quo_prime = vec_poly_prime[1,1..sca_q]
							// derive remainder polynomial
							vec_poly_rem       = (cols(vec_poly)      >sca_q) ? (J(1,cols(vec_poly_quo)      ,0),vec_poly[1,sca_q+1..cols(vec_poly)])             : (vec_poly)
							vec_poly_rem_prime = (cols(vec_poly_prime)>sca_q) ? (J(1,cols(vec_poly_quo_prime),0),vec_poly_prime[1,sca_q+1..cols(vec_poly_prime)]) : (vec_poly_prime)
							// derive multiplicative polynomial
							vec_poly_mul = J(1,vector_weights_k[1,l]+1,0)
							vec_poly_mul[1,1] = 1
							vec_poly_mul[1,vector_weights_k[1,l]+1] = 1
							vec_poly_mul_prime = vec_poly_mul[.,1],(vec_poly_mul[.,2..cols(vec_poly_mul)]:*vector_prime_k[1,l])
							// combine polynomials with primes
							vec_poly_prime = polymult(vec_poly_mul_prime, vec_poly_quo_prime)
							vec_poly_prime = polyadd(vec_poly_prime, vec_poly_rem_prime)
							// combine polynomials
							vec_poly = polymult(vec_poly_mul, vec_poly_quo)
							vec_poly = polyadd(vec_poly, vec_poly_rem)
							// search for MWCs			
							// one possible combination to obtain x seats
							if (rowmax(vec_poly) == 1) mat_poly_prime = vec_poly_prime
							// more than one possible combination to obtain x seats
							else {
								// one additional possibility
								if (all((vec_poly:-vec_poly_l_1):<= 1)) {
									// change from no possibility to one possibility
									mat_poly_prime = mat_poly_prime\J((rowmax(vec_poly)-rows(mat_poly_prime)),cols(mat_poly_prime),0)
									mat_a = selectindex(vec_poly:==1 :& (vec_poly:-vec_poly_l_1):==1)
									mat_poly_prime[1,mat_a] = vec_poly_prime[.,mat_a]			
									// one additional possibility
									mat_a = selectindex(vec_poly:>1 :& (vec_poly:-vec_poly_l_1):==1)
									mat_b = vec_poly[.,mat_a]
									mat_c = mat_a\mat_b
									for (col=1;col<=cols(mat_c);++col) {
										mat_poly_prime[mat_c[2,col],mat_c[1,col]] = vec_poly_prime[.,mat_c[1,col]]-vec_poly_prime_l_1[.,mat_c[1,col]]
									}	
								}
								// at least more than one additional possibility
								if (any((vec_poly:-vec_poly_l_1):> 1)) {
									// change from no possibility to one possibility
									mat_poly_prime = mat_poly_prime\J((rowmax(vec_poly)-rows(mat_poly_prime)),cols(mat_poly_prime),0)
									mat_a = selectindex(vec_poly:==1 :& (vec_poly:-vec_poly_l_1):==1)
									mat_poly_prime[1,mat_a] = vec_poly_prime[.,mat_a]
									// one additional possibility					
									mat_a = selectindex(vec_poly:>1 :& (vec_poly:-vec_poly_l_1):==1)
									mat_b = vec_poly[.,mat_a]
									mat_c = mat_a\mat_b
									for (col=1;col<=cols(mat_c);++col) {
										mat_poly_prime[mat_c[2,col],mat_c[1,col]] = vec_poly_prime[.,mat_c[1,col]]-vec_poly_prime_l_1[.,mat_c[1,col]]
									}
									// more than one additional possibility				
									mat_a = selectindex(vec_poly:>1 :& (vec_poly:-vec_poly_l_1):>1)
									mat_b = vec_poly[.,mat_a]
									mat_c = mat_a\mat_b
									mat_d = mat_c[1,.]:-vector_weights_k[1,l]
									for (col=1;col<=cols(mat_c);++col) {
										// change from more than one possibility by more than one possibility
										if (vec_poly_l_1[.,mat_c[1,col]]>0) {
											for (row=1;row<=mat_c[2,col]-vec_poly_l_1[.,mat_c[1,col]];++row) {
												mat_poly_prime[row+vec_poly_l_1[.,mat_c[1,col]],mat_c[1,col]] = mat_poly_prime[row,mat_d[1,col]]:*vector_prime_k[1,l]
											}
										}
										// change from no possibility by more than one possibility
										if (vec_poly_l_1[.,mat_c[1,col]]==0) {
											for (row=1;row<=mat_c[2,col];++row) {
												mat_poly_prime[row,mat_c[1,col]] = mat_poly_prime[row,mat_d[1,col]]:*vector_prime_k[1,l]
											}							
										}
									}	
								}
							}
							// update polynomial vector and polynomial vector with prime for l-1
							vec_poly_l_1 = vec_poly
							vec_poly_prime_l_1 = vec_poly_prime
						}
					}
					// relevant columns of polynomial vector and polynomial matrix with prime
					vec_poly = (cols(vec_poly)<sca_medianseat) ? J(rows(vec_poly),sca_q,0),vec_poly[.,sca_q+1..cols(vec_poly)] : J(rows(vec_poly),sca_q,0),vec_poly[.,sca_q+1..sca_medianseat]
					mat_poly_prime = (cols(vec_poly)<sca_medianseat) ? mat_poly_prime[.,sca_q+1..cols(vec_poly)] : mat_poly_prime[.,sca_q+1..sca_medianseat]
					// exclude players which are in no MWC	
					if (any(mat_poly_prime:>0)) {
						vec_mwc_prime = select(vec(mat_poly_prime),vec(mat_poly_prime):!=0)
						mat_results_mwc_help = J(rows(vec_mwc_prime),cols(vec_prime),0)
						//"--"
						//vec_mwc_prime
						//vec_prime
						//"--"
						for (player=1;player<=cols(vec_prime);++player) {
							//"player"
							//player
							//"--"
							if (any(mod(vec_mwc_prime,vec_prime[1,player]):==0)) {
							help = selectindex(mod(vec_mwc_prime,vec_prime[1,player]):==0)							
							//help
							mat_results_mwc_help[help,player] = J(rows(help),1,1)
							//mat_results_mwc_help
							}
						}
						mat_results_mwc_help[.,k] = J(rows(mat_results_mwc_help),1,1)
						if (k==1) mat_results_mwc = mat_results_mwc_help
						else mat_results_mwc = mat_results_mwc\mat_results_mwc_help
					}
				}
				if (sca_weight_0==1) {
					mat_mwcs_ind_help = uniqrows(mat_results_mwc)'
					mat_mwcs_ind = J(sca_players,cols(mat_mwcs_ind_help),0)
					mat_mwcs_ind[vec_players_weight_1,.] = mat_mwcs_ind_help
				}
				else mat_mwcs_ind = uniqrows(mat_results_mwc)'
				mat_mwcs_ind = mat_mwcs_ind[vec_order_orig,.]
				// identify MWCs				
				
				varnames = J(1,cols(mat_mwcs_ind),"")
				for (i=1; i<=cols(mat_mwcs_ind); i++) {
					varnames[i] = sprintf("%s%g", mwc, i)
					if (_st_varindex(varnames[i]) == .) {
						(void) _st_addvar("double", varnames[i])
					}
				}
				st_view(mwcs,.,varnames,touse)
				mwcs[.,.] = mat_mwcs_ind
			}
		// Generate MWC indicator-variables	end

		// Print output
			if (print!="noprint" & any(vec_options:!=0)) {
				statsdisplay = cols(stats)
				ll = strlen(statstring)+1
				printf("\n")
				printf("{txt} Quota (seats): %-9.0g \n", 
					sca_medianseat)
				printf("\n")
				printf("{txt} Player (weight) {c |} {col 4} %s\n", 
					statstring)
				printf(" {hline 16}{c +}{hline "+strofreal(ll)+"}\n")
				for (i=1;i<=sca_players;++i) {
					if (statsdisplay == 1) {
						display = "%7.0g  {c |}{res}  %9.3g"
						printf("{txt}%7.0g "+display+" \n",
							vec_players[i,.],mat_weights[i,.],stats[i,1..1])
					}
					if (statsdisplay == 2 & banzhaf != "nobanzhaf") {
						display = "%7.0g  {c |}{res}  %9.3g      %9.3g"
						printf("{txt}%7.0g "+display+" \n",
							vec_players[i,.],mat_weights[i,.],stats[i,1],stats[i,2])
					}
					if (statsdisplay == 2 & banzhaf == "nobanzhaf") {
						printf("{txt}%7.0g  %7.0g {c |} {res}"+display+" \n",
							vec_players[i,.],mat_weights[i,.],stats[i,1],stats[i,2])
					}
					if (statsdisplay == 3 & banzhaf != "nobanzhaf") {
						display = " %9.3g      %9.3g       %9.3g"
						printf("{txt}%7.0g %7.0g  {c |} {res}"+display+" \n",
							vec_players[i,.],mat_weights[i,.],stats[i,1],stats[i,2],stats[i,3])
					}
					if (statsdisplay == 3 & banzhaf == "nobanzhaf") {
						printf("{txt}%7.0g {c |} {res}"+display+" \n",
							vec_players[i,.],stats[i,1],stats[i,2],stats[i,3])
					}
					if (statsdisplay == 4) {
						printf("{txt}%7.0g %7.0g  {c |} {res}"+display+" \n",
							vec_players[i,.],mat_weights[i,.],stats[i,1],stats[i,2],stats[i,3],stats[i,4])
					}
				}
				printf("{txt} {hline 16}{c +}{hline "+strofreal(ll)+"}")
			}
		// Print output end
	}
end

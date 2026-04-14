function [Dq_alpha, Dq_beta, Dq_gamma, alpha_wM, beta_wM, gamma_wM, shannon] = hill_numbers(M,q)
% See Jost 2007
% M(i,j) biomass of species i at site j
%shannon is the Shannon index by site

assert( q > 0 && q ~=1.0 ) ;

B_site = sum(M,1) ;
norm_B = M*diag(1./(B_site+realmin)) ;
exponent = 1/(1-q) ;
Dq_alpha = mean(sum(norm_B.^q, 1)) .^exponent  ;
shannon = -sum(norm_B.*log(norm_B+realmin),1) ; 
alpha_wM = exp(mean(shannon))  ; %^1D_alpha in Hill numbers (Jost 2007)

norm_B = sum(M,2)/(sum(M(:))+realmin) ; %Metacommunity biomasses normalized
Dq_gamma = sum(norm_B.^q,1).^exponent ;

Dq_beta = Dq_gamma/(Dq_alpha+realmin) ;

gamma_wM = exp(-sum(norm_B.*log(norm_B+realmin),1)) ;%^1D_gamma in Hill numbers (Jost 2007)
beta_wM = gamma_wM/(alpha_wM+realmin) ;%^1D_beta in Hill numbers (Jost 2007)

end
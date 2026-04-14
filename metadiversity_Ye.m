function [alpha_wk, alpha_wM, alpha_uk, alpha_uM, ...
          beta_wM, beta_uM,...
          beta_Jwk,  beta_JwM,  beta_Juk, beta_JuM,...
          gamma_wM, gamma_uM,...
          M_uncut] = metadiversity_Ye(B, w, n_species, n_sites)
% Uses a formulation that is closer to Ye's

% B: matrix with n_species*n_sites columns. Each row represents biomasses
%    for each species at a site at a given time. The number of rows must be at
%    least w
% w: positive integer. Time window. Only the last w rows of B will be
%    considered for computing mean biomasses and diversity indices
% n_species, n_sites: self-explanatory

% alpha_wk: Alpha diversity. 1 x n_sites array of exp(shannon indices) per site
% alpha_wM: Alpha diversity. mean(alpha_wk).
% alpha_uk: Alpha diversity. 1 x n_sites array of richness per site
% alpha_uM: Alpha diversity. mean(alpha_uk).

% beta_wM: Beta diversity. gamma_wM/alpha_wM.
% beta_uM: Beta diversity. gamma_uM - alpha_uM.

% beta_Jwk: Beta diversity. 1 x n_sites array of 1 minus weighted Jaccard index.
% beta_JwM: Beta diversity. mean(beta_Jwk).
% beta_Juk: Beta diversity. 1 x n_sites array of 1 minus unweighted Jaccard index.
% beta_JuM: Beta diversity. mean(beta_Juk).

% gamma_wM: Gamma diversity. exp(Shannon index) for the metacommunity.
% gamma_uM: Gamma diversity. species richness for the metacommunity.

% M_uncut: n_species x n_sites matrix with the avg biomass over the last w
%    samples

assert(w>=1) ;
assert(size(B,2)==n_species*n_sites) ;
assert(size(B,1)>= w ) ;

X = B(end-w+1:end, :) ;
warning('Arbitrary extinction threshold of 1e-7 applied to local biomasses');

M = reshape(mean(X),n_species, n_sites) ;
M_uncut = M ;
M(M<1e-7) = 0 ;

B_site = sum(M,1) ;
norm_B = M*diag(1./(B_site+realmin)) ;
shannon = -sum(norm_B.*log(norm_B+realmin),1) ; %shannon is the Shannon index by site
alpha_wk = exp(shannon) ;
alpha_wM = exp(mean(shannon))  ; %^1D_alpha in Hill numbers (Jost 2007)

last_B = M ;
alpha_uk = sum(last_B>0,1) ;
alpha_uM = mean(alpha_uk) ; %^0D_alpha in Hill numbers (Jost 2007)

norm_B = sum(M,2)/(sum(M(:))+realmin) ; %Metacommunity biomasses normalized
gamma_wM = exp(-sum(norm_B.*log(norm_B+realmin))) ;%^1D_gamma in Hill numbers (Jost 2007)
gamma_uM = sum(norm_B>0) ; %^0D_gamma in Hill numbers (Jost 2007)

beta_wM = gamma_wM/(alpha_wM+realmin) ;%^1D_beta in Hill numbers (Jost 2007)
beta_uM = gamma_uM - alpha_uM ; %^0D_beta in Hill numbers (Jost 2007)

jaccard   = zeros(n_sites) ; %dissimilarities
jaccard_w = zeros(n_sites) ;

for i = 1:n_sites
    for j = (i+1):n_sites
        jaccard(i,j) = 1-sum(min(last_B(:,i)>0, last_B(:,j)>0)) / sum(max(last_B(:,i)>0, last_B(:,j)>0)) ;
        jaccard_w(i,j) = 1-sum(min(M(:,i), M(:,j))) / sum(max(M(:,i), M(:,j))) ;
        jaccard(j,i) = jaccard(i,j) ;
        jaccard_w(j,i) = jaccard_w(i,j) ;
    end
end

beta_Jwk  = sum(jaccard_w,2) / (n_sites-1) ;
beta_Juk = sum(jaccard,2) / (n_sites-1) ;
beta_JwM = mean(beta_Jwk) ;
beta_JuM = mean(beta_Juk) ;

end

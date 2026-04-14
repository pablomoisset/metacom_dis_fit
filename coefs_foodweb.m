function [Coef_A, ...
          Coef_H, ...
          Coef_C, ...
          Coef_r, ... 
          Coef_alpha, ...
          K_capacity, ...
          Coef_x,...
          Coef_e] = coefs_foodweb(A, M, foodweb_pars)
% Input:
%   A is the foodweb adjacency matrix. A(i,j)=1 iff j preys on i
%   TL(i) trophic position of species i. Basals have a zero TP.
%   M(i) body mass of species i (M is a column vector).   

%Attack rate (whithout Arrenhius correction)  
% j preys on i
Coef_A = foodweb_pars.a0 * (M.^foodweb_pars.ba * M'.^foodweb_pars.ca).*A ;

%Handling times (whithout Arrenhius correction)
% j preys on i
Coef_H = foodweb_pars.h0 * (M.^foodweb_pars.bh * M'.^foodweb_pars.ch).*A ;

% Interference. Column vector, because we assume it is predator dependent
% only
Coef_C = foodweb_pars.c0*ones(foodweb_pars.nspecies,1) ;
Coef_C(1:foodweb_pars.nbasals) = 0 ;

Coef_r = foodweb_pars.r0 * M.^foodweb_pars.br ;
Coef_r(foodweb_pars.nbasals+1:end) = 0 ;

% Competition only for plants. Must be symmetric
Coef_alpha = foodweb_pars.alpha ;

%Carrying capacity of plants. Column vector
K_capacity = foodweb_pars.K0 * M.^foodweb_pars.bK ;
K_capacity(foodweb_pars.nbasals+1:end) = 0 ;

Coef_x  = foodweb_pars.x0 * M.^foodweb_pars.bx ;
Coef_x(1:foodweb_pars.nbasals) = foodweb_pars.x0/foodweb_pars.x_cost_basals ; % Plants have an ad-hoc metabolic cost

Coef_e = [foodweb_pars.efficiency_on_plant*ones(foodweb_pars.nbasals,foodweb_pars.nspecies);...
          foodweb_pars.efficiency_on_animal*ones(foodweb_pars.nspecies-foodweb_pars.nbasals,foodweb_pars.nspecies)] .*A ;
end

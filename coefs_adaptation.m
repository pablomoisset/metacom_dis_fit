function [Coef_prio, Coef_OmegaR, Coef_kappaR, Coef_OmegaW, Coef_kappaW] = coefs_adaptation(p)

ones_nspecies = ones(size(p.W.M)) ;

Coef_prio   = p.adapt.pars.p*ones_nspecies ;      %hunger vs fear
Coef_OmegaR = p.adapt.pars.OmegaR*ones_nspecies ; %optimism
Coef_OmegaW = p.adapt.pars.OmegaW*ones_nspecies ;
Coef_kappaR = p.adapt.pars.kappaR*ones_nspecies ;
Coef_kappaW = p.adapt.pars.kappaW*ones_nspecies ;
end
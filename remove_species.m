function p_new = remove_species(p, victims)
%victims is a logical array
assert(numel(victims)==p.W.pars.nspecies) ;
p_new = p ;


p_new.W.adjacency = p.W.adjacency(~victims,~victims) ;
p_new.W.M = p.W.M(~victims) ;
p_new.W.A = p.W.A(~victims,~victims) ;
p_new.W.H = p.W.H(~victims,~victims) ;
p_new.W.C = p.W.C(~victims) ;
p_new.W.r = p.W.r(~victims) ;
p_new.W.alpha = p.W.alpha(~victims,~victims) ;
p_new.W.K = p.W.K(~victims) ;
p_new.W.x = p.W.x(~victims) ;
p_new.W.e = p.W.e(~victims,~victims) ;
p_new.W.trophic_level = p.W.trophic_level(~victims) ;
p_new.W.coherence = NaN; warning('Trophic coherence must be recomputed');
p_new.W.pars.nspecies = p.W.pars.nspecies - sum(victims) ;

p_new.adapt.prio   = p.adapt.prio(~victims) ; 
p_new.adapt.OmegaR = p.adapt.OmegaR(~victims) ; 
p_new.adapt.kappaR = p.adapt.kappaR(~victims) ; 
p_new.adapt.OmegaW = p.adapt.OmegaW(~victims) ; 
p_new.adapt.kappaW = p.adapt.kappaW(~victims) ; 

p_new.D.D_max = p.D.D_max(~victims) ;
p_new.D.Q0    = p.D.Q0(~victims);
p_new.D.U0    = p.D.U0(~victims);
p_new.D.beta  = p.D.beta(~victims);
idx = repmat(reshape(victims,1,p.W.pars.nspecies),1,p.aux.total_sites) ;
p_new.D.gamma = p.D.gamma(~idx,~idx);

p_new.B_init = p.B_init(~idx) ;
p_new.R_init = p.R_init(~idx) ;
p_new.W_init = p.W_init(~idx) ;
p_new.init   = p.init(~repmat(idx,1,3)) ;

p_new.aux = update_aux(p_new) ;
end
function aux = update_aux(p)

aux = struct() ;

aux.total_sites = p.L.pars.nsites+p.L.pars.species_pool ;

aux.inv_dist = kron(p.L.inv_dist, speye(p.W.pars.nspecies)) ;

aux.KL1 = kron(p.L.dist>0, speye(p.W.pars.nspecies)) ;   % Kronecker landscape1 (hat{1})

tmp_m = speye(aux.total_sites) ;
tmp_v = sparse(ones(aux.total_sites,1)) ;

if p.L.pars.species_pool
   tmp_m(end,end) = 0 ;
   tmp_v(end) = 0 ;
end

aux.Tr = kron(tmp_m, p.W.adjacency) ; %Block diagonal matrix with foodweb adjacency
aux.A  = kron(tmp_m, p.W.A) ; %Attack rates. nsites block diagonal. Blocks are nspecies^2
aux.H  = kron(tmp_m, p.W.H) ; %Handling times. nsites block diagonal. Blocks are nspecies^2
aux.C  = kron(tmp_v, p.W.C) ; %same as biomasses
aux.alpha = kron(tmp_m, p.W.alpha) ; %nsites block diagonal. Blocks are nspecies^2
aux.e  = kron(tmp_m, p.W.e) ; %nsites block diagonal. Blocks are nspecies^2

aux.r = kron(tmp_v, p.W.r) ; %nsites x 1 block matrix. Blocks are nspecies x 1
aux.K = kron(tmp_v, p.W.K).* ...
          kron(p.L.K_factor, sparse(ones(p.W.pars.nspecies,1))) ; %nsites x 1 block matrix. Blocks are nspecies x 1
aux.x = kron(tmp_v, p.W.x) ; %nsites x 1 block matrix. Blocks are nspecies x 1

%Arrhenius corrections and effect of risk on "x" (metabolism)
tmp = sparse(ones(p.W.pars.nspecies,1)) ;
aux.r = aux.r .* kron(arrhenius(p.L.T, p.W.pars.Er, p.W.pars.T0), tmp) ;

basal_idx = repmat(sum(p.W.adjacency).'==0, aux.total_sites, 1 ) ;
aux.x = aux.x .* kron(arrhenius(p.L.T, p.W.pars.Ex, p.W.pars.T0), tmp) ;
aux.x_risk = aux.x * p.W.pars.x_cost_animals ;
aux.x_risk(basal_idx) = aux.x(basal_idx) * p.W.pars.x_cost_basals ;


%aux.A, aux.H, aux.AH are nsites-block diagonal. Blocks are nspecies^2
aux.A = aux.A * diag(kron(arrhenius(p.L.T, p.W.pars.Ea, p.W.pars.T0),tmp)) ;
aux.H = aux.H * diag(kron(arrhenius(p.L.T, p.W.pars.Eh, p.W.pars.T0),tmp)) ;
aux.AH = aux.A .* aux.H ;
% aux.K, aux.C are (nsites * nspecies)x 1 vector
aux.K = aux.K .* kron(arrhenius(p.L.T, p.W.pars.EK, p.W.pars.T0),tmp) ;
aux.C = aux.C .* kron(arrhenius(p.L.T, p.W.pars.Ec, p.W.pars.T0),tmp) ;

aux.prio = kron(tmp_v, p.adapt.prio) ;

aux.OmegaR = kron(tmp_v, p.adapt.OmegaR) ;
aux.kappaR = kron(tmp_v, p.adapt.kappaR) ; %adaptation speed
aux.OmegaW = kron(tmp_v, p.adapt.OmegaW) ;
aux.kappaW = kron(tmp_v, p.adapt.kappaW) ; %adaptation speed

aux.D_max = sparse(repmat(p.D.D_max, p.aux.total_sites,1)) ;
aux.Q0    = kron(tmp_v, p.D.Q0) ;
aux.U0    = sparse(kron(diag(tmp_v)*p.L.dist>0, diag(p.D.U0))) ;

aux.beta  = diag(kron(tmp_v, p.D.beta))*aux.KL1 ;
if p.L.pars.species_pool %species pool to rest of landscape is passive
    aux.beta(end-p.W.pars.nspecies+1:end,:) = 0 ;
end
aux.Dd_ratio  = diag(sparse(aux.D_max))*aux.inv_dist ;

if ~p.use_sparse
    aux.inv_dist = full(aux.inv_dist) ;
    aux.KL1 = full(aux.KL1) ;
    aux.Tr = full(aux.Tr) ;
    aux.A = full(aux.A) ;
    aux.H = full(aux.H) ;
    aux.C = full(aux.C) ;
    aux.alpha = full(aux.alpha) ;
    aux.e = full(aux.e) ;
    aux.r = full(aux.r) ;
    aux.K = full(aux.K) ;
    aux.x = full(aux.x) ;
    aux.x_risk = full(aux.x_risk) ;
    aux.AH = full(aux.AH) ;
    aux.prio = full(aux.prio) ;
    aux.OmegaR = full(aux.OmegaR) ;
    aux.OmegaW = full(aux.OmegaW) ;
    aux.kappaR = full(aux.kappaR) ;
    aux.kappaW = full(aux.kappaW) ;
    aux.D_max = full(aux.D_max) ;
    aux.Q0 = full(aux.Q0) ;
    aux.U0 = full(aux.U0) ;
    aux.beta = full(aux.beta) ;
    aux.Dd_ratio = full(aux.Dd_ratio) ;
end

end
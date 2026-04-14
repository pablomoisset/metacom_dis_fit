function p_sim_ready = define_system(p)

% p.L.pars has:
%  seed: random seed, NaN or Inf means use current
%  nsites: number of sites
%  edges: number of links
%  excess: excess factor
%  centers: number of clusters
%  p_exponent: probability exponent
%  species_pool


lpar.seed = p.L.pars.seed ;
lpar.land_size = p.L.pars.nsites ;
lpar.edges = p.L.pars.edges ;
lpar.excess = p.L.pars.excess ; %Excess factor
lpar.centers = p.L.pars.centers ;
lpar.p_exponent = p.L.pars.p_exponent ;
lpar.non_intermitent = 0 ; % don't care value, feature not used
[dist, ~, ~ , xP, yP] = p_modular_landscape(lpar) ;

%warning('Forcing a vertical line as landscape') ;
%yP = linspace(0,500,p.L.pars.nsites).' ;
%xP = yP*0 ;
%dist = 500*diag(ones(p.L.pars.nsites-1,1),1) ;
%dist = (dist+dist.')/(p.L.pars.nsites-1) ;

p.aux.total_sites = p.L.pars.nsites+p.L.pars.species_pool ;

if p.L.pars.species_pool
   dist(p.aux.total_sites,:) = 1000.0 ;
   dist(:,p.aux.total_sites) = 0.0 ;
   xP(p.aux.total_sites) = 0 ;
   yP(p.aux.total_sites) = 0 ;
end

p.L.dist = dist ;
p.L.xP = xP ;
p.L.yP = yP ;
p.L.inv_dist = (dist>0)./(dist+realmin) ;

if p.use_sparse
    p.L.dist = sparse(p.L.dist) ;
    p.L.inv_dist = sparse(p.L.inv_dist) ;
end

%p.W.pars has:
%  seed
%  nspecies
%  nlinks
%  nbasals: number of basal species
%  temperature
%  efficiency

if isfinite(p.W.pars.seed)
    save_seed = rng  ;
    rng(p.W.pars.seed) ;
end
%p.W.pars.nbasals = 12 ;
%p.W.pars.nlinks = 2*p.W.pars.nspecies;
%p.W.pars.temperature = 0.2 ;
%[adjacency, TL, q]  = gppm_coherence( p.W.pars.nspecies, p.W.pars.nbasals,...
%                                      p.W.pars.nlinks, p.W.pars.temperature ) ;

warning('Forcing cascade/niche foodweb model. Ignores nbasals, nlinks and temperature')
warning('Limiting supremum trophic level to 6')
undef_foodweb = true ;
while undef_foodweb
    adjacency = gen_niche_model(p.W.pars.connectance, p.W.pars.nspecies, 1) ;
    adjacency = adjacency - diag(diag(adjacency)) ; %Remove cannibalism
    tmp1 = sum(adjacency,1)==0 ;
    p.W.pars.nbasals = sum(tmp1) ;
    if p.W.pars.nbasals == 0
        undef_foodweb = true ;
    else
        G = graph(adjacency+adjacency') ;
        [~, binsizes] = conncomp(G) ;
        undef_foodweb = numel(binsizes) > 1 ;
        if ~undef_foodweb
            tmp2 = sum(adjacency)~=0 ;
            adjacency = [adjacency(:,tmp1) adjacency(:,tmp2)] ; %Rearrange adjacency matrix so basals are at front
            adjacency = [adjacency(tmp1,:); adjacency(tmp2,:)] ;
            
            p.W.pars.nlinks = sum(adjacency(:)) ;
            TL = tropos(adjacency)' ;
            %TL = TP_shortestpath(adjacency)' ;
            undef_foodweb = ~all(isfinite(TL) & (TL<6) & (TL>=0)) ;
        end
    end
    q = 0 ; % Trophic coherence is set to zero. This number is not used 
end

warning('Added inter-plant competition, as Ye & Wang 2023 did') ;
p.W.pars.alpha = 0.2*rand(p.W.pars.nbasals) ; % Inter competition
p.W.pars.alpha = p.W.pars.alpha - diag(diag(p.W.pars.alpha)) + ...
                        p.W.pars.alpha_all * eye(p.W.pars.nbasals) ;            % Intra competition
p.W.pars.alpha(p.W.pars.nspecies, p.W.pars.nspecies) = 0 ;           % Extend the matrix with zeros

% Conversion efficiency
p.W.pars.e =  [p.W.pars.efficiency_on_plant*ones(p.W.pars.nbasals,1);... 
               p.W.pars.efficiency_on_animal*ones(p.W.pars.nspecies-p.W.pars.nbasals,1)] ;

%if isfinite(p.W.pars.seed)
%    rng(save_seed) ;
%end                          
                          
% Compute body masses
M = p.W.pars.m0 * (p.W.pars.Z.^(TL' + p.W.pars.omega*2.0*(rand(size(TL'))-0.5))) ;

if isfinite(p.W.pars.seed)
    rng(save_seed) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Coef_A, Coef_H, Coef_C,...
 Coef_r,...
 Coef_alpha, K_capacity,...
 Coef_x, Coef_e] = coefs_foodweb(adjacency,M,p.W.pars) ;
                          
p.W.adjacency = adjacency ;
p.W.M   = M ;      % Body weights
p.W.A   = Coef_A ; % Attack rates. Positive
p.W.H   = Coef_H ;
p.W.C   = Coef_C ;
p.W.r = Coef_r ;
p.W.alpha = Coef_alpha ;
p.W.K = K_capacity ;
p.W.x = Coef_x ;
p.W.e = Coef_e ;
p.W.trophic_level = TL ;
p.W.coherence = q ;
if p.use_sparse
    p.W.adjacency = sparse(p.W.adjacency) ;
    p.W.A = sparse(p.W.A) ;
    p.W.H = sparse(p.W.H) ;
    p.W.e = sparse(p.W.e) ;
    p.W.pars.alpha = sparse(p.W.pars.alpha) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Coef_D_max, Coef_Q0, Coef_U0 , Coef_beta, Coef_gamma] = coefs_migration(p) ;
p.D.D_max = Coef_D_max ;
p.D.Q0    = Coef_Q0 ;
p.D.U0    = Coef_U0 ;
p.D.beta  = Coef_beta ;
p.D.gamma = Coef_gamma ;
%warning('Should I sparsify more fields?')
if p.use_sparse
    p.D.gamma = sparse(p.D.gamma) ;
end

assert(all(Coef_gamma(:)<=1)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Coef_prio, Coef_OmegaR, Coef_kappaR, Coef_OmegaW, Coef_kappaW] = coefs_adaptation(p) ;

p.adapt.prio   = Coef_prio ;
p.adapt.OmegaR = Coef_OmegaR ;
p.adapt.kappaR = Coef_kappaR ;
p.adapt.OmegaW = Coef_OmegaW ;
p.adapt.kappaW = Coef_kappaW ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.aux = update_aux(p) ;

p.init = [p.B_init; p.B_init.*p.R_init; p.B_init.*p.W_init] ;

p_sim_ready = p ;
end

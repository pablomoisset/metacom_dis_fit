function [p,results] = run_cluster(...
   master_seed,...
   species_pool,... %boolean
   landscape_seed,...
   nsites,...
   edges,...
   excess,...
   centers,...
   p_exponent,...
   T_min,...
   T_max,...
   T_noise,...
   K_noise,... %factor
   foodweb_seed,...
   nspecies,...
   connectance,...
   competition,... %intra, for plants
   D0,...
   beta_mig,...
   gamma_0,...
   kappaR,...
   OmegaR,...
   Rmin,...
   x_cost_plants,...
   x_cost_animals,...
   priority,... %0: fear, 1:hunger
   kappaW,...
   OmegaW,...
   Wmin,...
   Wmax,...
   init_B_cont,...
   replicates,...
   tspan_transient_pars,... %triples [first last qty]
   tspan_final_pars,....    %triples [first last qty]
   name_prefix...
)

rng(master_seed+replicates(1)) ;

p = test_base(species_pool) ;

p.L.pars.seed = landscape_seed+replicates(1) ;
p.L.pars.nsites = nsites ;
p.L.pars.edges = edges ;
p.L.pars.excess = excess ;
p.L.pars.centers = centers ;
p.L.pars.p_exponent = p_exponent ;
%p.L.T = 293*ones(p.L.pars.nsites+p.L.pars.species_pool,1);
p.adapt.pars.kappaR = kappaR ; %0.25
p.adapt.pars.kappaW = kappaW ; %10
p.adapt.pars.epsilon=0.1;
    
p.W.pars.nspecies = nspecies ;
p.W.pars.connectance = connectance ;
p.W.pars.alpha_all = competition ;


p.R_init = repmat(ones(p.W.pars.nspecies,1), p.L.pars.nsites+p.L.pars.species_pool,1) ;
p.W_init = repmat(293*ones(p.W.pars.nspecies,1), p.L.pars.nsites+p.L.pars.species_pool,1) ;

p.adapt.pars.OmegaR = OmegaR ;
p.adapt.pars.Rmin = Rmin ;
p.adapt.pars.x_cost_plants = x_cost_plants ;
p.adapt.pars.x_cost_animals = x_cost_animals ;
p.adapt.pars.p = priority ;

p.adapt.pars.OmegaW = OmegaW ;
p.adapt.pars.Wmin = Wmin ;
p.adapt.pars.Wmax = Wmax ;

init_biomasses ;

%p = rand_W0(p, 10.0) ;

p.L.pars.T_max = T_max ;%310.0 ;
p.L.pars.T_min = T_min ;%280 ;
p.L.pars.T_noise = T_noise ;
p.L.T = 293*ones(p.L.pars.nsites+p.L.pars.species_pool, 1) ; %default temperature
p.L.K_factor = ones(p.L.pars.nsites+p.L.pars.species_pool, 1) ; % site factor for carrying capacity

p = rand_K_factor(p, K_noise) ; %1

p.D.pars.beta = beta_mig*ones(p.W.pars.nspecies,1);
p.D.pars.d0 = realmax*ones(p.W.pars.nspecies,1) ; %Numbers to make gamma=gamma_0 
%p.D.pars.M0 = realmin*p.W.pars.m0*p.W.pars.Z ; 


p.D.pars.gamma0 = gamma_0 ; %0.5
p.D.pars.D0 = D0 ; %1e-8 ;

init_biomasses ;

p = define_system(p) ;

p = rand_T(p) ;

tspan_transient = linspace(tspan_transient_pars(1), tspan_transient_pars(2),tspan_transient_pars(3));
tspan_final = linspace(tspan_final_pars(1), tspan_final_pars(2),tspan_final_pars(3));

%assert(all(p.D.beta==0)) ;

results = struct() ;

tmp = zeros(numel(replicates),1) ;
results.alpha_wM_list = tmp ;
results.alpha_uM_list = tmp ;
results.beta_wM_list  = tmp ;
results.beta_uM_list  = tmp ;
results.beta_JwM_list = tmp ;
results.beta_JuM_list = tmp ;
results.gamma_wM_list = tmp ;
results.gamma_uM_list = tmp ;
results.gamma_biomass_list = tmp ;

all_p = [] ;

p.tspan_transient_pars = tspan_transient_pars ;
p.tspan_final_pars = tspan_final_pars ;

original_p = p ;

for r = 1:numel(replicates)
    fprintf("Replicate=%d\n", replicates(r)) ;
    rng(master_seed+replicates(r))
    p = original_p ;
    p.L.pars.seed = landscape_seed+replicates(r) ;
    p.W.pars.seed = foodweb_seed+replicates(r) ;
    p = define_system(p) ;
    p = rand_K_factor(p, K_noise) ; %1
    p = rand_T(p) ;
   warning('W noise set to 0') ;
%    p = rand_W0(p, 10.0) ;
    p = rand_W0(p, 0.0) ;
    p.aux = update_aux(p) ;

    init_biomasses ;

    if ~isempty(name_prefix)
        p_save = p ;
        p_save.aux = [] ;
        p_save.init = [] ;
        all_p(r).p = p_save ;
    end
    
%    tic
    [t, X, B, R, W, Z, Y] = run_sim(p, tspan_transient, 5000);
%    toc
%    last_B = B(end,:) ;
%    last_X = X(end,:) ;
%    last_p = p ;
%    last_d=rhs_full_flows(0,last_X',last_p);
    
%    plot_tseries ;
    
    extinct = sum(reshape(B(end,:),p.W.pars.nspecies, p.L.pars.nsites),2)<1e-12 ;
    if sum(extinct)~=p.W.pars.nspecies
        p.init = X(end,:)' ;
        tmp = reshape(p.init,p.aux.total_sites*p.W.pars.nspecies,3) ;
        p.B_init = tmp(:,1) ;
        p.R_init = tmp(:,2)./(p.B_init+realmin) ;
        p.W_init = tmp(:,3)./(p.B_init+realmin) ;
        p = remove_species(p,extinct') ;
%        new_d=rhs_full_flows(0,p.init,p);
        %    tic
        [t, X, B, R, W, Z, Y] = run_sim(p, tspan_final, 5000);
        %    toc
 %           plot_tseries ;
    end
    w = 100;
    warning('Using a mean based metadiversity function') ;
    [alpha_wk, alpha_wM, alpha_uk, alpha_uM, ...
        beta_wM, beta_uM,...
        beta_Jwk,  beta_JwM,  beta_Juk, beta_JuM,...
        gamma_wM, gamma_uM,...
        M] = metadiversity_Ye(B, w, p.W.pars.nspecies, p.L.pars.nsites);
    
  %  mean_biomasses = M ;
  %  plot_metadiversity ;
  
    results.alpha_wM_list(r) =  alpha_wM ;
    results.alpha_uM_list(r) =  alpha_uM ;
    results.beta_wM_list(r)  =  beta_wM  ;
    results.beta_uM_list(r)  =  beta_uM  ;
    results.beta_JwM_list(r) =  beta_JwM ;
    results.beta_JuM_list(r) =  beta_JuM ;
    results.gamma_wM_list(r) =  gamma_wM ;
    results.gamma_uM_list(r) =  gamma_uM ;
    results.gamma_biomass_list(r) = sum(M(:)) ;
    results.M{r} = M ;
    results.extinct{r} = extinct ;
end

if ~isempty(name_prefix)
    save([name_prefix sprintf('_rep%d_%d.mat', min(replicates), max(replicates) )], 'all_p','results');
end


%MA=mean(alpha_wM_list,3);
%MBeta=mean(beta_JwM_list,3);
%MG=mean(gamma_wM_list,3);
%MBio=mean(gamma_biomass_list,3);
%MAu=mean(alpha_uM_list,3);

%sMA=std(alpha_wM_list,0,3);
%sMBeta=std(beta_JwM_list,0,3);
%sMG=std(gamma_wM_list,0,3);
%sMBio=std(gamma_biomass_list,0,3);
%sMAu=std(alpha_uM_list,0,3);
end

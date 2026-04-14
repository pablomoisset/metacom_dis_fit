
master_seed = 1 ;
species_pool = 'false' ;
landscape_seed = 1000 ;
nsites = 25%90 ;
%edges = 50%180 ;
excess = 1 ;%5 ;
centers = 3; %2 ;

p_exponent = 4 ;

%T_radius = 5 ;%15 ;
%T_min = 293-T_radius ;
%T_max = 293+T_radius ;
T_noise = 0.28 %0.5 ;
K_noise = 1.0 ;

foodweb_seed = 2000 ;
nspecies = 50 ;
connectance = 0.1 ;
competition = 1.0 ;
%D0 = 1e-10 ;
%beta_mig = 0.99;%3 niveles
gamma_0 = 0.0 ;
kappaR = 0.0 ;
OmegaR = 0 ;
Rmin = 0.5 ;
x_cost_plants = 0.5 ;
x_cost_animals = 2.0 ;
priority = 0.5 ; %0: fear, 1:hunger
%kappaW = 0 ;%10 ;
OmegaW = 0 ;
Wmin = 273 ;
Wmax = 313 ;
initial_B_cont = 1e-7 ;
replicates = '0:29' ;
tspan_transient_pars = '[0 5e8 1001]' ;
tspan_final_pars = '[0 5e8 1001]' ;   %triples [first last qty]
name_prefix = 'results_paper/R6_' ;


landscape_seed = 1 ;
foodweb_seed   = 1000 ;

for nsites = [25]
    if nsites == 3
        edges = 3 ;
    else
        edges = 2*nsites ;
    end
    for T_radius = 5% [0.5,5,15]
        T_min = 293-T_radius ;
        T_max = 293+T_radius ;
        for beta_mig =  [0 0.99]
            name = sprintf('ns%d_ne%d_Tr%g_beta%g_gamma%g_tnoise%g_Ref6.job',nsites,edges,T_radius,beta_mig,gamma_0,T_noise) ;
            create_jobs_helper
        end
    end
end


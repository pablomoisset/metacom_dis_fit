function p_no_sim = test_base(species_pool)
p = struct() ;
%rng(10)

p.W.pars.seed = 4 ;
p.W.pars.nspecies = 15 ;
p.W.pars.connectance = 0.15 ;
%p.W.pars.nlinks  = 67; %2*p.W.pars.nspecies ;
%p.W.pars.nbasals = 8 ;
%p.W.pars.temperature = 0.25 ;%0.001 ;
p.W.pars.efficiency_on_plant  = 0.545 ;
p.W.pars.efficiency_on_animal = 0.906 ;

p.L.pars.seed = 2 ;%50
p.L.pars.nsites = 10 ;% 20 ;
p.L.pars.edges = 20 ;%60 ;
p.L.pars.excess = 5 ;
p.L.pars.centers = 3 ;%3 ;
p.L.pars.p_exponent = 4 ;
p.L.pars.species_pool = species_pool ;

p.W.pars.m0 = 1e-2 ;        %Body mass of a basal (trophic level 0)
p.W.pars.Z  = 10 ;          %Predator/prey body mass ratio
p.W.pars.r0 = exp(-15.68) ; %Allometric constant for r
p.W.pars.x0 = exp(-16.54) ; %Allometric constant for x
p.W.pars.a0 = exp(-13.1)  ; %Allometric constant for a
p.W.pars.h0 = exp(9.66)   ; %Allometric constant for h
p.W.pars.K0 = 40 ;          %Allometric constant for K
p.W.pars.br = -0.25 ;       %Allometric exponent for r
p.W.pars.bx = -0.31 ;       %Allometric exponent for x
p.W.pars.ba =  0.25 ;       %Allometric exponent for a
p.W.pars.bh = -0.45 ;       %Allometric exponent for h
p.W.pars.bK =  0.28 ;       %Allometric exponent for K

p.W.pars.ca =  -0.8 ;       %Second allometric exponent for a
p.W.pars.ch =  0.47 ;       %Second allometric exponent for h

p.W.pars.alpha_all =  1 ; % Competition among plants

p.W.pars.c0 =  0.8 ; % Coefficient of saturation constant c_ij
p.W.pars.q = 1 ;

p.W.pars.Er = -0.84 ; %Activation energies
p.W.pars.Ex = -0.69 ;
p.W.pars.Ea = -0.38 ;
p.W.pars.Eh =  0.26 ;
p.W.pars.EK =  0.71 ;
p.W.pars.Ec = -0.65 ;

p.W.pars.T0 = 293 ; % About 20C in K


p.W.pars.x_cost_animals =  2.0 ;%      %1: no cost, inf: maximum cost
p.W.pars.x_cost_basals  =  0.5 ;%      %1: no cost, 0: maximum cost

p.W.pars.omega   = 0.1 ; %Noise for trophic level computation

p.D.pars.bD = 0.25 ; %Allometric exponent for maximum dispersal rate
%warning('Migration disabled')
p.D.pars.D0 = 1e-7; %Allometric coefficient for maximum dispersal rate
%p.D.pars.D0 = 0.0 ; %Allometric coefficient for maximum dispersal rate
%warning('p.D.pars.D0 = 0.0 ') ;

p.D.pars.Q0 = 1e-5 ; 
p.D.pars.U0 = 1e-6 ; 
p.D.pars.gamma0 = 0 ; 
p.D.pars.b_gamma = 0.5 ; 
p.D.pars.d0 = realmax*ones(p.W.pars.nspecies,1) ; %Numbers to make gamma=gamma_0 
p.D.pars.M0 = realmin*p.W.pars.m0*p.W.pars.Z ; 
p.D.pars.beta = 0*ones(p.W.pars.nspecies,1); %zero for passive dispersal
%p.D.pars.alpha = NaN ;

p.adapt.pars.p = 0.5 ;
p.adapt.pars.OmegaR = 0 ; 
p.adapt.pars.OmegaW = 0 ; 
p.adapt.pars.kappaR = 0.25; %warning('kappas = 0') 
p.adapt.pars.kappaW = 10; 
p.adapt.pars.Rmin = 0.5 ;
p.adapt.pars.Wmin = 270 ; 
p.adapt.pars.Wmax = 315 ; 

%p.adapt.pars.T = NaN ; 
p.adapt.pars.epsilon = 0.1 ; % For the clamping function

p.use_sparse = true ;

p.L.T = 293*ones(p.L.pars.nsites+p.L.pars.species_pool, 1) ; %default temperature
p.L.K_factor = ones(p.L.pars.nsites+p.L.pars.species_pool, 1) ; % site factor for carrying capacity

% Default initial conditions
if p.L.pars.species_pool
    p.B_init = 0.1*repmat(ones(p.W.pars.nspecies,1)/(p.L.pars.nsites+p.L.pars.species_pool-1), p.L.pars.nsites+p.L.pars.species_pool,1) ;
    p.B_init((end-p.W.pars.nspecies+1):end) = 1e-7*ones(p.W.pars.nspecies,1)/(p.D.pars.D0+realmin) ;
else
    p.B_init = 0.1*repmat(ones(p.W.pars.nspecies,1)/(p.L.pars.nsites+p.L.pars.species_pool), p.L.pars.nsites+p.L.pars.species_pool,1) ;    
end

p.R_init = repmat(ones(p.W.pars.nspecies,1), p.L.pars.nsites+p.L.pars.species_pool,1) ;
p.W_init = repmat(293*ones(p.W.pars.nspecies,1), p.L.pars.nsites+p.L.pars.species_pool,1) ;

p_no_sim = p ;

end

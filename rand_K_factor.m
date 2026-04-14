function p_new = rand_K_factor(p, noise)
    assert ( noise >= 0.0) ;
    p_new = p ;
    p_new.L.K_factor = 1.0 + noise*rand(p.L.pars.nsites+p.L.pars.species_pool, 1) ; % site factor for carrying capacity
end
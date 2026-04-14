function [Coef_D_max, Coef_Q0, Coef_U0 , Coef_beta , Coef_gamma] = coefs_migration(p)
% Coef_D_max is a nspecies x 1 vector
% Coef_Q0 is a nspecies x 1 vector
% Coef_U0 is a nspecies x 1 vector
% Coef_beta is a (nspecies) x 1 vector

% Coef_gamma is a block matrix (nsites x nsites) with diagonal blocks (nspecies x nspecies) 

ones_nspecies = ones(size(p.W.M)) ;
Coef_D_max = p.D.pars.D0 .* p.W.M.^p.D.pars.bD ;
Coef_Q0    = p.D.pars.Q0*ones_nspecies ;
Coef_U0    = p.D.pars.U0*ones_nspecies ;
Coef_beta = p.D.pars.beta ;

T1 = kron(p.L.dist, eye(p.W.pars.nspecies)) ;
T2 = kron(p.L.dist>0, diag(p.D.pars.d0)) ;
T3 = zeros(p.aux.total_sites*p.W.pars.nspecies) ;
nz = find(T2) ;
T3(nz) = 1./(T1(nz)+T2(nz)) ;
T4 = T2.*T3 ; %T4 holds d_i^0/(d_i^0+d_kl)


T5 = diag(p.W.M.^p.D.pars.b_gamma ./ (p.W.M.^p.D.pars.b_gamma+p.D.pars.M0.^p.D.pars.b_gamma) ) ;
T6 = kron(p.L.dist>0, T5) ;  %T5 holds M_i^b_gamma/(M_i^b_gamma+M_0^b_gamma)

Coef_gamma = p.D.pars.gamma0 * (T4.*T6) ;
end

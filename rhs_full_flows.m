function [dX, mig_in, mig_out] = rhs_full_flows(t, X, p)
%t
global rhs_counter rhs_mod
if mod(rhs_counter,rhs_mod)==0
    disp(t)
end
rhs_counter = rhs_counter + 1 ;
n = uint32(numel(X)) ;

%X = sparse(X) ;

B = sparse(X(1:idivide(n,3))) ;
Z = sparse(X(idivide(n,3)+1:2*idivide(n,3))) ;
Y = sparse(X(2*idivide(n,3)+1:end)) ;

B_hat = Z ;
R = Z./(B+realmin) ;
W = Y./(B+realmin) ;
fW = full(W) ;

%if any(W<288)&&t>2.072e7
    %W'
%    1;
%end

temps = temperatures(t,p) ;
T_tilde = kron(temps, ones(p.W.pars.nspecies, 1)) ;
K = p.aux.K ;
rho_temp = sparse(rho_temperature_adapt(fW,T_tilde)) ;

r_hat = p.aux.r .* rho_temp ;
x_hat = p.aux.x + (p.aux.x_risk-p.aux.x).*R ;
e_hat = p.aux.e * diag(rho_temp) ;

G = (p.aux.alpha*B)./(K+realmin) ;
G(K>0) = 1-G(K>0) ;

D_tilde = 1 + B_hat.*p.aux.C +  (B_hat'*p.aux.AH).' ; %denominators of functional responses

inv_D_tilde = diag(sparse(1./(D_tilde+realmin))) ;
F_noB = p.aux.A*inv_D_tilde ;

%Q_plus_tilde = p.aux.r.*G  - x_hat + R.*sum(diag(B_hat)*(F_noB.*p.aux.e), 1).' ;
Q_plus = r_hat.*G  - x_hat + R.*sum(diag(B_hat)*(F_noB.*e_hat), 1).' ;

Q_minus = R.*(F_noB*B_hat) ;

%if t>3.0025e7
%    t
%end
if any(p.aux.beta(:)~=0)
    %   Q_tilde = Q_plus_tilde - Q_minus ;
    % warning('Modificar esto y correr heatmap')
    Q_tilde = Q_plus - Q_minus ;
    
    f_hat = f(Q_tilde, p.aux.Q0) ;

    tmp = p.aux.KL1*diag(Q_tilde) ;
    g_hat = g(tmp - tmp.', p.aux.U0 ) ;
    
    D = (p.aux.beta.*(diag(f_hat)*(p.aux.KL1-p.D.gamma)+ p.D.gamma.*g_hat)+(p.aux.KL1-p.aux.beta)).*p.aux.Dd_ratio ;
else
    D = p.aux.KL1.*p.aux.Dd_ratio ; %Optimization for passive dispersal only
end
tmp = diag(B)*D ;

mig_in = sum(tmp,1)' ;
mig_out = sum(tmp,2) ; 

if p.L.pars.species_pool
   mig_out(end-p.W.pars.nspecies+1:end) = 0 ;
end

dB = (Q_plus-Q_minus).*B ;
small_B = B<1e-12 & dB<0 ;
dB(small_B) = -1e-8*B(small_B) ;
dB = dB + mig_in - mig_out ;


F = diag(B_hat)*F_noB ;

if any(p.adapt.kappaR>0)
    
    phi_R = 2*p.aux.prio.*(...
        sum((F-(diag(B_hat)*p.aux.A*diag((B_hat.*D_tilde.^-2).*p.aux.C))).*e_hat).' - ...
        (p.aux.x_risk-p.aux.x)) ...
        -2*(1-p.aux.prio).*((p.aux.A*inv_D_tilde - diag(B_hat)*(p.aux.A.*p.aux.AH)*(inv_D_tilde.^2))*B_hat) + p.aux.OmegaR; %nueva
    %    -2*(1-p.aux.prio).*((p.aux.A*inv_D_tilde -  p.aux.A*diag(B_hat).*p.aux.AH*(inv_D_tilde.^2))*B_hat) + p.aux.OmegaR; %  Estaba mala
    
%    candidate1 = p.aux.A*diag(B_hat).*p.aux.AH*(inv_D_tilde.^2);
 %   candidate2 = diag(B_hat)*(p.aux.A.*p.aux.AH)*(inv_D_tilde.^2) ;

    
    dZ = p.aux.kappaR.*d_clamp_v( (R-p.adapt.pars.Rmin)/(1-p.adapt.pars.Rmin),phi_R, p.adapt.pars.epsilon).*phi_R.*B + ...
        + dB.*R ;
    dZ = dZ + (diag(R)*D - D*diag(R)).'*B ;
else
    dZ = dB.*R + (diag(R)*D - D*diag(R)).'*B ; %Special case kappa_R==0
end

if any(p.adapt.kappaW>0)
    dr = drho_dW(fW,T_tilde) ;
    phi_W = dr .* (p.aux.r.*G.*(2*(G>0)-1) + R.*sum(F.*p.aux.e).') + p.aux.OmegaW ;
    
    dY = p.aux.kappaW.*d_clamp_v( (W-p.adapt.pars.Wmin)/(p.adapt.pars.Wmax-p.adapt.pars.Wmin), phi_W, p.adapt.pars.epsilon).*phi_W.*B + ...
        + W.*dB ;
    
    dY = dY + (diag(W)*D - D*diag(W)).'*B ;
else
    dY = W.*dB + (diag(W)*D - D*diag(W)).'*B ; %Special case Kappa_W==0
end

dX=full([dB; dZ; dY]) ;
end

function result = f(Q,Q0)
   result = Q.*(Q<0) ; %Q.*heaviside(-Q)
   result = result./((Q-Q0)+realmin) ;
end

function result = g(Q_dif, U0)
[i,j,v] = find(Q_dif+U0) ;
tmp = sparse(i,j,1./v,size(U0,1),size(U0,2),numel(v)) ;
result = (Q_dif.*tmp).*(Q_dif>0) ;
end


function s = d_clamp_v(T,dFit,epsilon)

s = T ; % Preallocate
idx = (dFit > 0 & T <= 1-epsilon) | (dFit <= 0 & T > epsilon) ;
s(idx) = 1 ;

idx = dFit > 0 & T > 1-epsilon ;
s(idx) = (1-T(idx))/epsilon ;

idx = dFit <= 0 & T <= epsilon ;
s(idx) = T(idx)/epsilon ;
end

function s = d_clamp(T,dFit,epsilon)

s = T ; % Preallocate

for i=1:numel(T)
    if dFit(i) > 0
        if T(i) <= 1-epsilon
            s(i) = 1 ;
        else
            s(i) = (1-T(i))/epsilon ;
        end
    else
        if T(i) > epsilon
            s(i) = 1 ;
        else
            s(i) = T(i)/epsilon ;
        end
    end
end


end

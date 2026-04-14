function [t, X, B, R, W, Z, Y] = run_sim(p, tspan, report_every)

init = p.init ;

global rhs_counter rhs_mod

rhs_counter = 0 ;
rhs_mod = report_every ;

rhs = @(t, X) rhs_full(t, X, p) ;
ev = @(t, X) extinction_events(t, X, p) ;

%options = odeset('RelTol',1e-7,'AbsTol',1e-16);
options = odeset('RelTol',1e-8,'AbsTol',1e-10); %Looks ok. Larger values lead to oscillations in R
%options = odeset('RelTol',1e-8,'AbsTol',1e-10);
options = odeset(options, 'Events', ev) ;
options = odeset(options, 'NonNegative', 1:numel(init)) ;

t = [] ;
X = [] ;
current_t = tspan(1) ;
disp('Integration begins.') ;

while current_t < tspan(end)
    [partial_t, partial_X,te,ye,ie] = ode45(rhs,tspan,init,options) ;
    current_t = partial_t(end) ;
    t = [t; partial_t] ;
    X = [X; partial_X] ;
    if ~isempty(ie)
       fprintf('At %d: Extinct', partial_t(end)) ;
       fprintf(' %d ', ie) ;
       fprintf('\n') ;
       init = ye(end,:)' ;
       for k=0:(p.aux.total_sites-1)
          init(ie + 3*k*p.W.pars.nspecies) = 0 ;                        % B=0
          init(ie + p.W.pars.nspecies + 3*k*p.W.pars.nspecies) = 0 ;    % Z=0
          init(ie + 2*p.W.pars.nspecies + 3*k*p.W.pars.nspecies) = 0 ;  % Y=0
       end
    end
    idx = find(current_t<tspan, 1, 'First') ;
    tspan = [current_t tspan(idx:end)] ;
end

n = uint32(size(X,2));
B = X(:,1:idivide(n,3)) ;
Z = X(:,(idivide(n,3)+1):(2*idivide(n,3))) ;
Y = X(:,(2*idivide(n,3)+1):end) ;

assert(all(isfinite(X(:)))) ;

if p.L.pars.species_pool
    B = B(:,1:end-p.W.pars.nspecies) ;
    Z = Z(:,1:end-p.W.pars.nspecies) ;
    Y = Y(:,1:end-p.W.pars.nspecies) ;
end

R = Z./(B+realmin) ;
W = Y./(B+realmin) ;

end
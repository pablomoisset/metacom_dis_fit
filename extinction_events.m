function [value,isterminal,direction] = extinction_events(~,X,p)
n = uint32(numel(X)) ;
B = X(1:idivide(n,3)) ;

BP = reshape(B,p.W.pars.nspecies,p.aux.total_sites) ;
% We incresed the extinction threshold to 1e-6 from 1e-3
value = sum(BP,2) - 1e-6 ;
direction = -ones(size(value)) ;
isterminal = ones(size(value)) ;

end
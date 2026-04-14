function rho_temp = rho_temperature_adapt(W,T_tilde)
% From Wagner, T., Schliep, E. M., North, J. S., Kundel, H., Custer, C. A., Ruzich, J. K., and Hansen, G. J. (2023).
% Predicting climate change impacts on poikilotherms using physiologically guided species abundance models.
% Proceedings of the National Academy of Sciences, 120(15):e2214199120

if any(size(W)~=size(T_tilde))
   assert(all(size(W)==size(T_tilde))) ;
end

rho_temp = zeros(size(W)) ;

idx = find(T_tilde<=W) ;

% test for revisor 1. fudge_factor should be 1 to return to the original setting
fudge_factor = 1/1.5 ;
rho_temp(idx) = exp(-(1-((T_tilde(idx)-273.15)./(W(idx)-273.15))).^2/(2/(4.1*fudge_factor))^2) ;

idx = find(W<T_tilde & (T_tilde<=W+10.2)) ;
rho_temp(idx) = 1 - ((T_tilde(idx)-W(idx))/10.2).^2 ;

end
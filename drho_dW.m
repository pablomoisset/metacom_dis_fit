function dr = drho_dW(W,T_tilde)
% From Wagner, T., Schliep, E. M., North, J. S., Kundel, H., Custer, C. A., Ruzich, J. K., and Hansen, G. J. (2023).
% Predicting climate change impacts on poikilotherms using physiologically guided species abundance models.
% Proceedings of the National Academy of Sciences, 120(15):e2214199120

assert(all(size(W)==size(T_tilde))) ;

dr = zeros(size(W)) ;

idx = find(T_tilde<=W) ;
W_idx = W(idx) ;
TT_idx = T_tilde(idx) ;
tmp = (TT_idx-273.15)./(W_idx-273.15) ;
dr(idx) = -exp(-(1-tmp).^2/(2/4.1)^2)*(4.1^2/2).*(W_idx-TT_idx).*tmp./(W_idx-273.15).^2 ;

idx = find(W<T_tilde & (T_tilde<=W+10.2)) ;
dr(idx) = (2/10.2^2)*(T_tilde(idx)-W(idx)) ;

end

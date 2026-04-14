function r = arrhenius(T, AE, T0)
   kT = 8.617E-5 ; %Boltzmann constant
   r = exp(AE.*(T0-T)./(kT*T*T0)) ;
end
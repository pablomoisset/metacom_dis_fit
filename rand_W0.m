function p_new = rand_W0(p, noise)

p_new = p ;

p_new.W_init = p.W_init + 2.0*(rand(size(p.W_init))-0.5)*noise ;
p_new.init = [p.B_init; p.B_init.*p.R_init; p.B_init.*p_new.W_init] ;

assert ( all(p_new.W_init > 0.0) ) ;
end

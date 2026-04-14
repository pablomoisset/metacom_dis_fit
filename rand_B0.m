function p_new = rand_B0(p, noise)
p_new = p ;

p_new.B_init = p.B_init .* (1.0 + 2.0*noise*(rand(size(p.B_init))-0.5)) ;
p_new.init = [p_new.B_init; p_new.B_init.*p.R_init; p_new.B_init.*p.W_init] ;

assert ( all(p_new.B_init >= 0.0) ) ;
end
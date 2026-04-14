%%{
function p_new = rand_T(p)
% p.L.pars.T_min
% p.L.pars.T_max
% p.L.pars.T_noise

assert(p.L.pars.T_noise>=0.0)

y_min = min(p.L.yP) ;
L = max(p.L.yP) - y_min ; %warning('Landscape extension fixed at 500') ;
D = p.L.pars.T_max - p.L.pars.T_min;
k = D/(L+100*realmin) ;

T = k*(p.L.yP - y_min) + p.L.pars.T_min ;

if L*p.L.pars.T_noise>0
    perm = randperm(p.L.pars.nsites) ;
    for first = perm
        p_choice = normpdf(p.L.yP, p.L.yP(first), L*p.L.pars.T_noise) ;
        second = find(mnrnd(1,p_choice/sum(p_choice))) ;
        tmp = T(first) ;
        T(first) = T(second) ;
        T(second) = tmp  ;
    end
end

p.L.T = T ;
p_new = p ;

end
%%}

%{
function p_new = rand_T(p)
% p.L.pars.T_min
% p.L.pars.T_max
% p.L.pars.T_noise

warning('Temps assigned on a cluster basis') ;

assert(p.L.pars.T_noise>=0.0)

centros = 1:p.L.pars.centers ;
if p.L.pars.nsites > 1
    DParches = squareform(pdist([p.L.xP p.L.yP],'euclidean')) ;
else
    DParches = 0 ;
end

[~, in_cluster] = min(DParches(1:end-p.L.pars.species_pool,centros),[],2) ;

idx = [] ;
for i = 1:p.L.pars.centers
    tmp = find(in_cluster==i) ;
    idx = [idx; tmp(randperm(numel(tmp)))] ;
end

if p.L.pars.nsites > 1
    T = (p.L.pars.T_max - p.L.pars.T_min)*(0:p.L.pars.nsites-1)/(p.L.pars.nsites-1) + p.L.pars.T_min ;
else
    T = (p.L.pars.T_max + p.L.pars.T_min)/2.0 ;
end

L = p.L.pars.nsites;
Y = 1:L ;
if p.L.pars.T_noise>0
    perm = randperm(p.L.pars.nsites) ;
    for first = perm
        p_choice = normpdf(Y, first, L*p.L.pars.T_noise) ;
        second = find(mnrnd(1,p_choice/sum(p_choice))) ;
        tmp = T(first) ;
        T(first) = T(second) ;
        T(second) = tmp  ;
    end
end

p.L.T(idx) = T ;
p_new = p ;

end

%}

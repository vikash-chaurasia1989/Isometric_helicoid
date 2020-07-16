function f = fun_p(p)

global tau 

[K,E]  = ellipke(p);


f = 1/p^2*(2*E/K - 1) - 4*tau^2;
end
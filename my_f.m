function vector = my_f(v,beta,gamma,delta)
vector = beta.*v+v.*(abs(v).^2*gamma+abs(v).^4*delta);
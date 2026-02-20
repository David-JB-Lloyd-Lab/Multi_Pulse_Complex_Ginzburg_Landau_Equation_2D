function Psi=my_df(V,w,beta,gamma,delta)

Psi=beta*w+2.0*gamma*abs(V).^2.0.*w+gamma*V.^2.0.*conj(w)+3.0*delta*abs(V).^4.0.*w ...
    +2.0*delta*abs(V).^2.0.*V.^2.0.*conj(w);

end
function [what,bhat]=solve_order_w_problem_2d_gmres(NN,grid_points,V_n,x2,y2,phi_g_n,phi_rx_n,phi_ry_n,psi_g_n,psi_rx_n,psi_ry_n,alpha,beta,gamma,delta,bhat);

N=2*(grid_points+1);
M=N;
tol=1E-8;

% Setup matrix to solve the real and imaginary parts of 
% w equation using finite differences

aR=real(alpha);
aI=imag(alpha);
bR=real(beta);
bI=imag(beta);
gR=real(gamma);
gI=imag(gamma);
dR=real(delta);
dI=imag(delta);

abbarR=aR*bR+aI*bI;
abbarI=-aR*bI+aI*bR;
agbarR=aR*gR+aI*gI;
agbarI=-aR*gI+aI*gR;
adbarR=aR*dR+aI*dI;
adbarI=-aR*dI+aI*dR;
absa=aR^2+aI^2;

hx=x2(2)-x2(1);
hy=y2(2)-y2(1);

Adiag=sparse(2*M*N,2*M*N);
b=zeros(2*M*N,1);
Ddiag=sparse(2*M*N,2*M*N);

icomps=zeros(1,2*(5*(M-2)*(N-2)+8*(N+M-4)+12+N*M));
jcomps=zeros(1,2*(5*(M-2)*(N-2)+8*(N+M-4)+12+N*M));
Acomps=zeros(1,2*(5*(M-2)*(N-2)+8*(N+M-4)+12+N*M));
Dcomps=zeros(1,2*(5*(M-2)*(N-2)+8*(N+M-4)+12+N*M));

total=5*(M-2)*(N-2)+8*(N+M-4)+12+N*M;

% constrct RHS of PDE

Vsum=zeros(2*grid_points+2,2*grid_points+2);

for j=1:NN
    Vsum(:,:)=Vsum(:,:)+V_n(:,:,j);
end

Phi=my_f(Vsum,beta,gamma,delta);
for j=1:NN
    Phi(:,:)=Phi(:,:)-my_f(V_n(:,:,j),beta,gamma,delta);
end

for k=1:N
    for j=1:M

    VsumR=real(Vsum(k,j));
    VsumI=imag(Vsum(k,j));
    absV=VsumR^2+VsumI^2;

    real_wr_term=abbarR+absV*agbarR+2.0*VsumR*(agbarR*VsumR+agbarI*VsumI) ...
        +absV^2*adbarR+4.0*VsumR*absV*(adbarR*VsumR+adbarI*VsumI);
    real_wi_term=abbarI+absV*agbarI+2.0*VsumI*(agbarR*VsumR+agbarI*VsumI) ...
        +absV^2*adbarI+4.0*VsumI*absV*(adbarR*VsumR+adbarI*VsumI);

    imag_wi_term=abbarR+absV*agbarR+2.0*VsumI*(-agbarI*VsumR+agbarR*VsumI) ...
        +absV^2*adbarR+4.0*VsumI*absV*(-adbarI*VsumR+adbarR*VsumI);
    imag_wr_term=-abbarI+absV*(-agbarI)+2.0*VsumR*(-agbarI*VsumR+agbarR*VsumI) ...
        +absV^2*(-adbarI)+4.0*VsumR*absV*(-adbarI*VsumR+adbarR*VsumI);

    if (k==1 && j==1)

        % Real equation - wr

        icomps(1)=(k-1)*M+j;
        icomps(2)=(k-1)*M+j;
        icomps(3)=(k-1)*M+j;

        jcomps(1)=(k-1)*M+j;
        jcomps(2)=(k-1)*M+j+1;
        jcomps(3)=(k-0)*M+j;

        Acomps(1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        Acomps(2)=1.0/hx^2.0*absa;
        Acomps(3)=1.0/hy^2.0*absa;

        Dcomps(1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        Dcomps(2)=1.0/hx^2.0*absa;
        Dcomps(3)=1.0/hy^2.0*absa;

        % Real equation - wi

        icomps(4)=(k-1)*M+j;

        jcomps(4)=(k-1)*M+j+N*M;

        Acomps(4)=real_wi_term;

        % Imag equation - wi

        icomps(total+1)=(k-1)*M+j+N*M;
        icomps(total+2)=(k-1)*M+j+N*M;
        icomps(total+3)=(k-1)*M+j+N*M;

        jcomps(total+1)=(k-1)*M+j+N*M;
        jcomps(total+2)=(k-1)*M+j+1+N*M;
        jcomps(total+3)=(k-0)*M+j+N*M;

        Acomps(total+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        Acomps(total+2)=1.0/hx^2.0*absa;
        Acomps(total+3)=1.0/hy^2.0*absa;

        Dcomps(total+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        Dcomps(total+2)=1.0/hx^2.0*absa;
        Dcomps(total+3)=1.0/hy^2.0*absa;

        % Imag equation - wr

        icomps(total+4)=(k-1)*M+j+M*N;

        jcomps(total+4)=(k-1)*M+j;

        Acomps(total+4)=imag_wr_term;


    elseif (k==1 && j==M)

        % Real equation - wr

        icomps(5*(M-2)+4+1)=(k-1)*M+j;
        icomps(5*(M-2)+4+2)=(k-1)*M+j;
        icomps(5*(M-2)+4+3)=(k-1)*M+j;

        jcomps(5*(M-2)+4+1)=(k-1)*M+j;
        jcomps(5*(M-2)+4+2)=(k-1)*M+j-1;
        jcomps(5*(M-2)+4+3)=(k-0)*M+j;

        % j,k
        Acomps(5*(M-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        % j-1,k
        Acomps(5*(M-2)+4+2)=1.0/hx^2.0*absa;
        % j,k+1
        Acomps(5*(M-2)+4+3)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(5*(M-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(5*(M-2)+4+2)=1.0/hx^2.0*absa;
        % j,k+1
        Dcomps(5*(M-2)+4+3)=1.0/hy^2.0*absa;


        % Real equation - wi

        icomps(5*(M-2)+4+4)=(k-1)*M+j;

        jcomps(5*(M-2)+4+4)=(k-1)*M+j+N*M;

        Acomps(5*(M-2)+4+4)=real_wi_term;

        % Imag equation - wi

        icomps(total+5*(M-2)+4+1)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+4+2)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+4+3)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+4+1)=(k-1)*M+j+N*M;
        jcomps(total+5*(M-2)+4+2)=(k-1)*M+j-1+N*M;
        jcomps(total+5*(M-2)+4+3)=(k-0)*M+j+N*M;

        % j,k
        Acomps(total+5*(M-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        % j-1,k
        Acomps(total+5*(M-2)+4+2)=1.0/hx^2.0*absa;
        % j,k+1
        Acomps(total+5*(M-2)+4+3)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(total+5*(M-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(total+5*(M-2)+4+2)=1.0/hx^2.0*absa;
        % j,k+1
        Dcomps(total+5*(M-2)+4+3)=1.0/hy^2.0*absa;
        
        % Imag equation - wr

        icomps(total+5*(M-2)+4+4)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+4+4)=(k-1)*M+j;

        Acomps(total+5*(M-2)+4+4)=imag_wr_term;

    elseif (j==1 && k==N)

        % Real equation - wr

        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+1)=(k-1)*M+j;
        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+2)=(k-1)*M+j;
        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+3)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+1)=(k-1)*M+j;
        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+2)=(k-1)*M+j+1;
        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+3)=(k-2)*M+j;

        % j,k
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        % j+1,k
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+2)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+3)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j+1,k
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+2)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+3)=1.0/hy^2.0*absa;
        
        % Real equation - wi

        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+4)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+4)=(k-1)*M+j+N*M;

        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+4)=real_wi_term;

        % Imag equation - wi

        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+1)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+2)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+3)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+1)=(k-1)*M+j+N*M;
        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+2)=(k-1)*M+j+1+N*M;
        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+3)=(k-2)*M+j+N*M;

        % j,k
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        % j+1,k
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+2)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+3)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j+1,k
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+2)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+3)=1.0/hy^2.0*absa;

        % Imag equation - wr

        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+4)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+4)=(k-1)*M+j;

        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+4)=imag_wr_term;

    elseif (j==M && k==N)

        % Real equation - wr

        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+1)=(k-1)*M+j;
        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+2)=(k-1)*M+j;
        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+3)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+1)=(k-1)*M+j;
        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+2)=(k-1)*M+j-1;
        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+3)=(k-2)*M+j;

        % j,k
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        % j-1,k
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+2)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+3)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+2)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+3)=1.0/hy^2.0*absa;

        % real equation - wi

        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+4)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+4)=(k-1)*M+j+N*M;

        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+4)=real_wi_term;

        % Imag equation - wi

        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+1)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+2)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+3)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+1)=(k-1)*M+j+N*M;
        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+2)=(k-1)*M+j-1+N*M;
        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+3)=(k-2)*M+j+N*M;

        % j,k
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        % j-1,k
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+2)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+3)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+2)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+3)=1.0/hy^2.0*absa;

        % imag equation - wr

        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+4)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+4)=(k-1)*M+j;

        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(M-2)+4+4)=imag_wr_term;

    elseif (k==1)

        % Real equation - wr

        icomps(5*(j-2)+4+1)=(k-1)*M+j;
        icomps(5*(j-2)+4+2)=(k-1)*M+j;
        icomps(5*(j-2)+4+3)=(k-1)*M+j;
        icomps(5*(j-2)+4+4)=(k-1)*M+j;

        jcomps(5*(j-2)+4+1)=(k-1)*M+j;
        jcomps(5*(j-2)+4+2)=(k-1)*M+j-1;
        jcomps(5*(j-2)+4+3)=(k-1)*M+j+1;
        jcomps(5*(j-2)+4+4)=(k-0)*M+j;

        % j,k
        Acomps(5*(j-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        % j-1,k
        Acomps(5*(j-2)+4+2)=1.0/hx^2.0*absa;
        % j+1,k
        Acomps(5*(j-2)+4+3)=1.0/hx^2.0*absa;
        % j,k+1
        Acomps(5*(j-2)+4+4)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(5*(j-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(5*(j-2)+4+2)=1.0/hx^2.0*absa;
        % j+1,k
        Dcomps(5*(j-2)+4+3)=1.0/hx^2.0*absa;
        % j,k+1
        Dcomps(5*(j-2)+4+4)=1.0/hy^2.0*absa;

        % Real equation - wi

        icomps(5*(j-2)+4+5)=(k-1)*M+j;

        jcomps(5*(j-2)+4+5)=(k-1)*M+j+N*M;

        Acomps(5*(j-2)+4+5)=real_wi_term;

        % Imag equation - wi

        icomps(total+5*(j-2)+4+1)=(k-1)*M+j+N*M;
        icomps(total+5*(j-2)+4+2)=(k-1)*M+j+N*M;
        icomps(total+5*(j-2)+4+3)=(k-1)*M+j+N*M;
        icomps(total+5*(j-2)+4+4)=(k-1)*M+j+N*M;

        jcomps(total+5*(j-2)+4+1)=(k-1)*M+j+N*M;
        jcomps(total+5*(j-2)+4+2)=(k-1)*M+j-1+N*M;
        jcomps(total+5*(j-2)+4+3)=(k-1)*M+j+1+N*M;
        jcomps(total+5*(j-2)+4+4)=(k-0)*M+j+N*M;

        % j,k
        Acomps(total+5*(j-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        % j-1,k
        Acomps(total+5*(j-2)+4+2)=1.0/hx^2.0*absa;
        % j+1,k
        Acomps(total+5*(j-2)+4+3)=1.0/hx^2.0*absa;
        % j,k+1
        Acomps(total+5*(j-2)+4+4)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(total+5*(j-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(total+5*(j-2)+4+2)=1.0/hx^2.0*absa;
        % j+1,k
        Dcomps(total+5*(j-2)+4+3)=1.0/hx^2.0*absa;
        % j,k+1
        Dcomps(total+5*(j-2)+4+4)=1.0/hy^2.0*absa;
        
        % Imag equation - wr

        icomps(total+5*(j-2)+4+5)=(k-1)*M+j+N*M;

        jcomps(total+5*(j-2)+4+5)=(k-1)*M+j;

        Acomps(total+5*(j-2)+4+5)=imag_wr_term;

    elseif (k==N)

        % Real equation - wr

        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+1)=(k-1)*M+j;
        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+2)=(k-1)*M+j;
        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+3)=(k-1)*M+j;
        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+4)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+1)=(k-1)*M+j;
        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+2)=(k-1)*M+j-1;
        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+3)=(k-1)*M+j+1;
        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+4)=(k-2)*M+j;

        % j,k
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        % j-1,k
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+2)=1.0/hx^2.0*absa;
        % j+1,k
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+3)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+4)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+2)=1.0/hx^2.0*absa;
        % j+1,k
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+3)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+4)=1.0/hy^2.0*absa;

        % Real equation - wi

        icomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+5)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+5)=(k-1)*M+j+N*M;

        Acomps(5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+5)=real_wi_term;

        % imag equation - wi

        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+1)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+2)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+3)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+4)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+1)=(k-1)*M+j+N*M;
        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+2)=(k-1)*M+j-1+N*M;
        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+3)=(k-1)*M+j+1+N*M;
        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+4)=(k-2)*M+j+N*M;

        % j,k
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        % j-1,k
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+2)=1.0/hx^2.0*absa;
        % j+1,k
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+3)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+4)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+2)=1.0/hx^2.0*absa;
        % j+1,k
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+3)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+4)=1.0/hy^2.0*absa;
        
        % imag equation - wr

        icomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+5)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+5)=(k-1)*M+j;

        Acomps(total+5*(M-2)+8+(N-2)*(10+6*(M-2))+5*(j-2)+4+5)=imag_wr_term;

    elseif (j==1)

        % Real equation - wr

        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+1)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+2)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+3)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+4)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+1)=(k-1)*M+j;
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+2)=(k-1)*M+j+1;
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+3)=(k-2)*M+j;
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+4)=(k-0)*M+j;

        % j,k
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        % j+1,k
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+2)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+3)=1.0/hy^2.0*absa;
        % j,k+1
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+4)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j+1,k
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+2)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+3)=1.0/hy^2.0*absa;
        % j,k+1
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+4)=1.0/hy^2.0*absa;

        % Real equation - wi

        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5)=(k-1)*M+j+N*M;

        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5)=real_wi_term;

        % Imag equation - wi

        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+1)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+2)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+3)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+4)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+1)=(k-1)*M+j+N*M;
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+2)=(k-1)*M+j+1+N*M;
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+3)=(k-2)*M+j+N*M;
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+4)=(k-0)*M+j+N*M;

        % j,k
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        % j+1,k
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+2)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+3)=1.0/hy^2.0*absa;
        % j,k+1
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+4)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j+1,k
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+2)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+3)=1.0/hy^2.0*absa;
        % j,k+1
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+4)=1.0/hy^2.0*absa;

        % Imag equation - wr

        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5)=(k-1)*M+j;

        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5)=imag_wr_term;

    elseif (j==M)

        % Real equation - wr
            
        icomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+1)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+2)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+3)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+4)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+1)=(k-1)*M+j;
        jcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+2)=(k-1)*M+j-1;
        jcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+3)=(k-2)*M+j;
        jcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+4)=(k-0)*M+j;

        % j,k
        Acomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        % j-1,k
        Acomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+2)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+3)=1.0/hy^2.0*absa;
        % j,k+1
        Acomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+4)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+2)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+3)=1.0/hy^2.0*absa;
        % j,k+1
        Dcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+4)=1.0/hy^2.0*absa;
        
        % Real equation - wi
            
        icomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+5)=(k-1)*M+j;

        jcomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+5)=(k-1)*M+j+N*M;

        Acomps(5*(M-2)+8+(k-1)*(10+6*(M-2))-5+5)=real_wi_term;

        % Imag equation - wi
            
        icomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+1)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+2)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+3)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+4)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+1)=(k-1)*M+j+N*M;
        jcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+2)=(k-1)*M+j-1+N*M;
        jcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+3)=(k-2)*M+j+N*M;
        jcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+4)=(k-0)*M+j+N*M;

        % j,k
        Acomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        % j-1,k
        Acomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+2)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+3)=1.0/hy^2.0*absa;
        % j,k+1
        Acomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+4)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+2)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+3)=1.0/hy^2.0*absa;
        % j,k+1
        Dcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+4)=1.0/hy^2.0*absa;
        
        % Imag equation - wr
            
        icomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+5)=(k-1)*M+j+N*M;

        jcomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+5)=(k-1)*M+j;
        
        Acomps(total+5*(M-2)+8+(k-1)*(10+6*(M-2))-5+5)=imag_wr_term;

    else

        % Real equation - wr

        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+1)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+2)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+3)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+4)=(k-1)*M+j;
        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+5)=(k-1)*M+j;
        
        % j,k
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+1)=(k-1)*M+j;
        % j-1,k
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+2)=(k-1)*M+j-1;
        % j+1,k
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+3)=(k-1)*M+j+1;
        % j,k-1
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+4)=(k-2)*M+j;
        % j,k+1
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+5)=(k+0)*M+j;

        % j,k
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+real_wr_term;
        % j-1,k
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+2)=1.0/hx^2.0*absa;
        % j+1,k
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+3)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+4)=1.0/hy^2.0*absa;
        % j,k+1
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+5)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+2)=1.0/hx^2.0*absa;
        % j+1,k
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+3)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+4)=1.0/hy^2.0*absa;
        % j,k+1
        Dcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+5)=1.0/hy^2.0*absa;
        
        % Real equation - wi

        icomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+6)=(k-1)*M+j;
        
        % j,k
        jcomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+6)=(k-1)*M+j+N*M;

        % j,k
        Acomps(5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+6)=real_wi_term;

        % Imag equation - wi

        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+1)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+2)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+3)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+4)=(k-1)*M+j+N*M;
        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+5)=(k-1)*M+j+N*M;
        
        % j,k
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+1)=(k-1)*M+j+N*M;
        % j-1,k
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+2)=(k-1)*M+j-1+N*M;
        % j+1,k
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+3)=(k-1)*M+j+1+N*M;
        % j,k-1
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+4)=(k-2)*M+j+N*M;
        % j,k+1
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+5)=(k+0)*M+j+N*M;

        % j,k
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa+imag_wi_term;
        % j-1,k
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+2)=1.0/hx^2.0*absa;
        % j+1,k
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+3)=1.0/hx^2.0*absa;
        % j,k-1
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+4)=1.0/hy^2.0*absa;
        % j,k+1
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+5)=1.0/hy^2.0*absa;

        % j,k
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+1)=-2.0*(1.0/hx^2+1.0/hy^2)*absa;
        % j-1,k
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+2)=1.0/hx^2.0*absa;
        % j+1,k
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+3)=1.0/hx^2.0*absa;
        % j,k-1
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+4)=1.0/hy^2.0*absa;
        % j,k+1
        Dcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+5)=1.0/hy^2.0*absa;
        
        % Imag equation - wr

        icomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+6)=(k-1)*M+j+N*M;
        
        % j,k
        jcomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+6)=(k-1)*M+j;

        % j,k
        Acomps(total+5*(M-2)+8+(k-2)*(10+6*(M-2))+5+6*(j-2)+6)=imag_wr_term;

    end
        % Real equation
        b((k-1)*M+j)=aR*real(-Phi(k,j))+aI*imag(-Phi(k,j));
        % Imag equation
        b((k-1)*M+j+M*N)=aR*imag(-Phi(k,j))-aI*real(-Phi(k,j));
    end
end

Adiag=sparse(icomps,jcomps,Acomps);
Ddiag=sparse(icomps,jcomps,Dcomps);

[z,output2,rels2,it2,resvec2]=gmres(Adiag,b,500,tol,[],Ddiag,[],bhat);

for k=1:N
    for j=1:M

        what(k,j)=z((k-1)*M+j)+1i*z((k-1)*M+j+N*M);

    end
end

bhat=b;

end
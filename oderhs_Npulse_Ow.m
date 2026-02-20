function dydt = oderhs_Npulse_Ow(t,y,alpha,beta,gamma,delta,N,L,V,phi_g,phi_rx,phi_ry,psi_g,psi_rx,psi_ry)

dydt = zeros(3*N,1);
F    = zeros(3*N,1);
H = zeros(3*N,1);
C = zeros(3*N,3*N);

rx=zeros(N,1);
ry=zeros(N,1);
g=zeros(N,1);

grid_points=150;

bhat=zeros(2*(2*(grid_points+1))^2,1);
V_n=zeros(2*grid_points+2,2*grid_points+2,N);

phi_g_n=zeros(2*grid_points+2,2*grid_points+2,N);
phi_rx_n=zeros(2*grid_points+2,2*grid_points+2,N);
phi_ry_n=zeros(2*grid_points+2,2*grid_points+2,N);

psi_g_n=zeros(2*grid_points+2,2*grid_points+2,N);
psi_rx_n=zeros(2*grid_points+2,2*grid_points+2,N);
psi_ry_n=zeros(2*grid_points+2,2*grid_points+2,N);

dpsi_g_x_n=zeros(2*grid_points+2,2*grid_points+2,N);
dpsi_g_y_n=zeros(2*grid_points+2,2*grid_points+2,N);

dpsi_rx_x_n=zeros(2*grid_points+2,2*grid_points+2,N);
dpsi_rx_y_n=zeros(2*grid_points+2,2*grid_points+2,N);

dpsi_ry_x_n=zeros(2*grid_points+2,2*grid_points+2,N);
dpsi_ry_y_n=zeros(2*grid_points+2,2*grid_points+2,N);

what=zeros(2*grid_points+2,2*grid_points+2);

for k=1:N
rx(k)= y(3*(k-1)+1);
ry(k)= y(3*(k-1)+2);
g(k)= y(3*(k-1)+3);
end

% Form maximum area upon which w can be calculated
dx=(min(L,L+min(rx))-max(-L,-L+max(rx)))/(2*grid_points+1);
x2 = [max(-L,-L+max(rx)):dx:min(L,L+min(rx))+dx/2];
dy=(min(L,L+min(ry))-max(-L,-L+max(ry)))/(2*grid_points+1);
y2 = [max(-L,-L+max(ry)):dy:min(L,L+min(ry))+dy/2]';

% Differentiate adjoint eigenmodes w.r.t x and y
dpsi_g_x=diffx(psi_g,1);
dpsi_g_y=diffy(psi_g,1);

dpsi_rx_x=diffx(psi_rx,1);
dpsi_rx_y=diffy(psi_rx,1);

dpsi_ry_x=diffx(psi_ry,1);
dpsi_ry_y=diffy(psi_ry,1);

for k=1:N

%    k

    % pulses numerically evalued on grid
    V_n(:,:,k)=exp(1i*g(k)).*V(x2-rx(k),y2-ry(k)); 

    % eigenmodes
    phi_g_n(:,:,k)=exp(1i*g(k)).*phi_g(x2-rx(k),y2-ry(k));
    phi_rx_n(:,:,k)=exp(1i*g(k)).*phi_rx(x2-rx(k),y2-ry(k));
    phi_ry_n(:,:,k)=exp(1i*g(k)).*phi_ry(x2-rx(k),y2-ry(k));

    % adjoint eigenmodes
    psi_g_n(:,:,k)=exp(1i*g(k)).*psi_g(x2-rx(k),y2-ry(k));
    psi_rx_n(:,:,k)=exp(1i*g(k)).*psi_rx(x2-rx(k),y2-ry(k));
    psi_ry_n(:,:,k)=exp(1i*g(k)).*psi_ry(x2-rx(k),y2-ry(k));

    % derivatives of adjoint eigenmodes
    dpsi_g_x_n(:,:,k)=exp(1i*g(k)).*dpsi_g_x(x2-rx(k),y2-ry(k));
    dpsi_g_y_n(:,:,k)=exp(1i*g(k)).*dpsi_g_y(x2-rx(k),y2-ry(k));

    dpsi_rx_x_n(:,:,k)=exp(1i*g(k)).*dpsi_rx_x(x2-rx(k),y2-ry(k));
    dpsi_rx_y_n(:,:,k)=exp(1i*g(k)).*dpsi_rx_y(x2-rx(k),y2-ry(k));

    dpsi_ry_x_n(:,:,k)=exp(1i*g(k)).*dpsi_ry_x(x2-rx(k),y2-ry(k));
    dpsi_ry_y_n(:,:,k)=exp(1i*g(k)).*dpsi_ry_y(x2-rx(k),y2-ry(k));

end    

% % Solve linear problem for what, and then remove multiples of homogeneous
% % solution


[what,bhat]=solve_order_w_problem_2d_gmres(N,grid_points,V_n,x2,y2,phi_g_n,phi_rx_n,phi_ry_n,psi_g_n,psi_rx_n,psi_ry_n, ...
    alpha,beta,gamma,delta,bhat);

% % Set up problem to subtract of multiple of homogeneous solution.
% 
for j=1:N
    for k=1:N
        if (j==k)
        C(3*j-2,3*k-2)=1.0;
        C(3*j-1,3*k-2)=0.0;
        C(3*j-0,3*k-2)=0.0;

        C(3*j-2,3*k-1)=0.0;
        C(3*j-1,3*k-1)=1.0;
        C(3*j-0,3*k-1)=0.0;

        C(3*j-2,3*k-0)=0.0;
        C(3*j-1,3*k-0)=0.0;        
        C(3*j-0,3*k-0)=1.0;
        else
        C(3*j-2,3*k-2)=dx*dy*real(sum(sum(phi_rx_n(:,:,k).*conj(psi_rx_n(:,:,j)))));
        C(3*j-1,3*k-2)=dx*dy*real(sum(sum(phi_rx_n(:,:,k).*conj(psi_ry_n(:,:,j)))));
        C(3*j-0,3*k-2)=dx*dy*real(sum(sum(phi_rx_n(:,:,k).*conj(psi_g_n(:,:,j)))));

        C(3*j-2,3*k-1)=dx*dy*real(sum(sum(phi_ry_n(:,:,k).*conj(psi_rx_n(:,:,j)))));
        C(3*j-1,3*k-1)=dx*dy*real(sum(sum(phi_ry_n(:,:,k).*conj(psi_ry_n(:,:,j)))));
        C(3*j-0,3*k-1)=dx*dy*real(sum(sum(phi_ry_n(:,:,k).*conj(psi_g_n(:,:,j)))));

        C(3*j-2,3*k-0)=dx*dy*real(sum(sum(phi_g_n(:,:,k).*conj(psi_rx_n(:,:,j)))));
        C(3*j-1,3*k-0)=dx*dy*real(sum(sum(phi_g_n(:,:,k).*conj(psi_ry_n(:,:,j)))));
        C(3*j-0,3*k-0)=dx*dy*real(sum(sum(phi_g_n(:,:,k).*conj(psi_g_n(:,:,j)))));
        end
    end
    F(3*j-2)=dx*dy*real(sum(sum(what(:,:).*conj(psi_rx_n(:,:,j)))));
    F(3*j-1)=dx*dy*real(sum(sum(what(:,:).*conj(psi_ry_n(:,:,j)))));
    F(3*j-0)=dx*dy*real(sum(sum(what(:,:).*conj(psi_g_n(:,:,j)))));
end

 H = C\F;

 w=what;
 for j=1:N
     w(:,:)=w(:,:)-(H(3*j-2)*phi_rx_n(:,:,j)+H(3*j-1)*phi_ry_n(:,:,j) ...
         +H(3*j-0)*phi_g_n(:,:,j));
 end

[Phi,G,Psi]=cap_functions(N,grid_points,V_n,w,beta,gamma,delta); 

% Extend matrix C to find RHS of ODE for time integration.

for j=1:N
    C(3*j-2,3*j-2)=C(3*j-2,3*j-2)+dx*dy*real(sum(sum(w(:,:).*conj(dpsi_rx_x_n(:,:,j)))));
    C(3*j-1,3*j-2)=C(3*j-1,3*j-2)+dx*dy*real(sum(sum(w(:,:).*conj(dpsi_rx_y_n(:,:,j)))));
    C(3*j-0,3*j-2)=C(3*j-0,3*j-2)+dx*dy*real(sum(sum(1i*w(:,:).*conj(psi_rx_n(:,:,j)))));
        
    C(3*j-2,3*j-1)=C(3*j-2,3*j-1)+dx*dy*real(sum(sum(w(:,:).*conj(dpsi_ry_x_n(:,:,j)))));
    C(3*j-1,3*j-1)=C(3*j-1,3*j-1)+dx*dy*real(sum(sum(w(:,:).*conj(dpsi_ry_y_n(:,:,j)))));
    C(3*j-0,3*j-1)=C(3*j-0,3*j-1)+dx*dy*real(sum(sum(1i*w(:,:).*conj(psi_ry_n(:,:,j)))));

    C(3*j-2,3*j-0)=C(3*j-2,3*j-0)+dx*dy*real(sum(sum(w(:,:).*conj(dpsi_g_x_n(:,:,j)))));
    C(3*j-1,3*j-0)=C(3*j-1,3*j-0)+dx*dy*real(sum(sum(w(:,:).*conj(dpsi_g_y_n(:,:,j)))));     
    C(3*j-0,3*j-0)=C(3*j-0,3*j-0)+dx*dy*real(sum(sum(1i*w(:,:).*conj(psi_g_n(:,:,j)))));

    F(3*j-2)=dx*dy*real(sum(sum((Phi(:,:)+G(:,:)+Psi(:,:,j)).*conj(psi_rx_n(:,:,j)))));
    F(3*j-1)=dx*dy*real(sum(sum((Phi(:,:)+G(:,:)+Psi(:,:,j)).*conj(psi_ry_n(:,:,j)))));
    F(3*j-0)=dx*dy*real(sum(sum((Phi(:,:)+G(:,:)+Psi(:,:,j)).*conj(psi_g_n(:,:,j)))));
end

H = C\F;

for j=1:3*N
    dydt(j,1)=H(j);
end

t

end 
 
 
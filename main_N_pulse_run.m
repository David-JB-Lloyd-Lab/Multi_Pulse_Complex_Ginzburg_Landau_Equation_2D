addpath('../chebfun')

load pulsedata2 % bR=-0.05
%load pulsedata3 % bR=-2.0

options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on','Events',@(t,y,flag) events_period2(t,y));

r0=2.56;
N=2;

tspan=[0 10^5]; 
V = ufull+1i.*vfull; 
gam = gR+1i.*gI; 
d = dR+1i.*dI; 
beta=bR+1i*bI;
alpha=aR+1i*aI;

y0=zeros(3*N,1);

% set intial condition (r_x1,r_y1,g_1,.....,r_xN,r_yN,g_N)
% Case N=2
y0=[-r0/2;0;pi/2+1E-12;r0/2;0;0];

[t,y] = ode15s(@(t,y) oderhs_Npulse_Ow(t,y,alpha,beta,gam,d,N,LL,V,phi0full,phi1cfull,phi1sfull,psi0full,psi1cfull,psi1sfull),tspan, y0,options);

ss=size(y);

rx=zeros(ss(1),N);
ry=zeros(ss(1),N);
g=zeros(ss(1),N);

for j=1:ss(1)
    for k=1:N
        rx(j,k)= y(j,3*(k-1)+1);
        ry(j,k)= y(j,3*(k-1)+2);
        g(j,k)= y(j,3*(k-1)+3);
    end
end

figure(2); plot(t(:),rx(:,1),'b-',t(:),rx(:,2),'r-')
figure(3); plot(t(:),ry(:,1),'b-',t(:),ry(:,2),'r-')
figure(4); plot(t(:),g(:,1),'b-',t(:),g(:,2),'r-')

figure(5);
plot(abs(rx(:,1)-rx(:,2)).*cos(g(:,1)-g(:,2)),abs(rx(:,1)-rx(:,2)).*sin(g(:,1)-g(:,2)),'b-','Linewidth',2)
hold on
plot(abs(rx(1,1)-rx(1,2)).*cos(g(1,1)-g(1,2)),abs(rx(1,1)-rx(1,2)).*sin(g(1,1)-g(1,2)),'ro','Linewidth',2)
plot(abs(rx(end,1)-rx(end,2)).*cos(g(end,1)-g(end,2)),abs(rx(end,1)-rx(end,2)).*sin(g(end,1)-g(end,2)),'gs','Linewidth',2)
hold off

axis ([-1.4 1.4 0.0 2.56])


drawnow

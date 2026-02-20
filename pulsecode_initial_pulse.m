clear
clc

addpath('../chebfun')

load initial_data_br02_bi7379_half
cheboppref.setDefaults('factory');
chebfunpref.setDefaults('factory');
cheboppref.setDefaults('plotting', 0.01,'bvpTol',1e-14);

bI_init=-37.97;
bI_init=-37.867;
omegar_init=real(sqrt(-(bR+1i*bI_init)/(aR+1i*aI)));
omegai_init=-imag(sqrt(-(bR+1i*bI_init)/(aR+1i*aI)));
aR = 1/2; aI = (1/2); bR = -0.02; gR = 1.8; gI = 1; dR = -0.05; dI = 0.05; % CGL params
LL= 10;dom = [0, LL]; x = chebfun('x', dom); N = chebop(0,LL);    % spatial grid                                                    % define Chebop and prefs
N = chebop(@(x,u,v,bI,omegar,omegai) ...
    [aR.*diff(u,2)-aI.*diff(v,2)+bR.*u-bI.*v+(u.^2+v.^2).*(gR.*u-gI.*v)+(u.^2+v.^2).^2.*(dR.*u-dI.*v)+c.*diff(u,1); ...
     aR.*diff(v,2)+aI.*diff(u,2)+bR.*v+bI.*u+(u.^2+v.^2).*(gR.*v+gI.*u)+(u.^2+v.^2).^2.*(dR.*v+dI.*u)+c.*diff(v,1)], [0, LL]);
N.bc = @(x,u,v,bI,omegar,omegai) [feval(diff(u,1),0);feval(diff(u,1),LL)+(omegar.*u(LL)+omegai.*v(LL));feval(diff(v,1),0);feval(diff(v,1),LL)-(-omegar.*v(LL)+omegai.*u(LL));sum(sech(x).*v);aR*omegar^2-aR*omegai^2+2*aI*omegar*omegai+bR;-2*aR*omegar*omegai+aI*omegar^2-aI*omegai^2+bI]; %sum(-sech(x).*tanh(x).*(u-sech(x)))-1
N.init = [ur3;ui3;bI_init;omegar_init;omegai_init];
tic, [u,v,bI,omegar,omegai] = N\0  ;t = toc; % solve Nonlinear system
bI
norm(N(x,u,v,bI,omegar,omegai),inf) % error

for two_dim=0.1:0.1:1.0

bI_init=bI;
omegar_init=omegar;
omegai_init=omegai;
ur3=u;
ui3=v;
    
aR = 1/2; aI = (1/2); gR = 1.8; gI = 1; dR = -0.05; dI = 0.05; % CGL params    
    
N = chebop(@(x,u,v,bI,omegar,omegai) ...
    [aR.*(x.*diff(u,2)+two_dim.*diff(u,1))-aI.*(x.*diff(v,2)+two_dim.*diff(v,1))+x.*bR.*u-x.*bI.*v+x.*(u.^2+v.^2).*(gR.*u-gI.*v)+x.*(u.^2+v.^2).^2.*(dR.*u-dI.*v)+x.*c.*diff(u,1); ...
     aR.*(x.*diff(v,2)+two_dim.*diff(v,1))+aI.*(x.*diff(u,2)+two_dim.*diff(u,1))+x.*bR.*v+x.*bI.*u+x.*(u.^2+v.^2).*(gR.*v+gI.*u)+x.*(u.^2+v.^2).^2.*(dR.*v+dI.*u)+x.*c.*diff(v,1)], [0, LL]);
N.bc = @(x,u,v,bI,omegar,omegai) [feval(diff(u,1),0);feval(diff(u,1),LL)+(omegar.*u(LL)+omegai.*v(LL));feval(diff(v,1),0);feval(diff(v,1),LL)-(-omegar.*v(LL)+omegai.*u(LL));sum(sech(x).*v);aR*omegar^2-aR*omegai^2+2*aI*omegar*omegai+bR;-2*aR*omegar*omegai+aI*omegar^2-aI*omegai^2+bI]; %sum(-sech(x).*tanh(x).*(u-sech(x)))-1
N.init = [ur3;ui3;bI_init;omegar_init;omegai_init];
tic, [u,v,bI,omegar,omegai] = N\0  ;t = toc; % solve Nonlinear system
bI
norm(N(x,u,v,bI,omegar,omegai),inf) % error
   
Nx=5000;
xx=[0:Nx]*LL/Nx;
u_temp=u(xx);
v_temp=v(xx);

ufull_temp = [u_temp(end:-1:2).';u_temp.'];
vfull_temp = [v_temp(end:-1:2).';v_temp.'];

%chebfunpref.setDefaults('splitting',true)

domfull = [-LL, LL]; xfull = chebfun('xfull', domfull);
ufull = chebfun(ufull_temp,domfull,'equi'); 
vfull = chebfun(vfull_temp,domfull,'equi');

figure(10); hold on
plot(xfull,ufull,'k-')
axis([-1.0 1.0 0 6])
hold off
figure(20); hold on
plot(xfull,vfull,'k-')
axis([-1.0 1.0 -1 2.5])
hold off

end

cheboppref.setDefaults('discretization',@ultraS);
m=1;
NN = chebop(@(x,u1,v1) ...
    [aR.*(x.^2.*diff(u1,2)+x.*diff(u1,1)-m^2.*u1)-aI.*(x.^2.*diff(v1,2)+x.*diff(v1,1)-m^2.*v1)+x.^2.*bR.*u1-x.^2.*bI.*v1+x.^2.*(u.^2+v.^2).*(gR.*u1-gI.*v1)+2.0*x.^2.*(u.*u1+v.*v1).*(gR.*u-gI.*v)+x.^2.*(u.^2+v.^2).^2.*(dR.*u1-dI.*v1)+4*x.^2.*(u.^2+v.^2).*(u.*u1+v.*v1).*(dR.*u-dI.*v); ...
     aI.*(x.^2.*diff(u1,2)+x.*diff(u1,1)-m^2.*u1)+aR.*(x.^2.*diff(v1,2)+x.*diff(v1,1)-m^2.*v1)+x.^2.*bI.*u1+x.^2.*bR.*v1+x.^2.*(u.^2+v.^2).*(gI.*u1+gR.*v1)+2.0*x.^2.*(u.*u1+v.*v1).*(gI.*u+gR.*v)+x.^2.*(u.^2+v.^2).^2.*(dI.*u1+dR.*v1)+4*x.^2.*(u.^2+v.^2).*(u.*u1+v.*v1).*(dI.*u+dR.*v)], [0, LL]);
NN.bc=@(x,u1,v1) [u1(0);u1(LL);v1(LL);v1(0)];
[Vodd,Dodd1]=eigs(NN,1,0);
eigs1odd = diag(Dodd1);
eigs1odd

NN = chebop(@(x,u1,v1) ...
    [aR.*(x.*diff(u1,2)+diff(u1,1))-aI.*(x.*diff(v1,2)+diff(v1,1))+x.*bR.*u1-x.*bI.*v1+x.*(u.^2+v.^2).*(gR.*u1-gI.*v1)+2.0*x.*(u.*u1+v.*v1).*(gR.*u-gI.*v)+x.*(u.^2+v.^2).^2.*(dR.*u1-dI.*v1)+4*x.*(u.^2+v.^2).*(u.*u1+v.*v1).*(dR.*u-dI.*v); ...
     aI.*(x.*diff(u1,2)+diff(u1,1))+aR.*(x.*diff(v1,2)+diff(v1,1))+x.*bI.*u1+x.*bR.*v1+x.*(u.^2+v.^2).*(gI.*u1+gR.*v1)+2.0*x.*(u.*u1+v.*v1).*(gI.*u+gR.*v)+x.*(u.^2+v.^2).^2.*(dI.*u1+dR.*v1)+4*x.*(u.^2+v.^2).*(u.*u1+v.*v1).*(dI.*u+dR.*v)], [0, LL]);
NN.bc=@(x,u1,v1) [feval(diff(u1),0);u1(LL);v1(LL);feval(diff(v1),0)];
[Veven,Deven1]=eigs(NN,1,0);
eigs1even = diag(Deven1);
eigs1even

% Now solve adjoint problem
NN = chebop(@(x,u1,v1) ...
    [aR.*(x.^2.*diff(u1,2)+x.*diff(u1,1)-m^2.*u1)+aI.*(x.^2.*diff(v1,2)+x.*diff(v1,1)-m^2.*v1)+x.^2.*bR.*u1+x.^2.*bI.*v1+x.^2.*(u.^2+v.^2).*(gR.*u1+gI.*v1)+2.0*x.^2.*u.*((gR.*u-gI.*v).*u1+(gI.*u+gR.*v).*v1)+x.^2.*(u.^2+v.^2).^2.*(dR.*u1+dI.*v1)+4*x.^2.*(u.^2+v.^2).*u.*((dR.*u-dI.*v).*u1+(dI.*u+dR.*v).*v1); ...
     -aI.*(x.^2.*diff(u1,2)+x.*diff(u1,1)-m^2.*u1)+aR.*(x.^2.*diff(v1,2)+x.*diff(v1,1)-m^2.*v1)-x.^2.*bI.*u1+x.^2.*bR.*v1+x.^2.*(u.^2+v.^2).*(-gI.*u1+gR.*v1)+2.0*x.^2.*v.*((gR.*u-gI.*v).*u1+(gI.*u+gR.*v).*v1)+x.^2.*(u.^2+v.^2).^2.*(-dI.*u1+dR.*v1)+4*x.^2.*(u.^2+v.^2).*v.*((dR.*u-dI.*v).*u1+(dI.*u+dR.*v).*v1)], [0, LL]);
NN.bc=@(x,u1,v1) [u1(0);u1(LL);v1(LL);v1(0)];
[Wodd,Dodd2] = eigs(NN,1,0);
eigs2odd = diag(Dodd2);
eigs2odd

NN = chebop(@(x,u1,v1) ...
    [aR.*(x.*diff(u1,2)+diff(u1,1))+aI.*(x.*diff(v1,2)+diff(v1,1))+x.*bR.*u1+x.*bI.*v1+x.*(u.^2+v.^2).*(gR.*u1+gI.*v1)+2.0*x.*u.*((gR.*u-gI.*v).*u1+(gI.*u+gR.*v).*v1)+x.*(u.^2+v.^2).^2.*(dR.*u1+dI.*v1)+4*x.*(u.^2+v.^2).*u.*((dR.*u-dI.*v).*u1+(dI.*u+dR.*v).*v1); ...
     -aI.*(x.*diff(u1,2)+diff(u1,1))+aR.*(x.*diff(v1,2)+diff(v1,1))-x.*bI.*u1+x.*bR.*v1+x.*(u.^2+v.^2).*(-gI.*u1+gR.*v1)+2.0*x.*v.*((gR.*u-gI.*v).*u1+(gI.*u+gR.*v).*v1)+x.*(u.^2+v.^2).^2.*(-dI.*u1+dR.*v1)+4*x.*(u.^2+v.^2).*v.*((dR.*u-dI.*v).*u1+(dI.*u+dR.*v).*v1)], [0, LL]);
NN.bc=@(x,u1,v1) [feval(diff(u1),0);u1(LL);v1(LL);feval(diff(v1),0)];
[Weven,Deven2] = eigs(NN,1,0);
eigs2even = diag(Deven2);
eigs2even

eig1 = Vodd(:,1); eig2 = Veven(:,1); 
eig1 = chebfun(eig1); eig2 = chebfun(eig2);
eig1r=eig1(:,1);eig1i=eig1(:,2);eig2r=eig2(:,1);eig2i=eig2(:,2);
phi1= eig1r+1i*eig1i;
phi0= eig2r+1i*eig2i;

A=sum(phi0*conj(phi0))/sum((u+1i*v)*(u-1i*v));
phi0=phi0/(sqrt(abs(A)))*sign(real(phi0(0)))*sign(-v(0));

ud=diff(u,1);
vd=diff(v,1);
A=sum(phi1*conj(phi1))/sum((ud+1i*vd)*(ud-1i*vd));
phi1=-phi1/sqrt(abs(A))*sign(real(phi1(0.1)))*sign(ud(0.1));

eig1 = Wodd(:,1); eig2 = Weven(:,1); 
eig1 = chebfun(eig1); eig2 = chebfun(eig2); 
eig1 = eig1(:,1)+1i*eig1(:,2); eig2 = eig2(:,1)+1i*eig2(:,2);
psi1=eig1;
psi0=eig2;

% Now convert the chebfuns on the polar domain [0,L] to the 2D domain
% [-L,L]x[-L,L]

LL=7;
Nx=1000;
Ny=1000;
xx=[0:Nx]*LL/Nx;
yy=[0:Ny]'*LL/Ny;

phi1c_temp=phi1(sqrt(xx.^2+yy.^2)).*cos(atan2(yy,xx));
phi1s_temp=phi1(sqrt(xx.^2+yy.^2)).*sin(atan2(yy,xx));
phi0_temp=phi0(sqrt(xx.^2+yy.^2)).*exp(0*1i*atan2(yy,xx));

psi1c_temp=psi1(sqrt(xx.^2+yy.^2)).*cos(atan2(yy,xx));
psi1s_temp=psi1(sqrt(xx.^2+yy.^2)).*sin(atan2(yy,xx));
psi0_temp=psi0(sqrt(xx.^2+yy.^2)).*exp(0*1i*atan2(yy,xx));

u_temp=u(sqrt(xx.^2+yy.^2)).*exp(0*1i*atan2(yy,xx));
v_temp=v(sqrt(xx.^2+yy.^2)).*exp(0*1i*atan2(yy,xx));

ufull_temp=[u_temp(end:-1:2,end:-1:2) u_temp(end:-1:2,1:end); ...
    u_temp(1:end,end:-1:2) u_temp(1:end,1:end)];
vfull_temp=[v_temp(end:-1:2,end:-1:2) v_temp(end:-1:2,1:end); ...
    v_temp(1:end,end:-1:2) v_temp(1:end,1:end)];

phi0full_temp=[phi0_temp(end:-1:2,end:-1:2) phi0_temp(end:-1:2,1:end); ...
    phi0_temp(1:end,end:-1:2) phi0_temp(1:end,1:end)];
phi1cfull_temp=[-phi1c_temp(end:-1:2,end:-1:2) phi1c_temp(end:-1:2,1:end); ...
    -phi1c_temp(1:end,end:-1:2) phi1c_temp(1:end,1:end)];
phi1sfull_temp=[-phi1s_temp(end:-1:2,end:-1:2) -phi1s_temp(end:-1:2,1:end); ...
    phi1s_temp(1:end,end:-1:2) phi1s_temp(1:end,1:end)];

psi0full_temp=[psi0_temp(end:-1:2,end:-1:2) psi0_temp(end:-1:2,1:end); ...
    psi0_temp(1:end,end:-1:2) psi0_temp(1:end,1:end)];
psi1cfull_temp=[-psi1c_temp(end:-1:2,end:-1:2) psi1c_temp(end:-1:2,1:end); ...
    -psi1c_temp(1:end,end:-1:2) psi1c_temp(1:end,1:end)];
psi1sfull_temp=[-psi1s_temp(end:-1:2,end:-1:2) -psi1s_temp(end:-1:2,1:end); ...
    psi1s_temp(1:end,end:-1:2) psi1s_temp(1:end,1:end)];

domfull = [-LL LL -LL LL]; xyfull = chebfun2('xyfull', domfull);
phi1cfull = chebfun2(phi1cfull_temp,domfull,'equi'); 
phi1sfull = chebfun2(phi1sfull_temp,domfull,'equi'); 
phi0full = chebfun2(phi0full_temp,domfull,'equi'); 

psi1cfull = chebfun2(psi1cfull_temp,domfull,'equi'); 
psi1sfull = chebfun2(psi1sfull_temp,domfull,'equi'); 
psi0full = chebfun2(psi0full_temp,domfull,'equi'); 

ufull = chebfun2(ufull_temp,domfull,'equi'); 
vfull = chebfun2(vfull_temp,domfull,'equi');

dx=(min(LL,LL)-max(-LL,-LL))/(2*1000+1);
xx2 = [max(-LL,-LL):dx:min(LL,LL)];
dy=(min(LL,LL)-max(-LL,-LL))/(2*1000+1);
yy2 = [max(-LL,-LL):dx:min(LL,LL)]';

phi1cfull_n=phi1cfull(xx2,yy2);
psi1cfull_n=psi1cfull(xx2,yy2);
phi1sfull_n=phi1sfull(xx2,yy2);
psi1sfull_n=psi1sfull(xx2,yy2);
phi0full_n=phi0full(xx2,yy2);
psi0full_n=psi0full(xx2,yy2);
ufull_n=ufull(xx2,yy2);
vfull_n=vfull(xx2,yy2);

A = dx*dy*real(sum(sum(phi1cfull_n.*conj(psi1cfull_n)))); psi1cfull = psi1cfull/A;
A = dx*dy*real(sum(sum(phi1sfull_n.*conj(psi1sfull_n)))); psi1sfull = psi1sfull/A;
A = dx*dy*real(sum(sum(phi0full_n.*conj(psi0full_n)))); psi0full = psi0full/A;

psi1cfull_n=psi1cfull(xx2,yy2);
psi1sfull_n=psi1sfull(xx2,yy2);
psi0full_n=psi0full(xx2,yy2);

dx*dy*sum(sum(phi1cfull_n.*conj(psi1cfull_n)))
dx*dy*sum(sum(phi1cfull_n.*conj(psi1sfull_n)))
dx*dy*sum(sum(phi1cfull_n.*conj(psi0full_n)))

dx*dy*sum(sum(phi1sfull_n.*conj(psi1cfull_n)))
dx*dy*sum(sum(phi1sfull_n.*conj(psi1sfull_n)))
dx*dy*sum(sum(phi1sfull_n.*conj(psi0full_n)))

dx*dy*sum(sum(phi0full_n.*conj(psi1cfull_n)))
dx*dy*sum(sum(phi0full_n.*conj(psi1sfull_n)))
dx*dy*sum(sum(phi0full_n.*conj(psi0full_n)))

save('pulsedata2.mat')



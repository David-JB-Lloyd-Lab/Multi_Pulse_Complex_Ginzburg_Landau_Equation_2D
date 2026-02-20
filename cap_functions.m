function [Phi,G,Psi]=cap_functions(N,grid_points,V_n,w,beta,gamma,delta)

Vsum=zeros(2*grid_points+2,2*grid_points+2);
Vsum_w=zeros(2*grid_points+2,2*grid_points+2);
Psi=zeros(2*grid_points+2,2*grid_points+2,N);

for j=1:N
    Vsum(:,:)=Vsum(:,:)+V_n(:,:,j);
end
Vsum_w(:,:)=Vsum(:,:)+w(:,:);

f_Vsum=my_f(Vsum,beta,gamma,delta);
df_Vsum=my_df(Vsum,w,beta,gamma,delta);

Phi=f_Vsum;
for j=1:N
    Phi(:,:)=Phi(:,:)-my_f(V_n(:,:,j),beta,gamma,delta);
    Psi(:,:,j)=df_Vsum(:,:)-my_df(V_n(:,:,j),w,beta,gamma,delta);
end

G(:,:)=my_f(Vsum_w,beta,gamma,delta)-df_Vsum(:,:)-f_Vsum(:,:);

end
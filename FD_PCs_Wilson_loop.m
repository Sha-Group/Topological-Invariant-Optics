% This MATLAB code caculates the Wilson loops for non-Hermitian systems based on the "right-right", "left-left", "right-left", "left-right" eigenstates.
% REF M. L. N. Chen, etal.¡°Comparative Study of Hermitian and Non-Hermitian
% Topological Dielectric Photonic Crystals,¡±Phys. Rev. A 104, 033501, 2021.  DOI: 10.1103/PhysRevA.104.033501
% Menglin L. N. Chen, the University of Hong Kong, Hong Kong
% email:menglin@connect.hku.hk

function FD_PCs;
clc;
clear;
% close all

flag=1;  % 1, "right-right"; 2, "left-left"; 3, "right-left"; 4, "left-right"

% initial parameter
aa=2.8;  % <3, nontrivial; =3,double Dirac cones; >3,trivial.
Px=2.8*6e-3;  % lattice constant along x (obl)
Py=Px;  % lattice constant along y (obl)
R=Px/aa;
a0=Px;

Nx=100;  % grids along x
Ny=round(Nx*Py/Px);  % grids along y
deltax=Px/Nx;
deltay=Py/Ny;
N=Nx*Ny;

% position of cylinder
xc1=R;yc1=0;
xc2=a0-R;
xc3=0;yc2=R;
xc4=a0-R;
xc5=a0;
xc6=0;yc3=(a0-R);
xc7=R;
xc8=a0;
xc9=R;yc4=a0;
xc10=a0-R;

% generalized coordinate
u1=[1,0,0];
u2=[1/2,sqrt(3)/2,0];
u3=[0,0,1];
g=[u1*u1',u1*u2',u1*u3';u2*u1',u2*u2',u2*u3';u3*u1',u3*u2',u3*u3'];
Area=abs(u1*cross(u2,u3)');

% material parameters
eps_arr_tensor=zeros(3,3);
mu_arr_tensor=zeros(3,3);
epr1=(12+0.05*i)*eye(3); % gain
epr2=(12-0.05*i)*eye(3); % loss
mur1=eye(3);
mur2=eye(3);

eps_arr=ones(1,N);
mu_arr=ones(1,N);

muinv11_arr=zeros(1,N);
muinv12_arr=zeros(1,N);
muinv21_arr=zeros(1,N);
muinv22_arr=zeros(1,N);

for m=1:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        flag1=shape_f((m-1)*deltax+0.5*deltax,(n-1)*deltay+0.5*deltay,xc1,xc4,xc6,xc8,xc9,yc1,yc2,yc3,yc4,g);
        flag2=shape_f((m-1)*deltax+0.5*deltax,(n-1)*deltay-0.5*deltay,xc1,xc4,xc6,xc8,xc9,yc1,yc2,yc3,yc4,g);
        flag3=shape_f((m-1)*deltax-0.5*deltax,(n-1)*deltay+0.5*deltay,xc1,xc4,xc6,xc8,xc9,yc1,yc2,yc3,yc4,g);
        flag4=shape_f((m-1)*deltax-0.5*deltax,(n-1)*deltay-0.5*deltay,xc1,xc4,xc6,xc8,xc9,yc1,yc2,yc3,yc4,g);
        
        flag5=shape_f2((m-1)*deltax+0.5*deltax,(n-1)*deltay+0.5*deltay,xc2,xc3,xc5,xc7,xc10,yc1,yc2,yc3,yc4,g);
        flag6=shape_f2((m-1)*deltax+0.5*deltax,(n-1)*deltay-0.5*deltay,xc2,xc3,xc5,xc7,xc10,yc1,yc2,yc3,yc4,g);
        flag7=shape_f2((m-1)*deltax-0.5*deltax,(n-1)*deltay+0.5*deltay,xc2,xc3,xc5,xc7,xc10,yc1,yc2,yc3,yc4,g);
        flag8=shape_f2((m-1)*deltax-0.5*deltax,(n-1)*deltay-0.5*deltay,xc2,xc3,xc5,xc7,xc10,yc1,yc2,yc3,yc4,g);
        
        
        eps_arr_tensor=(eps_arr(index)*eye(3)+(flag1+flag2+flag3+flag4)/4*(epr1-eye(3))+(flag5+flag6+flag7+flag8)/4*(epr2-eye(3)));
        eps_arr_tensor=(eps_arr_tensor*inv(g)*Area);
        eps_arr(index)=eps_arr_tensor(3,3);
        
        mu_arr_tensor=(mu_arr(index)*eye(3)+(flag1+flag2+flag3+flag4)/4*(mur1-eye(3))+(flag5+flag6+flag7+flag8)/4*(mur2-eye(3)));
        mu_arr_tensor=(mu_arr_tensor*inv(g)*Area);
        mu_arr_tensor=inv(mu_arr_tensor);
        
        muinv11_arr(index)=mu_arr_tensor(1,1);
        muinv12_arr(index)=mu_arr_tensor(1,2);
        muinv21_arr(index)=mu_arr_tensor(2,1);
        muinv22_arr(index)=mu_arr_tensor(2,2);
        
        line_eps(:,index)=[index];
        column_eps(:,index)=[index];
    end
end

mu11_inv=sparse(line_eps,column_eps,muinv11_arr,N,N);
mu12_inv=sparse(line_eps,column_eps,muinv12_arr,N,N);
mu21_inv=sparse(line_eps,column_eps,muinv21_arr,N,N);
mu22_inv=sparse(line_eps,column_eps,muinv22_arr,N,N);

epsV=sparse(line_eps,column_eps,eps_arr,N,N);

% %check the mesh
% for m=1:Nx
%     for n=1:Ny
%         X(m,n)=m*deltax+0.5*n*deltay;Y(m,n)=sqrt(3)/2*n*deltay;
%         index=(n-1)*Nx+m;
%         eps_arrPLOT(m,n)=eps_arr(index);
%     end
% end
%
% figure;
% pcolor(X,Y,imag(eps_arrPLOT));
% shading interp
% colormap jet
% colorbar
% axis equal


%%%%% wilson loop, assuming winding direction along x
N_Lband=1;
N_Hband=3;
N_max=N_Hband;

% winding step
Nkxg=5;
Pkxg=4*pi/sqrt(3)/Px;
kcxg=linspace(0,1,Nkxg)*Pkxg;

% path along y
Nkyg=21;
Pkyg=2*pi/sqrt(3)/Py;
kcyg=linspace(0,1,Nkyg)*Pkyg;

% projection on xoy coordinate
for m=1:Nkyg
    for n=1:Nkxg
        index=(m-1)*(Nkxg)+n;
        kcx(index)=kcxg(n)*sqrt(3)/2;
        kcy(index)=kcyg(m)-kcxg(n)*0.5;
    end
end

res1=zeros(1,Nkyg);
res2=zeros(1,Nkyg);
res3=zeros(1,Nkyg);

for m=1:Nkyg
    Nkyg-m
    kk0=(m-1)*(Nkxg)+1;
    
    Wp=eye(N_Hband-N_Lband+1);
    
    if flag==1 | flag==3
        [omega,field1]=eigs_PCs(Nx,Ny,Px,Py,kcx(kk0),kcy(kk0),N_max,mu11_inv,...
            mu12_inv,mu21_inv,mu22_inv,epsV,0);
    else
        [omega,field1]=eigs_PCs_h(Nx,Ny,Px,Py,kcx(kk0),kcy(kk0),N_max,mu11_inv,...
            mu12_inv,mu21_inv,mu22_inv,epsV,0);
    end
    
    k=[kcx(kk0),kcy(kk0),0];
    for mm=1:Nx
        for nn=1:Ny
            index=(nn-1)*Nx+mm;
            Lxg=(mm-1/2)*deltax;
            Lyg=(nn-1/2)*deltay;
            field1(index,:)=field1(index,:)*exp(j*k*u1'*Lxg)*exp(j*k*u2'*Lyg);
        end
    end
    
    kk=kk0+1;
    for n=2:Nkxg % this is the loop along kx
        
        if flag==1  % <1_R,2_R> <2_R,3_R><3_R,4_R><4_R,5_R>
            [omega,field2]=eigs_PCs(Nx,Ny,Px,Py,kcx(kk),kcy(kk),N_max,mu11_inv,...
                mu12_inv,mu21_inv,mu22_inv,epsV,0);
        elseif flag==2  % <1_L,2_L> <2_L,3_L><3_L,4_L><4_L,5_L>
            [omega,field2]=eigs_PCs_h(Nx,Ny,Px,Py,kcx(kk),kcy(kk),N_max,mu11_inv,...
                mu12_inv,mu21_inv,mu22_inv,epsV,0);
        elseif flag==3   % <1_R,2_L> <2_L,3_R><3_R,4_L><4_L,5_R>,Nkyg must be odd
            if (mod(n,2)==1)
                [omega,field2]=eigs_PCs(Nx,Ny,Px,Py,kcx(kk),kcy(kk),N_max,mu11_inv,...
                    mu12_inv,mu21_inv,mu22_inv,epsV,0);
            else
                [omega,field2]=eigs_PCs_h(Nx,Ny,Px,Py,kcx(kk),kcy(kk),N_max,mu11_inv,...
                    mu12_inv,mu21_inv,mu22_inv,epsV,0);
            end
        elseif flag==4  % <1_L,2_R> <2_R,3_L><3_L,4_R><4_R,5_L>,Nkyg must be odd
            if (mod(n,2)==0)
                [omega,field2]=eigs_PCs(Nx,Ny,Px,Py,kcx(kk),kcy(kk),N_max,mu11_inv,...
                    mu12_inv,mu21_inv,mu22_inv,epsV,0);
            else
                [omega,field2]=eigs_PCs_h(Nx,Ny,Px,Py,kcx(kk),kcy(kk),N_max,mu11_inv,...
                    mu12_inv,mu21_inv,mu22_inv,epsV,0);
            end
        end
        
        k=[kcx(kk),kcy(kk),0];
        for mm=1:Nx
            for nn=1:Ny
                index=(nn-1)*Nx+mm;
                Lxg=(mm-1/2)*deltax;
                Lyg=(nn-1/2)*deltay;
                field2(index,:)=field2(index,:)*exp(j*k*u1'*Lxg)*exp(j*k*u2'*Lyg);
            end
        end
        
        % gauge fixing
        if (n==Nkxg)
            if flag==1 | flag==3
                [omega,field1_add]=eigs_PCs(Nx,Ny,Px,Py,kcx(kk0),kcy(kk0),N_max,mu11_inv,...
                    mu12_inv,mu21_inv,mu22_inv,epsV,0);
                [omega,field2_add]=eigs_PCs(Nx,Ny,Px,Py,kcx(kk),kcy(kk),N_max,mu11_inv,...
                    mu12_inv,mu21_inv,mu22_inv,epsV,0);
            else
                [omega,field1_add]=eigs_PCs_h(Nx,Ny,Px,Py,kcx(kk0),kcy(kk0),N_max,mu11_inv,...
                    mu12_inv,mu21_inv,mu22_inv,epsV,0);
                [omega,field2_add]=eigs_PCs_h(Nx,Ny,Px,Py,kcx(kk),kcy(kk),N_max,mu11_inv,...
                    mu12_inv,mu21_inv,mu22_inv,epsV,0);
            end
            for mm=1:Nx
                for nn=1:Ny
                    index=(nn-1)*Nx+mm;
                    field2(index,:)=field2(index,:).*field1_add(index,:)./field2_add(index,:);
                end
            end
        end
        
        kk=kk+1;
        
        % field normalization
        for mm=1:N_max
            field1(:,mm)=field1(:,mm)/sqrt(field1(:,mm)'.*eps_arr*field1(:,mm));
            field2(:,mm)=field2(:,mm)/sqrt(field2(:,mm)'.*eps_arr*field2(:,mm));
        end
        
        Wloop=( field1(:,N_Lband:N_Hband)'.*repmat(eps_arr,N_Hband-N_Lband+1,1)*field2(:,N_Lband:N_Hband) );
        Wp=Wp*Wloop;
        field1=field2;
    end
    
    pp=angle(eigs(Wp))/pi; % phase of the eigenvalue of the wilson loop matrix
    pp=sort(pp);
    res1(1,m)=pp(1);
    res2(1,m)=pp(2);
    res3(1,m)=pp(3);
end

figure;
hold on
plot(res1(1,:),'ro');
plot(res2(1,:),'gx');
plot(res3(1,:),'bd');
axis([1 Nkyg -1 1])
xlabel('ky')
ylabel('W')
yticks([-1 0 1])


%% solve for the right eigenstates
function  [omega,fieldz]=eigs_PCs(Nx,Ny,Px,Py,kx,ky,N_max,mu11_inv,...
    mu12_inv,mu21_inv,mu22_inv,epsV,NN)

deltax=Px/Nx;
deltay=Py/Ny;

N=Nx*Ny;

% generalized coordinate
u1=[1,0,0];
u2=[1/2,sqrt(3)/2,0];
u3=[0,0,1];
g=[u1*u1',u1*u2',u1*u3';u2*u1',u2*u2',u2*u3';u3*u1',u3*u2',u3*u3'];
Area=abs(u1*cross(u2,u3)');

k=[kx,ky,0];

line_U1=zeros(2,N);
column_U1=zeros(2,N);
value_U1=zeros(2,N);
line_U2=zeros(2,N);
column_U2=zeros(2,N);
value_U2=zeros(2,N);
line_V1=zeros(2,N);
column_V1=zeros(2,N);
value_V1=zeros(2,N);
line_V2=zeros(2,N);
column_V2=zeros(2,N);
value_V2=zeros(2,N);

% center region_U1
for m=1:Nx-1
    for n=1:Ny
        index=(n-1)*Nx+m;
        line_U1(:,index)=[index,index];
        column_U1(:,index)=[index,index+1];
        value_U1(:,index)=[-1,1]/deltax;
    end
end
% right region_U1
for m=Nx:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        line_U1(:,index)=[index,index];
        column_U1(:,index)=[index,(n-1)*Nx+1];
        value_U1(:,index)=[-1,exp(-j*k*u1'*Px)]/deltax;
    end
end
% center region_U2
for m=1:Nx
    for n=1:Ny-1
        index=(n-1)*Nx+m;
        line_U2(:,index)=[index,index];
        column_U2(:,index)=[index,index+Nx];
        value_U2(:,index)=[-1,1]/deltay;
    end
end
% bottom region_U2
for m=1:Nx
    for n=Ny:Ny
        index=(n-1)*Nx+m;
        line_U2(:,index)=[index,index];
        column_U2(:,index)=[index,m];
        value_U2(:,index)=[-1,exp(-j*k*u2'*Py)]/deltay;
    end
end
% center region_V1
for m=2:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        line_V1(:,index)=[index,index];
        column_V1(:,index)=[index,index-1];
        value_V1(:,index)=[1,-1]/deltax;
    end
end
% left region_V1
for m=1:1
    for n=1:Ny
        index=(n-1)*Nx+m;
        line_V1(:,index)=[index,index];
        column_V1(:,index)=[index,(n-1)*Nx+Nx];
        value_V1(:,index)=[1,-exp(j*k*u1'*Px)]/deltax;
    end
end
% center region_V2
for m=1:Nx
    for n=2:Ny
        index=(n-1)*Nx+m;
        line_V2(:,index)=[index,index];
        column_V2(:,index)=[index,index-Nx];
        value_V2(:,index)=[1,-1]/deltay;
    end
end
% top region_V2
for m=1:Nx
    for n=1:1
        index=(n-1)*Nx+m;
        line_V2(:,index)=[index,index];
        column_V2(:,index)=[index,(Ny-1)*Nx+m];
        value_V2(:,index)=[1,-exp(j*k*u2'*Py)]/deltay;
    end
end

U1=sparse(line_U1,column_U1,value_U1,N,N);
U2=sparse(line_U2,column_U2,value_U2,N,N);
V1=sparse(line_V1,column_V1,value_V1,N,N);
V2=sparse(line_V2,column_V2,value_V2,N,N);

% TM polarization
A=(U1*(mu21_inv*V2-mu22_inv*V1)-U2*(mu11_inv*V2-mu12_inv*V1));
[V,D] = eigs(A,epsV,N_max,'sm');       % find eigenmodes and eigen-frequency

[omega,order]=sort(sqrt(abs(real(diag(D)))));  %  from small to large eigenvalues

fieldz=zeros(N,N_max);
for m=1:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        fieldz(index,:)=V(index,:);
    end
end

%% solve for the left eigenstates
function  [omega,fieldz]=eigs_PCs_h(Nx,Ny,Px,Py,kx,ky,N_max,mu11_inv,...
    mu12_inv,mu21_inv,mu22_inv,epsV,NN)

deltax=Px/Nx;
deltay=Py/Ny;

N=Nx*Ny;      % total grid numbers;

% generalized coordinate
u1=[1,0,0];
u2=[1/2,sqrt(3)/2,0];
u3=[0,0,1];
g=[u1*u1',u1*u2',u1*u3';u2*u1',u2*u2',u2*u3';u3*u1',u3*u2',u3*u3'];
Area=abs(u1*cross(u2,u3)');

k=[kx,ky,0];

line_U1=zeros(2,N);
column_U1=zeros(2,N);
value_U1=zeros(2,N);
line_U2=zeros(2,N);
column_U2=zeros(2,N);
value_U2=zeros(2,N);
line_V1=zeros(2,N);
column_V1=zeros(2,N);
value_V1=zeros(2,N);
line_V2=zeros(2,N);
column_V2=zeros(2,N);
value_V2=zeros(2,N);

% center region_U1
for m=1:Nx-1
    for n=1:Ny
        index=(n-1)*Nx+m;
        line_U1(:,index)=[index,index];
        column_U1(:,index)=[index,index+1];
        value_U1(:,index)=[-1,1]/deltax;
    end
end
% right region_U1
for m=Nx:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        line_U1(:,index)=[index,index];
        column_U1(:,index)=[index,(n-1)*Nx+1];
        value_U1(:,index)=[-1,exp(-j*k*u1'*Px)]/deltax;
    end
end
% center region_U2
for m=1:Nx
    for n=1:Ny-1
        index=(n-1)*Nx+m;
        line_U2(:,index)=[index,index];
        column_U2(:,index)=[index,index+Nx];
        value_U2(:,index)=[-1,1]/deltay;
    end
end
% bottom region_U2
for m=1:Nx
    for n=Ny:Ny
        index=(n-1)*Nx+m;
        line_U2(:,index)=[index,index];
        column_U2(:,index)=[index,m];
        value_U2(:,index)=[-1,exp(-j*k*u2'*Py)]/deltay;
    end
end
% center region_V1
for m=2:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        line_V1(:,index)=[index,index];
        column_V1(:,index)=[index,index-1];
        value_V1(:,index)=[1,-1]/deltax;
    end
end
% left region_V1
for m=1:1
    for n=1:Ny
        index=(n-1)*Nx+m;
        line_V1(:,index)=[index,index];
        column_V1(:,index)=[index,(n-1)*Nx+Nx];
        value_V1(:,index)=[1,-exp(j*k*u1'*Px)]/deltax;
    end
end
% center region_V2
for m=1:Nx
    for n=2:Ny
        index=(n-1)*Nx+m;
        line_V2(:,index)=[index,index];
        column_V2(:,index)=[index,index-Nx];
        value_V2(:,index)=[1,-1]/deltay;
    end
end
% top region_V2
for m=1:Nx
    for n=1:1
        index=(n-1)*Nx+m;
        line_V2(:,index)=[index,index];
        column_V2(:,index)=[index,(Ny-1)*Nx+m];
        value_V2(:,index)=[1,-exp(j*k*u2'*Py)]/deltay;
    end
end

U1=sparse(line_U1,column_U1,value_U1,N,N);
U2=sparse(line_U2,column_U2,value_U2,N,N);
V1=sparse(line_V1,column_V1,value_V1,N,N);
V2=sparse(line_V2,column_V2,value_V2,N,N);

% TM polarization
A=(U1*(mu21_inv*V2-mu22_inv*V1)-U2*(mu11_inv*V2-mu12_inv*V1));
[V,D] = eigs(A',epsV',N_max,'sm');       % find eigenmodes and eigen-frequency

[omega,order]=sort(sqrt(diag(D)));  %  from small to large eigenvalues

fieldz=zeros(N,N_max);
for m=1:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        fieldz(index,:)=V(index,:);
    end
end

%  model the unit cell
function flag=shape_f(x,y,xc1,xc4,xc6,xc8,xc9,yc1,yc2,yc3,yc4,g)
flag=0;
radius1=2.0e-3;
radius2=2.0e-3;
distance1=[x-xc1,y-yc1,0]*g*[x-xc1;y-yc1;0];
distance4=[x-xc4,y-yc2,0]*g*[x-xc4;y-yc2;0];
distance6=[x-xc6,y-yc3,0]*g*[x-xc6;y-yc3;0];
distance8=[x-xc8,y-yc3,0]*g*[x-xc8;y-yc3;0];
distance9=[x-xc9,y-yc4,0]*g*[x-xc9;y-yc4;0];
if (distance1<radius1^2 | (distance4)<radius2^2|...
        (distance6)<radius2^2| (distance8)<radius2^2 | ...
        (distance9)<radius2^2)
    flag=1;
end
function flag=shape_f2(x,y,xc2,xc3,xc5,xc7,xc10,yc1,yc2,yc3,yc4,g)
flag=0;
radius1=2.0e-3;
radius2=2.0e-3;
distance2=[x-xc2,y-yc1,0]*g*[x-xc2;y-yc1;0];
distance3=[x-xc3,y-yc2,0]*g*[x-xc3;y-yc2;0];
distance5=[x-xc5,y-yc2,0]*g*[x-xc5;y-yc2;0];
distance7=[x-xc7,y-yc3,0]*g*[x-xc7;y-yc3;0];
distance10=[x-xc10,y-yc4,0]*g*[x-xc10;y-yc4;0];

if ((distance2)<radius2^2 | (distance3)<radius2^2| ...
        (distance5)<radius2^2| (distance7)<radius2^2 | ...
        (distance10)<radius2^2)
    flag=1;
end



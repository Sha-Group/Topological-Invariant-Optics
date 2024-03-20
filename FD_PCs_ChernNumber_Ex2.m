% compute band structure of 2D PCs with TM polarization
% dielectric cylinders in air
% v1.0 Wei Sha@ZJU   Email: weisha@zju.edu.cn
% v2.0 Ran Zhao@AHU  Email: rzhao@ahu.edu.cn
% reference: First-principle calculation of Chern number in gyrotropic photonic crystals
% https://doi.org/10.1364/OE.380077

function FD_PCs;
clc;
clear;
Px=1;     % lattice constant along x
Py=1;     % lattice constant along y
N_sam=10; % Brillouin sampling points ***
N_max=3;  % maximum band number ***

Polarflag=0; %极化方式，0，TM极化；1，TE极化

shifted=0.0;
Nkx=4;Nky=4;;
Pkx=(2*pi+shifted)/Px;
Pky=(2*pi+shifted)/Py;
deltakx=Pkx/Nkx;
deltaky=Pky/Nky;
NK=Nkx*Nky;
kcx=zeros(NK,4);
kcy=zeros(NK,4);
% eps_arr=ones(1,N)
% deltax=0.01;%0.025; % space step along x ****
% deltay=0.01;%0.025;% space step along y ****
% omega=1,f=1;
% [omega,f]=chern_PCs(Px,Py,kx(1),ky(1),N_max,0);
fieldTot=zeros(N_max,1);
ChernNumber=zeros(N_max,1);
fieldtemp=zeros(N_max,1);
N_Lband=1;
N_Hband=2;
fieldTotComposite4=0;
ChernNumberComposite=0;
fieldtempComposite=0;
fieldTotComposite12=0;fieldTotComposite13=0;fieldTotComposite14=0;
%定义画图，离散布里渊区域所需的变量；
%==========meep colorbar========
% mincolor    = [0 0 1]; % red
% mediancolor = [1 1 1]; % white
% maxcolor    = [1 0 0]; % blue
% color1=zeros(NK,3);
% ColorMapSize = 64;
% int1 = zeros(ColorMapSize,3);
% int2 = zeros(ColorMapSize,3);
% for k=1:3
%     int1(:,k) = linspace(mincolor(k), mediancolor(k), ColorMapSize);
%     int2(:,k) = linspace(mediancolor(k), maxcolor(k), ColorMapSize);
% end
% meep = [int1(1:end-1,:); int2];
% p=patch(kcx,kcy,color1);
% set(p,'EdgeColor',[1,1,1]/2,'LineWidth',1)
% colormap(meep);
% axis equal
% colorbar('hide');


tic;
for m=1:Nkx
    for n=1:Nky

        index=(n-1)*Nkx+m;                 
%         kcx(index,:)=[(m-1)*deltakx-0.5*Pkx,(m)*deltakx-0.5*Pkx,(m)*deltakx-0.5*Pkx,(m-1)*deltakx-0.5*Pkx];
%         kcy(index,:)=[(n-1)*deltaky-0.5*Pky,(n-1)*deltaky-0.5*Pky,(n)*deltaky-0.5*Pky,(n)*deltaky-0.5*Pky];   
        kcx(index,:)=[(m-1)*deltakx,(m)*deltakx,(m)*deltakx,(m-1)*deltakx];
        kcy(index,:)=[(n-1)*deltaky,(n-1)*deltaky,(n)*deltaky,(n)*deltaky];         
        
       [omega,field1,eps_arr]=eigs_PCs(Px,Py,kcx(index,1),kcy(index,1),N_max,0,Polarflag);
       [omega,field2,eps_arr]=eigs_PCs(Px,Py,kcx(index,2),kcy(index,2),N_max,0,Polarflag);       
       [omega,field3,eps_arr]=eigs_PCs(Px,Py,kcx(index,3),kcy(index,3),N_max,0,Polarflag);       
       [omega,field4,eps_arr]=eigs_PCs(Px,Py,kcx(index,4),kcy(index,4),N_max,0,Polarflag); 
       for mm=1:N_max
                 fieldtemp(mm)=field1(:,mm)'.*eps_arr*field2(:,mm)/abs( field1(:,mm)'.*eps_arr*field2(:,mm)  )*...
                 field2(:,mm)'.*eps_arr*field3(:,mm)/abs( field2(:,mm)'.*eps_arr*field3(:,mm)  )*...
                 field3(:,mm)'.*eps_arr*field4(:,mm)/abs( field3(:,mm)'.*eps_arr*field4(:,mm)  )*...
                 field4(:,mm)'.*eps_arr*field1(:,mm)/abs( field4(:,mm)'.*eps_arr*field1(:,mm)  );
           fieldTot(mm)=fieldTot(mm)+imag(log(fieldtemp(mm))); 
       end
       color1(index,1)=imag(log(fieldtemp(2)));
   fieldt1=( field1(:,N_Lband:N_Hband)'.*repmat(eps_arr,N_Hband-N_Lband+1,1)*field2(:,N_Lband:N_Hband) )*...
           ( field2(:,N_Lband:N_Hband)'.*repmat(eps_arr,N_Hband-N_Lband+1,1)*field3(:,N_Lband:N_Hband) )*...
           ( field3(:,N_Lband:N_Hband)'.*repmat(eps_arr,N_Hband-N_Lband+1,1)*field4(:,N_Lband:N_Hband) )*...
           ( field4(:,N_Lband:N_Hband)'.*repmat(eps_arr,N_Hband-N_Lband+1,1)*field1(:,N_Lband:N_Hband) );
       fieldtempComposite=det(fieldt1)/abs(det(fieldt1));
      fieldTotComposite4=fieldTotComposite4+imag(log(fieldtempComposite)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fieldt1=( field1(:,1:2)'.*repmat(eps_arr,2-1+1,1)*field2(:,1:2) )*...
           ( field2(:,1:2)'.*repmat(eps_arr,2-1+1,1)*field3(:,1:2) )*...
           ( field3(:,1:2)'.*repmat(eps_arr,2-1+1,1)*field4(:,1:2) )*...
           ( field4(:,1:2)'.*repmat(eps_arr,2-1+1,1)*field1(:,1:2) );
       fieldtempComposite=det(fieldt1)/abs(det(fieldt1));
      fieldTotComposite12=fieldTotComposite12+imag(log(fieldtempComposite)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fieldt1=( field1(:,1:3)'.*repmat(eps_arr,3-1+1,1)*field2(:,1:3) )*...
           ( field2(:,1:3)'.*repmat(eps_arr,3-1+1,1)*field3(:,1:3) )*...
           ( field3(:,1:3)'.*repmat(eps_arr,3-1+1,1)*field4(:,1:3) )*...
           ( field4(:,1:3)'.*repmat(eps_arr,3-1+1,1)*field1(:,1:3) );
       fieldtempComposite=det(fieldt1)/abs(det(fieldt1));
      fieldTotComposite13=fieldTotComposite13+imag(log(fieldtempComposite)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fieldt1=( field1(:,1:10)'.*repmat(eps_arr,10-1+1,1)*field2(:,1:10) )*...
%            ( field2(:,1:10)'.*repmat(eps_arr,10-1+1,1)*field3(:,1:10) )*...
%            ( field3(:,1:10)'.*repmat(eps_arr,10-1+1,1)*field4(:,1:10) )*...
%            ( field4(:,1:10)'.*repmat(eps_arr,10-1+1,1)*field1(:,1:10) );
%        fieldtempComposite=det(fieldt1)/abs(det(fieldt1));
%       fieldTotComposite14=fieldTotComposite14+imag(log(fieldtempComposite)); 
      
    end
end

ChernNumber(:)=fieldTot(:)/2/pi
ChernNumberComposite4=fieldTotComposite4/2/pi
ChernNumberComposite12=fieldTotComposite12/2/pi
ChernNumberComposite13=fieldTotComposite13/2/pi
% ChernNumberComposite14=fieldTotComposite14/2/pi
toc;






% sampling at the edge of Brillouin zone
kx1=linspace(0,pi/Px,N_sam); ky1=linspace(0,0,N_sam);
kx2=linspace(pi/Px,pi/Px,N_sam); ky2=linspace(0,pi/Py,N_sam);
radius=sqrt((pi/Px)^2+(pi/Py)^2);
radius_sam=linspace(radius,0,round(N_sam*sqrt(2)));
kx3=radius_sam*cos(pi/4); ky3=radius_sam*sin(pi/4);

kx=[kx1(1:end),kx2(2:end),kx3(2:end-1),kx1(1)];
ky=[ky1(1:end),ky2(2:end),ky3(2:end-1),ky1(1)];
omega_f=zeros(N_max,length(kx));

% key points at Brillouin zone for drawing eigenmodes
Gamma=1;
Tao=length(kx1);
M=length(kx1)+length(kx2)-1;


% calculate eigenvalue
for m=1:length(kx)
    [omega_f(:,m),field1,eps_arr]=eigs_PCs(Px,Py,kx(m),ky(m),N_max,0,Polarflag);
end

% draw band diagram
figure(1)
for m=1:N_max
    plot(omega_f(m,:)*Px/2/pi,'rx-')
    hold on
    axis([1,length(kx),min(min(omega_f))*Px/2/pi,max(max(omega_f))*Px/2/pi]);
    xlabel('Bloch Wavenumber')
    ylabel('Normalized Eigen-frequency')
end

% draw eigenmode (first four eigenmodes at M point)
[omega_r,field1,eps_arr]=eigs_PCs(Px,Py,kx(M),ky(M),N_max,1,Polarflag);

% solve the eigen-equation with a fixed kx and ky
% N_max is the maximum eigenvalues of interest
% Px Py are lattice constants along x and y directions
% flag=1 for drawing eigenmodes
%((eps_arr.').*V(:,1))'*V(:,4)

% function [omega]=eigs_PCs(Px,Py,kx,ky,N_max,flag)
 function [omega,fieldz,eps_arr]=eigs_PCs(Px,Py,kx,ky,N_max,flag,Polarflag)
deltax=0.02;%0.025; % space step along x ****
deltay=0.02;%0.025;% space step along y ****
% epr1=15;       % dielectric ****
% mur1=14;
epr1=[15,0,0;0,15,0;0,0,15];       % dielectric ****
% mur1=[1,0,0;0,1,0;0,0,1];
% epr1=[14,-12.4*j,0;12.4*j,14,0;0,0,15];      % dielectric ****
mur1=[14,12.4*j,0;-12.4*j,14,0;0,0,15];

kappa1=12.4;
Nx=Px/deltax; % grids along x
Ny=Py/deltay; % grids along y
N=Nx*Ny;      % total grid numbers;
x_center=(Nx)/2*deltax; %  x center
y_center=(Ny)/2*deltay; %  y center
line_eps=zeros(1,N);
column_eps=zeros(1,N);

eps_arr=ones(1,N);  % permittivity of air (initial condition)
epsinv11_arr=zeros(1,N);
epsinv12_arr=zeros(1,N);
epsinv21_arr=zeros(1,N);
epsinv22_arr=zeros(1,N);
epsinv33_arr=zeros(1,N);

mu_arr=ones(1,N);  %   of air (initial condition)
muinv11_arr=zeros(1,N);
muinv12_arr=zeros(1,N);
muinv21_arr=zeros(1,N);
muinv22_arr=zeros(1,N);
muinv33_arr=zeros(1,N);

kappa_arr=zeros(1,N);
eps_arr_tensor=zeros(3,3);
mu_arr_tensor=zeros(3,3);

% averaged scheme for permittivity
for m=1:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        flag1=shape_f((m-1)*deltax,(n-1)*deltay,x_center,y_center);
        flag2=shape_f((m-1)*deltax,(n)*deltay,x_center,y_center);
        flag3=shape_f((m)*deltax,(n-1)*deltay,x_center,y_center);
        flag4=shape_f((m)*deltax,(n)*deltay,x_center,y_center);
        
        flag1u=shape_f((m-1)*deltax,(n-1)*deltay,x_center,y_center);
        flag2u=shape_f((m-1)*deltax,(n)*deltay,x_center,y_center);
        flag3u=shape_f((m)*deltax,(n-1)*deltay,x_center,y_center);
        flag4u=shape_f((m-1)*deltax,(n-1)*deltay,x_center,y_center);      
%         eps_arr(index)=eps_arr(index)+(flag1+flag2+flag3+flag4)/4*(epr1(3,3)-1);
%         eps_arr_tensor=inv(  eps_arr(index)*eye(3)+(flag1+flag2+flag3+flag4)/4*(epr1-eye(3))      );
       eps_arr_tensor=(  eps_arr(index)*eye(3)+(flag1+flag2+flag3+flag4)/4*(epr1-eye(3))      ); 
       eps_arr(index)=eps_arr_tensor(3,3);       
       eps_arr_tensor=inv(eps_arr_tensor);       
       
        epsinv11_arr(index)=eps_arr_tensor(1,1);       
        epsinv12_arr(index)=eps_arr_tensor(1,2);   
        epsinv21_arr(index)=eps_arr_tensor(2,1);       
        epsinv22_arr(index)=eps_arr_tensor(2,2);  
        epsinv33_arr(index)=eps_arr_tensor(3,3);        
        
        
%         mu_arr(index)=mu_arr(index)+(flag1u+flag2u+flag3u+flag4u)/4*(mur1(3,3)-1);  
%         mu_arr_tensor=inv(  mu_arr(index)*eye(3)+(flag1u+flag2u+flag3u+flag4u)/4*(mur1-eye(3))      );        
        mu_arr_tensor=(  mu_arr(index)*eye(3)+(flag1+flag2+flag3+flag4)/4*(mur1-eye(3))     );           
        mu_arr(index)=mu_arr_tensor(3,3); 
        mu_arr_tensor=inv(mu_arr_tensor);
        
        muinv11_arr(index)=mu_arr_tensor(1,1);       
        muinv12_arr(index)=mu_arr_tensor(1,2);   
        muinv21_arr(index)=mu_arr_tensor(2,1);       
        muinv22_arr(index)=mu_arr_tensor(2,2);  
        muinv33_arr(index)=mu_arr_tensor(3,3);          
        
        kappa_arr(index)=(flag1u+flag2u+flag3u+flag4u)/4*kappa1;
        
        line_eps(:,index)=[index];
        column_eps(:,index)=[index];        
    end
end

% sparse matrix structure (line number, column number, and values)
%  PREFORMATTED
%  TEXT
% 

% value=zeros(5,N);
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
        value_U1(:,index)=[-1,exp(-j*kx*Px)]/deltax;
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
        value_U2(:,index)=[-1,exp(-j*ky*Py)]/deltay;
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
        value_V1(:,index)=[1,-exp(j*kx*Px)]/deltax;
    end
end


%center region_V2
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
        value_V2(:,index)=[1,-exp(j*ky*Py)]/deltay;
    end
end




U1=sparse(line_U1,column_U1,value_U1,N,N);
U2=sparse(line_U2,column_U2,value_U2,N,N);
V1=sparse(line_V1,column_V1,value_V1,N,N);
V2=sparse(line_V2,column_V2,value_V2,N,N);
epsV=sparse(line_eps,column_eps,eps_arr,N,N);
muV=sparse(line_eps,column_eps,mu_arr,N,N);
% epsV=sparse(line_eps,column_eps,eps_arr,N,N);
% mu11_inv=sparse(line_eps,column_eps,mu_arr./(mu_arr.^2-kappa_arr.^2),N,N);
% mu12_inv=sparse(line_eps,column_eps,-j*kappa_arr./(mu_arr.^2-kappa_arr.^2),N,N);


eps11_inv=sparse(line_eps,column_eps,epsinv11_arr,N,N);
eps12_inv=sparse(line_eps,column_eps,epsinv12_arr,N,N);
eps21_inv=sparse(line_eps,column_eps,epsinv21_arr,N,N);
eps22_inv=sparse(line_eps,column_eps,epsinv22_arr,N,N);

mu11_inv=sparse(line_eps,column_eps,muinv11_arr,N,N);
mu12_inv=sparse(line_eps,column_eps,muinv12_arr,N,N);
mu21_inv=sparse(line_eps,column_eps,muinv21_arr,N,N);
mu22_inv=sparse(line_eps,column_eps,muinv22_arr,N,N);
mu33_inv=sparse(line_eps,column_eps,muinv33_arr,N,N);

if(Polarflag==0)
    %Solve the eigen-mode of the TM mode wave.
    A=(U1*(mu21_inv*V2-mu22_inv*V1)-U2*(mu11_inv*V2-mu12_inv*V1));
    [V,D] = eigs(A,epsV,N_max,'sm');       % find eigenmodes and eigen-frequency
else
    %Solve the eigen-mode of the TE mode wave.
    A=-(V1*(eps21_inv*U2-eps22_inv*U1)-V2*(eps11_inv*U2-eps12_inv*U1));
    [V,D] = eigs(A,muV,N_max,'sm');       % find eigenmodes and eigen-frequency
    eps_arr=mu_arr;
end

% diag(D);


[omega,order]=sort(sqrt(abs(real(diag(D)))));  %  from small to large eigenvalues



%测试不同特征值的正交性%%%%%%%
% V(:,:)'.*repmat(eps_arr,N_max,1)*V(:,:)



fieldz=zeros(N,N_max);

for m=1:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        
        Lx=(m-1/2)*deltax;
        Ly=(n-1/2)*deltay;
        fieldz(index,:)=V(index,:)*exp(j*kx*Lx)*exp(j*ky*Ly);
%           fieldz(index,:)=V(index,:);           
    end
end


% draw eigenmodes
field_2d=zeros(Nx,Ny);
if (flag==1) % flag control
    figure(2)
    for mm=1:N_max
        field_1d=V(:,order(mm));  % eigenmodes from low to high order
        for m=1:Nx
            for n=1:Ny
                index=(n-1)*Nx+m;
                field_2d(m,n)=field_1d(index);
            end
        end
        subplot(2,2,mm)
        pcolor(real(field_2d));
        shading interp
        colorbar
        text_um=num2str(mm);
        title(text_um)
        axis equal
        axis([1 Nx 1 Ny])
    end
end



%  model the unit cell
function flag=shape_f(x,y,xc,yc)
flag=0;
radius=0.11;  % radius of circle in unit cells ***
if (((x-xc)^2+(y-yc)^2)<radius^2)
    flag=1;
end

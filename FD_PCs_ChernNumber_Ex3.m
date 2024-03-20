% compute band structure of 2D PCs with TM/TE polarization
% set polarflag=0 TM polarization; set polarflag=1 TE polarization;
% dielectric cylinders in air
% v1.0 Wei Sha@ZJU   Email: weisha@zju.edu.cn
% v2.0 Ran Zhao@AHU  Email: rzhao@ahu.edu.cn
% reference: First-principle calculation of Chern number in gyrotropic photonic crystals
% https://doi.org/10.1364/OE.380077

function FD_PCs; 
Pxg=1;     % lattice constant along xg
Pyg=1;     % lattice constant along yg
N_sam=10;  % Brillouin sampling points ***
N_max=2;   % maximum band number ***

polarflag=0; %极化方式，0，TM极化；1，TE极化

shifted=0.0;
Nkxg=4;Nkyg=4;
Pkxg=(4*pi)/Pxg/sqrt(3);
Pkyg=(4*pi)/Pyg/sqrt(3);
deltakxg=Pkxg/Nkxg;
deltakyg=Pkyg/Nkyg;
NK=Nkxg*Nkyg;
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
fieldTotComposite=0;
ChernNumberComposite=0;
fieldtempComposite=0;



tic;
for m=1:Nkxg
    for n=1:Nkyg

        index=(n-1)*Nkxg+m;                    
        
        kcxg=[(m-1)*deltakxg,(m)*deltakxg,(m)*deltakxg,(m-1)*deltakxg];
        kcyg=[(n-1)*deltakyg,(n-1)*deltakyg,(n)*deltakyg,(n)*deltakyg];   
        kcx(index,:)=kcxg*sqrt(3)/2;
        kcy(index,:)=kcyg-kcxg*0.5;        
%         kcx(index,:)=[(m-1)*deltakx,(m)*deltakx,(m)*deltakx,(m-1)*deltakx];
%         kcy(index,:)=[(n-1)*deltaky,(n-1)*deltaky,(n)*deltaky,(n)*deltaky];         
        
       [omega,field1,eps_arr]=eigs_PCs(Pxg,Pyg,kcx(index,1),kcy(index,1),N_max,0,polarflag);
       [omega,field2,eps_arr]=eigs_PCs(Pxg,Pyg,kcx(index,2),kcy(index,2),N_max,0,polarflag);       
       [omega,field3,eps_arr]=eigs_PCs(Pxg,Pyg,kcx(index,3),kcy(index,3),N_max,0,polarflag);       
       [omega,field4,eps_arr]=eigs_PCs(Pxg,Pyg,kcx(index,4),kcy(index,4),N_max,0,polarflag); 
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

      fieldTotComposite=fieldTotComposite+imag(log(fieldtempComposite)); 
  
    end
end

ChernNumber(:)=fieldTot(:)/2/pi
ChernNumberComposite=fieldTotComposite/2/pi
toc;





% sampling at the edge of Brillouin zone
kx1=linspace(0,0,round(N_sam*sqrt(3))); ky1=linspace(0,2*pi/Pyg/sqrt(3),round(N_sam*sqrt(3)));
kx2=linspace(0,2*pi/Pxg/3,N_sam); ky2=linspace(2*pi/Pyg/sqrt(3),2*pi/Pyg/sqrt(3),N_sam);

radius=sqrt((2*pi/Pxg/3)^2+(2*pi/Pyg/sqrt(3))^2);
radius_sam=linspace(radius,0,round(N_sam*2));
kx3=radius_sam*cos(pi/3); ky3=radius_sam*sin(pi/3);

kx=[kx1(1:end),kx2(2:end),kx3(2:end-1),kx1(1)];
ky=[ky1(1:end),ky2(2:end),ky3(2:end-1),ky1(1)];
omega_f=zeros(N_max,length(kx));

% key points at Brillouin zone for drawing eigenmodes
Gamma=1;
Tao=length(kx1);
M=length(kx1)+length(kx2)-1;

% calculate eigenvalue
for m=1:length(kx)
    [omega_f(:,m),field1,eps_arr]=eigs_PCs(Pxg,Pyg,kx(m),ky(m),N_max,0,polarflag);
end

% draw band diagram
figure(1)
for m=1:N_max
    plot(omega_f(m,:)*Pxg/2/pi,'rx-')
    hold on
    axis([1,length(kx),min(min(omega_f))*Pxg/2/pi,max(max(omega_f))*Pxg/2/pi]);
    xlabel('Bloch Wavenumber')
    ylabel('Normalized Eigen-frequency')
end

% draw eigenmode (first four eigenmodes at M point)
[omega_r,field1,eps_arr]=eigs_PCs(Pxg,Pyg,kx(M),ky(M),N_max,1,polarflag);

% solve the eigen-equation with a fixed kx and ky
% N_max is the maximum eigenvalues of interest
% Px Py are lattice constants along x and y directions
% flag=1 for drawing eigenmodes
%((eps_arr.').*V(:,1))'*V(:,4)
function  [omega,fieldz,eps_arr]=eigs_PCs(Pxg,Pyg,kx,ky,N_max,flag,polarflag)

deltaxg=0.02*Pxg;%0.025; % space step along u1 ****
deltayg=0.02*Pyg;%0.025;% space step along u2 ****

epr1=[15,0,0;0,15,0;0,0,15];       % dielectric ****
%mur1=[1,0,0;0,1,0;0,0,1]; 
%epr1=[14,12.4*j,0;-12.4*j,14,0;0,0,15];
% mur1=[13.82,10.931*j,0;-10.931*j,13.82,0;0,0,13.82];
 mur1=[0.8736,-0.6671*j,0;0.6671*j,0.8736,0;0,0,0.8736];
kappa1=-0.6671;
Nxg=Pxg/deltaxg; % grids along u1 ****
Nyg=Pyg/deltayg; % grids along u2 ****
N=Nxg*Nyg;      % total grid numbers;
x_center=(Nxg)/2*deltaxg; %  x center
y_center=(Nyg)/2*deltayg; %  y center
line_eps=zeros(1,N);
column_eps=zeros(1,N);

eps_arr=ones(1,N);  % permittivity of air (initial condition)
mu_arr=ones(1,N);  %   of air (initial condition)
muinv11_arr=zeros(1,N);
muinv12_arr=zeros(1,N);
muinv21_arr=zeros(1,N);
muinv22_arr=zeros(1,N);
kappa_arr=zeros(1,N);
% averaged scheme for permittivity

%广义坐标的三个分量；
u2=[1/2,sqrt(3)/2,0];
u1=[1,0,0];
u3=[0,0,1];
g=[u1*u1',u1*u2',u1*u3';u2*u1',u2*u2',u2*u3';u3*u1',u3*u2',u3*u3'];
% g=inv(g);
Area=abs(u1*cross(u2,u3)');

k=[kx,ky,0];
eps_arr_tensor=zeros(3,3);
mu_arr_tensor=zeros(3,3);
for m=1:Nxg
    for n=1:Nyg
        index=(n-1)*Nxg+m;
        flag1=shape_f((m-1)*deltaxg,(n-1)*deltayg,x_center,y_center,g);
        flag2=shape_f((m-1)*deltaxg,(n)*deltayg,x_center,y_center,g);
        flag3=shape_f((m)*deltaxg,(n-1)*deltayg,x_center,y_center,g);
        flag4=shape_f((m)*deltaxg,(n)*deltayg,x_center,y_center,g);
%         eps_arr(index)=eps_arr(index)+(flag1+flag2+flag3+flag4)/4*(epr1(3,3)-1);
%         eps_arr_tensor=inv(eps_arr(index)*inv(g)*Area);   
       eps_arr_tensor=(  eps_arr(index)*eye(3)+(flag1+flag2+flag3+flag4)/4*(epr1-eye(3))      );  

       eps_arr_tensor=(eps_arr_tensor*inv(g)*Area);  
       eps_arr(index)=eps_arr_tensor(3,3);       
       eps_arr_tensor=inv(eps_arr_tensor);
       
        epsinv11_arr(index)=eps_arr_tensor(1,1);       
        epsinv12_arr(index)=eps_arr_tensor(1,2);   
        epsinv21_arr(index)=eps_arr_tensor(2,1);       
        epsinv22_arr(index)=eps_arr_tensor(2,2);  
        epsinv33_arr(index)=eps_arr_tensor(3,3);          
        
%         eps_arr(index)=eps_arr_tensor(3,3);        
%         mu_arr(index)=mu_arr(index)+(flag1+flag2+flag3+flag4)/4*(mur1(3,3)-1);    
        mu_arr_tensor=(  mu_arr(index)*eye(3)+(flag1+flag2+flag3+flag4)/4*(mur1-eye(3))     );           
        mu_arr_tensor=(mu_arr_tensor*inv(g)*Area);
        mu_arr(index)=mu_arr_tensor(3,3); 
        mu_arr_tensor=inv(mu_arr_tensor);
        
        muinv11_arr(index)=mu_arr_tensor(1,1);       
        muinv12_arr(index)=mu_arr_tensor(1,2);   
        muinv21_arr(index)=mu_arr_tensor(2,1);       
        muinv22_arr(index)=mu_arr_tensor(2,2);  
        muinv33_arr(index)=mu_arr_tensor(3,3);   
        
        kappa_arr(index)=(flag1+flag2+flag3+flag4)/4*kappa1;
        
        line_eps(:,index)=[index];
        column_eps(:,index)=[index];        
    end
end
% eps_arr=eps_arr.*Area;



% sparse matrix structure (line number, column number, and values)
% five-point difference
%%
% 
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
for m=1:Nxg-1
    for n=1:Nyg
        index=(n-1)*Nxg+m;
        line_U1(:,index)=[index,index];
        column_U1(:,index)=[index,index+1];
        value_U1(:,index)=[-1,1]/deltaxg;
    end
end
% right region_U1
for m=Nxg:Nxg
    for n=1:Nyg
        index=(n-1)*Nxg+m;
        line_U1(:,index)=[index,index];
        column_U1(:,index)=[index,(n-1)*Nxg+1];
%         value_U1(:,index)=[-1,exp(-j*kx*Px)]/deltax;
        value_U1(:,index)=[-1,exp(-j*k*u1'*Pxg)]/deltaxg;        
        
    end
end

% center region_U2
for m=1:Nxg
    for n=1:Nyg-1
        index=(n-1)*Nxg+m;
        line_U2(:,index)=[index,index];
        column_U2(:,index)=[index,index+Nxg];
        value_U2(:,index)=[-1,1]/deltayg;
    end
end

% bottom region_U2
for m=1:Nxg
    for n=Nyg:Nyg
        index=(n-1)*Nxg+m;
        line_U2(:,index)=[index,index];
        column_U2(:,index)=[index,m];
        value_U2(:,index)=[-1,exp(-j*k*u2'*Pyg)]/deltayg;
    end
end




% center region_V1
for m=2:Nxg
    for n=1:Nyg
        index=(n-1)*Nxg+m;
        line_V1(:,index)=[index,index];
        column_V1(:,index)=[index,index-1];
        value_V1(:,index)=[1,-1]/deltaxg;
    end
end
% left region_V1
for m=1:1
    for n=1:Nyg
        index=(n-1)*Nxg+m;
        line_V1(:,index)=[index,index];
        column_V1(:,index)=[index,(n-1)*Nxg+Nxg];
        value_V1(:,index)=[1,-exp(j*k*u1'*Pxg)]/deltaxg;
    end
end


%center region_V2
for m=1:Nxg
    for n=2:Nyg
        index=(n-1)*Nxg+m;
        line_V2(:,index)=[index,index];
        column_V2(:,index)=[index,index-Nxg];
        value_V2(:,index)=[1,-1]/deltayg;
    end
end
% top region_V2
for m=1:Nxg
    for n=1:1
        index=(n-1)*Nxg+m;
        line_V2(:,index)=[index,index];
        column_V2(:,index)=[index,(Nyg-1)*Nxg+m];
        value_V2(:,index)=[1,-exp(j*k*u2'*Pyg)]/deltayg;
    end
end




U1=sparse(line_U1,column_U1,value_U1,N,N);
U2=sparse(line_U2,column_U2,value_U2,N,N);
V1=sparse(line_V1,column_V1,value_V1,N,N);
V2=sparse(line_V2,column_V2,value_V2,N,N);
epsV=sparse(line_eps,column_eps,eps_arr,N,N);
muV=sparse(line_eps,column_eps,mu_arr,N,N);
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

if(polarflag==0)
    %TM polarization
    A=(U1*(mu21_inv*V2-mu22_inv*V1)-U2*(mu11_inv*V2-mu12_inv*V1));
    [V,D] = eigs(A,epsV,N_max,'sm');       % find eigenmodes and eigen-frequency
elseif(polarflag==1)
    %TE polarization
    A=-(V1*(eps21_inv*U2-eps22_inv*U1)-V2*(eps11_inv*U2-eps12_inv*U1));
    [V,D] = eigs(A,muV,N_max,'sm');       % find eigenmodes and eigen-frequency    
    eps_arr=mu_arr;    
end



[omega,order]=sort(sqrt(abs(real(diag(D)))));  %  from small to large eigenvalues



fieldz=zeros(N,N_max);

for m=1:Nxg
    for n=1:Nyg
        index=(n-1)*Nxg+m;
        
        Lxg=(m-1/2)*deltaxg;
        Lyg=(n-1/2)*deltayg;
        fieldz(index,:)=V(index,:)*exp(j*k*u1'*Lxg)*exp(j*k*u2'*Lyg);
%           fieldz(index,:)=V(index,:);           
    end
end



% draw eigenmodes

% X=g(1,1)*(1:Nxg)'+g(1,2)*(1:Nyg)';
% Y=sqrt(3)/2*(1:Nyg)';
X=zeros(Nxg,Nyg);Y=zeros(Nxg,Nyg);
field_2d=zeros(Nxg,Nyg);
if (flag==1) % flag control
    figure(2)
    for mm=1:N_max
        field_1d=V(:,order(mm));  % eigenmodes from low to high order
        for m=1:Nxg
            for n=1:Nyg
                X(m,n)=m+0.5*n;Y(m,n)=sqrt(3)/2*n;
                index=(n-1)*Nxg+m;
                field_2d(m,n)=field_1d(index);
            end
        end
        subplot(2,2,mm)
        pcolor(X,Y,real(field_2d));
        shading interp
        colorbar
        text_um=num2str(mm);
        title(text_um)
        axis equal
%         axis([1 Nxg 1 Nyg])
    end
end

%  model the unit cell
function flag=shape_f(x,y,xc,yc,g)
Pxg=1.0;Pyg=1.0;
flag=0;
radius=0.2*Pxg;  % radius of circle in unit cells ***
xc=1/3*Pxg;yc=1/3*Pyg;
distance=[x-xc,y-yc,0]*g*[x-xc;y-yc;0];
if ((distance)<radius^2)
    flag=1;
end
xc=2/3*Pxg;yc=2/3*Pyg;
distance=[x-xc,y-yc,0]*g*[x-xc;y-yc;0];
if ((distance)<radius^2)
    flag=1;
end
% xc=0;yc=1;
% distance=[x-xc,y-yc,0]*g*[x-xc;y-yc;0];
% if ((distance)<radius^2)
%     flag=0;
% end
% xc=1;yc=1;
% distance=[x-xc,y-yc,0]*g*[x-xc;y-yc;0];
% if ((distance)<radius^2)
%     flag=0;
% end

function [sol_int] = fun_HZ_int(in_m)
%fun_HZ_int calculates the integration coefficients for the Herve-Zaoui Problem
%   INPUT: is a matrix called in_m containing the bluk moduli k, shear
%   moduli mu, and the radii of the phases (R1<R2<R3,...,Rn=NaN)
%   e.g.: recycled concrete with ITZ
%k1=25;mu1=26;R1=1; %original agg
%k2=20;mu2=12;R2=1.2; %old cement paste
%k3=6;mu3=3;R3=1.3; %new ITZ
%k4=7;mu4=4;R4=1.9; %new paste
%in_m=[k1, mu1, R1;...
%    k2, mu2, R2;...
%    k3, mu3, R3;...
%    k4, mu4, R4];
%   OUTPUT: sol_int...matrix with integration coefficients for all phases
%           note that Avol_mat=sol_int(:,5) and Adev_mat=sol_int(:,7)
%           are the matrices with infinite vol/dev part of strain
%           concentration tensors
%


%% START
% check input
assert(sum(sum(in_m<0))==0,'ERROR: negative input -> check input');
assert(issorted(in_m(:,3)),'ERROR: radii not increasing -> check input')

%last line is the matrix around the n-layered inclusion
n=length(in_m(:,1))-1; %number of layers

% add Poisson's ratios to last column
nu_column=zeros(n+1,1);
for i=1:(n+1)
    [E,nu]=fun_Enu_from_kmu(in_m(i,1),in_m(i,2));
    nu_column(i)=nu;
end
in_m=[in_m,nu_column];

%% isotropic pressure
% transition matrices
N_cell=cell(n,1);
Q_cell=cell(n,1);
Qi=eye(2);
for i=1:n
    
    Ji=[in_m(i,3),1/in_m(i,3)^2;...
        3*in_m(i,1),-4*in_m(i,2)/in_m(i,3)^3];
    Jip1=[in_m(i,3),1/in_m(i,3)^2;...
        3*in_m(i+1,1),-4*in_m(i+1,2)/in_m(i,3)^3];
    
    Ni=inv(Jip1)*Ji;
    N_cell{i}=Ni;
    Qi=Ni*Qi;
    Q_cell{i}=Qi;
end

% integration constants
F_mat=zeros(n+1,1);
G_mat=zeros(n+1,1);
F_mat(1)=1/Q_cell{n}(1,1);
G_mat(1)=0;
for i=2:(n+1)
    F_mat(i)=Q_cell{i-1}(1,1)/Q_cell{n}(1,1);
    G_mat(i)=Q_cell{i-1}(2,1)/Q_cell{n}(1,1);
end
%disp('isotropic constants [F1,F2...Fn+1;G1,G2,...Gn+1]')
%disp([F_mat';G_mat'])

%% simple shear
% transition matrices
M_cell=cell(n,1);
P_cell=cell(n,1);
Pi=eye(4);
for i=1:n
    nu=in_m(i,4);
    Li=[in_m(i,3),-6*nu/(1-2*nu)*in_m(i,3)^3,3/in_m(i,3)^4,(5-4*nu)/(1-2*nu)*1/in_m(i,3)^2;...
        in_m(i,3),-(7-4*nu)/(1-2*nu)*in_m(i,3)^3,-2/in_m(i,3)^4,2/in_m(i,3)^2;...
        in_m(i,2),3*nu/(1-2*nu)*in_m(i,2)*in_m(i,3)^2,-12/in_m(i,3)^5*in_m(i,2),2*(nu-5)/(1-2*nu)*in_m(i,2)/in_m(i,3)^3;...
        in_m(i,2),-(7+2*nu)/(1-2*nu)*in_m(i,2)*in_m(i,3)^2,8/in_m(i,3)^5*in_m(i,2),2*(1+nu)/(1-2*nu)*in_m(i,2)/in_m(i,3)^3];
    nu=in_m(i+1,4);
    Lip1=[in_m(i,3),-6*nu/(1-2*nu)*in_m(i,3)^3,3/in_m(i,3)^4,(5-4*nu)/(1-2*nu)*1/in_m(i,3)^2;...
        in_m(i,3),-(7-4*nu)/(1-2*nu)*in_m(i,3)^3,-2/in_m(i,3)^4,2/in_m(i,3)^2;...
        in_m(i+1,2),3*nu/(1-2*nu)*in_m(i+1,2)*in_m(i,3)^2,-12/in_m(i,3)^5*in_m(i+1,2),2*(nu-5)/(1-2*nu)*in_m(i+1,2)/in_m(i,3)^3;...
        in_m(i+1,2),-(7+2*nu)/(1-2*nu)*in_m(i+1,2)*in_m(i,3)^2,8/in_m(i,3)^5*in_m(i+1,2),2*(1+nu)/(1-2*nu)*in_m(i+1,2)/in_m(i,3)^3];
    Mi=inv(Lip1)*Li;
    M_cell{i}=Mi;
    Pi=Mi*Pi;
    P_cell{i}=Pi;
end

% integration constants
Adev_mat=zeros(n+1,1);
A_mat=zeros(n+1,1);
B_mat=zeros(n+1,1);
C_mat=zeros(n+1,1);
D_mat=zeros(n+1,1);
A_mat(1)=P_cell{n}(2,2)/(P_cell{n}(1,1)*P_cell{n}(2,2)-P_cell{n}(1,2)*P_cell{n}(2,1));
B_mat(1)=-P_cell{n}(2,1)/(P_cell{n}(1,1)*P_cell{n}(2,2)-P_cell{n}(1,2)*P_cell{n}(2,1));
for i=2:(n+1)
    W_mat=P_cell{i-1}./(P_cell{n}(1,1)*P_cell{n}(2,2)-P_cell{n}(1,2)*P_cell{n}(2,1))*[P_cell{n}(2,2);-P_cell{n}(2,1);0;0];
    A_mat(i)=W_mat(1);
    B_mat(i)=W_mat(2);
    C_mat(i)=W_mat(3);
    D_mat(i)=W_mat(4);
end

% infinite strain concentration tensors
Adev_mat(1)=1/A_mat(n+1)*(A_mat(1)-21/5*(in_m(1,3)^5)/((1-2*in_m(1,4))*(in_m(1,3)^3))*B_mat(1));
for i=2:(n+1)
    Adev_mat(i)=1/A_mat(n+1)*(A_mat(i)-21/5*(in_m(i,3)^5-in_m(i-1,3)^5)/((1-2*in_m(i,4))*(in_m(i,3)^3-in_m(i-1,3)^3))*B_mat(i));
end

%disp('shear constants [A1,A2...An+1;B1,B2,...Bn+1;C1,C2,...Cn+1;D1,D2,...Dn+1]')
%disp([A_mat';B_mat';C_mat';D_mat'])
sol_int=[A_mat,B_mat,C_mat,D_mat,F_mat,G_mat,Adev_mat];

end


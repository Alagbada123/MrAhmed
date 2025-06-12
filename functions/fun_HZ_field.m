function [eps_cart_cell,si_cart_cell,r_list,angle_list] = fun_HZ_field(in_m,Einf,sol_int,steps,steps_angle)
%fun_HZ_field Solution for strain and stress fields for Herve_Zaoui Problem
%   INPUT:  in_m...matrix with bulk moduli, shear moduli, and radii
%           Einf...strain at infinite boundary, Cartesian components in
%           compressed notation
%           sol_int...solution of the integration (output of fun_HZ_int)
%           steps... discretization of radial coordinate in each phase
%           steps_angle... discretization of position angles
%   OUTPUT: eps_cart_cell{aziit,zeniit} and si_cart_cell{aziit,zeniit} are
%   cell arrays referring to different position angles between zero and 90
%   for both, azimuth and zenith angle, respectively. The elements in this
%   array are matrices whereby the column(i) represent the Cartesian
%   strain/stress components (in compressed notation) at radial coordinate
%   r(i).

%% START
% check input
assert(sum(sum(in_m<0))==0,'ERROR: negative input -> check input');
assert(issorted(in_m(:,3)),'ERROR: radii not increasing -> check input')

%last line is the matrix around the n-layered inclusion
n=length(in_m(:,1))-1; %number of layers

% add Poisson's ratios to last column
nu_column=zeros(n+1,1);
for i=1:(n+1)
    [~,nu]=fun_Enu_from_kmu(in_m(i,1),in_m(i,2));
    nu_column(i)=nu;
end
in_m=[in_m,nu_column];

% integration coefficients
    A_mat=sol_int(:,1);
    B_mat=sol_int(:,2);
    C_mat=sol_int(:,3);
    D_mat=sol_int(:,4);
    F_mat=sol_int(:,5);
    G_mat=sol_int(:,6);
    

%% displacements, strains and stresses
% discretization and pointer
r=linspace(0,in_m(1,3),steps);
p=ones(1,n*steps);
for i=1:n
    if i<n
        r=[r,linspace(in_m(i,3),in_m(i+1,3),steps)];
    else
        r=[r,linspace(in_m(i,3),2*in_m(i,3),steps)];
    end
end
for i=1:(n+1)
    p(((i-1)*steps+1):(i*steps))=i;
end

% 1) ISOTROPIC PART
eps_iso_sph=zeros(6,(n+1)*steps);
for ri=1:(n+1)*steps
    eps_iso_sph(:,ri)=[F_mat(p(ri))-2*G_mat(p(ri))/(r(ri)^3);F_mat(p(ri))+G_mat(p(ri))/(r(ri)^3);F_mat(p(ri))+G_mat(p(ri))/(r(ri)^3);0;0;0];
end

% 2) SHEAR PART
epsrr_shear=zeros(1,(n+1)*steps);
epsrt_shear=zeros(1,(n+1)*steps);
epsrp_shear=zeros(1,(n+1)*steps);
epstt_shear=zeros(1,(n+1)*steps);
epstp_shear=zeros(1,(n+1)*steps);
epspp_shear=zeros(1,(n+1)*steps);
eps_shear_sph=zeros(6,(n+1)*steps);
eps_shear_cart=zeros(6,(n+1)*steps);

% position angles
phi_glob_list=zeros(steps_angle);
theta_glob_list=zeros(steps_angle);
theta_loc_cell=repmat({zeros(steps_angle,steps_angle)},1,5);
phi_loc_cell=repmat({zeros(steps_angle,steps_angle)},1,5);
eps_shear_sph_cell=cell(5,steps_angle,steps_angle);

% Transformation matrices
Q2cell=cell(1,5);
Q2cell{1}=[1/sqrt(2),1/sqrt(2),0;-1/sqrt(2),1/sqrt(2),0;0,0,1];
Q2cell{2}=[0,0,1;1,0,0;0,1,0]*[1/sqrt(2),0,-1/sqrt(2);0,1,0;1/sqrt(2),0,1/sqrt(2)];
Q2cell{3}=[0,1/sqrt(2),1/sqrt(2);0,-1/sqrt(2),1/sqrt(2);1,0,0];
Q2cell{4}=[1,0,0;0,1,0;0,0,1];
Q2cell{5}=[0,0,1;0,1,0;-1,0,0];
Q4cell=cell(1,5);
for i=1:5
Q4cell{i}=[Q2cell{i}(1,1)^2,Q2cell{i}(1,2)^2,Q2cell{i}(1,3)^2,Q2cell{i}(1,2)*Q2cell{i}(1,3)*sqrt(2),Q2cell{i}(1,1)*Q2cell{i}(1,3)*sqrt(2),Q2cell{i}(1,1)*Q2cell{i}(1,2)*sqrt(2);...
                Q2cell{i}(2,1)^2,Q2cell{i}(2,2)^2,Q2cell{i}(2,3)^2,Q2cell{i}(2,2)*Q2cell{i}(2,3)*sqrt(2),Q2cell{i}(2,1)*Q2cell{i}(2,3)*sqrt(2),Q2cell{i}(2,1)*Q2cell{i}(2,2)*sqrt(2);...
                Q2cell{i}(3,1)^2,Q2cell{i}(3,2)^2,Q2cell{i}(3,3)^2,Q2cell{i}(3,2)*Q2cell{i}(3,3)*sqrt(2),Q2cell{i}(3,1)*Q2cell{i}(3,3)*sqrt(2),Q2cell{i}(3,1)*Q2cell{i}(3,2)*sqrt(2);...
                sqrt(2)*Q2cell{i}(3,1)*Q2cell{i}(2,1),sqrt(2)*Q2cell{i}(3,2)*Q2cell{i}(2,2),sqrt(2)*Q2cell{i}(3,3)*Q2cell{i}(2,3),Q2cell{i}(3,2)*Q2cell{i}(2,3)+Q2cell{i}(3,3)*Q2cell{i}(2,2),Q2cell{i}(3,1)*Q2cell{i}(2,3)+Q2cell{i}(3,3)*Q2cell{i}(2,1),Q2cell{i}(3,1)*Q2cell{i}(2,2)+Q2cell{i}(3,2)*Q2cell{i}(2,1);...
                sqrt(2)*Q2cell{i}(3,1)*Q2cell{i}(1,1),sqrt(2)*Q2cell{i}(3,2)*Q2cell{i}(1,2),sqrt(2)*Q2cell{i}(3,3)*Q2cell{i}(1,3),Q2cell{i}(3,2)*Q2cell{i}(1,3)+Q2cell{i}(3,3)*Q2cell{i}(1,2),Q2cell{i}(3,1)*Q2cell{i}(1,3)+Q2cell{i}(3,3)*Q2cell{i}(1,1),Q2cell{i}(3,1)*Q2cell{i}(1,2)+Q2cell{i}(3,2)*Q2cell{i}(1,1);...
                sqrt(2)*Q2cell{i}(2,1)*Q2cell{i}(1,1),sqrt(2)*Q2cell{i}(2,2)*Q2cell{i}(1,2),sqrt(2)*Q2cell{i}(2,3)*Q2cell{i}(1,3),Q2cell{i}(2,2)*Q2cell{i}(1,3)+Q2cell{i}(2,3)*Q2cell{i}(1,2),Q2cell{i}(2,1)*Q2cell{i}(1,3)+Q2cell{i}(2,3)*Q2cell{i}(1,1),Q2cell{i}(2,1)*Q2cell{i}(1,2)+Q2cell{i}(2,2)*Q2cell{i}(1,1)];
end


for i=1:5
    for aziit=1:steps_angle
        phi_glob=(aziit-1)/(steps_angle-1)*90*pi/180; %global
        phi_glob_list(aziit)=phi_glob;
        for zeniit=1:steps_angle
            theta_glob=(zeniit-1)/(steps_angle-1)*90*pi/180;
            theta_glob_list(zeniit)=theta_glob;
            
            % r vector in global and local Cartesian base
            r_glob=[sin(theta_glob)*cos(phi_glob),sin(theta_glob)*sin(phi_glob),cos(theta_glob)]';
            r_loc=Q2cell{i}*r_glob;
            
            % Position angles in local definition
            r_loc_proj=[r_loc(1),r_loc(2),0]';
            r_loc_proj=r_loc_proj/norm(r_loc_proj);
            theta_loc=acos([0,0,1]*r_loc);
            phi_loc=atan2(r_loc_proj(2),r_loc_proj(1)); %atan=GK/AK
            phi_loc(isnan(phi_loc))=0; % if NaN then it does not matter, I use zero then
            theta_loc_cell{i}(aziit,zeniit)=theta_loc;
            phi_loc_cell{i}(aziit,zeniit)=phi_loc;
            
            % transformation matrix from loc Cart to loc sph
            Q4loc=fun_Q4_standard(phi_loc,theta_loc);
            
            % transformation matrix from glob Cart to glob sph
            Q4glob=fun_Q4_standard(phi_glob,theta_glob);

            for ri=1:(n+1)*steps
                epsrr_shear(ri)=(A_mat(p(ri))-18*in_m(p(ri),4)*B_mat(p(ri))*r(ri)^2/(1-2*in_m(p(ri),4))-12*C_mat(p(ri))/r(ri)^5-(2*(5-4*in_m(p(ri),4)))*D_mat(p(ri))/((1-2*in_m(p(ri),4))*r(ri)^3))*sin(theta_loc)^2*cos(2*phi_loc);
                epsrt_shear(ri)=sin(theta_loc)*cos(2*phi_loc)*cos(theta_loc)*(2*B_mat(p(ri))*in_m(p(ri),4)*r(ri)^7+7*B_mat(p(ri))*r(ri)^7+2*A_mat(p(ri))*in_m(p(ri),4)*r(ri)^5-A_mat(p(ri))*r(ri)^5-2*D_mat(p(ri))*in_m(p(ri),4)*r(ri)^2-2*D_mat(p(ri))*r(ri)^2+16*C_mat(p(ri))*in_m(p(ri),4)-8*C_mat(p(ri)))/((-1+2*in_m(p(ri),4))*r(ri)^5);
                epsrp_shear(ri)=-sin(theta_loc)*sin(2*phi_loc)*(2*B_mat(p(ri))*in_m(p(ri),4)*r(ri)^7+7*B_mat(p(ri))*r(ri)^7+2*A_mat(p(ri))*in_m(p(ri),4)*r(ri)^5-A_mat(p(ri))*r(ri)^5-2*D_mat(p(ri))*in_m(p(ri),4)*r(ri)^2-2*D_mat(p(ri))*r(ri)^2+16*C_mat(p(ri))*in_m(p(ri),4)-8*C_mat(p(ri)))/((-1+2*in_m(p(ri),4))*r(ri)^5);
                epstt_shear(ri)=cos(2*phi_loc)*(-14*cos(theta_loc)^2*B_mat(p(ri))*in_m(p(ri),4)*r(ri)^7+14*cos(theta_loc)^2*B_mat(p(ri))*r(ri)^7+2*cos(theta_loc)^2*A_mat(p(ri))*in_m(p(ri),4)*r(ri)^5+10*B_mat(p(ri))*r(ri)^7*in_m(p(ri),4)-cos(theta_loc)^2*A_mat(p(ri))*r(ri)^5-7*B_mat(p(ri))*r(ri)^7+4*cos(theta_loc)^2*D_mat(p(ri))*in_m(p(ri),4)*r(ri)^2+cos(theta_loc)^2*D_mat(p(ri))*r(ri)^2-14*cos(theta_loc)^2*C_mat(p(ri))*in_m(p(ri),4)+7*cos(theta_loc)^2*C_mat(p(ri))-3*D_mat(p(ri))*r(ri)^2+10*C_mat(p(ri))*in_m(p(ri),4)-5*C_mat(p(ri)))/(r(ri)^5*(-1+2*in_m(p(ri),4)));
                epstp_shear(ri)=-(-4*B_mat(p(ri))*in_m(p(ri),4)*r(ri)^7+7*B_mat(p(ri))*r(ri)^7+2*A_mat(p(ri))*in_m(p(ri),4)*r(ri)^5-A_mat(p(ri))*r(ri)^5+4*D_mat(p(ri))*in_m(p(ri),4)*r(ri)^2-2*D_mat(p(ri))*r(ri)^2-4*C_mat(p(ri))*in_m(p(ri),4)+2*C_mat(p(ri)))*sin(2*phi_loc)*cos(theta_loc)/((-1+2*in_m(p(ri),4))*r(ri)^5);
                epspp_shear(ri)=-cos(2*phi_loc)*(10*cos(theta_loc)^2*B_mat(p(ri))*in_m(p(ri),4)*r(ri)^7-7*cos(theta_loc)^2*B_mat(p(ri))*r(ri)^7-14*B_mat(p(ri))*r(ri)^7*in_m(p(ri),4)+14*B_mat(p(ri))*r(ri)^7+2*A_mat(p(ri))*r(ri)^5*in_m(p(ri),4)-A_mat(p(ri))*r(ri)^5-3*cos(theta_loc)^2*D_mat(p(ri))*r(ri)^2+10*cos(theta_loc)^2*C_mat(p(ri))*in_m(p(ri),4)+4*D_mat(p(ri))*r(ri)^2*in_m(p(ri),4)-5*cos(theta_loc)^2*C_mat(p(ri))+D_mat(p(ri))*r(ri)^2-14*C_mat(p(ri))*in_m(p(ri),4)+7*C_mat(p(ri)))/(r(ri)^5*(-1+2*in_m(p(ri),4)));
                eps_shear_sph(:,ri)=[epsrr_shear(ri);epstt_shear(ri);epspp_shear(ri);epstp_shear(ri)*sqrt(2);epsrp_shear(ri)*sqrt(2);epsrt_shear(ri)*sqrt(2)]; % in lokal spherical base                
            end
            eps_shear_cart=transpose(Q4loc)*eps_shear_sph;  % in lokal Cartesian base
            eps_shear_cart=transpose(Q4cell{i})*eps_shear_cart;  % in global Cartesian base
            eps_shear_sph=Q4glob*eps_shear_cart;  % in global spherical base
            eps_shear_sph_cell{i,aziit,zeniit}=eps_shear_sph; 
        end
    end
end
r_list=r;
angle_list=phi_glob_list;

%% superimpose all 6 strain states
% split up into 5 loading cases
TraceE=(Einf(1)+Einf(2)+Einf(3));
Einf_iso=TraceE/3*[1,1,1,0,0,0]';
einf=Einf-Einf_iso;
gamma_vect=[1/sqrt(2)*einf(6),1/sqrt(2)*einf(5),1/sqrt(2)*einf(4),einf(1),einf(3)];

% initialize cells
eps_sph_cell=cell(steps_angle,steps_angle);
eps_cart_cell=cell(steps_angle,steps_angle);
si_sph_cell=cell(steps_angle,steps_angle);
si_cart_cell=cell(steps_angle,steps_angle);

for aziit=1:steps_angle
    phi=(aziit-1)/(steps_angle-1)*90*pi/180;
    for zeniit=1:steps_angle
        theta=(zeniit-1)/(steps_angle-1)*90*pi/180;
        
        % stress/strain in global spherical coordinates
        eps_sph_cell{aziit,zeniit}=TraceE/3*eps_iso_sph+eps_shear_sph_cell{1,aziit,zeniit}*gamma_vect(1)+eps_shear_sph_cell{2,aziit,zeniit}*gamma_vect(2)+eps_shear_sph_cell{3,aziit,zeniit}*gamma_vect(3)+eps_shear_sph_cell{4,aziit,zeniit}*gamma_vect(4)+eps_shear_sph_cell{5,aziit,zeniit}*gamma_vect(5);
        for ri=1:(n+1)*steps
            si_sph_cell{aziit,zeniit}(:,ri)=fun_Cfromkmu(in_m(p(ri),1),in_m(p(ri),2))*eps_sph_cell{aziit,zeniit}(:,ri);
        end
        
        % stress/strain in global Cartesian coordinates
        Q4=fun_Q4_standard(phi,theta); %global position angles
        Q4t=transp(Q4);
        eps_cart_cell{aziit,zeniit}=Q4t*eps_sph_cell{aziit,zeniit};
        si_cart_cell{aziit,zeniit}=Q4t*si_sph_cell{aziit,zeniit};
    end
end
end


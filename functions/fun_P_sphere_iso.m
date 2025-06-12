% function for computation of P-tensor for spherical inclusions in 
% isotropic matrix material
%
%-----------------------------------------------------------------------
% INPUT: matrix stiffness tensor C0 (isotropic)
%
% OUTPUT: P-tensor
%
%-----------------------------------------------------------------------
%
function[Psph]=fun_P_sphere_iso(C0)

% computation of elastic constants of matrix material
D0=inv(C0);
E0=1/D0(1,1);
mu0=0.5*C0(4,4);
k0=E0*mu0/(3*(3*mu0-E0));

alpha0=3*k0/(3*k0+4*mu0);
beta0=6*(k0+2*mu0)/(5*(3*k0+4*mu0));

I=[1 0 0 0 0 0;...
   0 1 0 0 0 0;...
   0 0 1 0 0 0;...
   0 0 0 1 0 0;...
   0 0 0 0 1 0;...
   0 0 0 0 0 1];

K=[1/3 1/3 1/3 0 0 0;...
   1/3 1/3 1/3 0 0 0;...
   1/3 1/3 1/3 0 0 0;...
   0   0   0   0 0 0;...
   0   0   0   0 0 0;...
   0   0   0   0 0 0];

J=I-K;

Ssph=alpha0*K+beta0*J;
Psph=Ssph*inv(C0);

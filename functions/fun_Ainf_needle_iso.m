% function for computation of Ainf for cylindric isotropic inclusions
% orientated in all space directions, embedded in isotropic matrix material
%
%-----------------------------------------------------------------------
% INPUT: inclusion stiffness tensor Cinc (isotropic)
%        infinite matrix stiffness tensor Cinf (isotropic)
%
% OUTPUT: P-tensor
%
%-----------------------------------------------------------------------
%
function[Ainfneedle]=fun_Ainf_needle_iso(Cinc,Cinf)

% computation of elastic constants of matrix material
muinc=0.5*Cinc(6,6);
kinc=Cinc(1,1)-4/3*0.5*Cinc(6,6);

muinf=0.5*Cinf(6,6);
kinf=Cinf(1,1)-4/3*0.5*Cinf(6,6);

Avolinf_needle = (1/3)*(3*kinf+muinc+3*muinf)/(3*kinc+muinc+3*muinf);
Adevinf_needle = (1/10)*(9*kinc*muinc^2*kinf+84*kinf*muinf^3+64*muinf^4+21*kinc*muinc^2*muinf+81*kinc*kinf*muinf^2+120*kinf*muinc*muinf^2+63*kinc*muinf^3+184*muinc*muinf^3+90*kinf*muinf*kinc*muinc+36*kinf*muinf*muinc^2+156*kinc*muinc*muinf^2+72*muinc^2*muinf^2)/((muinf+muinc)*(3*kinf*muinc+7*muinf*muinc+3*kinf*muinf+muinf^2)*(3*kinc+muinc+3*muinf));

I=[1 0 0 0 0 0;...
   0 1 0 0 0 0;...
   0 0 1 0 0 0;...
   0 0 0 1 0 0;...
   0 0 0 0 1 0;...
   0 0 0 0 0 1];

J=[1/3 1/3 1/3 0 0 0;...
   1/3 1/3 1/3 0 0 0;...
   1/3 1/3 1/3 0 0 0;...
   0   0   0   0 0 0;...
   0   0   0   0 0 0;...
   0   0   0   0 0 0];

K=I-J;

Ainfneedle=3*Avolinf_needle*J+2*Adevinf_needle*K;
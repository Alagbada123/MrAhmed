% function for computation of P-tensor for ellipsoidal inclusions in 
% transversely isotropic matrix material with axis 3 as axis of rotational 
% symmetry
%
%-----------------------------------------------------------------------
% INPUT: 
% C0       ... (trans. iso.) matrix stiffness tensor C0
% aspRat   ... aspect ratio (a1/a2)
% slendRat ... slenderness ratio (a1/a3)
%
% OUTPUT: P-tensor
% Pellips ... P-tensor
%
%-----------------------------------------------------------------------
%
function[Pellips]=fun_P_ellipsoid_transiso(C0,aspRat,slendRat)

% determination of (independent) components of stiffness tensor C0
C1111=C0(1,1);
C3333=C0(3,3);
C1122=C0(1,2);
C2233=C0(2,3);
C2323=C0(4,4)/2;


% computation of tensor components by integration over unit sphere
P1111=1/(16*pi)*dblquad(@fun_Ptrielli_1111,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);
P1122=1/(16*pi)*dblquad(@fun_Ptrielli_1122,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);
P1133=1/(16*pi)*dblquad(@fun_Ptrielli_1133,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);
P2222=1/(16*pi)*dblquad(@fun_Ptrielli_2222,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);
P2233=1/(16*pi)*dblquad(@fun_Ptrielli_2233,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);
P3333=1/(16*pi)*dblquad(@fun_Ptrielli_3333,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);
P2323=1/(16*pi)*dblquad(@fun_Ptrielli_2323,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);
P1313=1/(16*pi)*dblquad(@fun_Ptrielli_1313,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);
P1212=1/(16*pi)*dblquad(@fun_Ptrielli_1212,0,2*pi,0,pi,1e-6,@quadl,C1111,C1122,C2233,C3333,C2323,aspRat,slendRat);


% compilation of P-tensor in compressed notation
Pellips=zeros(6,6);
Pellips(1,1)=P1111;
Pellips(2,2)=P2222;
Pellips(3,3)=P3333;
Pellips(1,2)=P1122;
Pellips(2,1)=Pellips(1,2);
Pellips(1,3)=P1133;
Pellips(3,1)=Pellips(1,3);
Pellips(2,3)=P2233;
Pellips(3,2)=Pellips(2,3);

Pellips(4,4)=2*P2323;
Pellips(5,5)=2*P1313;
Pellips(6,6)=2*P1212;

% function for computation of transformation matrix for transformation of
% 4th order tensors,
% Definition of base vectors according to Physics
%
%-----------------------------------------------------------------------
% INPUT: azi  ... azimuth angle, reference = e_x
%        zeni ... azimuth angle, reference = e_z
% OUTPUT: Q4-Matrix
%
%-----------------------------------------------------------------------
%
function[Q4_standard]=fun_Q4_standard(azi,zeni)
Q21=cos(azi)*cos(zeni); %Definition Standard!
Q22=sin(azi)*cos(zeni);
Q23=-sin(zeni);
Q31=-sin(azi);
Q32=cos(azi);
Q33=0;
Q11=cos(azi)*sin(zeni);
Q12=sin(azi)*sin(zeni);
Q13=cos(zeni);

Q4_standard=[Q11^2,Q12^2,Q13^2,Q12*Q13*sqrt(2),Q11*Q13*sqrt(2),Q11*Q12*sqrt(2);...
    Q21^2,Q22^2,Q23^2,Q22*Q23*sqrt(2),Q21*Q23*sqrt(2),Q21*Q22*sqrt(2);...
    Q31^2,Q32^2,Q33^2,Q32*Q33*sqrt(2),Q31*Q33*sqrt(2),Q31*Q32*sqrt(2);...
    sqrt(2)*Q31*Q21,sqrt(2)*Q32*Q22,sqrt(2)*Q33*Q23,Q32*Q23+Q33*Q22,Q31*Q23+Q33*Q21,Q31*Q22+Q32*Q21;...
    sqrt(2)*Q31*Q11,sqrt(2)*Q32*Q12,sqrt(2)*Q33*Q13,Q32*Q13+Q33*Q12,Q31*Q13+Q33*Q11,Q31*Q12+Q32*Q11;...
    sqrt(2)*Q21*Q11,sqrt(2)*Q22*Q12,sqrt(2)*Q23*Q13,Q22*Q13+Q23*Q12,Q21*Q13+Q23*Q11,Q21*Q12+Q22*Q11];
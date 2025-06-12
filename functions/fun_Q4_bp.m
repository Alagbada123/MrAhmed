% function for computation of transformation matrix for transformation of
% 4th order tensors,
% Definition of base vectors according to Bernhard Pichler, thus for
% azi=zeni=0, the e_x,e_y,e_z base vectors are equal to the base vectors
% e_zeni, e_azi, e_r
% This defnition deviates from standard definition, but allows to directly
% use this matrix for stroud integration
%
%-----------------------------------------------------------------------
% INPUT: azi  ... azimuth angle, reference = e_x
%        zeni ... azimuth angle, reference = e_z
% OUTPUT: Q4-Matrix
%
%-----------------------------------------------------------------------
%
function[Q4_bp]=fun_Q4_bp(azi,zeni)
Q11=cos(azi)*cos(zeni); %Definition Bernhard!
Q12=sin(azi)*cos(zeni);
Q13=-sin(zeni);
Q21=-sin(azi);
Q22=cos(azi);
Q23=0;
Q31=cos(azi)*sin(zeni);
Q32=sin(azi)*sin(zeni);
Q33=cos(zeni);

Q4_bp=[Q11^2,Q12^2,Q13^2,Q12*Q13*sqrt(2),Q11*Q13*sqrt(2),Q11*Q12*sqrt(2);...
    Q21^2,Q22^2,Q23^2,Q22*Q23*sqrt(2),Q21*Q23*sqrt(2),Q21*Q22*sqrt(2);...
    Q31^2,Q32^2,Q33^2,Q32*Q33*sqrt(2),Q31*Q33*sqrt(2),Q31*Q32*sqrt(2);...
    sqrt(2)*Q31*Q21,sqrt(2)*Q32*Q22,sqrt(2)*Q33*Q23,Q32*Q23+Q33*Q22,Q31*Q23+Q33*Q21,Q31*Q22+Q32*Q21;...
    sqrt(2)*Q31*Q11,sqrt(2)*Q32*Q12,sqrt(2)*Q33*Q13,Q32*Q13+Q33*Q12,Q31*Q13+Q33*Q11,Q31*Q12+Q32*Q11;...
    sqrt(2)*Q21*Q11,sqrt(2)*Q22*Q12,sqrt(2)*Q23*Q13,Q22*Q13+Q23*Q12,Q21*Q13+Q23*Q11,Q21*Q12+Q22*Q11];
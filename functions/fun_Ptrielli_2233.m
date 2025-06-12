% P2233 for ellipsoidal inclusions 
% with aspect ratio asp=a1/a2 and
% slenderness ratio slend=a1/a3
% in transversely isotropic matrix 
%----------------------------------
function[z]=fun_Ptrielli_2233(phi,theta,C1111,C1122,C2233,C3333,C2323,asp,slend);

z = (4..*asp.^2.*slend.^2.*cos(theta).^2.*sin(phi).^2.*sin(theta).^3.* ((-1..*C2233 - 1..*C2323).*C2323.*slend.^2.*cos(theta).^2 + (-0.5.*C1111.*C2233 + 0.5.*C1122.*C2233 - 0.5.*C1111.*C2323 + 0.5.*C1122.*C2323).*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2))./ (C2323.^2.*C3333.*slend.^6.*cos(theta).^6 - 2..*C2323.*(0.5.*C2233.^2 + 1..*C2233.*C2323 - 0.75.*C1111.*C3333 + 0.25.*C1122.*C3333).*slend.^4.*cos(theta).^4.*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2 + slend.^2.*cos(theta).^2.* ((C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*cos(phi).^4 + 1..*asp.^2.* (C1122.*C2233.*(1..*C2233 + 2..*C2323) + 1..*C1111.^2.*C3333 + C1111.*(-1..*C2233.^2 - 2..*C2233.*C2323 + 2..*C2323.^2 - 1..*C1122.*C3333)).* cos(phi).^2.*sin(phi).^2 + asp.^4.*(C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*sin(phi).^4).*sin(theta).^4 + 0.5.*C1111.*(1..*C1111 - 1..*C1122).*C2323.*(cos(phi).^6 + 3..*asp.^2.*cos(phi).^4.*sin(phi).^2 + 3..*asp.^4.*cos(phi).^2.*sin(phi).^4 + asp.^6.*sin(phi).^6).*sin(theta).^6);

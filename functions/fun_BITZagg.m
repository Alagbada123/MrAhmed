function [ BITZsand_sph ] = fun_BITZagg(kagg,kITZ,muagg,muITZ)
%fun_BITZagg Stress concentration tensor to ITZ
%   calculates the spherical components of the 4th order stress concentration tensor
%   BITZagg in the theta,phi,r-base for downscaling aggregate stresses to
%   ITZ stresses; all phases are isotropic, spherical aggregates, perfect
%   bond-related continutity conditions
                    deltaITZsand= kagg*muagg*(3*kITZ+4*muITZ);
                    
                    BITZsandrrrr= 1.;
                    BITZsandtttt= muITZ*(3*kagg*kITZ+2*kagg*muITZ+2*kITZ*muagg)/deltaITZsand;
                    BITZsandpppp= BITZsandtttt;
                    BITZsandttpp= 2*muITZ*(kITZ*muagg-kagg*muITZ)/deltaITZsand;
                    BITZsandpptt= BITZsandttpp;
                    BITZsandttrr= (3*kagg*kITZ*(muagg-muITZ)-2*muagg*muITZ*(kagg-kITZ))/deltaITZsand;
                    BITZsandpprr= BITZsandttrr;
                    BITZsandrtrt= 0.5;
                    BITZsandrprp= 0.5;
                    BITZsandtptp= muITZ/(2*muagg);
                    
                    BITZsand_sph= [BITZsandtttt, BITZsandttpp, BITZsandttrr, 0., 0., 0.;...
                        BITZsandpptt, BITZsandpppp, BITZsandpprr, 0., 0., 0.;...
                        0., 0., BITZsandrrrr, 0., 0., 0.;...
                        0., 0., 0., 2*BITZsandrprp, 0., 0.;...
                        0., 0., 0., 0., 2*BITZsandrtrt, 0.;...
                        0., 0., 0., 0., 0., 2*BITZsandtptp];


end


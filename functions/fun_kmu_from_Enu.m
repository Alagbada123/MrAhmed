function [k,mu] = fun_kmu_from_Enu(E,nu)
%calculate bulk (k) and %shear modulus (mu) from Young's modulus (E) 
%and Poisson's ratio (nu)
k=E/3/(1-2*nu); 
mu=E/2/(1+nu);

end


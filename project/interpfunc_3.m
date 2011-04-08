function [Phi,dPhi_dxi,dPhi_deta]=interpfunc_3(xi,eta)

n=length(xi);

for i=1:n
    Phi(i,1) = 1-xi(i)-eta(i);
    Phi(i,3) = xi(i);
    Phi(i,2) = eta(i);
    
    dPhi_dxi(i,1) = -1;
    dPhi_dxi(i,3) =  1;
    dPhi_dxi(i,2) =  0;
    
    dPhi_deta(i,1) = -1;
    dPhi_deta(i,3) =  0;
    dPhi_deta(i,2) =  1;
end
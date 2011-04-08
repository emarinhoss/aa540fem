function [Phi,dPhi_dxi,dPhi_deta]=interpfunc_4(xi,eta)

n=length(xi);

for i=1:n
    Phi(i,1) = 0.25*(1-xi(i))*(1-eta(i));
    Phi(i,2) = 0.25*(1+xi(i))*(1-eta(i));
    Phi(i,3) = 0.25*(1+xi(i))*(1+eta(i));
    Phi(i,4) = 0.25*(1-xi(i))*(1+eta(i));
    
    dPhi_dxi(i,1) = 0.25*(1-eta(i));
    dPhi_dxi(i,2) = 0.25*(1-eta(i));
    dPhi_dxi(i,3) = 0.25*(1+eta(i));
    dPhi_dxi(i,4) = 0.25*(1+eta(i));
    
    dPhi_deta(i,1) = 0.25*(1-xi(i));
    dPhi_deta(i,2) = 0.25*(1+xi(i));
    dPhi_deta(i,3) = 0.25*(1+xi(i));
    dPhi_deta(i,4) = 0.25*(1-xi(i));
end
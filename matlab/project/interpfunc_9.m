function [Phi,dPhi_dxi,dPhi_deta]=interpfunc_9(xi,eta)

n=length(xi);

for i=1:n
    Phi(i,1) = 0.25*(xi(i)^2-xi(i))*(eta(i)^2-eta(i));
    Phi(i,2) = 0.50*(1-xi(i)^2)*(eta(i)^2-eta(i));
    Phi(i,3) = 0.25*(xi(i)^2+xi(i))*(eta(i)^2-eta(i));
    Phi(i,4) = 0.50*(xi(i)^2-xi(i))*(1-eta(i)^2);
    Phi(i,5) = (1-xi(i)^2)*(1-eta(i)^2);
    Phi(i,6) = 0.50*(xi(i)^2+xi(i))*(1-eta(i)^2);
    Phi(i,7) = 0.25*(xi(i)^2-xi(i))*(eta(i)^2+eta(i));
    Phi(i,8) = 0.50*(1-xi(i)^2)*(eta(i)^2+eta(i));
    Phi(i,9) = 0.25*(xi(i)^2+xi(i))*(eta(i)^2+eta(i));
    
    dPhi_dxi(i,1) = 0.25*(xi(i)^2-xi(i))*(2*eta(i)-1);
    dPhi_dxi(i,2) = 0.50*(-2*xi(i))*(eta(i)^2-eta(i));
    dPhi_dxi(i,3) = 0.25*(2*xi(i)+1)*(eta(i)^2-eta(i));
    dPhi_dxi(i,4) = 0.50*(2*xi(i)-1)*(1-eta(i)^2);
    dPhi_dxi(i,5) = (-2*xi(i))*(1-eta(i)^2);
    dPhi_dxi(i,6) = 0.50*(2*xi(i)+1)*(1-eta(i)^2);
    dPhi_dxi(i,7) = 0.25*(2*xi(i)-1)*(eta(i)^2+eta(i));
    dPhi_dxi(i,8) = 0.50*(-2*xi(i))*(eta(i)^2+eta(i));
    dPhi_dxi(i,9) = 0.25*(2*xi(i)+1)*(eta(i)^2+eta(i));
    
    dPhi_deta(i,1) = 0.25*(xi(i)^2-xi(i))*(2*eta(i)-1);
    dPhi_deta(i,2) = 0.50*(1-xi(i)^2)*(2*eta(i)-1);
    dPhi_deta(i,3) = 0.25*(xi(i)^2+xi(i))*(2*eta(i)-1);
    dPhi_deta(i,4) = 0.50*(xi(i)^2-xi(i))*(-2*eta(i));
    dPhi_deta(i,5) = (1-xi(i)^2)*(-2*eta(i));
    dPhi_deta(i,6) = 0.50*(xi(i)^2+xi(i))*(-2*eta(i));
    dPhi_deta(i,7) = 0.25*(xi(i)^2-xi(i))*(2*eta(i)+1);
    dPhi_deta(i,8) = 0.50*(1-xi(i)^2)*(2*eta(i)+1);
    dPhi_deta(i,9) = 0.25*(xi(i)^2+xi(i))*(2*eta(i)+1);
end
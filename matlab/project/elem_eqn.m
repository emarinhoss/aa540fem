function [Ke,fe]=elem_eqn(x,y,Phi,dPhi_dxi,dPhi_deta,w)

n = length(x);

Ke=zeros(n,n); fe=zeros(1,n);
% add information about interpolation functions
% and natural coordinates
X=0; Y=0;
dx_dxi = 0;
dx_deta= 0;
dy_dxi = 0;
dy_deta= 0;
for i=1:n
    X = X + x(i)*Phi(i);
    Y = Y + y(i)*Phi(i);
    
    dx_dxi = dx_dxi + x(i)*dPhi_dxi(i);
    dx_deta= dx_deta+ x(i)*dPhi_deta(i);
    dy_dxi = dy_dxi + y(i)*dPhi_dxi(i);
    dy_deta= dy_deta+ y(i)*dPhi_deta(i);
end

% determinant of the jacobian
hs = dx_dxi*dy_deta-dx_deta*dy_dxi;

% add information about anisotropic conductivity
% tensor and forcing vector
[k11,k12,k21,k22,f]=conductivity_and_forcing(X,Y);

dxi_dx =  1/hs*dy_deta; dxi_dy = -1/hs*dx_deta;
deta_dx= -1/hs*dy_dxi;  deta_dy=  1/hs*dx_dxi;

for i=1:n
    for j=1:n
        F1 = k11*(dxi_dx*dPhi_dxi(i)+deta_dx*dPhi_deta(i))...
                *(dxi_dx*dPhi_dxi(j)+deta_dx*dPhi_deta(j));
            
        F2 = k12*(dxi_dx*dPhi_dxi(i)+deta_dx*dPhi_deta(i)) ...
                *(dxi_dy*dPhi_dxi(j)+deta_dy*dPhi_deta(j));
            
        F3 = k21*(dxi_dy*dPhi_dxi(i)+deta_dy*dPhi_deta(i)) ...
                *(dxi_dx*dPhi_dxi(j)+deta_dx*dPhi_deta(j));
            
        F4 = k22*(dxi_dy*dPhi_dxi(i)+deta_dy*dPhi_deta(i)) ...
                *(dxi_dy*dPhi_dxi(j)+deta_dy*dPhi_deta(j));
            
        Ke(i,j) = w*(F1+F2+F3+F4)*hs;
        
    end
    
    fe(i) = w*f*Phi(i)*hs;
    
end
function [K,f]=dirichlet(K,f,bc,val)

n = length(bc);
NG= length(f);

for k=1:n
    m=bc(k);
    for i=1:NG
        f(i) = f(i) - K(i,m)*val;
        K(m,i) = 0;
        K(i,m) = 0;
    end
    K(m,m) = 1.0;
    f(m) = val;
end
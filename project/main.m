%%%%%%%
%
%
%
%
%
%%%%%%%

clear all; clc

%% User inputs
a = 4;
b = 6;

%%%%% element along x and y must be equal %%%%%
elems = 50;

%%%%% Element type %%%%
elem_type = 3;
% 1 --> 3-node linear
% 2 --> 4-node linear
% 3 --> 9-node quadratic

%%%%% Boundary Condition types %%%%
% 0 --> Dirichlet BC's
% 1 --> Neumann BC's
top_bc_type   = 0;
left_bc_type  = 1;
right_bc_type = 1;
bottom_bc_type= 0;

%%%%%% Values at the boundaries %%%%%
top_bc_val   = 100;
left_bc_val  = 0;
right_bc_val = 0;
bottom_bc_val= 0;

%%%%% Order of Gauss-Legendre quadrature %%%%
% values from 1-4, any number outside that range will default to 1
order=1;
%%%%%%%%%%%%%%%%%%

%% Select Geometry depending on element types to be used
switch elem_type
    case 1
        %number of interpolation functions
        N = 3;
        
        %
        [bcm,x,y,X,Y,bc_nodes] = geometry(a,b,elems,1);
        
        %number of nodes along the x and y directions
        nodes = elems/2+1;
        
        % get values to perform quadrature integrations
        % for triangular element only using first order.
        [xi,wi] = gauss_legendre_quad(1); eta=xi;
        
        % Interpolation and geometry calculations
        [Phi,dPhi_dxi,dPhi_deta]=interpfunc_3(xi,eta);
        
        % size of connectivity matrix
        last = max(size(bcm));
        
    case 2
        %number of interpolation functions
        N = 4;
        
        %get geometry information
        [bcm,x,y,X,Y,bc_nodes] = geometry(a,b,elems,2);
        
        %number of nodes along the x and y directions
        nodes = elems+1;
        
        % get values of xi and eta to perform quadrature integrations
        % by specifying the order of accuracy
        [xi,wi] = gauss_legendre_quad(order); eta=xi;
        
        % Interpolation and geometry calculations
        [Phi,dPhi_dxi,dPhi_deta]=interpfunc_4(xi,eta);
        
        % size of connectivity matrix
        last = max(size(bcm));
        
    case 3
        %number of interpolation functions
        N = 9;
        
        %get geometry information
        [bcm,x,y,X,Y,bc_nodes] = geometry(a,b,elems,3);
        
        %number of nodes along the x and y directions
        nodes = 2*elems+1;
        
        % get values of xi and eta to perform quadrature integrations
        % by specifying the order of accuracy
        [xi,wi] = gauss_legendre_quad(order); eta=xi;
        
        % Interpolation and geometry calculations
        [Phi,dPhi_dxi,dPhi_deta]=interpfunc_9(xi,eta);
        
        % size of connectivity matrix
        last = max(size(bcm));
    otherwise
        disp('Number entered corresponds to no element type!')
        return
end

disp('Finished generating Grid.')

%% Build Global Stiffness matrix

% alocate memory for global stiffness matrix
K = zeros(nodes*nodes,nodes*nodes);
% alocate memory for global forcing vector
F = zeros(nodes*nodes,1);

xcoords_e = zeros(1,N);
ycoords_e = zeros(1,N);
for e=1:last
    for k=1:N
        xcoords_e(k)=x(bcm(e,k));
        ycoords_e(k)=y(bcm(e,k));
    end
    
    if (length(xi) > 1)
        Ke = zeros(N,N);
        for j=1:length(xi)
            [Kxi,fe]=elem_eqn(xcoords_e, ycoords_e ...
                             ,Phi(j,:),dPhi_dxi(j,:),dPhi_deta(j,:),wi(j));
                    
            Ke = Ke + Kxi;
        end
    else
        [Ke,fe]=elem_eqn(xcoords_e, ycoords_e ...
                        ,Phi,dPhi_dxi,dPhi_deta,wi);
    end
                    
    for k=1:N
        m = bcm(e,k);
        for l=1:N
            n = bcm(e,l);
            K(m,n) = K(m,n) + Ke(k,l);
        end
        F(m) = F(m) + fe(k);
    end
end

disp('Finished Assembling Global Stifness Matrix.')
%% Apply boundary conditions
bc_type=[top_bc_type right_bc_type left_bc_type bottom_bc_type];
bc_val =[top_bc_val right_bc_val left_bc_val bottom_bc_val]; 
for j=1:4
     if bc_type(j)==0
         [K,F] = dirichlet(K,F,bc_nodes(j,:),bc_val(j));
%      else
%          
     end
end

% ss=reshape(F,nodes,nodes);
% surf(ss)
% figure(2), spy(K)
% break

disp('Boundary Conditions Applied.')

%% Solve for nodal Temperatures
disp('Solving Equations ...')
T = K\F;
%% Plot the solution
sol=reshape(T,nodes,nodes);
contour(X,Y,sol,20), axis equal
colorbar
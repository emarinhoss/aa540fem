function [b,x,y,X,Y,bcs]=geometry(a,b,elems,type)

%%% inputs
%
% a - domain size in the x-direction
% b - domain size in the y-direction
% elems - number of nodes in one row
%         in the x and y directions
%
%%% returns
%
% b - bolean connectivity matrix
% x - nodal x coordinates
% y - nodal y coordinates

switch type
    case 1
        nodes = elems/2+1;
        xx=linspace(0,a,nodes);
        yy=linspace(0,b,nodes);
        [X,Y]=meshgrid(xx,yy);
        x = reshape(X,1,nodes*nodes);
        y = reshape(Y,1,nodes*nodes);
        
        %% 3-node elements
        %b = zeros(elems*elems/4,3);
        
        b(1,1:3) = [1 2 (1+nodes)];
        b(2,1:3) = [(nodes+2) (nodes+1) 2];
        
        for i=3:elems
            b(i,:) = b(i-2,:) + 1;
        end
        
        for i=1:elems-1
            b((elems*i+1):(elems*(i+1)),:) = b((elems*(i-1)+1):(elems*i),:) + nodes;
        end
        
        bcs = zeros(4,nodes);
        % top bc nodes
        bcs(1,:) = nodes:nodes:length(x);
        % right bc nodes
        bcs(2,:) = length(x)-nodes+1:1:length(x);
        % left bc nodes
        bcs(3,:) = 1:1:nodes;
        % bottom bc nodes
        bcs(4,:) = 1:nodes:length(x)-nodes+1;
        
    case 2
        nodes = elems+1;
        xx=linspace(0,a,nodes);
        yy=linspace(0,b,nodes);
        [X,Y]=meshgrid(xx,yy);
        x = reshape(X,1,nodes*nodes);
        y = reshape(Y,1,nodes*nodes);
        
        %% 4-node elements
        
        b(1,1:4) = [1 (nodes+1) (nodes+2) 2];

        
        for i=2:elems
            b(i,:) = b(i-1,:) + 1;
        end
        
        for i=1:elems-1
            b((elems*i+1):(elems*(i+1)),:) = b((elems*(i-1)+1):(elems*i),:) + nodes;
        end
        
        bcs = zeros(4,nodes);
        % top bc nodes
        bcs(1,:) = nodes:nodes:length(x);
        % right bc nodes
        bcs(2,:) = length(x)-nodes+1:1:length(x);
        % left bc nodes
        bcs(3,:) = 1:1:nodes;
        % bottom bc nodes
        bcs(4,:) = 1:nodes:length(x)-nodes+1;
        
    case 3
        nodes = 2*elems+1;
        xx=linspace(0,a,nodes);
        yy=linspace(0,b,nodes);
        [X,Y]=meshgrid(xx,yy);
        x = reshape(X,1,nodes*nodes);
        y = reshape(Y,1,nodes*nodes);
        
        %% 9-node elements
        
        b(1,1:9) = [1 2 3 202 203 204 403 404 405];

        
        for i=2:elems
            b(i,:) = b(i-1,:) + 1;
        end
        
        for i=1:elems-1
            b((elems*i+1):(elems*(i+1)),:) = b((elems*(i-1)+1):(elems*i),:) + 2*nodes;
        end
        
        bcs = zeros(4,nodes);
        % top bc nodes
        bcs(1,:) = (nodes-1)*nodes+1:1:nodes*nodes;
        % right bc nodes
        bcs(2,:) = nodes:nodes:nodes*nodes;
        % left bc nodes
        bcs(3,:) = 1:nodes:(nodes-1)*nodes+1;
        % bottom bc nodes
        bcs(4,:) = 1:1:nodes;
end
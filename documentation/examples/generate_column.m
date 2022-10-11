%% Create a structured mesh for the bi-clamped beam
% X = nodal coordinates
% T = element connectivity and data
% row,col = for displaying design using imagesc (as a matrix)
function [X,T,row,col,solids,voids,F,freedofs] = generate_column(sizex,sizey,helem,doplot)
%% Define grid parameters 
dx = helem;             % element size in x-direction
dy = helem;             % element size in y-direction
l_load = dx*4;            % distribution length for point load 
%% Create nodal grid for FEA and optimization
[Xnode,Ynode] = meshgrid(0:dx:sizex,0:dy:sizey);
nodelist(:,1) = Xnode(:);
nodelist(:,2) = Ynode(:);
xnodeout = [];
ynodeout = [];
xynodeout = intersect(xnodeout,ynodeout);
nodelist_clean = nodelist;
nodelist_clean(xynodeout,:) = [];
X = nodelist_clean; % For output
% Plot nodes as points
if (doplot)
    figure(1);
    plot(nodelist_clean(:,1),nodelist_clean(:,2),'o');
    hold on;
    axis equal
    axis tight
end
%% Create element grid for FEA and optimization
[Xelem,Yelem] = meshgrid(dx/2:dx:sizex-dx/2,dy/2:dy:sizey-dy/2);
elemlist(:,1) = Xelem(:);
elemlist(:,2) = Yelem(:);
xelemout = [];
yelemout = []; 
xyelemout = intersect(xelemout,yelemout);
elemlist_clean = elemlist;
elemlist_clean(xyelemout,:) = [];
% Create element connectivity
nelem = size(elemlist_clean,1);
T = zeros(nelem,12);
% T(1:4) = [node1 node2 node3 node4]
% T(5) = 0/1/2 (design / void / solid)
% T(6) = 0/1 (domain / padding region)
% T(7:8) = [x_centroid y_centroid]
% T(9:12) = [left bottom right top] BCs for PDE filter
% This is a rectangular mesh - easy to create based on X and Y
for ex = 1:size(Xelem,2)
    for ey = 1:size(Yelem,1)
        e = (ex-1)*size(Xelem,1)+ey;
        n1 = (ex-1)*size(Xnode,1)+ey;
        n2 = n1+size(Xnode,1);
        n3 = n2+1;
        n4 = n1+1;
        T(e,1:4) = [n1 n2 n3 n4];
    end
end
for e = 1:nelem
    % Centroid coordinates
    x_cent = elemlist_clean(e,1);
    y_cent = elemlist_clean(e,2);
    T(e,7:8) = [x_cent y_cent];
    % Assign voids
    % No voids in the current representation of the column
    % Assign solids (for loaded and supported regions)
%     if (x_cent>sizex-helem && y_cent>sizey/2-l_load/2 && y_cent<sizey/2+l_load/2); T(e,5:6) = [2 0]; end
% 	if (x_cent<helem && y_cent>sizey/2-l_load/2 && y_cent<sizey/2+l_load/2); T(e,5:6) = [2 0]; end
    % Assign non-Neumann borders (for augmented PDE filter)
    % 0 = Neumann (supports, load)
    % 1 = Robin BCs with l_s=l_o (free faces)
    % 2 = xi_corner (Dirichlet-like)
    % Default is Neumann [left bottom right top]
    if (x_cent<helem); T(e,9) = 1; end % left face
    if (x_cent<helem && y_cent>sizey/2-l_load/2 && y_cent<sizey/2+l_load/2); T(e,9) = 0; end % left face, near support
    if (y_cent<helem); T(e,10) = 1; end % bottom face
    if (x_cent>sizex-helem); T(e,11) = 1; end % right face
    if (x_cent>sizex-helem && y_cent>sizey/2-l_load/2 && y_cent<sizey/2+l_load/2); T(e,11) = 0; end % right face, near load
    if (y_cent>sizey-helem); T(e,12) = 1; end % top face  
end
% Plot elements
if (doplot)
    figure(2);
    hold on;
    axis equal
    axis tight
    for e = 1:nelem
        if (T(e,5) == 0) % Regular element
            plot(T(e,7),T(e,8),'ro'); 
        end
        if (T(e,5) == 2) % Solid element
            plot(T(e,7),T(e,8),'k*'); 
        end
        if (T(e,5) == 1) % Void element
            plot(T(e,7),T(e,8),'m+'); 
        end
    end
end
% Plot boundary elements for filter BCs
if (doplot)
    figure(3);
    hold on;
    axis equal
    axis tight
    for e = 1:nelem
        if (T(e,9) == 0) 
            plot(T(e,7),T(e,8),'m+'); 
        end
        if (T(e,9) == 1) 
            plot(T(e,7),T(e,8),'ro'); 
        end
        if (T(e,9) == 2) 
            plot(T(e,7),T(e,8),'k*'); 
        end
        if (T(e,10) == 0) 
            plot(T(e,7),T(e,8),'m+'); 
        end
        if (T(e,10) == 1) 
            plot(T(e,7),T(e,8),'ro'); 
        end
        if (T(e,10) == 2) 
            plot(T(e,7),T(e,8),'k*'); 
        end
        if (T(e,11) == 0) 
            plot(T(e,7),T(e,8),'m+'); 
        end
        if (T(e,11) == 1) 
            plot(T(e,7),T(e,8),'ro'); 
        end
        if (T(e,11) == 2) 
            plot(T(e,7),T(e,8),'k*'); 
        end
        if (T(e,12) == 0) 
            plot(T(e,7),T(e,8),'m+'); 
        end
        if (T(e,12) == 1) 
            plot(T(e,7),T(e,8),'ro'); 
        end
        if (T(e,12) == 2) 
            plot(T(e,7),T(e,8),'k*'); 
        end
    end
end
%% Create matrix representation of topology
xmin = dx/2; 
ymin = dy/2; ymax = sizey-dy/2;
nrows = (ymax-ymin)/dy + 1;
col = zeros(nelem,1);
row = col;
for e = 1:nelem
    col(e) = (T(e,7)-xmin)/dx + 1;
    row(e) = nrows - ((T(e,8)-ymin)/dy);
end
%% Define loads and supports 
solids = find(T(:,5)==2);   % find solid elements
voids = find(T(:,5)==1);   % find void elements
% loadedelem = solids;
% nloadedelem = size(loadedelem,1);
loadnode1 = find(X(:,1)==sizex); % Find nodes with x=sizex
loadnode2 = find(X(:,2)==sizey/2); % Find nodes with y=sizey/2
loadnode = intersect(loadnode1,loadnode2);
loadeddof = 2*loadnode-1;
loadmag = -1e4;
F = sparse(loadeddof,1,loadmag,2*size(nodelist_clean,1),1);
% Supports
[supnodes1,~] = find(X(:,1)==0); % Find nodes with x=0
[supnodes2,~] = find(X(:,2)==sizey/2); % Find nodes with y=sizey/2
supnodes = intersect(supnodes1,supnodes2);
supdofs1 = [2*supnodes-1,2*supnodes];
[supnodes1,~] = find(X(:,1)==sizex); % Find nodes with x=sizex
[supnodes2,~] = find(X(:,2)==sizey/2); % Find nodes with y=sizey/2
supnodes = intersect(supnodes1,supnodes2);
supdofs2 = 2*supnodes;
supdofs = union(supdofs1,supdofs2);
alldofs = 1:2*size(nodelist_clean,1);
freedofs = setdiff(alldofs,supdofs)';
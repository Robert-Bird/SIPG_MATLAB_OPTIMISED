clear;
nD = 2;                                                                    % number of dimensions
v = 0.2;                                                                   % Poisson's ratio
E = 1e3;                                                                   % Youngs Modulus
pr = 100;                                                                  % Pressure

nen  = 4;                                                                  % Number of element nodes
nel_half = 29;
nel=nel_half*2;                                                            % Number of elements along x or y axis
nels = nel^2;                                                              % Number of elements
nelblo_size = 30000;                                                       % Number of elements in a block loop
num_faces = 4;                                                             % Number of element faces

if nen == 4
    [coord,etpl_face,etpl,ed] = DGmeshGenPlateHole(nel,nel);               % Hole 4 noded quad mesh
elseif nen == 8
    [coord,etpl_face,etpl,ed] = DGmeshGenPlateHole_8node(nel,nel);         % Hole 8 noded quad mesh
end

k = Linear2D(nen,nels,num_faces,nD,etpl,etpl_face,ed,coord,E,v,nelblo_size); % Optimised stiffness matrix generation
k=k+k';

% Example Boundary conditions and problem: hole in an infinite plate
[f,BC,d,c]=F_BC(nel,nel,etpl,ed,coord,pr,nen,nD);                          % Boundary conditions
k(BC,:)=[]; k(:,BC)=[]; f(BC)=[]; c(BC)=1;                                 % Applying boundary conditions
d(c~=1)=k\f;                                                               % Solving the linear system

stressplotter(E,v,coord,etpl,nD,nen,d,nels)                                % Stress plotter
if nels<500
    displacement_plotter(d,etpl,coord,nen,nel,nel)                         % Displacement plotter
end









% Author: Dimitrios Nentidis, Mechanical Engineer in the making (4th year)
% 6th of May 2023

% the following script comes with some limitations, although with a few
% alterations it can be expanded to work for more cases
% limitations:
% - assuming a rectangular cross section, to alter this directly input Iz,Iy
% - assuming a uniform cross section in the entirety of the beam length
% - equal spacing of the nodes, this cannot be altered


% Define material properties
answer1 = inputdlg({"Young's modulus (Pa)", "Poisson's ratio", 'density (kg/m^3)'},...
    'Material Properties', [1 51],{'2.1e11', '0.3', '7800'}); 
% material properties
E   = str2double(answer1{1});
nu  = str2double(answer1{2});
rho = str2double(answer1{3});
G = E/(2*(1+nu));  % Torsional Rigidity (Pa)

% Define Geometry
answer2 = inputdlg({"Beam Length (meters)", "Cross Section Width (meters)", 'Cross Section Height (meters'},...
    'Geometrical properties', [1 51],{'0.5', '0.003', '0.005'}); 
% geometry
Length = str2double(answer2{1});
width  = str2double(answer2{2});
height = str2double(answer2{3});

% cross-sectional properties 
A = height*width;           % area (m^2)
Iy = height*(width^3)/12;   % Iy area moment of inertia (m^4)
Iz = (height^3)*width/12;   % Iz area moment of inertia (m^4)
J = Iz + Iy;                % J polar moment of inertia (m^4)


answer3 = inputdlg({'Number of Nodes'}, ...
    'How many nodes', [1 51], {'8'});
nodes = str2double(answer3{1});
nel = nodes - 1;
L = Length/nel; % element length, assuming equal spacing

%% Kstiff for beam with 6dofs per node, 2 nodes


Kel= [E*A/L       0        0        0        0        0        -E*A/L        0        0        0        0        0;
        0    12*E*Iz/L^3   0        0        0   6*E*Iz/L^2       0   -12*E*Iy/L^3    0        0        0   6*E*Iz/L^2;
        0         0  12*E*Iy/L^3    0  -6*E*Iy/L^2    0           0          0  -12*E*Iy/L^3   0   -6*E*Iy/L^2   0;
        0         0        0      G*J/L      0        0           0          0        0      -G*J/L     0        0;
        0         0  -6*E*Iy/L^2    0     4*E*Iy/L    0           0          0    6*E*Iy/L^2   0     2*E*Iy/L    0;
        0     6*E*Iz/L^2   0        0        0    4*E*Iz/L        0    -6*E*Iz/L^2    0        0        0      2*E*Iz/L;
     -E*A/L       0        0        0        0        0         E*A/L        0        0        0        0        0;
        0   -12*E*Iz/L^3   0        0        0  -6*E*Iz/L^2       0    12*E*Iy/L^3    0        0        0  -6*E*Iz/L^2;
        0         0 -12*E*Iy/L^3    0   6*E*Iy/L^2    0           0          0   12*E*Iy/L^3   0    6*E*Iy/L^2   0;
        0         0        0     -G*J/L      0        0           0          0        0       G*J/L     0        0;
        0         0  -6*E*Iy/L^2    0     2*E*Iy/L    0           0          0    6*E*Iy/L^2   0     4*E*Iy/L    0;
        0     6*E*Iz/L^2   0        0        0    2*E*Iz/L        0     -6*E*Iz/L^2   0        0        0     4*E*Iz/L];

%% for nodes

beam_dofs = 6*nodes;                       % 6 dofs per node, 3 translational, 3 rotational
Kel_b = zeros(beam_dofs, beam_dofs, nel);  % this allows each element to position in the beam

for el=1:7
   for i=1:12
       for j=1:12
           n = (el-1)*6 + i;
           m = (el-1)*6 + j;
           Kel_b(n, m, el) = Kel(i, j);   % the actual position inside the beam
       end
   end
end
Kpb = zeros(beam_dofs, beam_dofs);        % initializing the beam Kstiff
for el=1:7
    Kpb = Kpb + Kel_b(:,:, el);           % adding every positioned element
end

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

%% mass matrix for 6 dofs per node, 2 node element
%M = zeros(150, 150);

% rod
Mrod = [rho*A*L/3  rho*A*L/6;   
        rho*A*L/6  rho*A*L/3]; 
% torsion
Mtrs = [rho*J*L/3  rho*J*L/6;  % [th_x_1-th_x_1 th_x_2-th_x_1;
        rho*J*L/6  rho*J*L/3]; %  th_x_1-th_x_2 th_x_2-th_x_2];

% bending X-Z or X-Y plane
Mbnd = [156*rho*A*L/420      22*L*rho*A*L/420     54*rho*A*L/420    -13*L*rho*A*L/420;           
        22*L*rho*A*L/420     4*L^2*rho*A*L/420    13*L*rho*A*L/420  -3*L^2*rho*A*L/420;           
        54*rho*A*L/420       13*L*rho*A*L/420     156*rho*A*L/420   -22*L*rho*A*L/420;           
        -13*L*rho*A*L/420   -3*L^2*rho*A*L/420   -22*L*rho*A*L/420   4*L^2*rho*A*L/420];



Mel= [Mrod(1,1)   0         0         0         0         0       Mrod(1,2)     0         0         0         0         0;
        0      Mbnd(1,1)    0         0         0      Mbnd(1,2)     0       Mbnd(1,3)    0         0         0      Mbnd(1,4);
        0         0      Mbnd(1,1)    0      Mbnd(1,2)    0          0          0      Mbnd(1,3)    0      Mbnd(1, 4)   0;
        0         0        0       Mtrs(1,1)    0         0          0          0         0      Mtrs(1,2)    0         0;
        0         0     Mbnd(2,1)     0      Mbnd(2,2)    0          0          0      Mbnd(2,3)    0      Mbnd(2,4)    0;
        0      Mbnd(2,1)   0          0         0      Mbnd(2,2)     0       Mbnd(2,3)    0         0         0      Mbnd(2,4);
      Mrod(2,1)   0        0          0         0         0       Mrod(2,2)     0         0         0         0         0;
        0      Mbnd(3,1)   0          0         0      Mbnd(3,2)     0       Mbnd(3,3)    0         0         0      Mbnd(3,4);
        0         0     Mbnd(3,1)     0      Mbnd(3,2)    0          0          0      Mbnd(3,3)    0      Mbnd(3,4)    0;
        0         0        0       Mtrs(2,1)    0         0          0          0         0      Mtrs(2,2)    0         0;
        0         0     Mbnd(4,1)     0      Mbnd(4,2)    0          0          0      Mbnd(4,3)    0      Mbnd(4,4)    0;
        0      Mbnd(4,1)   0          0         0      Mbnd(4,2)     0       Mbnd(4,3)    0         0         0      Mbnd(4,4)];
  
%% for 8 node

Mel_b = zeros(6*nodes, 6*nodes, nodes-1);
for el=1:nodes-1
   for i=1:12
       for j=1:12
           n = (el-1)*6 + i;
           m = (el-1)*6 + j;
           Mel_b(n, m, el) = Mel(i, j);
       end
   end
end
Mpb = zeros(6*nodes, 6*nodes);
for el=1:nodes-1
    Mpb = Mpb + Mel_b(:,:, el);
end

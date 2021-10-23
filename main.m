% Displacement controlled Newton-Raphson algorithm
% This code calculates the nonlinear force-displacement curve (assuming material nonlinearity, geometric linearity) 
% using displacement controlled Newton-Raphson algorithm.
% The input files are (nodes.xlsx, restraints.xlsx,
% elements_information.xlsx) defining the geometry of the system, the
% restraints of the nodes and the properties of the elements. The outputs
% are .txt files for the displacement, the corresponding force and the residual force for all iterations.
%Author: Nicole Widmer
%Last update: 23/10/2021


close all
clear all
clc

% Nodes
nodes = readmatrix("nodes.txt"); % m       %define coordinates of every node

% Restraints
rests = readmatrix("restraints.txt");   %restraints at nodes (1: DOF restraint, 0: DOF free)

% Material stiffness
E=200e9; % Pa
fy = 500e6; % Pa
ey = fy/E;
alpha = 0.02;   % post-yielding stiffness: alpha*E

% Element information
elems = readmatrix("elements_information.txt");  % for every element: Node1, Node2, area, pre-yielding stiffness

% Load pattern
Pn = [0;0;4615;7993;0;0]*10^3; %applied external force

% Define steps
steps = 3;          %number of loadsteps
iter_lim = 20;      %maximum number of iterations per loadstep
U_max = 0.150; %m   %final displacement
incr = 0.01; %m     %displacement increment applied per loadstep

n = size(elems, 1); % Number of elements
ndof = size(rests); % Degrees of freedom
u = zeros(size(Pn)); % Define zero initial displacements

%initialize variables
U2 = 0;
dU1 = 0;
incr_step = 0;

% Define convergence limit
conv_limit = norm(Pn)*10^(-6);

% Define inital load factor
lam = 0;
counter = 0;

% Construct bm that defines the dof to be controlled
bm = [0 1];

% Modify the stopping condition
while norm(U2)<norm(U_max-incr)
    incr_step = incr_step+1;

    for i = 1:iter_lim
        counter = counter+1

        % Define displacement increment
        if i == 1
            du = [dU1; -incr];
            R = [0;0];
        else
            du = [du(1); 0];
            R = R(:, i-1);
        end
        
        % Calculate the global and reduced stiffness matrix
        [K_g, K_r, T, L] = get_material_stiffness_matrix(nodes, elems, rests); 

        % Solve the equation -R = dl*P - K_r*du for the
        % unkown displacements (dU1 in this case) and labda increments
        dl = (-bm*inv(K_r)*R + bm*du)/(bm*inv(K_r)*Pn(find(~rests)));% Calculate the increment in load factor (lambda)
        du = inv(K_r)*(dl*Pn(find(~rests))+R); % Calculate the increment in displacements
        u(find(~rests)) = u(find(~rests)) + du; % Update displacements
        u_red = u(find(~rests));
        U2 = u_red(2);
        lam = lam + dl;  % Update load factors

        % Update the stiffness based on the strain limit
        for m = 1:n
            eps(m) = get_strain(elems, u, L, T, m);

            if abs(eps(m)) >= ey
                elems(m, 4) = alpha * E;
            else
                elems(m, 4) =  E;
            end 

        end

        % Calculate the external load corresponding to the displacement u
        P_red = lam * Pn(find(~rests)); 

        % Obtain the internal resistance forces
        Pr = get_nodal_forces(elems, u, L, T, fy, ey, n, ndof, rests);
        Pr_red = Pr(find(~rests));
        R(:, i) = P_red - Pr_red % Caluclate the residual
        monitor_R(counter, :) = R(:,i);

        % Check whether convergence is achieved
        if norm(R(:,i))<conv_limit
            dU1 = u_red(1);
            break;
        end
    end

    U_step(:,incr_step)=u;
    Pr_step(:,incr_step)= Pr;
    
end

%create output files
writematrix(U_step,'U')
writematrix(Pr_step,'F')
writematrix(monitor_R,'R')

% This file is used to create the potential files.
%clear all; close all;

%% MoS2 potential calculated by DFT

% Above the sulfur
z_S = [1.5,1.75,2,2.25,2.5,2.75,3,3.5,4,4.25,4.5,4.75,5,5.5,6];
V_S = [3.9991,1.894,0.8359,0.3437,0.1259,0.0339,-0.0028,-0.0157,-0.0105,-0.0081,-0.0063,-0.0046,-0.003,-0.0009,-0.0005]*1000;

% Above the molybdinum
z_Mo = [1.5,1.75,2,2.25,2.5,2.75,3,3.5,4,4.25,4.5,4.75,5,5.5,6];
V_Mo = [0.5259,0.2621,0.1146,0.0374,-0.0004,-0.0168,-0.0223,-0.0182,-0.0105,-0.0082,-0.0065,-0.0049,-0.0032,-0.001,-0.0007]*1000;

% Above the hollow site
z_H = [1.5,1.75,2,2.25,2.5,2.75,3,3.5,4,4.25,4.5,4.75,5,5.5,6];
V_H = [0.5255,0.2812,0.1325,0.0493,0.0066,-0.0132,-0.0206,-0.018,-0.0105,-0.0082,-0.0066,-0.0049,-0.0033,-0.001,-0.0008]*1000;

figure
hold on
plot(z_Mo, V_Mo, 'DisplayName', 'Mo')
plot(z_S, V_S, 'DisplayName', 'S')
plot(z_H, V_H, 'DisplayName', 'Hollow')
hold off
ylim([-25, 100])
legend()
xlabel('z/A')
ylabel('V/meV')
grid on
title('DFT potential')

%% Generate a potential

% the constants used for generating potential (Corrugated Morse potential)
D = 22.3; % Well depth
alpha = 1.1; % Well width
beta = 0.1;
gamma = 0.2;
%lattice constant
a = 2.84;

[V, V0, V1, param] = mos2_potential(D, alpha, beta, gamma, a);

V_tmp = V0(1,1,:);
figure
plot(param.z, V_tmp(:));
xlabel('z/A')
ylabel('V/meV')
ylim([-20, 100])
grid on
title('Uncorrugated part of the potential')

if true
    V_tmp = V1(1,1,:);
    figure
    plot(param.z, V_tmp(:))
    xlabel('z/A')
    ylabel('V_1/meV')
    ylim([-20,0])
    grid on
    title('Corrugated part of the potential')
end

%% Plot the potential
how_many_cell = param.how_many_cell;
gridp = param.gridp;

figure
contourf(z,  linspace(0, a*how_many_cell, gridp*how_many_cell), reshape(V(1,:,:), [gridp*how_many_cell,100]), linspace(-20,100,24))
xlabel('z/A')
ylabel('x/A')
colorbar
xlim([-1,4])
title('Potential in z, used in simulation')

if false
    figure
    contourf(z,  linspace(0, 4*a, 4*gridp), repmat(reshape(V(1,:,:), [gridp,100]),[4,1]), linspace(-20,100,24))
    xlabel('z/A')
    ylabel('x/A')
    colorbar
    xlim([-1,4])
    title('Potential in z, extended')
end

inds0 = zeros(gridp*how_many_cell,gridp*how_many_cell);

for ia1=1:length(param.i1)
    for ia2=1:length(param.i2)
        for indz=1:length(z)-1
            if V(ia1,ia2,indz) > 0 && V(ia1,ia2,indz+1) < 0
                m = (V(ia1,ia2,indz) - V(ia1, ia2, indz+1))/(indz - (indz+1));
                c = V(ia1, ia2,indz) - m*indz;
                inds0(ia1,ia2) = -c/m;
            end
        end
    end
end
pot_height = z(1) + inds0*(z(2)-z(1));

figure
surf(linspace(0, a*how_many_cell, gridp*how_many_cell), linspace(0, a*how_many_cell, gridp*how_many_cell), pot_height)
xlabel('x/A')
ylabel('y/A')
title('Equipotential V=0, used in simulation')

if false
    figure
    surf(linspace(0, 4*a, 4*gridp), linspace(0, 4*a, 4*gridp), repmat(pot_height, [4,4]))
    xlabel('x/A')
    ylabel('y/A')
    title('Equipotential V=0, extended')
end


%% Optimise parameters for the DFT calculation



%% MultiScat format
% now we have got the matrix V, it's time to create the files used by
% Multiscat.
potStructArray.V = V;
% this creates a struct called "potStructArray with a field V,
% and this field is defined by the matrix V in your MATLAB s Workspace.
Multiscat.PreparePotentialFiles(potStructArray);
% this function creates the pot10001.in file, when executing this line,
% make sure that Multiscat.m is in your current folder.
Multiscat.prepareFourierLabels(V);
% this function creates the FourierLabels.in file, which is also essential
% for Multiscat to get the information about potential.
potStructArray.a1=a1*how_many_cell; potStructArray.a2 = a2*how_many_cell;
% make sure you have the base vectors a1 and a2 in your Workspace.
% The orientation of a1 or a2 does not matter, you just need to get
%the length and angle between them right. As an example, for
% LiF (001) surface, a1=[2.84, 0] and a2=[0, 2.84].
% Note that the length unit is angstrom.
potStructArray.zmin = z(1);
potStructArray.zmax = z(end);
potStructArray.zPoints = length(z);
% where we have an array called z that stores the values of z being
% considered. So zmin and zmax are just the first and the last element
% of z, respectively. The length of z is the number of z points used.
confStruct = Multiscat.createConfigStruct(potStructArray);
Multiscat.prepareConfigFile(confStruct);
%the above two lines create Multiscat.conf

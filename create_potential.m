% This file is used to create the potential files.
%clear all; close all;

%% Generate a potential

% z points, note that the points we take should be evenly spaced
z = linspace(-2,6,100);
%lattice constant
a = 2.84;
% base vectors in real space
a1=[a, 0];
a2=[0, a];
% reciprocal base lattice vectors
[b1, b2] = Reciprocal(a1, a2);
% grid points parallel to the surface in one dimension
gridp = 16;
% point grids in a1 and a2 directions
how_many_cell = 1;
i1 = 1:gridp*how_many_cell;
i2 = 1:gridp*how_many_cell;

% potential matrix
V = zeros(length(i1),length(i2),length(z));

% the constants used for generating potential (Corrugated Morse potential)
D = 7.63; % Well depth
alpha = 1.1; % Well width
z0 = 1.0; % z offset

% now assign values to V0 - the Morse potential without the corrugation
V0 = D*exp(2*alpha*(z0-z))-2*D*exp(alpha*(z0-z));
figure
plot(z, V0);
xlabel('z/A')
ylabel('V/meV')
ylim([-20, 100])
grid on
title('Uncorrugated part of the potential')
V0 = repmat(reshape(V0,1,1,[]),size(V,1),size(V,2),1);

% repeat matrix V0 to prepare to add it to V
V = V + V0;
for ia1=1:length(i1)
    for ia2=1:length(i2)
        X(ia1,ia2) = (a1(1)*ia1+a2(1)*ia2)./gridp; %#ok<SAGROW> %the x grid points
        Y(ia1,ia2) = (a1(2)*ia1+a2(2)*ia2)./gridp; %#ok<SAGROW> %the y grid points
    end
end

% the constant used to create V1: the ratio of the corrugated part of the
% potential to the uncorrugated
beta = 0.10;
V1 = -2*beta*D*exp(2*alpha*(z0-z));

if false
    figure
    plot(z, V1)
    xlabel('z/A')
    ylabel('V_1/meV')
    ylim([-20,0])
    grid on
    title('Corrugated part of the potential')
end

dimple = false;

% Add the extra corrugated part according to the corrugation Q
Q = zeros(size(X));
for ia1=1:length(i1)
    for ia2=1:length(i2)
        Q(ia1,ia2) = cos(2*pi*X(ia1,ia2)/a)+cos(2*pi*Y(ia1,ia2)/a);
        if dimple
            Q(ia1, ia2) = Q(ia1, ia2) + Gaussian2D(X(ia1,ia2), Y(ia1,ia2), [4.14, 4.14], 1, 8);
        end
        for indz=1:length(z)
            V(ia1,ia2,indz) = V(ia1,ia2,indz)+V1(indz)*Q(ia1,ia2);
            %add V1*Q to V
        end
    end
end


%% Plot the potential

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

for ia1=1:length(i1)
    for ia2=1:length(i2)
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
%this function calculates the reciprocal lattice vectors b1 and b2
%from the real space basis vectors a1 and a2

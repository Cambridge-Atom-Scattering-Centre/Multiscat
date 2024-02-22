function  equipotential_plot(V, a, z, pot_level)

if nargin == 3
    pot_level = 0;
end

how_many_cell = 2;
gridp = 64;

i1 = 1:gridp*how_many_cell;
i2 = 1:gridp*how_many_cell;

inds0 = zeros(gridp*how_many_cell,gridp*how_many_cell);

for ia1=1:length(i1)
    for ia2=1:length(i2)
        for indz=1:length(z)-1
            if V(ia1,ia2,indz) > pot_level && V(ia1,ia2,indz+1) < pot_level
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
zlabel('z/A')
title(['Equipotential V=' num2str(pot_level) 'meV, used in simulation'])

end
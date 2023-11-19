% Calculates the reciprocal lattice vectors of a surface given the real
% space lattice vectors.
function [b1, b2] = Reciprocal(a1, a2)
    a3 = [0 0 1];
    crsprod = a3*(cross([a1 0], [a2 0]))';
    % the volume spanned by a1, a2, and a3
    b1 = 2*pi*cross([a2 0],a3)/crsprod;
    b1 = b1(1:2);
    b2 = 2*pi*cross(a3,[a1 0])/crsprod;
    b2 = b2(1:2);
end
function [V, V0, V1, param] = mos2_potential(D, alpha, beta, gamma, a)
    % z points, note that the points we take should be evenly spaced
    zrange = [-2,6];
    z = linspace(zrange(1),zrange(2),100);

    % base vectors in real space
    a1=[a, 0];
    a2=[0, a];
    % reciprocal base lattice vectors
    [b1, b2] = Reciprocal(a1, a2);
    % grid points parallel to the surface in one dimension
    gridp = 32;
    % point grids in a1 and a2 directions
    how_many_cell = 2;
    i1 = 1:gridp*how_many_cell;
    i2 = 1:gridp*how_many_cell;

    % potential matrix
    V = zeros(length(i1),length(i2),length(z));
    
    % the constants used for generating potential (Corrugated Morse potential)
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
    V1 = -2*beta*D*exp(2*alpha*(z0-z));

    % Add the extra corrugated part according to the corrugation Q
    Q = zeros(size(X));
    Q2 = zeros(size(X));
    for ia1=1:length(i1)
        for ia2=1:length(i2)
            % First order Fourier component
            Q(ia1,ia2) = cos(2*pi*X(ia1,ia2)/a) + cos(2*pi*(X(ia1,ia2)*0.5 + Y(ia1,ia2)*sqrt(3)/2)/a) + cos(2*pi*(X(ia1,ia2)*0.5 - Y(ia1,ia2)*sqrt(3)/2)/a);
            % Second order Fourier component
            x0 = a*tand(30)/2;
            y0 = a/2;
            x = X(ia1,ia2) - x0;
            y = Y(ia1,ia2) - y0;
            Q2(ia1,ia2) = cos(2*pi*x/a) + cos(2*pi*(x*0.5 + y*sqrt(3)/2)/a) + cos(2*pi*(x*0.5 - y*sqrt(3)/2)/a);

            for indz=1:length(z)
                V(ia1,ia2,indz) = V(ia1,ia2,indz) - V1(indz)*Q(ia1,ia2) - gamma*V1(indz)*Q2(ia1,ia2);
                %add V1*Q to V
            end
        end
    end
    param.z = z;
    param.gridp = gridp;
    param.how_many_cell = how_many_cell;
    param.X = X;
    param.Y = Y;
    param.z0 = z0;
    param.i1 = i1;
    param.i2 = i2;
end
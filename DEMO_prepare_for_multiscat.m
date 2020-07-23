global data_file

a=4.52; % ice lattice constant(?)
HeMass = 3;
rotMat = [1 0;0 1]; % No rotation
fVbasis = @structureAnalysis.Vbasis_atoms;

[ki, Ei, Ti] = beamprops('energy',8,HeMass);
[a1, a2, b1, b2] = structureAnalysis.unitCellVecs(1,a,rotMat);
G = structureAnalysis.spanVecs(b1,b2,-8:8,-8:8);

a1UnitCells = 5;
a2UnitCells = 5;
pointsInUnitCell = 32;
zPoints = 550;
zmin = 0;
zmax = 20;
rBasisZmin=2;

if zPoints>550, disp('WARNING: zPoints<=500, if you need more points, increase NZFIXED_MAX in multiscat.inc'); end

%% Generate flavours of the Potential


% Create flavours for the basis (of the lattice unitcell)

rBasisCase=[1]; % see rBasis_Cls.rBasis_Func

% Depending on the rBasisCase, different fields are expected for
% the structure basisArgs. (See rBasis_Cls.rBasis_Func)

clear basisArgs
zref = [1.0];
for j=1:length(zref)
    basisArgs(j).zref=zref(j);
end

paramIndMat = CustomFuncs.allPerm(1:length(basisArgs));
if iscell(paramIndMat)
    tmp = paramIndMat;
    paramIndMat = cell2mat(tmp);
    clear tmp
end

%%

for i=1:size(paramIndMat,2)
    rBasis = rBasis_Cls.rBasis_Func(rBasisCase,a1,a2,rBasisZmin,basisArgs(paramIndMat(1,i)));
    [X,Y,Z,V] = structureAnalysis.constructRealSpace(a1UnitCells, a2UnitCells, a1, a2, pointsInUnitCell, ...
        zPoints, zmin, zmax, rBasis, fVbasis);
    potStructArray(i).V=V;
    potStructArray(i).a1UnitCells = a1UnitCells;
    potStructArray(i).a2UnitCells = a2UnitCells;
    potStructArray(i).pointsInUnitCell = pointsInUnitCell;
    potStructArray(i).zPoints = zPoints;
    potStructArray(i).zmin = zmin;
    potStructArray(i).zmax = zmax;
    potStructArray(i).basisArgs = basisArgs((paramIndMat(1,i)));
    potStructArray(i).rBasisCase = rBasisCase;
    potStructArray(i).rBasisZmin = rBasisZmin;
    potStructArray(i).rBasis = rBasis;
    potStructArray(i).fVbasis = fVbasis;
    potStructArray(i).a1=a1;
    potStructArray(i).a2=a2;
end


realSpaceStructArray.X=X; realSpaceStructArray.Y=Y; realSpaceStructArray.Z=Z;
realSpaceStructArray.a1=a1; realSpaceStructArray.a2=a2;

%% Work out scattering conditions

% Multiscat returns the intensity of a 2D G channels, for
% {alpha,phi,E_beam} scattering conditions, were alpha is
% the incident angle and phi is the azimuthatl angle. At
% the HeSE "language", these are {gamma,alpha,E_beam}, we
% use those in the current script. We now define the
% experimental equivalent conditions. For example, we
% could choose a fixed total scattering angle. That would
% result many executions of Multiscat.

thetaRange=[-85 85]; constTheta = [];
gammaRange=[]; constGamma = 0;
alphaRange=[359 31];

[gamma, theta, alpha, G2b_measured] = structureAnalysis.findScatteringGeometries(ki, gammaRange,alphaRange,thetaRange,G, constGamma, constTheta);

%hold on; plot(G2b_measured(:,1),G2b_measured(:,2),'o');
%axis([-10 10 -10 10])

clear calcStrct

calcStrct.thetaRange=thetaRange;
calcStrct.constTheta=constTheta;
calcStrct.gammaRange=gammaRange;
calcStrct.constGamma=constGamma;
calcStrct.alphaRange=alphaRange;
calcStrct.HeMass = HeMass;
calcStrct.a1=a1;
calcStrct.a2=a2;
calcStrct.ki=ki;
calcStrct.Ei=Ei;
calcStrct.G=G;
calcStrct.gamma=gamma;
calcStrct.theta=theta;
calcStrct.alpha=alpha;
calcStrct.G2b_measured=G2b_measured;

%% Prepare the files for Multiscat

define_datafile_mfile('mltsct',1);
if exist(data_file,'file')
    save(data_file,'potStructArray', 'realSpaceStructArray','calcStrct', '-append')
else
    save(data_file,'potStructArray', 'realSpaceStructArray','calcStrct')
end
load(data_file)
confStruct = Multiscat.createConfigStruct(data_file,potStructArray, calcStrct);
Multiscat.PreparePotentialFiles(potStructArray);
Multiscat.prepareFourierLabels(potStructArray(1).V);
Multiscat.prepareConfigFile(confStruct);



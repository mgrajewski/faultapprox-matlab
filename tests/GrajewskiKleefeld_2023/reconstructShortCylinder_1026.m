% Reconstruction of the short cylinder in "Detecting and approximating
% decision boundaries in low dimensional spaces", section 4.3

% Author: Matthias Grajewski (grajewski@fh-aachen.de) and Andreas Kleefeld
% (a.kleefeld@fz-juelich.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)

% txt-file containing the far field data
filename = 'farFieldData1026/farFieldShortCylinder.txt';
resultfile = 'raw_results/pointsShortCylinder.txt';
normalsfile = 'raw_results/normalsShortCylinder.txt';
statsfile = 'raw_results/shortCylinder1026.csv';

global ExtendedStats
ExtendedStats = Statistics;

MyProb = ProblemDescr;

% factorKirsch contains Kirsch's factorization algorithm
MyProb.Testfunc = @factorKirsch;

% No vtu-output. Instead, we write the resulting point set along with
% the normal vectors in txt-Files in order to triangulate this set using
% open3D in a separate python script.
MyProb.OutputFileVTU = '';

% Enable extended statistics in order to store intermediate results like
% number of factorization calls.
MyProb.extendedStats = true;
MyProb.verboseMode = true;

% define domain [-2, 2]^3
MyProb.Xmin = [-2 -2 -2];
MyProb.Xmax = [2 2 2];

% parameter settings for fault approximation
MySettings = FaultApproxParameters;
MySettings.NumPointsLocal = 7;

% epsilon_gap
MySettings.maxDistForSurfacePoints = 0.25;

% maximal number of sweeps in adaptive refinement
MySettings.maxiterAdapt = 15;

% max tolerable error during refinement
MySettings.errMax = 0.01;

% point accuracy in normal direction
MySettings.abstolBisection = 0.001;


% Compute the initial point set consisting of 200 Halton-distributed points
% on [0,1]^3.
PointSet = CreateHaltonSet(200, 3, 1);

% Map the Halton points to to [Xmin, Xmax].
for idim = 1: 3
    PointSet(:,idim) = (MyProb.Xmax(idim) - MyProb.Xmin(idim))*PointSet(:,idim) + MyProb.Xmin(idim);
end

% Prepare computation of inverse problems (taken from Andreas Kleefeld's
% sample script).
global M
global d
global V
global sigma
global wavenumber

fprintf('Reading the data...\n')
X=load(filename,'-ASCII');          % contains M incident directions, 
                                    % M measuring points,
                                    % and M x M far-field data (3D)

d=X(:,1:3);                         % incident directions
far=X(:,7)+1i*X(:,8);               % far-field data
M=sqrt(size(d,1));                  % size of far-field matrix (M x M)
fprintf('Matrix Size is %d x %d\n',M,M)

wavenumber=2;                       % wave number 
A=reshape(far,M,M);                 % far-field matrix 
fprintf('Creating SVD...\n')
[~,sigmam,V]=svd(A);                % SVD needed for the 
sigma=diag(sigmam);                 % factorization method

% compute approximation to object boundary
fprintf('compute approximation to object boundary\n')
Subdomains = faultApprox(PointSet, MyProb, MySettings);

% Save the points as csv-file.
points = Subdomains{1}{2}{2};
writematrix(points, resultfile);

% Save the normals as csv-file.
normals = Subdomains{1}{2}{3};
writematrix(normals, normalsfile);

% save the statistical results
line = [ 'ncalls_', ExtendedStats.pos_in_code{2}, ', ', int2str(ExtendedStats.ncalls{2})];
writelines(line, statsfile, WriteMode='overwrite')
for i =3:size(ExtendedStats.pos_in_code,2)
    line = [ 'ncalls_', ExtendedStats.pos_in_code{i}, ', ', int2str(ExtendedStats.ncalls{i} - ExtendedStats.ncalls{i-1})];
    writelines(line, statsfile, WriteMode='append')
end

for i =3:size(ExtendedStats.pos_in_code,2)
    line = [ 'S12_', ExtendedStats.pos_in_code{i}, ', ', int2str(ExtendedStats.nPointsSurf{i-2}{1,2})];
    writelines(line, statsfile, WriteMode='append')
end
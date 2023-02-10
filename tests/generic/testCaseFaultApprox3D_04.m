% 3D test case 04

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.05;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 7;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'results/testFaultApprox3D_04.vtu';
MyProb.Testfunc = @testFuncFaultApprox3D_04;
MyProb.Xmin = [0 0 0] ;
MyProb.Xmax = [1 1 1];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;

npointsperSide = 10;

PointSet = 0:1/(npointsperSide-1):1;
[X, Y, Z] = ndgrid(PointSet);
X = reshape(X, [1 (npointsperSide)^3]);
Y = reshape(Y, [1 (npointsperSide)^3]);
Z = reshape(Z, [1 (npointsperSide)^3]);

PointSet = [X' Y' Z'];

disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox3D_04    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);
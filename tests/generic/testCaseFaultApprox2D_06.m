% 2D test case 06

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.05;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 5;
eps = 1e-8;
MySettings.errMin = 1e-4;
MySettings.errMax = 0.002;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'results/testFaultApprox2D_06.vtu';
MyProb.Testfunc = @testFuncFaultApprox2D_06;
MyProb.Xmin = [-eps -eps];
MyProb.Xmax = [1+eps 1+eps];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;

npointsperSide = 10;

PointSet = 0:1/(npointsperSide-1):1;
[X, Y] = ndgrid(PointSet);
X = reshape(X, [1 (npointsperSide)^2]);
Y = reshape(Y, [1 (npointsperSide)^2]);

PointSet = [X' Y'];

disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox2D_06    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);
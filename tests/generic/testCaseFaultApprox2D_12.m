% 2D test case 12

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.02;
MySettings.abstolBisection = 0.0005;
MySettings.NumPointsLocal = 7;
MySettings.errMax = 0.002;
MySettings.errMin = 1e-4;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'results/testFaultApprox2D_12.vtu';
MyProb.Testfunc = @testFuncFaultApprox2D_12;
MyProb.Xmin = [0 0];
MyProb.Xmax = [1 1];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;

npointsperSide = 20;

PointSet = 0:1/(npointsperSide-1):1;
[X, Y] = ndgrid(PointSet);
X = reshape(X, [1 (npointsperSide)^2]);
Y = reshape(Y, [1 (npointsperSide)^2]);

PointSet = [X' Y'];

disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox2D_12    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);
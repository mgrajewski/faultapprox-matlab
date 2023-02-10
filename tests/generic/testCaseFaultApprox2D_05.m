% 2D test case 05

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.05;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 7;
MySettings.errMin = 1e-4;
MySettings.errMax = 0.002;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'results/testFaultApprox2D_05.vtu';
MyProb.Testfunc = @testFuncFaultApprox2D_05;
MyProb.Xmin = [0 0];
MyProb.Xmax = [1 1];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;

PointSet = CreateHaltonSet(100,2,1);

disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox2D_05    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);
% 2D test case 04

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.05;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 5;
MySettings.errMax = 0.002;
MySettings.errMin = 1e-4;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'results/testFaultApprox2D_04.vtu';
MyProb.Testfunc = @testFuncFaultApprox2D_04;
MyProb.Xmin = [0 0];
MyProb.Xmax = [1 1];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;

PointSet = CreateHaltonSet(100,2,1);


disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox2D_04    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);
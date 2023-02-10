% 2D test case 01

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.1;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 10;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'results/testFaultApprox2D_01.vtu';
MyProb.Testfunc = @testFuncFaultApprox2D_01;
MyProb.Xmin = [0 0];
MyProb.Xmax = [1 1];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;
PointSet = CreateHaltonSet(20, 2, 1);

disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox2D_01    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);
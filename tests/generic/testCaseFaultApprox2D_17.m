% 2D test case 17

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.1;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 5;
MySettings.errMax = 0.002;
MySettings.errMin = 1e-4;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'results/testFaultApprox2D_17.vtu';
MyProb.Testfunc = @testFuncFaultApprox2D_17;
MyProb.Xmin = [0 0];
MyProb.Xmax = [1 1];

PointSet = CreateHaltonSet(30,2,1);

disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox2D_17    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);
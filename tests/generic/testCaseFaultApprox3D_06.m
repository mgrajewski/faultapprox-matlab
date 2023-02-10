% 3D test case 06

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)

MySettings = FaultApproxParameters;
MySettings.NumPointsLocal = 7;
MySettings.maxDistForSurfacePoints = 0.15;
MySettings.abstolBisection = 0.001;
MySettings.maxiterAdapt = 3;
MySettings.errMax = 0.005;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'results/testFaultApprox3D_06.vtu';
MyProb.Testfunc = @testFuncFaultApprox3D_06;
MyProb.Xmin = [-1.5 -1.5 -1.5];
MyProb.Xmax = [1.5 1.5 1.5];
MyProb.verboseMode = true;

PointSet = CreateHaltonSet(500,3,1);

for idim = 1: 3
    PointSet(:,idim) = (MyProb.Xmax(idim) - MyProb.Xmin(idim))*PointSet(:,idim) + MyProb.Xmin(idim);
end

disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox3D_06    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);

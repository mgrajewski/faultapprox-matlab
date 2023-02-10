% 3D test case 08

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
disp(' __________________________________');
disp('/                                  \');
disp('|    test case FaultApprox3D_08    |');
disp('\__________________________________/');

resultfile = 'results/points_TestCaseFaultApprox3D_08.txt';
normalsfile = 'results/normals_TestCaseFaultApprox3D_08.txt';

MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.1;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 10;
MySettings.errMax = 0.005;
MySettings.maxiterAdapt = 3;
MySettings.maxTrialsForFillingGaps = 5;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = '';
MyProb.Testfunc = @testFuncFaultApprox3D_08;
MyProb.Xmin = [0 0 0] ;
MyProb.Xmax = [1 1 1];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;

PointSet = CreateHaltonSet(200,3,1);

for idim = 1: 3
    PointSet(:,idim) = (MyProb.Xmax(idim) - MyProb.Xmin(idim))*PointSet(:,idim) + MyProb.Xmin(idim);
end

Subdomains = faultApprox(PointSet, MyProb, MySettings);

for iclass = 1:3
    for jclass = iclass+1:3
        points = Subdomains{iclass}{jclass}{2};
        writematrix(points, [resultfile, '_', int2str(iclass), '_', int2str(jclass), '.txt']);
    
        normals = Subdomains{iclass}{jclass}{3};
        writematrix(normals, [normalsfile, '_', int2str(iclass), '_', int2str(jclass), '.txt']);
    end
end


% 3D test case 07

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
resultfile = 'results/points_TestCaseFaultApprox3D_07.txt';
normalsfile = 'results/normals_TestCaseFaultApprox3D_07.txt';

MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.05;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 10;
MySettings.errMax = 0.005;
MySettings.maxiterAdapt = 15;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = '';
MyProb.Testfunc = @testFuncFaultApprox3D_07;
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
disp('|    test case FaultApprox3D_07    |');
disp('\__________________________________/');

Subdomains = faultApprox(PointSet, MyProb, MySettings);

for iclass = 1:3
    for jclass = iclass+1:3
        points = Subdomains{iclass}{jclass}{2};
        writematrix(points, [resultfile, '_', int2str(iclass), '_', int2str(jclass), '.txt']);
    
        normals = Subdomains{iclass}{jclass}{3};
        writematrix(normals, [normalsfile, '_', int2str(iclass), '_', int2str(jclass), '.txt']);
    end
end

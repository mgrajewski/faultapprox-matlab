% 2D test case in "Detecting and approximating decision boundaries in low
% dimensional spaces", section 2.1, Test problem 2.3

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
statsfile = 'testCasePaper2D_01_stats.csv';

MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.05;
MySettings.abstolBisection = 0.001;
MySettings.errMax = 0.001;
MySettings.errMin = 1e-4;
MySettings.NumPointsLocal = 10;

MyProb = ProblemDescr;
MyProb.OutputFileVTU = 'testCase2D_01.vtu';
MyProb.Testfunc = @testFunc2D_01;

% domain Omega
MyProb.Xmin = [0 0];
MyProb.Xmax = [1 1];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;

% compute initial point set
PointSet = CreateHaltonSet(50, 2, 1);

% map to [Xmin, Xmax]
for idim = 1: 2
    PointSet(:,idim) = (MyProb.Xmax(idim) - MyProb.Xmin(idim))*PointSet(:,idim) + MyProb.Xmin(idim);
end

Subdomains = faultApprox(PointSet, MyProb, MySettings);

writelines('testCase 2D_01', statsfile, WriteMode='overwrite')
% save the statistical results
line = [ 'ncalls_', ExtendedStats.pos_in_code{1}, ', ', int2str(ExtendedStats.ncalls{1})];
writelines(line, statsfile, WriteMode='append')
line = [ 'ncalls_', ExtendedStats.pos_in_code{2}, ', ', int2str(ExtendedStats.ncalls{2})];
writelines(line, statsfile, WriteMode='append')

for ipos =3:size(ExtendedStats.pos_in_code,2)
    line = [ 'ncalls_', ExtendedStats.pos_in_code{ipos}, ', ', int2str(ExtendedStats.ncalls{ipos}-ExtendedStats.ncalls{ipos-1})];
    writelines(line, statsfile, WriteMode='append')
end

for ipos =3:size(ExtendedStats.pos_in_code,2)
    for iclass = 1:3
        for jclass = iclass+1:3
            line = [ 'S', int2str(iclass), int2str(jclass), '_', ExtendedStats.pos_in_code{ipos}, ', ', int2str(ExtendedStats.nPointsSurf{ipos-2}{iclass,jclass})];
            writelines(line, statsfile, WriteMode='append')
        end
    end

end
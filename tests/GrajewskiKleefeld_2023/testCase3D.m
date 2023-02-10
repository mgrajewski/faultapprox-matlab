% 3D test case in "Detecting and approximating decision boundaries in low
% dimensional spaces", section 2.2

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
statsfile = 'testCasePaper3D_stats.csv';
resultfile = 'raw_results/testCasePaper3D_points';
normalsfile = 'raw_results/testCasePaper3D_normals';


MySettings = FaultApproxParameters;
MySettings.maxDistForSurfacePoints = 0.05;
MySettings.abstolBisection = 0.001;
MySettings.NumPointsLocal = 10;
MySettings.errMax = 0.002;
MySettings.errMax = 0.001;
MySettings.maxiterAdapt = 4;
MySettings.maxTrialsForFillingGaps = 15;


MyProb = ProblemDescr;
MyProb.OutputFileVTU = '';
MyProb.Testfunc = @testFunc3D;
MyProb.Xmin = [0 0 0] ;
MyProb.Xmax = [1 1 1];
MyProb.verboseMode = true;
MyProb.extendedStats = true;

global ExtendedStats;
ExtendedStats = Statistics;

PointSet = CreateHaltonSet(200,3,1);

Subdomains = faultApproximation(PointSet, MyProb, MySettings);

% write final point sets
for iclass = 1:3
    for jclass = iclass+1:3
        points = Subdomains{iclass}{jclass}{2};
        writematrix(points, [resultfile, '_', int2str(iclass), '_', int2str(jclass), '.txt']);
    
        normals = Subdomains{iclass}{jclass}{3};
        writematrix(normals, [normalsfile, '_', int2str(iclass), '_', int2str(jclass), '.txt']);
    end
end

writelines('testCase 3D', statsfile, WriteMode='overwrite')
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

            title = ['"x", "y", "z"'];
            actualFileName = [resultfile, '_', ExtendedStats.pos_in_code{ipos}, '_', int2str(iclass), '_', int2str(jclass), '.txt'];
            writelines(title, actualFileName, WriteMode='overwrite')
            % save the intermediate point clouds
            writematrix(ExtendedStats.PointSetsSurf{ipos-2}{iclass, jclass}{1}, actualFileName, WriteMode='append');
        end
    end

end
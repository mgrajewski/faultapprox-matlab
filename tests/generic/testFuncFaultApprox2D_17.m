% test function for test case 2D_17: 3 subdomains meeting in a point inside
% the domain, straight lines

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_17(PointSet, ProblemDescr)

    % test case: line segments
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,2) > 1.1 -PointSet(:,1)) = 2;
    ClassOfPoints(PointSet(:,2) > 0.5) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
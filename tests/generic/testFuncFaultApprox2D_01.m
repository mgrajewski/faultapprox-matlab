% function for test case 2D_01: 2 subdomains, one straight line in between

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_01(PointSet, ProblemDescr)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(atan(5*PointSet(:,1)+ 5*PointSet(:,2)) < 0.5) = 2;
    ClassOfPoints(PointSet(:,2) > PointSet(:,1)*0.7+0.5) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
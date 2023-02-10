% function for test case 2D_02: 2 subdomains, subdomain 2: rectangle
% attached to left domain boundary

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_02(PointSet, ProblemDescr)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,1) > 0.55) = 2;
    ClassOfPoints(PointSet(:,2) > 0.78) = 2;
    ClassOfPoints(PointSet(:,2) < 0.33) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
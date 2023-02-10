% test function for test case 2D_11: 2 subdomains, subdomain 2 consists of
% one small half-circle attached to the top side of [0,1]^2.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_11(PointSet, ProblemDescr)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints((PointSet(:,2)-1.0).*(PointSet(:,2)-1.0) + (PointSet(:,1)-0.5).*(PointSet(:,1)-0.5) < 0.1) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
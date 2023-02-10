% function for test case 2D_03: 4 subdomains, subdomain 3: square in the
% middle of the domain, straight lines only

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_03(PointSet, ProblemDescr)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,2) <= PointSet(:,1) & PointSet(:,2)>= 1 - PointSet(:,1)) = 2;
    ClassOfPoints(PointSet(:,2) >= PointSet(:,1) & PointSet(:,2)>= 1 - PointSet(:,1)) = 4;
    ClassOfPoints(PointSet(:,1) >= 0.3 & PointSet(:,1) < 0.7 & PointSet(:,2) > 0.3 & PointSet(:,2) < 0.7) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
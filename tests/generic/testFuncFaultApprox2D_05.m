% function for test case 2D_05: 3 subdomains, subdomain 1 consists of two
% components, straight boundaries

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_05(PointSet, ProblemDescr)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,1) < PointSet(:,2) - 0.4) = 2;
    ClassOfPoints(PointSet(:,1) > PointSet(:,2) - 0.2) = 2;
    ClassOfPoints(PointSet(:,1) > PointSet(:,2) + 0.3) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
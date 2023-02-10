% function for test case 2D_04: 3 subdomains, subdomain 2 consists of two
% components, curved boundaries

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_04(PointSet, ProblemDescr)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,2) > 0.7+ 0.1*(sin(12*PointSet(:,1)))) = 2;
    ClassOfPoints(PointSet(:,2) < 0.3+ 0.1*(sin(12*PointSet(:,1)))) = 2;
    ClassOfPoints((PointSet(:,2)-0.5).*(PointSet(:,2)-0.5) + (PointSet(:,1)-1).*(PointSet(:,1)-1) >= 0.2) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
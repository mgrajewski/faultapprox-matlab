% test function for test case 2D_13: 2 subdomains, subdomains 2 consists of
% two components, subdomain 1 of three

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_13(PointSet, ProblemDescr)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints((PointSet(:,2)-0.72).*(PointSet(:,2)-0.72) + (PointSet(:,1)-0.0).*(PointSet(:,1)-0.0)<0.03 & ...
                  (PointSet(:,2)-0.72).*(PointSet(:,2)-0.72) + (PointSet(:,1)-0.0).*(PointSet(:,1)-0.0)> 0.01) = 2;
    
    ClassOfPoints((PointSet(:,2)-0.28).*(PointSet(:,2)-0.28) + (PointSet(:,1)-0.0).*(PointSet(:,1)-0.0)<0.03 & ...
                  (PointSet(:,2)-0.28).*(PointSet(:,2)-0.28) + (PointSet(:,1)-0.0).*(PointSet(:,1)-0.0)> 0.01) = 2;

    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
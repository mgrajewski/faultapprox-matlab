% function for test case 3D_02: 2 subdomains, the subdomain 2 consists of
% two facets inside the domain

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox3D_02(PointSet, ProblemDescr)

    % 3D example: fault surface consisting of two planes
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(-PointSet(:,1) + 2*PointSet(:,2) < 0.5 & 4*PointSet(:,1) + PointSet(:,2) < 3.2) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,3) < ProblemDescr.Xmin(3)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;
    ClassOfPoints(PointSet(:,3) > ProblemDescr.Xmax(3)) = -1;

end

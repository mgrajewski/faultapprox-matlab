% function for test case 2D_06: 2 subdomains, subdomain 3 is a circle
% around (0.5, 0.5) with radius 1/sqrt{5}.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_06(PointSet, ProblemDescr)

    % simple 2D example with closed boundary segment (circle)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints((PointSet(:,1) - 0.5).*(PointSet(:,1) - 0.5) + (PointSet(:,2) - 0.5).*(PointSet(:,2) - 0.5) < 0.2) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
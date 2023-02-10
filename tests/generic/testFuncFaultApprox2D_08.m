% function for test case 2D_08: 2 subdomains, which consist of two
% components each.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_08(PointSet, ProblemDescr)

    % 2D-test case fo subdomains consisting of several components
    
    ClassOfPoints = 2*ones(size(PointSet,1),1);
    
    % Ball with radius 0.15 around (1,1)
    ClassOfPoints((PointSet(:,1)-1.0).*(PointSet(:,1)-1.0) + (PointSet(:,2)-1.0).*(PointSet(:,2)-1.0) < 0.0225) = 1;

    % second component of subdomain I is a kind of stripe
    ClassOfPoints((PointSet(:,2) < 0.5*PointSet(:,1)-0.1) & (PointSet(:,2) > 0.5*PointSet(:,1)-0.3)) = 1;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
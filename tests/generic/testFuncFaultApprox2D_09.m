% test function for test case 2D_09: 3 subdomains, subdomain 1 consists of
% two vertical stripes above 0.7 and below 0.3; subdomain 2 consists of
% two triangles starting in (0.85, 0.7) and (0.85, 0.3) at the right
% domain boundary

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_09(PointSet, ProblemDescr)

    % test case for subdomains consisting of several components
    ClassOfPoints = 3*ones(size(PointSet,1),1);
    
    % subdomain 1 consists of two vertical stripes above 0.7 and below 0.3
    ClassOfPoints(PointSet(:,2) > 0.7) = 1;
    ClassOfPoints(PointSet(:,2) < 0.3) = 1;

    % subdomain 2 consists of two triangles starting in (0.85, 0.7) and
    % (0.85, 0.3) at the right domain boundary
    ClassOfPoints((PointSet(:,2) < 0.5*PointSet(:,1)+ 0.275) & (PointSet(:,2) > -0.5*PointSet(:,1) + 1.125)) = 2;
    ClassOfPoints((PointSet(:,2) < 0.5*PointSet(:,1)- 0.125) & (PointSet(:,2) > -0.5*PointSet(:,1) + 0.725)) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
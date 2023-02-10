% function for test case 2D_07: 2 subdomains, subdomain 2 consists of two
% separated half-circles attached to the right side of the domain [0,1]^2.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_07(PointSet, ProblemDescr)
    % 2D-test case for subdomains consisting of sevferal components
    ClassOfPoints = ones(size(PointSet,1),1);
    
    % Ball with radius 0.15 around (1,0.7)
    ClassOfPoints((PointSet(:,1)-1.0).*(PointSet(:,1)-1.0) + (PointSet(:,2)-0.7).*(PointSet(:,2)-0.7) < 0.0225) = 2;

    % Ball with radius 0.15 around (1,0.3)
    ClassOfPoints((PointSet(:,1)-1.0).*(PointSet(:,1)-1.0) + (PointSet(:,2)-0.3).*(PointSet(:,2)-0.3) < 0.0225) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
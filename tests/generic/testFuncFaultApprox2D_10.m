% test function for test case 2D_10: 2 subdomains, subdomain 2 consists of
% three small half-circles attached to the bottom and left side of the
% domain boundary

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)MySettings = FaultApproxParameters;
function ClassOfPoints = testFuncFaultApprox2D_10(PointSet, ProblemDescr)
    % 2D-test case for subdomains consisting of several components

    ClassOfPoints = ones(size(PointSet,1),1);
    
    % Ball with radius 0.15 around (0,0.7)
    ClassOfPoints((PointSet(:,1)).*(PointSet(:,1)) + (PointSet(:,2)-0.7).*(PointSet(:,2)-0.7) < 0.0225) = 2;

    % Ball with radius 0.15 around (0,0.3)
    ClassOfPoints((PointSet(:,1)).*(PointSet(:,1)) + (PointSet(:,2)-0.3).*(PointSet(:,2)-0.3) < 0.0225) = 2;

    % Ball with radius 0.2 around (0.5,0.0)
    ClassOfPoints((PointSet(:,1)-0.5).*(PointSet(:,1)-0.5) + PointSet(:,2).*PointSet(:,2) < 0.04) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
% function for test case 3D_05: sphere with radius 1/3 around
% [0.5, 0.5, 0.5] (inside: 2, outside: 1)

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox3D_05(PointSet, ProblemDescr)

    % 3D example: sphere with radius 1/3 around [0.5, 0.5, 0.5] (inside:
    % 2, outside: 1)
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints((PointSet(:,1)-0.5).*(PointSet(:,1)-0.5) + ...
                  (PointSet(:,2)-0.5).*(PointSet(:,2)-0.5) + ...
                  (PointSet(:,3)-0.5).*(PointSet(:,3)-0.5) < 0.1111111) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,3) < ProblemDescr.Xmin(3)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;
    ClassOfPoints(PointSet(:,3) > ProblemDescr.Xmax(3)) = -1;
end
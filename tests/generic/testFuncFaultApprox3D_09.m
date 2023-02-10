% function for test case 3D_09: 3 subdomains, parts of two
% concentric spheres around (0,0,0)

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox3D_09(PointSet, ProblemDescr)

    % 3D example: part of two concentric spheres
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,1).*PointSet(:,1) + PointSet(:,2).*PointSet(:,2) + PointSet(:,3).*PointSet(:,3) < 0.25) = 2;
    ClassOfPoints(PointSet(:,1).*PointSet(:,1) + PointSet(:,2).*PointSet(:,2) + PointSet(:,3).*PointSet(:,3) < 0.12) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,3) < ProblemDescr.Xmin(3)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;
    ClassOfPoints(PointSet(:,3) > ProblemDescr.Xmax(3)) = -1;

end
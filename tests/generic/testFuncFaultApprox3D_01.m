% function for test case 3D_01: 2 subdomains on [0,1]^3, subdomain 2 is
% part of a sphere around (0,0,0) with radius 1/sqrt{5}.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox3D_01(PointSet, ProblemDescr)

    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,1).*PointSet(:,1) + PointSet(:,2).*PointSet(:,2) + ...
                  PointSet(:,3).*PointSet(:,3) < 0.2) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
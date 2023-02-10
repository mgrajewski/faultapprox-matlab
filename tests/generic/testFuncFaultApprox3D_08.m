% function for test case 3D_08: 3 subdomains, subdomain 3 is a
% sphere with radius 0.3 around (0.5, 0.5, 0.5); subdomains 1 and 2
% split the unit cube in the middle

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox3D_08(PointSet, ProblemDescr)

    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,1) > 0.5) = 2;
    ClassOfPoints( (PointSet(:,1) -0.5).^2 + (PointSet(:,2)).^2 + ...
                   (PointSet(:,3) -0.5).^2 < 0.09) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,3) < ProblemDescr.Xmin(3)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;
    ClassOfPoints(PointSet(:,3) > ProblemDescr.Xmax(3)) = -1;

end
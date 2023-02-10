% function for test case 3D_06: 2 subdomains, subdomain 2 has the
% shape of a diamond (section of the l_1-unit cube and the l_oo unit
% cube (inside: 2, outside: 1)

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox3D_06(PointSet, ProblemDescr)

    % 3D example: diamond (section of the l_1-unit cube and the l_oo unit
    % cube
    % (inside: 2, outside: 1)
    ClassOfPoints = ones(size(PointSet,1),1);
    
    % get angles from cartesian coordinates
    
    ClassOfPoints(vecnorm(PointSet', Inf) <= 1 & vecnorm(PointSet', 1) <= 1.5,:)  = 2;
        
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,3) < ProblemDescr.Xmin(3)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;
    ClassOfPoints(PointSet(:,3) > ProblemDescr.Xmax(3)) = -1;
end
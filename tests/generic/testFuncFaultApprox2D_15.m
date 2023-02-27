% test function for test case 2D_15: 3 subdomains meeting in a point inside
% the domain, curved lines, acute angle in one boundary

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab) 
function ClassOfPoints = testFuncFaultApprox2D_15(PointSet, ProblemDescr)

    % simple 2D example: part of a sphere
    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,2) > 0.4 & PointSet(:,1) + PointSet(:,2) <1.2) = 2;
    %ClassOfPoints(PointSet(:,2) > 1/4*sin(6*PointSet(:,1))+0.6 & PointSet(:,1) + PointSet(:,2) >=1.2) = 3;
    ClassOfPoints(PointSet(:,2) > 1/6*sin(12*PointSet(:,1))+0.7 & PointSet(:,1) + PointSet(:,2) >=1.2) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
% test function for test case 2D_12: 4 subdomains, subdomains 1,2 and 4
% consists two components; straight lines

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_12(PointSet, ProblemDescr)
    ClassOfPoints = 3*ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,2) > 0.7 | PointSet(:,2) < 0.3) = 2;

    ClassOfPoints(PointSet(:,2) > 0.7 & PointSet(:,1) < 0.2 & PointSet(:,2) < 0.9) = 1;
    ClassOfPoints(PointSet(:,2) < 0.3 & PointSet(:,1) < 0.2 & PointSet(:,2) > 0.1) = 1;

    ClassOfPoints(PointSet(:,2) >= 0.9 & PointSet(:,2) > PointSet(:,1) + 0.7) = 4;
    ClassOfPoints(PointSet(:,2) <= 0.1 & PointSet(:,2) < -PointSet(:,1) + 0.3) = 4;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
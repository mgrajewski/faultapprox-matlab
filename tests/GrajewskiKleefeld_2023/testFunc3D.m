% function for the 3D example in "Detecting and
% approximating decision boundaries in low dimensional spaces", section 2.2

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = TestFunc3D(PointSet, ProblemDescr)

    ClassOfPoints = ones(size(PointSet,1),1);
    ClassOfPoints(PointSet(:,2) + 0.1* PointSet(:,3) > 0.7 + 0.1*sin(10*PointSet(:,1).^1.5)+ 0.05*sin(5*PointSet(:,3).^1.5)) = 3;
    ClassOfPoints( (PointSet(:,1) -1.0).^6 + (PointSet(:,2) - 0.5).^6 + (PointSet(:,3) -0.5).^6 < 0.002) = 2;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,3) < ProblemDescr.Xmin(3)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;
    ClassOfPoints(PointSet(:,3) > ProblemDescr.Xmax(3)) = -1;

end
% Function for Test problem 4.2 in "Detecting and approximating decision
% boundaries in low dimensional spaces", section 4.1.
% It is based on Example 1 from Allasia, Giampietro et al., Efficient
% approximation algorithms. Part I: approximation of unknown fault lines
% from scattered data, Dolomites Research Notes On Approximation,
% Vol. 3, (2010), p. 7-38

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFunc2D_4_2(PointSet, ProblemDescr)

    ClassOfPoints = ones(size(PointSet,1),1);

    ClassOfPoints(PointSet(:,1) > 0.5) = 2;
    ClassOfPoints(PointSet(:,1) > 0.6) = 3;
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
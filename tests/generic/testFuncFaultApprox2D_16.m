% test function for test case 2D_16: 6 subdomains meeting in a point inside
% the domain, curved lines

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = testFuncFaultApprox2D_16(PointSet, ProblemDescr)

    % test case: many segments meeting in the middle
    center = [0.5, 0.5];

    nclasses = 6;

    aux = PointSet - center;
    args = atan2(aux(:,1), aux(:,2));
    dist = vecnorm(aux,2,2);

    ClassOfPoints = ones(size(PointSet,1),1);
    for iclass = 1: nclasses
        ClassOfPoints(args + 1.2*dist.^0.5> -pi+ (iclass-1)/nclasses*2*pi & args +1.2*dist.^0.5 <= -pi+ (iclass)/nclasses*2*pi) = iclass;
    end

    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end

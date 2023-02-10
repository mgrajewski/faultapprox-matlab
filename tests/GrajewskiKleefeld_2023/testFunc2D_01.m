function ClassOfPoints = testFunc2D_01(PointSet, ProblemDescr)
    ClassOfPoints = ones(size(PointSet,1),1);
    %ClassOfPoints(PointSet(:,2) > 0.6 + 0.1*(sin(10*pi*PointSet(:,1).*PointSet(:,1)))) = 2;
    ClassOfPoints(PointSet(:,2) > 0.6 + 0.1*(sin(10*pi*PointSet(:,1).^1.5))) = 2;
    ClassOfPoints((PointSet(:, 2) - 0.5).^6 + (PointSet(:, 1) - 1.0).^6 <= 0.005) = 3;
    
    
    ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
    ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
    ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
    ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;

end
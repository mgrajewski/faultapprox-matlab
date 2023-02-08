% This function is an auxiliary function for adaptive refinement of a
% set of triplets PointsInPlane representing a part of a fault
% surface. We assume that a Delaunay triangulation exists with
% connectivity list ConnectivityList. Moreover, we assume that an RBF
% interpolation of PointsInPlane is given by the vector coeffs and
% scale.
% We want to estimate the maximum error of linear interpolation on the
% triangle with index itri.
% For estimating the maximum error on itri, we rely on
% Waldron, Shayne: The error in linear interpolation at the vertices of a
% simplex, SIAM J. Numer. Analysis, Vol. 35, pp. 1191-1200, 1998
%
% In this work, formula (4.5) of theorem 4.1 states
% || f - I_f||_oo <= 0.5(R^2-d^2) |f|_2,oo,T
%
% Here, R denotes the radius of the circumcircle and d the distance of the
% center c of the circumcircle to the triangle T = ConnectivityList(itri,:)
% with vertices PointsInPlane(ConnectivityList(itri,:)).
%
% We estimate |f|_2,oo,T by estimating the second derivative of $f$
%
% Input:

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function errMax = estimateMaxErr3D(PointsInPlane, ConnectivityList, itri, coeffs, scale)

    iidx = ConnectivityList(itri,:);
    
    numPoints = size(PointsInPlane, 1);
 
    ScaleVec = scale*ones(1, numPoints);
    
    triPoints = PointsInPlane(iidx, :);
    
    triPoints = triPoints - triPoints(2,:);

    [curvature, circleCenter] = curvCenterFromThreePoints(triPoints);
    if (curvature > 1e-8)
        radius = 1.0/curvature;    
    else
        radius = 1e10;
        circleCenter = a;
    end
    
    % The center of the circumcircle is inside the triangle.
    if inpolygon(circleCenter(1), circleCenter(2), triPoints(:,1), triPoints(:,2))
        d = 0;

    % The center of the circumcircle is outside the triangle.
    else
        % compute the distance to the triangle as minimum of the distances
        % to the three lines
        edge1 = triPoints(2,:) - triPoints(1,:);
        edge2 = triPoints(3,:) - triPoints(1,:);
        edge3 = triPoints(3,:) - triPoints(2,:);
        
        p1 = ((circleCenter-triPoints(1,:))*edge1')/(edge1*edge1')*edge1;
        d1 = norm(circleCenter- p1 - triPoints(1,:));

        p2 = ((circleCenter-triPoints(1,:))*edge2')/(edge2*edge2')*edge2;
        d2 = norm(circleCenter- p2 - triPoints(1,:));

        p3 = ((circleCenter-triPoints(2,:))*edge3')/(edge3*edge3')*edge3;
        d3 = norm(circleCenter- p3 - triPoints(2,:));

        d = min([d1, d2, d3]);
    end
      
    % compute second derivative
    phiddot = zeros(4, 2,2);
    bary = 1/3*(PointsInPlane(iidx(1), :) + PointsInPlane(iidx(2), :)+ PointsInPlane(iidx(3), :));
    evalPoints = [PointsInPlane(iidx, :); bary];
    for i = 1: numPoints
        phiddot = phiddot + coeffs(i)* Gaussian_second_der(evalPoints, PointsInPlane(i,:), ...
                                      ScaleVec(i));
    end
         
    normPhiddot = zeros(4,1);
    for i = 1: 4
        normPhiddot(i) = max(max(abs(phiddot(i,:,:))));
    end

    normPhiddot = max(normPhiddot);
    errMax = (radius*radius- d*d)*normPhiddot;
end
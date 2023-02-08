% This function approximates the unit normal vector to a surface given by
% the points in Points.
% We compute the optimally fitting plane or line in the sense that the sum
% of the squared distances between the points and the plane is minimal; we
% use its normal vector as approximation.
% It is well known (see e.g. Shakarji, M.: Least-Squares Fitting Algorithms
% of the NIST Algorithm Testing System, J. Res. Natl. Inst. Stand. Technol.
% 103, 633 (1998), https://nvlpubs.nist.gov/nistpubs/jres/103/6/j36sha.pdf)
% that the optimal fitting plane contains the midpoint.
% Therefore, we shift the point such that 0 is the new midpoint beforehand.
% The normal vector is the right singular vector with respect to the
% smallest singular value of the point set. As the singular values are
% provided in descending order, this is the "last" right singular vector.
% For an explanation and the algorithm itself, we refer to the
% aforementioned article.
%
% Input:
% - PointSet: set of points to fit a line/plane to
%
% Output:
% - normal: approximate normal vector

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function normal = computeNormalVec(PointSet)
    [numPointsLocal, ndim] = size(PointSet);
    if ( numPointsLocal >= ndim)

        % mean value of coordinates
        xmean = 1/size(PointSet, 1)*ones(1, size(PointSet, 1))*PointSet;

        % shift point set such that its mean is zero
        PointSet = PointSet - xmean;
                   
        [~,~,Q] = svd(PointSet);
        normal = Q(:,ndim);
    else
        warning('Not enough points provided for computing a normal vector.')
        normal = zeros(1, ndim);
    end
end
% This function computes the maximal angle between consecutive segments of
% a polygonal line given by the points in PointSet.
%
% Input:
% - PointSet: set of points which defines a polygonal line
%
% Output:
% - maxAngle: maximal angle between to consecutive line segments (measured
% in rad)

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function maxAngle = computeMaximalAngle(PointSet)

    epstol = 1e-10;

    % number of points
    npoints = size(PointSet,1);

    % segments of the polygonal line
    segs = PointSet(1:npoints-1,:) - PointSet(2:npoints,:);

    % square of the norm of the segments
    normSegs = (segs.*segs)*[1;1];

    if any(normSegs < epstol*epstol)
        error(['Almost duplicate points detected. ' ...
               'This preventes reliably computing the maximal angle.'])
    end

    % We compute the angle by its definition
    %angle(a,b) = arccos(<a,b>/(||a|| ||b||))
    denominator = sqrt(normSegs(1:npoints-2).*normSegs(2:npoints-1));
    
    nominator = zeros(npoints-2,1);
    for i = 1:npoints-2
        nominator(i) = segs(i,:)*segs(i+1,:)';
    end
    
    % We employ the minimum, as cosine = 1 for parallel segments and
    % cosine = -1 for angles approaching pi.
    % It may happen that due to rouding errors, nomintor/denominator <
    % -1.0. Therefore, we safeguard the result with the maximum.
    maxAngle = acos(max(min(nominator./denominator), -1.0));
end
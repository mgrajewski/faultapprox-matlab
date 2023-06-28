% This function computes npoints equidistant points on a polygonal line
% given by PoyLine.
%
% Input:
% - PolyLine: ordered point set, which defines a polygonal line
% - number of points to compute
%
% Output:
% - PolyPoints, vector containing the equidistant points

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function PolyPoints = getPointsFromPoly(PolyLine, npoints)

    % compute distance vector and arclength
    % works in 2D only
    if (size(PolyLine,2) ~= 2)
        warning('Function getPointsFromPoly works for 2D only')
    else
        numPointsPoly = size(PolyLine,1);
        LineSegs = PolyLine(2:numPointsPoly,:) - PolyLine(1:numPointsPoly-1,:);
        PolyPoints = zeros(npoints,2);

        ArcLength = zeros(numPointsPoly,1);
        for i = 2: numPointsPoly
            ArcLength(i) = ArcLength(i-1) + norm(LineSegs(i-1,:),2);
        end

        arcLengthTotal = ArcLength(numPointsPoly);
        for ipoint = 1: npoints
            arcLengthNow = arcLengthTotal/(npoints-1)*(ipoint-1);

            istartIdx = find(ArcLength <= arcLengthNow, 1, 'last');
            arcLengthStart = ArcLength(istartIdx);

            if (istartIdx < numPointsPoly)
                arcLengthEnd = ArcLength(istartIdx+1);
                parVal = (arcLengthNow - arcLengthStart)/(arcLengthEnd - arcLengthStart);
                PolyPoints(ipoint,:) = PolyLine(istartIdx,:) + parVal*(PolyLine(istartIdx+1,:)-PolyLine(istartIdx,:));

            % very last point: just take this one
            else
                PolyPoints(ipoint,:) = PolyLine(numPointsPoly,:);
            end
            
        end

    end
end
% This function computes the distances of the points in PointSetTest to a
% polygonal line. It works for two dimensions only.
% For any point in PointSetTest, we compute the distance to the polygon
% defined in PolyLine. This distance is either the distance to one of the
% points in PolyLine, or the length of the orthogonal projection onto one
% of the line segments. 
% The orthogonal projection on the line defined by
% PolyLine[i], PolyLine[i+1] is, setting
% line_seg = PolyLine[i+1] - PolyLine[i]:
%
%    p(test point) = PolyLine[i] + alpha * line_seg ,
%    alpha = <test point, line_seg> /norm(line_seg)^2
%
% If 0 <= alpha <= 1, the projection is on line_seg, and we compute the 
% distance to the line segments in this case.
%
% Input:
% - PolyLine: ordered point set, which defines a polygonal line
% - PointSetTest: set of points to test on another polygonal line
%
% Output:
% - DistVec, vector of euclidean distances to PolyLine

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function DistVec = distOfPolyLines(PolyLine, PointSetTest)

    epsLenSq = 1e-20;
    epsDistSq = 1-08;
    onLineSegment = false;

    % works in 2D only
    if (size(PolyLine,2) ~= 2)
        warning('Function distOfPolyLines works for 2D only')
    else
        numPointsPoly = size(PolyLine,1);
        numPointsTest = size(PointSetTest,1);
        LineSegs = PolyLine(2:numPointsPoly,:) - PolyLine(1:numPointsPoly-1,:);
        DistVec = zeros(1, numPointsTest);


        for ipoint = 1: numPointsTest

            DistVec(ipoint) = 1e10;

            for iseg = 1: numPointsPoly-1

                % square of the length of line segment
                segLengthSq = LineSegs(iseg,:)*LineSegs(iseg,:)';
                % The projection-based test is reliable only for segments
                % longer than eps_len. If the segment is too small, we
                % just skip this line segment und continue with the next
                % one.
                if segLengthSq > epsLenSq
                    alpha = (PointSetTest(ipoint,:) - PolyLine(iseg,:))*LineSegs(iseg,:)'/segLengthSq;
                    
                    % If so, the projected test point is on the current
                    % line segment.
                    if (0 <= alpha && alpha <= 1)
                        dist = PointSetTest(ipoint,:) - PolyLine(iseg,:) - alpha*LineSegs(iseg,:);

                        % distance to line segment
                        DistVec(ipoint) = min(DistVec(ipoint), norm(dist,2));
                        
                        break
                    end
                else
                    continue
                end
            end


            % test the points on the polygonal line directly
            for jpoint = 1: numPointsPoly
                DistVec(ipoint) = min(DistVec(ipoint), norm(PointSetTest(ipoint,:) - PolyLine(jpoint,:), 2));
            end
        end
    end
end
% This function checks if any point in PointSetTest (in the following 
% called test point) is on the polygonal line defined by PolyLine. It works
% for two dimensions only.
% If the test point is on the polygonal line, it is on one of the line
% segments (or very close to one). This is the case, if the orthogonal
% projection of the test point on the corresponding line is on the line
% segment and if the distance to the corresponding line is zero or very 
% small. The orthogonal projection on the line defined by
% ref_point[i], ref_point[i+1] is, setting
% line_seg = ref_point[i+1] - ref_point[i]:
%
%    p(test point) = ref_point[i] + alpha * line_seg ,
%    alpha = <test point, line_seg> /norm(line_seg)^2
%
% If 0 <= alpha <= 1, the projection is on line_seg, and we compute the 
% distance to the line segments in this case.
%
% Input:
% - PolyLine: ordered point set, which defines the polygonal line
% - PointSetTest: set of points to test to be on the polygonal line given
%   by PolyLine
%
% Output:
% true, if any of the points in PointSetTest is on the polygonal line given
% by PolyLine; false otherwise

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function onLineSegment = isOnPolyLine(PolyLine, PointSetTest)

    epsLenSq = 1e-20;
    epsDistSq = 1-08;
    onLineSegment = false;

    % works in 2D only
    if (size(PolyLine,2) ~= 2)
        warning('Function IsOnPolyLine works for 2D only')
    else
        numPointsPoly = size(PolyLine,1);
        numPointsTest = size(PointSetTest,1);
        LineSegs = PolyLine(2:numPointsPoly,:) - PolyLine(1:numPointsPoly-1,:);
        
        for ipoint = 1: numPointsTest
            for iseg = 1: numPointsPoly-1

                % square of the length of line segment
                segLengthSq = LineSegs(iseg,:)*LineSegs(iseg,:)';
                % If a segment is shorter than 1e-10, we can not seriously
                % decide if a point is close to this line segment. In such 
                % case, we just skip this line segment und continue with
                % the next one.
                if segLengthSq > epsLenSq
                    alpha = (PointSetTest(ipoint,:) - PolyLine(iseg,:))*LineSegs(iseg,:)'/segLengthSq;
                    
                    if (0 <= alpha && alpha <= 1)
                        dist = PointSetTest(ipoint,:) - PolyLine(iseg,:) - alpha*LineSegs(iseg,:);

                        % square of the distance.
                        dist_sq = dist*dist';
                        
                        if (dist_sq < epsDistSq)
                            onLineSegment = true;
                            return
                        end
                        break
                    end
                else
                    continue
                end
            end
        end
    end
end
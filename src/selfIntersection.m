% This function tests, if a polygonal line given by the points in Points,
% which are assumed to be ordered, intersects itself.
% This can be a hint for failed sorting.
% We test any line segment with any following line segment if they
% intersect. This can be efficiently done using determinants. For details,
% we refer to the documentation in polyLinesIntersect.m
%
%Input:
% - Points: ordered point set defining the polygonal line
%
% Output:
% - self_intersect: true, if the polygonal line intersects itself

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function self_intersect = selfIntersection(Points)

    self_intersect = false;

    % works in 2D only
    if (size(Points,2) ~= 2)
        warning('This function works for 2D only')
    else
    
        epsLoc = 1e-10;

        numPoints = size(Points,1);

        if (numPoints < 2)
            return
        end

        LineSegs = Points(1:numPoints-1,:) - Points(2:numPoints,:);

        % This is actually the l1-norm. We scale the line segments for
        % numerical stability.
        NormLineSegs = abs(LineSegs)*[1;1];
        LineSegs = LineSegs./NormLineSegs;

        for iseg = 1: numPoints-1

            % subsequent segments never intersect, thus start the loop at iseg+2
            for jseg = iseg+2:numPoints-1
                % if greater 0: all points of the other segment on the same
                % side
                aux1 = (LineSegs(iseg,1)*(Points(iseg,2)-Points(jseg,2)) - LineSegs(iseg,2)*(Points(iseg,1)-Points(jseg,1))) * ...
                      (LineSegs(iseg,1)*(Points(iseg,2)-Points(jseg+1,2)) - LineSegs(iseg,2)*(Points(iseg,1)-Points(jseg+1,1)));

                aux2 = (LineSegs(jseg,1)*(Points(jseg,2)-Points(iseg,2)) - LineSegs(jseg,2)*(Points(jseg,1)-Points(iseg,1))) * ...
                      (LineSegs(jseg,1)*(Points(jseg,2)-Points(iseg+1,2)) - LineSegs(jseg,2)*(Points(jseg,1)-Points(iseg+1,1)));

                if(aux1 < -epsLoc && aux2 < -epsLoc)
                    self_intersect = true;
                    return
                end
            end
        end
    end
end
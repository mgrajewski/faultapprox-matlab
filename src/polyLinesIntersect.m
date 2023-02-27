% This function tests, if two polygonal lines given by the points in
% Points1 and Points2, which are assumed to be ordered, intersect.
% We test any line segment from the first line with any each line segment 
% of the second one if they intersect.
% This can be efficiently done using determinants.
%                  i+1
%                  x
%                 /
%                /
%     j x-------/------------------x j+1
%              /
%             /
%            x
%           i
% The determinant in 2D determines whether three points are sorted
% clockwise or counterclockwise. If det((i,i+1), (i,j)) and
% det((i,i+1), (j+1,i)) are both positive of both negative, then both
% points j and j+1 are both left or both right from the line segment from 
% i to i+1. Then, the line segment from i to i+1 cannot intersect the line
% segment from j to j+1. However, if the signs differ, this does not
% necessarily mean that the two line segments intersect:
%     j x-----------------------x j+1
%             i+1 x
%                /
%               /
%            i x
% So, if j and j+1 are on different sides with respect to the line segment
% from i to i+1, then we repeat our test for points i and i+1 with respect
% to the line segment from j to j+1.
% This function works in two dimensions only.
%
% Input:
% - Poly1, Poly2: ordered point sets defining the two polygonal lines
%
% Output:
% - doIntersect: true, if the two lines intersect

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function doIntersect = polyLinesIntersect(Poly1, Poly2)

    doIntersect = false;

    % works in 2D only
    if (size(Poly1,2) ~= 2 || size(Poly2,2) ~= 2)
        warning('This function works for 2D only')
    else
    
        epsLoc = 1e-10;

        numPoints1 = size(Poly1,1);
        numPoints2 = size(Poly2,1);

        if (numPoints1 < 2 || numPoints2 < 2)
            return
        end

        LineSegs1 = Poly1(1:numPoints1-1,:) - Poly1(2:numPoints1,:);
        LineSegs2 = Poly2(1:numPoints2-1,:) - Poly2(2:numPoints2,:);
        
        % This is actually the l1-norm. We scale the line segments for
        % numerical stability.
        normLineSegs1 = abs(LineSegs1) * [1;1];
        LineSegs1 = LineSegs1./normLineSegs1;

        for iseg = 1: numPoints1-1
            for jseg = 1:numPoints2-1
                % If greater 0: all points of the other segment on the same
                % side.
                aux1 = (LineSegs1(iseg,1)*(Poly1(iseg,2) - Poly2(jseg,2)) - ...
                        LineSegs1(iseg,2)*(Poly1(iseg,1) - Poly2(jseg,1))) * ...
                       (LineSegs1(iseg,1)*(Poly1(iseg,2) - Poly2(jseg+1,2)) - ...
                        LineSegs1(iseg,2)*(Poly1(iseg,1) - Poly2(jseg+1,1)));

                % Exclude the geometric situation shown above: test, if i
                % and i+1 are on different sides with respect to the line
                % segment from j to j+1.
                if (aux1 < -epsLoc)
                    aux2 = (LineSegs2(jseg,1)*(Poly2(jseg,2) - Poly1(iseg,2)) - ...
                            LineSegs2(jseg,2)*(Poly2(jseg,1) - Poly1(iseg,1))) * ...
                           (LineSegs2(jseg,1)*(Poly2(jseg,2) - Poly1(iseg+1,2)) - ...
                            LineSegs2(jseg,2)*(Poly2(jseg,1) - Poly1(iseg+1,1)));
    
                    if(aux2 < -epsLoc)
                        doIntersect = true;
                        return
                    end
                end
            end
        end
    end
end
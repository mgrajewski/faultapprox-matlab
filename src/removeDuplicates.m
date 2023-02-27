% This function removes duplicate triplets in PointSetI/PointSetJ. Triplets
% closer than epstol are considered identical.
%
% Input: 
% - PointsIclass, PointsJclass: set of triplets describing a fault line
% - epsLoc: tolerance for assuming identity
%
% Output:
% - PointsIclass, PointsJclass: set of triplets describing a fault line,
%   minus duplicate triplets

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointSetI, PointSetJ] = removeDuplicates(PointSetI, PointSetJ, epsLoc)
numAux = size(PointSetI, 1);

    % Inner loop is pointless if the point with index i has been
    % removed before.
    for i = 1: numAux
        if (PointSetI(i,1) ~= -42)
            for j = i+1: numAux
               if(abs(norm(PointSetI(i,:) - PointSetI(j,:))) < epsLoc)
                    PointSetI(j,:) = -42;
                    PointSetJ(j,:) = -42;
               end
            end
        end
    end
    PointSetI = PointSetI(PointSetI(:,1) ~= -42,:);
    PointSetJ = PointSetJ(PointSetJ(:,1) ~= -42,:);
end
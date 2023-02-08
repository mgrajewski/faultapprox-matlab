% This function is a placeholder for a classification procedure.
% We assume that we classify points according to some function. For
% maximal flexibility, this function is provided as a function pointer
% in ProblemDescr. Its behaviour can be controlled by a structure
% FunctionsParameters which is part of the list of function parameters.
% We do not specify what this structure contains nor how its content
% acts on function behaviour but delegate this to the user. We however
% specify the parameter list of our function: 
% f(PointSet, FunctionParameters).
%
% Input:
% - PointSet: set of points to classify
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
%
% Output:
% - ClassOfPoints: integer array containing the class values of the points
%   in Points

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function ClassOfPoints = computeClassification(PointSet, ProblemDescr)

    global ncalls
        
    % there is a test function given
    if (~strcmp(ProblemDescr.Testfunc, ''))
        
        if (size(PointSet,1) > 0)
            
            % update number of function evaluations
            ncalls(1) = ncalls(1) + 1;
            ncalls(2) = ncalls(2) + size(PointSet,1);
            
            ClassOfPoints = ProblemDescr.Testfunc(PointSet, ProblemDescr);
        else
            ClassOfPoints = zeros(0,1);
        end
        
    % There is no test function given: skip the computation.
    else
        error('No function given in ProblemDescr')
    end
    if (size(PointSet, 2) >= 2)
        ClassOfPoints(PointSet(:,1) < ProblemDescr.Xmin(1)) = -1;
        ClassOfPoints(PointSet(:,2) < ProblemDescr.Xmin(2)) = -1;
        ClassOfPoints(PointSet(:,1) > ProblemDescr.Xmax(1)) = -1;
        ClassOfPoints(PointSet(:,2) > ProblemDescr.Xmax(2)) = -1;
    end
    
    if (size(PointSet,2) == 3)
        ClassOfPoints(PointSet(:,3) < ProblemDescr.Xmin(3)) = -1;
        ClassOfPoints(PointSet(:,3) > ProblemDescr.Xmax(3)) = -1;
    end
end
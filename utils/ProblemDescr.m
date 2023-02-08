% The reason for collecting all test problem-related data in a class is
% reproducibility. In many situations, one forgets to document all
% parameters and test settings properly such that after some time, this
% information is lost and one cannot reproduce the test and its
% results. In such a case, the test and all related effort was
% pointless.
% I decided to separate problem-related data and approximation-related
% data in order to perform parameter studies conveniently. In this
% case, all problem-related data remain unchanged whereas the
% approximation-related data may change. 

%hierarchy for input data:
% 1) data specified by PointSet and PointVec
% 2) if these are not present, consider input files (Excel file first,
%    then csv file)
% 3) if these are not present, consider test function as a pointer to a
% function which is then used to create test data

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
classdef ProblemDescr
    properties
        
        % string for saving some comments
        comments = '';
        
        % name for Excel input file (can be helpful to read in data from
        % actual measurements or simulations)
        InputFileExcel = '';

        % name for CSV input file (same reason)
        InputFileCSV = '';
        
        % function pointer to test function (NULL if not provided)
        Testfunc = '';
        
        % lower boundaries of the domain (necessary if a function pointer
        % only is provided). For explicitly given point sets, Xmin is
        % ignored.
        Xmin = -42;
        
        % upper boundaries of the domain (necessary if a function pointer
        % only is provided). For explicitly given point sets, Xmax is
        % ignored.
        Xmax = -42;
        
        % domain for RBF approximation
        % if no VisPolygon is given, take either the hypercube [Xmin, Xmax]
        % if no point set is provided or the bounding box of the point set
        % if such is provided
        DomainPolygon = -42;
        
        % file name of VTU output file (no output is written if not
        % specified)
        OutputFileVTU = '';

        % name of data field inside VTU
        NameOfDataInVTU = 'func';

        % vector with coordinates of the data points (format: npoints x
        % ndim)
        PointVec = -42;
        
        % function values at the data points (format: npoints x 1)
        PointSet = -42;
        
        % scaling values at the data points
        ScaleVec = -42
        
        % description structure for external models (for maximal
        % flexibility, we just incorporate a structure)
        FunctionParameters = '';

        % if true, more internal information is logged
        extendedStats = false;

        % if true, more internal information is displayed during the run
        verboseMode;
        
    end
end
% Export2VTU exports pointwise data to a vtu-file which can be read in by
% e.g. Paraview.
% The routine works in 2D or 3D. In 2D, triangular, quadrilateral and
% (partially) polygonal meshes are supported, in 3D tetrahedral and
% hexehadral meshes for volumes and triangular meshes for surfaces.
% The mesh may be unstructured. Mixed meshes are supported.
% 
% Input:
% - PointSet(numberOfDataPoints, ndim): numberOfDataPoints stands for
%   the number of data points and ndim \in {2,3} for the dimension of the
%   data.
% - triang(numberOfCells, maxPointsPerCell+1): represents the triangulation
%   corresponding to PointSet. The indices of the points are 0-based; the
%   last column in this array contains the number of points of the current
%   cell.
% - PointData(numberOfDataPoints, numberOfDataSets): numberOfDataSets
%   represents the number of data sets in the vtu-file
% - PointDataNames = {'NameOfSet1'; 'NameOfSet2', ...}: cell array
%   containing the names of the data sets which appear in the vtu file
% - filename: the file name of the vtu file. If the ending .vtu is missing,
%   it is automatically appended.
% - VecData(numberOfDataPoints, 3, numberOfDataSets): data sets for
%   vector-valued quantities (optional). Paraview requires 3D vectors even
%   if the vector field is 2D.
% - VecDataNames = {'NameOfSet1'; 'NameOfSet2', ...}: names of the
%   optional vector-valued quantities (optional)
% - filename: name of the vtu-file to write
%
% Output:
% <no direct output>; vtu-file written on working directory
%

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function Export2VTU(PointSet, triang, PointData, PointDataNames, ...
                    VecData, VecDataNames, filename)

    if (nargin == 5)
        filename = VecData;
        numberOfVecDataSets = 0;

    elseif (nargin < 5)
        error('Export2VTU: At least 5 arguments needed.')

    elseif (nargin == 7)
        % determine the number of vector data sets to write
        dimVecDataSet = size(size(VecData),2);

        if (dimVecDataSet ==2)
            numberOfVecDataSets = 1;
        else
            numberOfVecDataSets = size(VecData,dimVecDataSet);
        end

        cellVecDataNames = cellstr(VecDataNames);

        aux = size(cellVecDataNames,1);

        if (aux ~= numberOfVecDataSets)
            error(strcat('The number of vector data sets (', ...
                  strcat(int2str(numberOfVecDataSets), ...
                ') does not match the number of names for data sets.')))
        end
    end

    % determine the dimension of data space
    [~,ndim] = size(PointSet);

    % determine the number of point data sets to write
    [~, numberOfDataSets] = size(PointData);

    cellPointDataNames = cellstr(PointDataNames);

    [aux, ~] = size(cellPointDataNames);

    if (aux ~= numberOfDataSets)
        error(strcat('The number of point data sets (', ...
              strcat(int2str(numberOfDataSets), ...
              ') does not match the number of names for data sets.')))
    end

    % number of vertices
    numberOfDataPoints = length(PointData(:,1));

    % number of cells in the triangulation
    [numberOfCells, aux] = size(triang);

    % points per cell in triangulation
    PointsPerCell = triang(:,aux);
    maxPointsPerCell = aux-1;

    % open file for writing and append file extension if missing
    if ~strcmp(filename(end-3:end), '.vtu')
        filename = [filename '.vtu'];
    end
    file = fopen(filename, 'wt');

    if (ndim == 2)
        zcomponent = zeros(numberOfDataPoints,1);
    else
        zcomponent = PointSet(:,3);
    end

    % write header of the vtu file
    fprintf(file, '<?xml version="1.0"?>\n');
    aux = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n';
    fprintf(file, aux);

    fprintf(file, '  <UnstructuredGrid>\n');
    fprintf(file, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n', ...
            numberOfDataPoints, numberOfCells);

    % data points
    fprintf(file, '      <Points>\n');
    fprintf(file, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');

    auxArr = [PointSet(:,1) PointSet(:,2) zcomponent]';
    fprintf(file, '          %.4e %.4e %.4e\n', auxArr);

    fprintf(file, '        </DataArray>\n');
    fprintf(file, '      </Points>\n');

    % cell definitions
    fprintf(file, '      <Cells>\n');
    fprintf(file, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');

    Polygons = cell(maxPointsPerCell,1);
    Offsets = zeros(numberOfCells,1);
    PolyTypes = zeros(numberOfCells,1);
    iidxStart = 1;
    absOffset = 0;
    type = 0;

    formatSpecStart = '           %d';
    for ipointsPerCell = 1: maxPointsPerCell
        Polygons{ipointsPerCell} = triang(PointsPerCell == ipointsPerCell, 1:ipointsPerCell);
        numberCellsCurrentType = size(Polygons{ipointsPerCell},1);
        iidxEnd = iidxStart + numberCellsCurrentType - 1;
        Offsets(iidxStart:iidxEnd) = absOffset + (1:numberCellsCurrentType)*ipointsPerCell;
        absOffset = Offsets(max(iidxEnd,1));

        formatSpec = strcat(formatSpecStart, '\n');
        if (ndim ==2)
            % VTK_Vertex
            if(ipointsPerCell == 1)
                type = 1;
            % VTK_Line
            elseif(ipointsPerCell == 2)
                type = 3;
            % VTK_Triangle
            elseif(ipointsPerCell == 3)
                type = 5;
            % VTK_Quad
            elseif(ipointsPerCell == 4)
                type = 9;
            % VTK_Polygon
            else
                type = 7;
            end
        else
            % VTK_Vertex
            if(ipointsPerCell == 1)
                type = 1;
            % VTK_Line
            elseif(ipointsPerCell == 2)
                type = 3;
            % VTK_Triangle (for surface meshes)
            elseif(ipointsPerCell == 3)
                type = 5;
            % VTK_Tetra
            elseif (ipointsPerCell == 4)
                type = 10;
            % VTK_Hexahedron
            elseif (ipointsPerCell == 8)
                type = 12;
            % VTK_Polygon
            else
                type = 7;
            end
        end

        PolyTypes(iidxStart:iidxEnd) = type*ones(numberCellsCurrentType, 1);
        iidxStart = iidxEnd+1;

        fprintf(file, formatSpec, Polygons{ipointsPerCell}');
        formatSpecStart = strcat(formatSpecStart, ' %d');
    end

    fprintf(file, '        </DataArray>\n');

    % offsets
    fprintf(file, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(file,'           %d\n', Offsets);
    fprintf(file, '        </DataArray>\n');

    % cell types
    fprintf(file, '        <DataArray type="UInt8" Name="types" format="ascii">\n');
    fprintf(file,'           %i\n', PolyTypes);
    fprintf(file, '        </DataArray>\n');

    fprintf(file, '      </Cells>\n');

    % point data
    fprintf(file, '      <PointData>\n'); % def of std value
    for i = 1:numberOfDataSets
        fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n', char(cellPointDataNames(i)));
        fprintf(file, '          %.6e\n', PointData(:,i));
        fprintf(file, '        </DataArray>\n');
    end

    for i = 1:numberOfVecDataSets
        fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="ascii">\n', char(cellVecDataNames(i)));
            fprintf(file, '          %.4e %.4e %.4e\n', VecData(:,:,i)');
        fprintf(file, '        </DataArray>\n');
    end

    fprintf(file, '      </PointData>\n');

    % footer section
    fprintf(file, '    </Piece>\n');
    fprintf(file, '  </UnstructuredGrid>\n');
    fprintf(file, '</VTKFile>\n');

    fclose(file);
end
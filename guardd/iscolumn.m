%ISCOLUMN True if input is a column vector.
%   ISCOLUMN(V) returns logical 1 (true) if SIZE(V) returns [n 1] 
%   with a nonnegative integer value n, and logical 0 (false) otherwise.
%
%   See also ISSCALAR, ISVECTOR, ISROW, ISMATRIX, ISNUMERIC, ISLOGICAL, ISCHAR, 
%            ISEMPTY, SIZE.

%   Copyright 2009-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/07/23 16:04:52 $
%   Built-in function.
%
% ALSO WRITTEN AS EXPLICIT FUNCTION BY IAN KLECKNER
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License

function isColumn = iscolumn(V)
    % Get size of V
    [nRow nCol] = size(V);
    
    % Check if size is [n 1]
    isColumn = ( nCol == 1 && nRow > 0 );
end

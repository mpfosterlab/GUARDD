function multiline_string = cellArrayToMultilineString( cell_array_string  )
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/05/22 Start code
%
% Convert cell array to multi-line string
%
%  This is useful for reading / displaying text to a MATLAB GUI edit box
%  with multi-line input

% Convert each line of text into a cell array
Nlines = size(cell_array_string,1);
multiline_string = '';
for l = 1:Nlines
    line = cell_array_string{l};
    line = line{1};
    multiline_string = strvcat( multiline_string, line );
end

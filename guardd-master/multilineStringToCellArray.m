function cell_array_string = multilineStringToCellArray( multiline_string )
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/05/22 Start code
%
% Convert multi-line string into a cell array
%  with each array element as a distinct line of text
%
%  This is useful for reading / displaying text to a MATLAB GUI edit box
%  with multi-line input

% Convert each line of text into a cell array
Nlines      = size(multiline_string,1);
cell_array_string = cell(Nlines,1);
for l = 1:Nlines
    %if( Nlines > 1 )
    %    cell_array_string{l} = multiline_string{l};
    %else
    %    cell_array_string{l} = multiline_string;
    %end
    % Remove trailing spaces with deblank
    cell_array_string{l} = deblank( multiline_string(l,:) );
end
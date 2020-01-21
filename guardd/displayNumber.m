function s = displayNumber( x, varargin )
% Return a string corresponding to a formatted number OR a "?" (if NaN)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/05/01 Start coding
% 
% OPTIONS
%  Can also supply format for number via varargin 
%  Ex: formattedNumber = displayNumber( unformattedNumber, '%0.3f' )
%

if( isnan(x) || isempty(x) )
    s = '?';
else
    % If the user wants to format the number (e.g., %0.1f, %0.2f, %03d)
    if( nargin == 2 )
        eval( sprintf('s=sprintf(\''%s\'',x);',varargin{1}) );
    else
        s = sprintf('%0.1f',x);
    end
end

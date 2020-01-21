classdef Note < handle
    % Note class: stores cell array for multiline notes, timestamp,  name
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/05/22 Start coding
    %
        
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        
        name        = 'DefaultName';    % Name of note
        text        = {};               % Contents of note
        time        = NaN;              % Time at creation of note
        
    end
    
    methods
        
        function obj = Note(varargin)
            %% Default constructor (optional arguments for text, name)
            obj.time = now;
            
            if( nargin >= 1 )
                % Set the text
                obj.setText( varargin{1} );
            end
            
            if( nargin >= 2 )
                obj.name = varargin{2};
            end
        end        
        
        function text_oneliner = getTextOneLiner(obj)
            %% Format text for one line output
            text_oneliner = '';
            Nlines = max(size(group.note.text));
            for l = 1:Nlines
                if( l == 1 )
                    text_oneliner = sprintf('%s', group.note.text{l});
                else
                    text_oneliner = sprintf('%s [CR] %s', text_oneliner, group.note.text{l});
                end
            end
        end
        
        function setName(obj, name)
            %% Set the name
            obj.name = name;
        end
        
        function setText(obj, text)
            %% Set the text (make sure it ends up as a CELL of STRINGS
            obj.text = {};
            
            if( iscell(text) )
                Nlines = size(text,1);
                for l = 1:Nlines
                    line = text{l};
                    
                    % If it is a cell in a cell, then extract the string
                    % inside
                    if( iscell(line) )
                        line = line{1};
                    end
                    
                    obj.text{l} = line;
                end                        
                
            else
                obj.text{1} = text;
            end
        end
        
        function setNameAndText(obj, name, text)
            %% Set the name and the text
            obj.name = name;
            
            obj.setText( text );
        end            
    end
end
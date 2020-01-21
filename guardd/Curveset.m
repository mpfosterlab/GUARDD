classdef Curveset < handle
    % Curveset class: for a NMR probe (index+atom), which includes 1+ curve
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/01/17 Start coding
    %
    % Properties and methods for a single RD curveset
    
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        
        name        = 'CurvesetName';   % Name of curve set        
        index       = 0;        % Residue number
        atom        = '?';      % Name of atom (N, H\alpha, C\delta_1, etc.)
        residue     = '???';    % Name of residue (Ile, Leu, Arg, etc.)        
        
        notes       = '';       % Notes about the group (e.g., from fitting)
        
        Nc          = 0;        % Number of curves in the set
        curves      = {};       % Array of curves in the set
        
        % Cell matrix for parameters to describe this class via table
        %  Used in generateParamTable()        
        param_info = {...
            'name:(String) Name of curve (e.g., Leu 22\delta_1)', ...      
            'index:(Number) Residue number (e.g., 10, 26, 55)', ...  
            'atom:(String) Name of atom (e.g., N, H\alpha, C\delta_1)', ...  
            'residue:(String) Name of residue (e.g., Ile, Leu, Arg)' ...
                      };
    end
    
    methods
        function addCurve( obj, curve )
            %% Add the new curve to the end of the list
            % Since Curveset is a subclass of handle, the added curve
            % points to its place in memory, and subsequent changes to
            % either source curve or obj.curves{obj.Nc+1} are linked
            obj.curves{obj.Nc+1} = curve;
            obj.Nc = obj.Nc+1;
        end
        
        function curve_in_set = containsCurve( obj, curve )
            %% (Boolean) The supplied curve is in the set
                        
            % Check each curve in this set for a match
            curve_in_set = false;
            for c = 1:obj.Nc
                % If curve is found, then return early
                if( curve == obj.curves{c} )
                    curve_in_set = true;
                    return
                end
            end
        end
        
        function curveset_copy = copy( obj )
            %% Create a new curveset with all properties of the current curveset            
            % DO NOT COPY the curves (these are basic data, and linked to dataset)
            % User can now change which curves are part of set
            curveset_copy              = Curveset;
            
            curveset_copy.index        = obj.index;
            curveset_copy.atom         = obj.atom;
            curveset_copy.residue      = obj.residue;
            curveset_copy.name         = sprintf('Cp(%s)',obj.name);            
            curveset_copy.notes        = obj.notes;
            
            % Curves are not NEW, but still linked to original
            %  If original curves are changed, so are those in the copy curveset
            curveset_copy.Nc            = obj.Nc;      
            curveset_copy.curves        = obj.curves;
            
            curveset_copy.param_info    = obj.param_info;
        end
        
        function isequal = eq(cs1, cs2)
            %% (Boolean) The two curves equal in all properties (Used as "cs1==cs2")
            
            isequal = cs1.index == cs2.index && ...
                      strcmp(cs1.atom, cs2.atom) && ...
                      strcmp(cs1.residue, cs2.residue) && ...
                      strcmp(cs1.name, cs2.name) && ...
                      cs1.Nc == cs2.Nc;
                  
            % Compare each curve in the set too, if necessary
            if( isequal )
                for c = 1:length(cs1.curves)
                    % If curves are not equal, then return early
                    if( ~(cs1.curves{c} == cs2.curves{c}) )
                        isequal = false;
                        return
                    end
                end
            end            
        end
        
        function curve_array = getCurves(obj)           
            %% Create a curve for each RD curve in the dataset
            curve_array = {};            
        end   
        
                
        function [isMQ, cFirstMQ] = isMQ( obj )
            %% (Boolean) The curveset contains a MQ curve
            
            isMQ        = false;
            cFirstMQ    = NaN;
            for c = 1:obj.Nc
                % If there is one MQ curve, then the whole set is MQ, and return early
                if( ~obj.curves{c}.SQX )
                    isMQ = true;
                    cFirstMQ = c;
                    return
                end
            end
        end
        
        function outputSpecs(obj, FILE, varagin)
            %% Output groups to a readable file
            % FILE is an already opened handle, or =1 for command window
            % First optional argument is TAB, which specifies a string to
            % use before each newline (e.g., '  ' or '    ')
            
            if( nargin == 3 )
                while( iscell(varagin) )
                    varagin = varagin{1};
                end
                TAB = varagin;
            else
                TAB = '';
            end
            
            fprintf(FILE, '\n%s%s\t%s', TAB,     'Name', obj.name);
            fprintf(FILE, '\n%s%s\t%d', TAB,     'Index', obj.index);
            fprintf(FILE, '\n%s%s\t%s', TAB,     'Atom', obj.atom);
            fprintf(FILE, '\n%s%s\t%s', TAB,     'Residue', obj.residue);
            fprintf(FILE, '\n%s%s\t%d', TAB,     'NumCurvesets', obj.Nc);
            
            for c = 1:obj.Nc
                fprintf(FILE, '\n\n%sCurve %d/%d', TAB, c, obj.Nc);
                
                curve = obj.curves{c};
                curve.outputSpecs(FILE, '      ');
            end
        end
        
        function removeCurve( obj, curve )
            %% Find the curve, then remove it (but not from the dataset)
            for c = 1:obj.Nc
                if( curve == obj.curves{c} )                    
                    % Here, "c" is the element to be eliminated
                    %  E.g., Remove eement #3 from a set of 5 elements
                    %  Set #3=#4, #4=#5, and #5=[], then Ntot=Ntot-1
                    for c1 = c:obj.Nc-1
                        obj.curves{c1} = obj.curves{c1+1};
                    end
                    
                    % Set final element to null
                    obj.curves(obj.Nc) = [];
                    
                    % Update the number of elements
                    obj.Nc = obj.Nc-1;
                    return                    
                end
            end            
        end
        
        function saveDataFromTable( obj, table_data )
            %% Set parameters using values from table            
            % Save each property, row by row
            %  Column 1 contains property name
            %  Column 2 contains value entered into the table
            for row = 1:size(table_data,1)
                propety_name = table_data{row,1};
                eval( sprintf('obj.%s = table_data{row,2};', propety_name) );                
            end
        end
        
        function setAssignment( obj, index, atom, residue )
            %% Set the specs for the assignment of NMR signal (and default name)
            obj.index = index;
            obj.atom = atom;
            obj.residue = residue;
            
            % E.g., Leu 12\delta_2
            obj.name = sprintf('%s %d%s', residue, index, atom);
            
            % Set if for each curve as well
            for c = 1:obj.Nc
               obj.curves{c}.setAssignment( index, atom, residue ); 
            end
        end
        
        function setAtom_DoNotChangeCurves( obj, atom )
            %% Set the atom and name for the curveset but do not change the curves within
            % This was used to swap stereoassignments 1<->2 in TRAP data
            obj.atom = atom;
            obj.name = sprintf('%s %d%s', obj.residue, obj.index, obj.atom);            
        end
        
         function obj = setName( obj, name )
            %% Set a non-default name
            obj.name = name;
         end
         
         function obj = setNotes( obj, notes )
            %% Set the notes for the group
            obj.notes = notes;
        end 
    end
end


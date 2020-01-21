classdef Curve < handle
    % curve class: for a single RD curve
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/01/17 Start coding
    %
    % Properties and methods for a single RD curve
    
    properties (SetAccess = private)        
        name        = 'CurveName'; % Name of curve (Leu 12\delta_2)
        index       = 0;        % Residue number
        atom        = '?';      % Name of atom (N, H\alpha, C\delta_1, etc.)
        residue     = '???';    % Name of residue (Ile, Leu, Arg, etc.)        
        notes       = '';       % Notes about the group (e.g., from fitting)
        
        Nobs        = 0;        % Number of observations (length(vcpmg))
        vcpmg       = [];       % Array of vcpmg values (Hz)
        R2eff       = [];       % Array of R2eff values (Hz)
        eR2eff      = [];       % Array of errors in R2eff (Hz)        
               
        B0          = 0;        % Magnetic field strength (MHz)
        Temp        = 0;        % Temperature (K)
        Temp_C      = 0;        % Temperature (C)
        TCPMG       = 0;        % Total CPMG time in pulse sequence (sec)
        SQX         = false;    % (boolean) data set is SQ for nucleus X
        dataset     = [];       % Data set number associated with this curve
        
        AX_String       = 'AX';     % 13C or 15N
        gammaX_relative = NaN;      % Defined using AX_String
        
        % Cell matrix for parameters to describe this class via table
        %  Used in generateParamTable()        
        param_info = {...
            'name:(String) Name of curve (e.g., Leu 22\delta_1)', ...
            'index:(Number) Residue number (e.g., 10, 26, 55)', ...  
            'atom:(String) Name of atom (e.g., N, H\alpha, C\delta_1)', ...  
            'residue:(String) Name of residue (e.g., Ile, Leu, Arg)' ...
            ...            
            'B0:(MHZ) Magnetic field strength', ...
            'Temp_C:(C) Temperature', ...
            'TCPMG:(sec) Total CPMG time in pulse sequence', ...
            'SQX:(Boolean) Data set is SQ for nucleus X' ...
                      };
    end
    
    methods        
        
        function ccopy = copy( obj )
            %% Create a new curve with all properties of the current curve            
            ccopy           = Curve;
            ccopy.index     = obj.index;
            ccopy.atom      = obj.atom;
            ccopy.residue   = obj.residue;
            ccopy.name      = obj.name;
            ccopy.notes     = obj.notes;
            
            ccopy.Nobs      = obj.Nobs;
            ccopy.vcpmg     = obj.vcpmg;
            ccopy.R2eff     = obj.R2eff;
            ccopy.eR2eff    = obj.eR2eff;
            
            ccopy.B0        = obj.B0;
            ccopy.Temp      = obj.Temp;
            ccopy.Temp_C    = obj.Temp_C;
            ccopy.TCPMG     = obj.TCPMG;
            ccopy.SQX       = obj.SQX;
            ccopy.dataset   = obj.dataset;    
            
            ccopy.AX_String         = obj.AX_String;
            ccopy.gammaX_relative   = obj.gammaX_relative;
            ccopy.param_info        = obj.param_info;
        end
        
        function isequal = eq(curve1, curve2)
            %% (Boolean) The two curves equal in all properties (Used as "curve1==curve2")
            
            isequal = curve1.index == curve2.index && ...
                      strcmp(curve1.atom, curve2.atom) && ...
                      strcmp(curve1.residue, curve2.residue) && ...
                      curve1.Nobs == curve2.Nobs && ...
                      sum(curve1.vcpmg == curve2.vcpmg)==length(curve1.vcpmg) && ...
                      sum(curve1.R2eff == curve2.R2eff)==length(curve1.R2eff) && ...
                      sum(curve1.eR2eff == curve2.eR2eff)==length(curve1.eR2eff) && ...
                      curve1.B0 == curve2.B0 && ...
                      curve1.Temp == curve2.Temp && ...
                      curve1.TCPMG == curve2.TCPMG && ...
                      curve1.SQX == curve2.SQX;
        end
        
        function exportData(obj, FILE)
            %% Export curves in dataset to a file
            % FILE = already opened file handle (1 => command window)
                        
            fprintf(FILE, '\n%s\t%d', 'INDEX',   obj.index);
            fprintf(FILE, '\n%s\t%s', 'ATOM',    obj.atom);
            fprintf(FILE, '\n%s\t%s', 'RESIDUE',    obj.residue);
            
            fprintf(FILE, '\n%s\t%s\t%s\t%s', 'OBS', 'VCPMG', 'R2', 'ERROR');
            for o = 1:obj.Nobs
                fprintf(FILE, '\n%d\t%f\t%f\t%f', o, obj.vcpmg(o), obj.R2eff(o), obj.eR2eff(o));
            end            
            fprintf(FILE, '\n%s', 'ADDDATA');
        end
        
        function group_array = getGroups()
            % Return all the groups this curve is associated with
            
            group_array = {};
            
        end
        
        function specs_string = getSpecsString(obj)
            %% Return a string of specifications
            if( obj.SQX )
                SQ_string = 'SQ';
            else
                SQ_string = 'MQ';
            end
            specs_string = sprintf('%d-%s-%dC', round(obj.B0), SQ_string, obj.Temp_C);            
        end
        
        function fromSameProbe = isFromSameProbe(curve1, curve2)
            %% (Boolean) The curves result from the same NMR probe (same index and atom)
            fromSameProbe = curve1.index == curve2.index && ...
                            strcmp(curve1.atom, curve2.atom);
            
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
            fprintf(FILE, '\n%s%s\t%f\t%s', TAB, 'B0',    obj.B0,    'MHz');
            fprintf(FILE, '\n%s%s\t%f\t%s', TAB, 'Temp',    obj.Temp_C,    'C');
            fprintf(FILE, '\n%s%s\t%f\t%s', TAB, 'TCPMG',    obj.TCPMG,    'sec');
            fprintf(FILE, '\n%s%s\t%d',     TAB, 'SQX',    obj.SQX');
            fprintf(FILE, '\n%s%s\t%s',     TAB, 'AX',    obj.AX_String);            
            fprintf(FILE, '\n%s%s\t%d',     TAB, 'Nobs', obj.Nobs);
            
            % Output all the observations
            fprintf(FILE, '\n%s%s\t%s\t%s\t%s',     TAB, 'OBS', 'vCPMG(Hz)', 'R2Eff(Hz)', 'eR2Eff(Hz)');
            for o = 1:obj.Nobs
                fprintf(FILE, '\n%s%d\t%f\t%f\t%f',     TAB, o, obj.vcpmg(o), obj.R2eff(o), obj.eR2eff(o));
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
            obj.updateName();
            % Convert Celcius to Kelvin
            obj.Temp = obj.Temp_C + 273;
        end
        
        function setAssignment( obj, index, atom, residue )
            %% Set the specs for the assignment of NMR signal
            obj.index = index;
            obj.atom = atom;
            obj.residue = residue;
            % E.g., Leu 12\delta_2
            obj.updateName();
        end
        
        function obj = setConditions( obj, B0, Temp, TCPMG, SQX, AX_String )            
            %% Set specs for experimental conditions
            obj.B0       = B0;
            obj.Temp     = Temp;
            obj.Temp_C   = Temp-273;
            obj.TCPMG    = TCPMG;
            obj.SQX      = SQX;
            
            obj.AX_String = AX_String;
            obj.gammaX_relative = Session.getGammaRelative(AX_String);           
        end
        
        function obj = setData( obj, vcpmg, R2eff, eR2eff )
            %% Set data vectors (option to set error vector too)
            data_are_valid = true;
            
            % Data vectors must be same length
            if( length(vcpmg) == length(R2eff) )
                
                obj.vcpmg   = vcpmg;
                obj.R2eff   = R2eff;
                obj.Nobs    = length(vcpmg);
                           
                % TODO % Permit function call without errors specified
                %if( narargin == 3 )
                    if( length(eR2eff) == length(vcpmg) )
                        obj.eR2eff  = eR2eff;
                    else
                        data_are_valid = false;
                    end
                %end                
                
            else
                data_are_valid = false;
            end
            
            if( ~data_are_valid )
                error('Must have same number of points in vCPMG and R2eff (and optional eR2eff)');
            end
        end
        
        function obj = setDataset( obj, dataset )
            %% Store a link to the dataset associated with this curve
            obj.dataset = dataset;            
        end
        
        function obj = setNotes( obj, notes )
            %% Set the notes for the group
            obj.notes = notes;
        end 
                
        function updateName( obj )
            %% Update the name of the residue based on residue, index, atom
            obj.name = sprintf('%s %d%s', obj.residue, obj.index, obj.atom);
        end
        
    end
end


classdef ParamDisplay < handle
    % ParamDisplay class: holds information for display of the parameters
    %  How many subplots, and what should be displayed on X and Y
    %  Used in in Cluster Display and sequence mapping
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/05/10 Start coding
    % 2011/05/20 Add histogram
    % 2011/09/05 Add PhiexX
    
    properties (SetAccess = private)
        %% SetAccess = private properties can be read via ".", but not set
        
        % The function of this class is to store an array of X,Y plots for
        % displaying certain RD parameters
        %
        % There are certain PLOT TYPES which specify X and Y parameters
        %  Arrhenius    -> X = Eab  Y = Eba
        %  van't Hoff   -> X = dH   Y = dS
        %  kA,kB        -> X = kA   Y = kB, each at some temperature
        %  Custom       -> X = (whatev)     Y = (whatev)
        
        % For the custom type, X and Y must be selected from below
        % These are PARAMETER TYPES
        % dwH
        % dwX
        % Pa      -> pick T
        % kex     -> pick T
        % R20     -> pick T -> pick B0
        % Rex     -> pick T -> pick B0
        % kA      -> pick T
        % kB      -> pick T
        % K       -> pick T
        % dH
        % dS
        % Ea_AB
        % Ea_BA
        param_type_string = {...
                            'Residue', ...
                            'dwH_ppm',  ...
                            'dwX_ppm',  ...
                            'Pa',   ...
                            'kex',  ...
                            'R20',  ...
                            'Rex',  ...
                            'alpha', ...
                            'PhiexX', ...
                            'kA',   ....
                            'kB',   ...
                            'K',    ...
                            'dH',   ...
                            'dS',   ...
                            'Ea(A->B)', ...
                            'Ea(B->A)'  ...
                            };
                        
        % Each SUBPLOT therefore specifies (1) PLOT TYPE, (2) PARAMETER
        % TYPE for X, (3) PARAMETER TYPE for Y
        
        % Subplots are displayed in rows and columns
        Nsubplots           = 0;
        subplot_rows        = 0;
        subplot_cols        = 0;
        
        % Holds the plot type (number referring to one of the strings above)
        %  for each of the Nsubplots
        % What type of plot number to make? (1=> Arrhenius, 2=>van't Hoff, 3=>kA,kB, etc.)
        plot_type_num_Array     = [];
        
        % For the X-axis, what is the parameter number desired? (1st, 2nd, 3rd, etc.)
        param_num_X_Array    = [];
        param_num_Y_Array    = [];
        
        % For the X-axis parameter, what is the temperature number desired? (1st, 2nd, 3rd, etc.)        
        Temp_num_X_Array     = [];
        Temp_num_Y_Array     = []; 
        
        % For the X-axis parameter, what is the B0 field strength number desired? (1st, 2nd, 3rd, etc.)
        B0_num_X_Array       = [];
        B0_num_Y_Array       = [];    
        
        % (Bool) display a histogram for this plot
        displayHistogram_Array = [];        

        % Parent session that owns this group
        parentSession = [];
    %end
    
    %properties (SetAccess = private, GetAccess = private)
        % This should be calculated for each accession
        plot_type_string = {};        
        
        % Actual B0 values in the session (unique, and sorted)
        B0_values   = [];
        NB0         = 0;
        
        % Actual tempererature values (unique and sorted)
        Temp_values = [];
        NTemp       = 0;
    end
    
    methods
        function obj = ParamDisplay( parentSession )
            %% Default constructor requires parent session
            obj.parentSession = parentSession;
        end
        
        function plot_type_string = getPlotTypeString( obj )
            %% Get the list of plot types (depends on number of temperatures )            
            obj.updateParameters();
            
            plot_type_string = obj.plot_type_string;                      
        end
        
        function uniqueB0_String = getUniqueB0String( obj )
            %% Return string array of unique temperature from session
            obj.updateParameters();

            uniqueB0_String{obj.NB0} = '';
            for b = 1:obj.NB0
                uniqueB0_String{b} = sprintf('%0.1f MHz', obj.B0_values(b));
            end
        end
        
        function uniqueTemp_String = getUniqueTempertureString( obj )
            %% Return string array of unique temperature from session
            obj.updateParameters();

            uniqueTemp_String{obj.NTemp} = '';
            for t = 1:obj.NTemp
                uniqueTemp_String{t} = sprintf('%dC', obj.Temp_values(t)-273);
            end
        end
        
        function value = getValue( obj, subplot_selected, feature_string, XorY_string )
            %% Set the subplot type number for the subplot selected
            if( subplot_selected > obj.Nsubplots )
                error('Invalid subplot specified');
            end
                        
            switch( upper(feature_string) )
                case 'PLOT'
                    value = obj.plot_type_num_Array(subplot_selected);
                     
                case 'PARAM'
                    switch( upper(XorY_string) )
                        case 'X'
                            index   = obj.param_num_X_Array(subplot_selected);
                            value   = obj.param_type_string{index};
                        case 'Y'
                            index   = obj.param_num_Y_Array(subplot_selected);
                            value   = obj.param_type_string{index};
                        otherwise
                            error('For XorY_string, use X or Y only');
                    end                    
                    
                case 'TEMP'
                    switch( upper(XorY_string) )
                        case 'X'                            
                            index   = obj.Temp_num_X_Array(subplot_selected);
                            value   = obj.Temp_values(index);
                        case 'Y'
                            index   = obj.Temp_num_Y_Array(subplot_selected);
                            value   = obj.Temp_values(index);
                        otherwise
                            error('For XorY_string, use X or Y only');
                    end
                    
                case 'B0'
                    switch( upper(XorY_string) )
                        case 'X'
                            index   = obj.B0_num_X_Array(subplot_selected);
                            value   = obj.B0_values(index);
                        case 'Y'
                            index   = obj.B0_num_Y_Array(subplot_selected);
                            value   = obj.B0_values(index);
                        otherwise
                            error('For XorY_string, use X or Y only');
                    end
                    
                case {'HIST', 'HISTO', 'HISTOGRAM'}
                    value = obj.displayHistogram_Array(subplot_selected);
                    
                otherwise
                    error('Invalid feature_string provided. Check code for options');
            end
        end
        
        function setPlottingSpecs( obj, subplot_selected, feature_string, XorY_string, number )
            %% Set the subplot type number for the subplot selected
            if( subplot_selected > obj.Nsubplots )
                error('Invalid subplot specified');
            end
                        
            switch( upper(feature_string) )
                case 'PLOT'
                    if( number > length(obj.plot_type_string) )
                        error('Invalid plot type number specified');
                    end
                    
                    obj.plot_type_num_Array(subplot_selected)  = number;
                     
                case 'PARAM'
                    switch( upper(XorY_string) )
                        case 'X'
                            obj.param_num_X_Array(subplot_selected)     = number;
                        case 'Y'
                            obj.param_num_Y_Array(subplot_selected)     = number;
                        otherwise
                            error('For XorY_string, use X or Y only');
                    end                    
                    
                case 'TEMP'
                    switch( upper(XorY_string) )
                        case 'X'
                            obj.Temp_num_X_Array(subplot_selected)     = number;
                        case 'Y'
                            obj.Temp_num_Y_Array(subplot_selected)     = number;
                        otherwise
                            error('For XorY_string, use X or Y only');
                    end
                    
                case 'B0'
                    switch( upper(XorY_string) )
                        case 'X'
                            obj.B0_num_X_Array(subplot_selected)     = number;
                        case 'Y'
                            obj.B0_num_Y_Array(subplot_selected)     = number;
                        otherwise
                            error('For XorY_string, use X or Y only');
                    end
                    
                case {'HIST', 'HISTO', 'HISTOGRAM'}
                    obj.displayHistogram_Array(subplot_selected)     = number;
                    
                otherwise
                    error('Invalid feature_string provided. Check code for options');
            end
        end
        
        function updateNumSubplots(obj, subplot_rows, subplot_cols )
            %% Set the number of rows and columns in the subplot
            obj.subplot_rows    = subplot_rows;
            obj.subplot_cols    = subplot_cols;
            
            % The new number of subplots will affect the ones that are
            % currently displayed
            Nsubplots_New       = subplot_rows * subplot_cols;
            
            if( Nsubplots_New > obj.Nsubplots )
                % If there needs to be a new subplot, then add one more
                % element to each array                
                obj.plot_type_num_Array(end+1:Nsubplots_New)     = 1;
        
                obj.param_num_X_Array(end+1:Nsubplots_New)      = 1;
                obj.param_num_Y_Array(end+1:Nsubplots_New)      = 1;

                obj.B0_num_X_Array(end+1:Nsubplots_New)         = 1;
                obj.B0_num_Y_Array(end+1:Nsubplots_New)         = 1;

                obj.Temp_num_X_Array(end+1:Nsubplots_New)       = 1;
                obj.Temp_num_Y_Array(end+1:Nsubplots_New)       = 1;
                
                obj.displayHistogram_Array(end+1:Nsubplots_New) = false;
                
            else
                % Remove some subplots
                obj.plot_type_num_Array     = obj.plot_type_num_Array(1:Nsubplots_New);
        
                obj.param_num_X_Array       = obj.param_num_X_Array(1:Nsubplots_New);
                obj.param_num_Y_Array       = obj.param_num_Y_Array(1:Nsubplots_New);

                obj.B0_num_X_Array          = obj.B0_num_X_Array(1:Nsubplots_New);
                obj.B0_num_Y_Array          = obj.B0_num_Y_Array(1:Nsubplots_New);

                obj.Temp_num_X_Array        = obj.Temp_num_X_Array(1:Nsubplots_New);
                obj.Temp_num_Y_Array        = obj.Temp_num_Y_Array(1:Nsubplots_New);
                
                if( isempty(obj.displayHistogram_Array) )
                    obj.displayHistogram_Array(1:Nsubplots_New) = true;
                else
                    obj.displayHistogram_Array  = obj.displayHistogram_Array(1:Nsubplots_New);
                end
            end
            
            obj.Nsubplots = Nsubplots_New;
        end
        
        function updateParameters(obj)
            %% Update the parameters as obtained from parent session
            obj.B0_values   = sort( unique( obj.parentSession.getB0Array() ) );
            obj.NB0  = length( obj.B0_values );
            
            obj.Temp_values   = sort( unique( obj.parentSession.getTemperatureArray() ) );
            obj.NTemp  = length( obj.Temp_values );
            
            % Make the plot_type_string
            obj.plot_type_string    = {};
            obj.plot_type_string{1} = 'vant Hoff [dH vs dS]';
            obj.plot_type_string{2} = 'Arrhenius [Ea(AB) vs. Ea(BA)]';
            for t = 1:obj.NTemp
                obj.plot_type_string{t+2} = sprintf('Kinetic rate (%dC) [kA vs kB]',obj.Temp_values(t)-273);
            end
            obj.plot_type_string{obj.NTemp+3} = 'Custom';
        end        
    end
end
classdef Session < handle
    % Session class: The highest-order organizational structure of data (contains Datasets, and Groups)
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/01/28 Start coding
    % 2011/05/01 Load data from script file
    % 2011/06/07 Add settings GUI and param_info
    % 2011/09/11 MC randomization mode A or B (MC_RANDOMIZE_MODE)
    % 2017/04/25 Added 19F (IK)
    %
    % Properties and methods for data belonging to the entiring program
    
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        % TODO % Add molecule class to wrap this up (e.g., to compare RD
        % curves from multiple molecules / sessions)
        %
        % TODO % Copy curves, curvessets ...to -> Group, Curveset (RC menu)
        % TODO % Copy curve, curveset, group
        
        name        = '?';      % Name of session
        description = '?';        
        
        Nds         = 0;        % Number of data sets loaded
        datasets    = {};       % Array of datasets currently loaded
        
        Ng          = 0;        % Number of groups created
        groups      = {};       % Array of groups created
        
        Ndc             = 0;    % Number of DisplayClusters for visualizing results
        displayClusters = {};   % Array of DisplayClusters for visualizing results
        
        notes           = {};   % Notes on the session, using the Note class (all groups)
        %note_names      = {};   % Name of each note
        %note_times      = [];   % Time of each note
        Nnotes          = 0;
        
        sequence_array = {};    % Assignments for each residue in the molecule
        
        START_DATE              = '';   % String for when the Session was created
        
        %% Display structures        
        % Holds parameter display information for the Cluster Display GUI
        paramDisplay = [];
                
        %% Settings        
        
        % Cell matrix for parameters to describe this class via table
        %  Used in generateParamTable()
        param_info = {...
            'outputDir:(String) Name of output directory', ...
            ...
            'Nmc:(Numeric) Number of iterations in Monte Carlo error analysis [100]', ...
            'MC_RANDOMIZE_MODE:(String) MC randomization using (A) fitted residuals, or (B) experimental errors', ...
            'MIN_F_ERROR:(Numeric) Minimum fractional error enforced in data [0.02]', ...
            ...
            'OUTPUT_DEBUG_UPDATE_FIT_PARAMS:(Boolean) Assign independent parameters and dependent scaling facotrs for group fit' ...
            'OUTPUT_DEBUG_ERRORS:(Boolean) Monte Carlo error estimation' ...
            'OUTPUT_DEBUG_CALCULATE_ALPHA:(Boolean) Calculating exchange timescale alpha' ...
            'OUTPUT_DEBUG_GENERATE_GROUPS:(Boolean) Generating NEW groups via index/atom (Session class)' ...
            'OUTPUT_DEBUG_RATE_ANALYSIS:(Boolean) Rate analysis (Arrhenius and vant Hoff)' ...            
            'OUTPUT_DEBUG_LOAD_SCRIPT:(Boolean) Executing script file (Session class)' ...
            'OUTPUT_DEBUG_ADD_DATA:(Boolean) Data are added to a dataset (Dataset class)' ...            
            'OUTPUT_DEBUG_ENFORCE_MIN_ERROR:(Boolean) Enforcing minimum error in data (Dataset class)' ...
            'OUTPUT_DEBUG_READ_NLIN_FILE:(Boolean) Reading NLIN.TAB file (Dataset class)' ...
            'OUTPUT_DEBUG_SET_SEQUENCE:(Boolean) Setting residue names to each index (Dataset class)' ...
            'OUTPUT_DEBUG_CALCULATE_R2EFF:(Boolean) Calculating R2Eff from Intensity and vCPMG (Dataset class)' ...            
            'DRIVE_USER_INSANE:(Boolean) Make the user frustrated beyond his/her mental capacity' ...
            ...
            'FONTSIZE_EXTRA_SMALL:(Numeric) Font size [10]', ...
            'FONTSIZE_SMALL:(Numeric) Font size [12]', ...
            'FONTSIZE_MEDIUM:(Numeric) Font size [14]', ...
            'FONTSIZE_LARGE:(Numeric) Font size [16]', ...
            'FONTNAME:(String) Name of font (Arial, Helvetica)', ...
            ...
            'MARKERSIZE:(Numeric) Size for marker in plots', ...
            'LINEWIDTH:(Numeric) Default linewidth in plots', ...
            'LINEWIDTH_SMALL:(Numeric) Default small linewidth in plots', ...
            ...
            'F_TEXT_LABEL:(Numeric) Analysis section: Fraction of plot range to move text label (next to data point) [0.01]', ...
            'F_PLOT_RANGE:(Numeric) Analysis section: The fraction of maximum data to which plot limits are extended [0.02]' ...
                      };
                  
        % Output directory
        outputDir                   = '';
        
        %% Plotting
        % Fonts and plots
        FONTSIZE_EXTRA_SMALL   = 10;
        FONTSIZE_SMALL         = 12;
        FONTSIZE_MEDIUM        = 14;
        FONTSIZE_LARGE         = 16;
        FONTNAME               = 'Arial';

        MARKERSIZE             = 8;
        LINEWIDTH              = 2;
        LINEWIDTH_SMALL        = 1;
        
        % Colors and symbols used for plotting
        plotColorOrder = [...
             0         0    1.0000; ...
             0    0.5000         0; ...
        1.0000         0         0; ...
             0    0.7500    0.7500; ...
        0.7500         0    0.7500; ...
        0.7500    0.7500         0; ...
        0.2500    0.2500    0.2500 ...
                ];

        plotSymbolOrder = {'o', 's', 'd', '^', 'v', '<', '>'};

        %% Calculations
        % Statistical significance level for model testing [0.01]
        ALPHA                  = 0.01;

        % Number of Monte-carlo simulation steps for error estimates [100]
        Nmc                  = 100;
        
        % MC randomization using (A) fitted residuals, or (B) experimental errors
        MC_RANDOMIZE_MODE = 'A';

        % Minimum fractional error enforced in data [0.01]
        MIN_F_ERROR            = 1e-5;

        % Scales the maximum random adjustment in initial guess in generating
        % initial fitting conditions (both grid search and error est.) [0.1] (12/10/09)
        F_RAND_MC_GUESS      = 0.1;

        % Top fraction of grid points used to display SSE(p)-space [0.5]
        % Can adjust such that top fraction is >100 points
        chi2TopFraction   = 0.5;
        
        % Parameters to display in chi2map
        pChi2Map             = [1 2 3 4 5 6];

        % Top fraction of fits used for fitting SSE(p) in bins [0.5]
        % Can adjust such that top fraction is >100 points
        SSE_FIT_TOP_FRACTION   = 0.5;

        % Number of parameter bins used for fitting SSE(p) [20]
        SSE_FIT_N_BINS         = 20;

        % Analysis section: The fraction of maximum data to which plot limits are extended [0.02]
        F_PLOT_RANGE           = 0.02;

        % Analysis section: Fraction of plot range to move text label (next to data point) [0.01]
        F_TEXT_LABEL           = 0.01;

        %% Debugging output flags
        
        % Assign independent parameters and dependent scaling facotrs for group fit
        OUTPUT_DEBUG_UPDATE_FIT_PARAMS = true;
        
        % During Monte Carlo error estimation
        OUTPUT_DEBUG_ERRORS        = false;    
        
        % Calculating exchange timescale alpha
        OUTPUT_DEBUG_CALCULATE_ALPHA         = false;    
        
        % Output figure during rate analysis
        OUTPUT_DEBUG_RATE_ANALYSIS = false;    
        
        % Data are added to a dataset (Dataset class)
        OUTPUT_DEBUG_ADD_DATA = false;
        
        % Calculating R2Eff from Intensity and vCPMG (Dataset class)
        OUTPUT_DEBUG_CALCULATE_R2EFF = false;
        
        % Enforcing minimum error in data (Dataset class)
        OUTPUT_DEBUG_ENFORCE_MIN_ERROR = false;
        
        % Reading NLIN.TAB file (Dataset class)
        OUTPUT_DEBUG_READ_NLIN_FILE = false;
        
        % Setting residue names to each index (Dataset class)
        OUTPUT_DEBUG_SET_SEQUENCE = false;
        
        % Generating NEW groups via index/atom (Session class)
        OUTPUT_DEBUG_GENERATE_GROUPS = false;
        
        % Executing script file (Session class)
        OUTPUT_DEBUG_LOAD_SCRIPT = false;
        
        % In case the user wants to be frustrated beyond his/her mental capacity
        DRIVE_USER_INSANE = false;
    end
    
    properties( Constant )
        % Universal constants
        R                           = 1.9858775;    % Gas constnat cal/K/mol
        gamma_H                     = 42.576;       % 1H MHz/T
        
        
        AX_name_array               = {'1H', '13C', '15N', '19F'};
        AX_gamma_relative_array     = [1.000, 0.25143, 0.10136, 0.9407];
        
        
        Np                  = 5;    % Number of parameters
        Nr                  = 8;    % Number of results (Np + Rex, alpha, PhiexX)
        
        gamma_X             = 999;     % Should be changed later
        param_name          = { 'dwH', 'dwX', 'Pa', 'kex', 'R20', ...
                                'Rex', 'alpha', 'PhiexX', 'dH', 'Eab' };
        param_name_plot     = { '\Delta\omega_H', '\Delta\omega_X', 'P_A', 'k_e_x', 'R_2^0', ...
                                'R_{ex}', '\alpha', '\Phi_{ex}^{X}', '\DeltaH', 'E_{A}^{A\rightarrowB}'};
        AX_name_plot        = '^{A}X';

        % Parameter settings (p = 1-7)
        % dwH, dwX, PA, kex, R20, Rex, alpha        
        param_units_natural        = {'rad/s', 'rad/s', '', '/s', 'Hz', ...
                                    'Hz', '', 'Hz^2', 'cal/mol', 'cal/mol'};
        param_units_display        = {'Hz', 'Hz', '%', '/s', 'Hz', ...
                                      'Hz', '', 'Hz^2', 'kcal/mol', 'kcal/mol'};
        convert_units_to_display   = [ 1/(2*3.1416) 1/(2*3.1416) 100 1 1 ...
                                       1 1 1 1/1e3 1/1e3];

        % Use same limits for each dataset (SQ will have max dwH=0 set later)
        min_dwH_ppm     = 0;
        min_dwX_ppm     = 0;
        min_pa          = 0.5;
        min_kex         = 1;
        min_R20         = 0;

        max_dwH_ppm     = 5;
        max_dwX_ppm     = 10;
        max_pa          = 1.0;
        max_kex         = 15000;
        max_R20         = 1000;

        % Grid search ranges (note this changes for the nucleus type)
        grid_p_init    = [ 0.0 0.1 0.500 10.0 ];
        grid_p_final   = [ 1.0 5.0 0.999 6000 ];
        grid_p_steps   = [ 4   4   4     15   ];
        
        % Initial guess for enthalpy and activation energy
        dH_Init         = -10e3;
        Ea_Init         = 1e3;
    end
    
    properties (SetAccess = private, GetAccess = private)
        % Cannot read nor write via ".", except within class functions              
        
        % TODO % Log the sequence_filename once it is used
        sequence_filename   = '';   % Name of the sequence file for residue_array

    end
    
    methods
        function obj = Session()
            %% Constructor
            obj.START_DATE  = datestr(now, 'yyyy.mm.dd-HH-MM');
            obj.outputDir   = ['output-', datestr(now, 'yyyy.mm.dd-HH-MM')];     
            
            obj.paramDisplay = ParamDisplay(obj);
        end
                
        function addDataset( obj, dataset )
            %% Add the new dataset to the end of the list
            % Since Group is a subclass of handle, the added curveset
            % points to its place in memory, and subsequent changes to
            % either source curveset or obj.curvesets{obj.Ncs+1} are linked
            
            % Set the sequence array, if possible
            if( ~isempty(obj.sequence_array) )
                dataset.setSequence( obj.sequence_array );
            end
            obj.datasets{obj.Nds+1} = dataset;            
            obj.Nds = obj.Nds+1;
            
            % Make sure this session is the parent for the dataset
            dataset.setParentSession(obj);
        end
        
        function addDisplayCluster( obj, displayCluster )
            %% Add the new displaycluster to the end of the list
            obj.displayClusters{obj.Ndc+1} = displayCluster;            
            obj.Ndc = obj.Ndc+1;
        end
        
        function addGroup( obj, group )
            %% Add the new group to the end of the list
            % Since Group is a subclass of handle, the added curveset
            % points to its place in memory, and subsequent changes to
            % either source groups or obj.groups{obj.Ncs+1} are linked
            obj.groups{obj.Ng+1} = group;
            obj.Ng = obj.Ng+1;
        end
        
        function addNote( obj, text, name )
            %% Add a new note to the session
            Nnotes              = obj.Nnotes+1;            
            obj.notes{Nnotes}   = Note( text, name );            
            obj.Nnotes          = Nnotes;
        end
        
        function curve_in_session = containsCurve( obj, curve )
            %% (Boolean) The supplied group is in the session
            % TODO % (LOW PRIORITY) containsCurve
        end
        
        function curveset_in_session = containsCurveset( obj, curveset )
            %% (Boolean) The supplied group is in the session
            % TODO % (LOW PRIORITY) containsCurveset
        end
        
        function dataset_in_session = containsDataset( obj, dataset )            
            %% (Boolean) The supplied dataset is in the session
            % TODO % (LOW PRIORITY) containsDataset                        
        end
        
        function group_in_session = containsGroup( obj, group )
            %% (Boolean) The supplied group is in the session
            % TODO % (LOW PRIORITY) containsGroup
        end
        
        function matrixDisplay = convertMatrixToDisplayUnits(obj, matrixNatural)
            %% Convert units from NATURAL (rad/sec, fraction) -> DISPLAY (Hz, %)
            Nctot           = size(matrixNatural,1);
            Nptot           = size(matrixNatural,2);
            convert_units   = obj.convert_units_to_display(1:Nptot);
            
            % In case there are more units to convert (e.g., resultsMatrix
            % includes columns 6 and 7 for Rex and alpha)
            if(Nptot > length(convert_units))
                convert_units(end+1:Nptot) = 1;
            end
            matrixDisplay = matrixNatural .* (ones(Nctot,1) * convert_units);
        end
        
        function matrixNatural = convertMatrixToNaturalUnits(obj, matrixDisplay)
            %% Convert units from NATURAL (rad/sec, fraction) -> DISPLAY (Hz, %)
            Nctot           = size(matrixDisplay,1);
            Nptot           = size(matrixDisplay,2);
            convert_units   = obj.convert_units_to_display(1:Nptot);
            
            % In case there are more units to convert (e.g., resultsMatrix
            % includes columns 6 and 7 for Rex and alpha)
            if(Nptot > length(convert_units))
                convert_units(end+1:Nptot) = 1;
            end
            matrixNatural = matrixDisplay ./ (ones(Nctot,1) * convert_units);
        end
        
        function valueDisplay = convertParameterToDisplayUnits(obj, valueNatural, paramNameOrNumber )
            %% Convert units from NATURAL (rad/sec, fraction) -> DISPLAY (Hz, %)
            % valueDisplay is the display value
            % paramNameOrNumber is the paramter number (1-5) OR the parameter name
            
            if( ischar(paramNameOrNumber) )
                isParam     = strcmpi( paramNameOrNumber, obj.param_name );
                p           = find(isParam);
                
                if( isempty(p) )
                    error('Invalid parameter name specified');
                end
                
            elseif( isnumeric(paramNameOrNumber) )
                p           = paramNameOrNumber;
                
                if( p > length(obj.convert_units_to_display) )
                    error('Invalid parameter number specified');
                end
            end
            
            % Make the conversion
            valueDisplay = valueNatural .* obj.convert_units_to_display(p);
        end
        
        function valueNatural = convertParameterToNaturalUnits(obj, valueDisplay, paramNameOrNumber)
            %% Convert units from NATURAL (rad/sec, fraction) -> DISPLAY (Hz, %)
            % valueDisplay is the display value
            % paramNameOrNumber is the paramter number (1-5) OR the parameter name
            
            if( ischar(paramNameOrNumber) )
                isParam     = strcmpi( paramNameOrNumber, obj.param_name );
                p           = find(isParam);
                
                if( isempty(p) )
                    error('Invalid parameter name specified');
                end
                
            elseif( isnumeric(paramNameOrNumber) )
                p           = paramNameOrNumber;
                
                if( p > length(obj.convert_units_to_display) )
                    error('Invalid parameter number specified');
                end
            end
            
            % Make the conversion
            valueNatural = valueDisplay ./ obj.convert_units_to_display(p);
        end
        
        function contains_ok_fit = displayClustersContainOKFit(obj)
            %% (Bool) At least one display cluster contains an OK fit
            contains_ok_fit = false;
            
            for dc = 1:obj.Ndc
                if( obj.displayClusters{dc}.containsOKFit() )
                    % Found an OK fit! Return early
                    contains_ok_fit = true;
                    return
                end
            end
        end
        
        function display_note = displayNote(obj, n)
            %% Return multiline string for showing note in GUI edit box
            %multiline_note = cellArrayToMultilineString( obj.notes{n} );
            
            note_session    = obj.notes{n}.text;
            Nlines          = size(note_session,1);
            line            = 0;
            display_note    = cell(Nlines, 1);
            for l = 1:Nlines
                line = line+1;
                display_note{line} = note_session{l};
            end
        end
        
        function exportDatasets(obj, FILE)
            %% Export datasets to a file
            % FILE = already opened file handle (1 => command window)
            for ds = 1:obj.Nds                
                fprintf(FILE, '\n\n# Dataset %d/%d', ds, obj.Nds);
                
                dataset = obj.datasets{ds};
                dataset.exportData(FILE);
            end            
        end
        
               
        function isequal = eq(s1, s2)
            % TODO % (LOW PRIORITY) test if two sessions are equal
        end        
        
        function generateGroups( obj, varargin )
            %% Generate minimal set of NEW groups to partition curves via NMR probe (index/atom)            
            
            OUTPUT_DEBUG = obj.OUTPUT_DEBUG_GENERATE_GROUPS;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);                    
            end
            
            % Get the initial number of groups
            Ng_initial = obj.Ng;            
            
            SPECIFY_INDEX   = false;
            SPECIFY_ATOM    = false;
            SPECIFY_RESIDUE = false;
            SPECIFY_TEMP    = false;
            SPECIFY_B0      = false;
            SPECIFY_SQX     = false;
            NAME_SUFFIX     = '';
            
            if( OUTPUT_DEBUG )
                varargin
                nargin
            end
            
            % Go through each pair of input arguments
            for v = 1:2:nargin-1
                property    = varargin{v};
                value       = varargin{v+1};
                
                if( OUTPUT_DEBUG )
                    fprintf('\n\tWorking on %s = %s (or %d)', property, value, value);
                end
                
                switch( upper(property) )
                    case 'INDEX'
                        SPECIFY_INDEX   = true;
                        INDEX           = value;                               
                        NAME_SUFFIX     = sprintf('%sI%d ', NAME_SUFFIX, INDEX);
                        
                    case 'ATOM'
                        SPECIFY_ATOM    = true;
                        ATOM            = value;
                        NAME_SUFFIX     = sprintf('%sA%s ', NAME_SUFFIX, ATOM);
                        
                    case 'RESIDUE'
                        SPECIFY_RESIDUE = true;
                        RESIDUE         = value;
                        NAME_SUFFIX     = sprintf('%sR%s ', NAME_SUFFIX, RESIDUE);
                                                
                    case 'TEMP'
                        SPECIFY_TEMP    = true;
                        TEMP            = value;      
                        NAME_SUFFIX     = sprintf('%s%dC ', NAME_SUFFIX, TEMP-273);
                        
                    case 'TEMPC'
                        SPECIFY_TEMP    = true;
                        TEMP            = value+273;      
                        NAME_SUFFIX     = sprintf('%s%dC ', NAME_SUFFIX, TEMP-273);
                        
                    case 'B0'
                        SPECIFY_B0      = true;
                        B0              = value;
                        NAME_SUFFIX     = sprintf('%s%dMHz ', NAME_SUFFIX, B0);
                        
                    case 'SQX'
                        SPECIFY_SQX     = true;
                        SQX             = value;
                        if( SQX )
                            NAME_SUFFIX     = sprintf('%sSQ ', NAME_SUFFIX);
                        else
                            NAME_SUFFIX     = sprintf('%sMQ ', NAME_SUFFIX);
                        end
                        
                    otherwise
                        error('Invalid property specified: "%s"', property);
                end
            end
            
            % Trim trailing space, if needed
            if( ~isempty(NAME_SUFFIX) )
                NAME_SUFFIX = NAME_SUFFIX(1:end-1);
            end
                        
            % Iterate through each curve in each dataset and create a new
            % group if one with its index does not exist already
            % ! Each curve must be assigned a group
            for ds = 1:obj.Nds
                dataset = obj.datasets{ds};
                for c = 1:dataset.Nc
                    curve = dataset.curves{c};
                    
                    % Do not use this curve if it does not match the
                    % constraits that are imposed
                    %  continue => Go to next iteration of for loop
                    if( SPECIFY_INDEX && curve.index ~= INDEX )
                        continue;
                    end
                    
                    if( SPECIFY_ATOM && ~strcmp(curve.atom,ATOM) )
                        continue;
                    end
                    
                    if( SPECIFY_RESIDUE && ~strcmp(curve.residue, RESIDUE) )
                        continue;
                    end             
                    
                    if( SPECIFY_TEMP && curve.Temp ~= TEMP )
                        continue;
                    end
                    
                    if( SPECIFY_B0 && curve.B0 ~= B0 )
                        continue;
                    end
                    
                    if( SPECIFY_SQX && curve.SQX ~= SQX )
                        continue;
                    end
                                        
                    % Check if any NEW groups contains this probe
                    curveset_for_probe = [];
                    for g = Ng_initial+1:obj.Ng
                        group = obj.groups{g};
                        
                        % Check if the probe already has a group
                        curveset_for_probe = group.getCurvesetForProbe(curve.index, curve.atom);
                        
                        % Add the probe to the existing group's curveset
                        if( ~isempty(curveset_for_probe) )
                            curveset_for_probe.addCurve(curve);
                            break
                        end
                    end
                    
                    % If no curveset was found in any group
                    %  then add a new group and curveset
                    if( isempty(curveset_for_probe) )
                        % Add a new group
                        new_group = Group(obj);
                        
                        if( ~isempty(NAME_SUFFIX) )
                            NAME = sprintf('%s %s', curve.name, NAME_SUFFIX);
                        else
                            NAME = curve.name;
                        end
                        new_group.setName(NAME);
                        new_group.setIndex(curve.index);
                        obj.addGroup(new_group)
                        
                        % Add a new curveset
                        new_curveset = Curveset;
                        new_curveset.setAssignment(curve.index, curve.atom, curve.residue);
                        new_group.addCurveset(new_curveset);
                        
                        % Add the curve
                        new_curveset.addCurve(curve);
                    end
                end
            end
        end
        
        function generateCurvesetsForGroup( obj, g )
            %% Generate minimal set of NEW curvesets to partition curves via NMR probe (index/atom)            
            % Get the initial number of groups
            % g = number (or class) for the group for curvesets to be added to
            
            % Get the group to add the curvesets to
            if( isnumeric(g) )
                group = obj.groups{g};
            elseif( isa(g,'Group') )
                group = g;
            else
                error('Must supply "g" as group number or Group instance');
            end
            
            % Iterate through each curve in each dataset and create a new
            % group if one with its index does not exist already
            % ! Each curve must be assigned a group
            for ds = 1:obj.Nds
                dataset = obj.datasets{ds};
                for c = 1:dataset.Nc
                    curve = dataset.curves{c};
                    
                    % Check if the NMR probe already has a curveset in this
                    % group
                    curveset_for_probe = group.getCurvesetForProbe(curve.index, curve.atom);

                    % If the curveset for this probe already exists...
                    if( ~isempty(curveset_for_probe) )
                        % ...Add the probe to the existing group's curveset
                        curveset_for_probe.addCurve(curve);
                     
                    % If no curveset was found in any group...
                    else
                                                                
                        % ...create a new curveset for this group
                        new_curveset = Curveset;
                        new_curveset.setAssignment(curve.index, curve.atom, curve.residue);
                        group.addCurveset(new_curveset);
                        
                        % Add the curve
                        new_curveset.addCurve(curve);
                    end
                end
            end
        end
        
        function B0_Array = getB0Array( obj )
            %% Get array of B0 field strengths from all groups
            
            B0_Array = [];
            b = 1;
            
            for g = 1:obj.Ng
                group = obj.groups{g};
                for cs = 1:group.Ncs
                    curveset = group.curvesets{cs};
                    for c = 1:curveset.Nc
                        curve = curveset.curves{c};
                        
                        % Check / replace limits
                        B0_Array(b) = curve.B0;
                        b = b+1;
                    end
                end
            end
        end
        
        function curve_array = getCurves(obj)
            %% Return array of pointers to each curve in session
            % TODO % Do NOT create copy of curve using "Curve"
            
            curves = {};
        end
                
        function curveset_array = getCurvesets(obj)
            %% Return array of pointers to each curveset in session
            % TODO % Do NOT create copy of curve using "Curveset"
            curvesets = {};
        end
        
        function [vcpmg_min, vcpmg_max, R2eff_min, R2eff_max] = getDataLimits( obj )
            %% Return vcpmg and R2eff limits on data by parsing each curve
            % TODO % Unnecessary? Use group.getDataLimits() instead?
            
            vcpmg_min   = Inf;        % Smallest vcpmg value in all curves
            vcpmg_max   = -Inf;        % Largest vcpmg value in all curves
            R2eff_min   = Inf;        % Smallest R2eff value in all curves
            R2eff_max   = -Inf;        % Largest R2eff value in all curves
            
            for g = 1:obj.Ng
                group = obj.groups{g};
                for cs = 1:group.Ncs
                    curveset = group.curvesets{cs};
                    for c = 1:curveset.Nc
                        curve = curveset.curves{c};
                        
                        % Check / replace limits
                        vcpmg_min = min( [vcpmg_min curve.vcpmg] );
                        vcpmg_max = max( [vcpmg_max curve.vcpmg] );
                        R2eff_min = min( [R2eff_min curve.R2eff] );
                        R2eff_max = max( [R2eff_max curve.R2eff] );
                    end
                end
            end
        end
        
        function dislaycluster_string_array = getDisplayClusterStringArray(obj)
            %% Return cell array of displaycluster strings
            dislaycluster_string_array = {''};            
            for dc = 1:obj.Ndc
                dislaycluster_string_array{dc} = sprintf('(%d) %s',  dc, obj.displayClusters{dc}.name);                
            end
        end
        
        function group = getGroup( obj, index, atom )
            %% Return the group(s) with these indices and atoms
            % TODO % (LOW PRIORITY) containsGroup
        end
        
        function Nctot = getNumCurves( obj )
            %% Get total number of curves among all groups
            Nctot = 0;
            for g = 1:obj.Ng
                Nctot = Nctot + obj.groups{g}.getNumCurves();
            end
        end
        
        function Ng_DisplayClusters = getNumGroupsInDisplayClusters(obj)
            %% Return number of groups in all dislay clusters
            Ng_DisplayClusters = 0;
            for dc = 1:obj.Ndc
                Ng_DisplayClusters = Ng_DisplayClusters + obj.displayClusters{dc}.Ng;
            end
        end
                
        function [title, axis_label, DATA_MIN] = getPlotLabels( obj, param_name, Temp, B0 )
            %% Return plot title and axis label for a given parameter name
            DATA_MIN            = NaN;
            
            % Check for parameter requirements
            [needsTemp, needsB0] = obj.getParamRequirements( param_name );
            
            switch( upper(param_name) )   
                case 'RESIDUE'
                    DATA_MIN    = 0;
                    title       = sprintf('Residue Number');
                    axis_label  = sprintf('Residue Number');
                    
                case 'DWH'
                    DATA_MIN    = 0;
                    title       = sprintf('^1H Chemical shift difference');
                    axis_label  = sprintf('|\\Delta\\omega_H| (Hz)');   
                    
                case {'DWH_PPM', 'DWHPPM'}
                    DATA_MIN    = 0;
                    title       = sprintf('^1H Chemical shift difference');
                    axis_label  = sprintf('|\\Delta\\omega_H| (ppm)');   
                    
                case 'DWX'
                    DATA_MIN    = 0;
                    title    = sprintf('%s Chemical shift difference', obj.AX_name_plot );
                    axis_label   = sprintf('|%s| (Hz)', obj.param_name_plot{2});
                    
                case {'DWX_PPM', 'DWXPPM'}
                    DATA_MIN    = 0;
                    title    = sprintf('%s Chemical shift difference', obj.AX_name_plot );
                    axis_label   = sprintf('|%s| (ppm)', obj.param_name_plot{2});
                    
                case 'PA'
                    DATA_MIN    = NaN;
                    title       = sprintf('Major state population at %d^oC', Temp-273);
                    axis_label  = sprintf('P_A (%%)');
                    
                case 'KEX'
                    DATA_MIN    = 0;
                    title    = sprintf('Total exchange rate at %d^oC', Temp-273);
                    axis_label   = sprintf('k_{ex} (/s)');

                case 'R20'
                    DATA_MIN    = 0;
                    title    = sprintf('Intrinsic R_2 relaxation rate at %d^oC at %0.0f MHz', Temp-273, B0);
                    axis_label   = sprintf('R_2^0 (Hz)');

                case 'REX'
                    DATA_MIN    = 0;
                    title    = sprintf('Exchange broadening at %d^oC at %0.0f MHz', Temp-273, B0);
                    axis_label   = sprintf('R_{ex} (Hz)');

                case 'ALPHA'
                    DATA_MIN    = 0;
                    title    = sprintf('Exch. timescale \\alpha at %d^oC', Temp-273);
                    axis_label   = sprintf('\\alpha (-)');
                    
                case 'PHIEXX'
                    DATA_MIN    = 0;
                    title    = sprintf('Exch. factor at %d^oC at %0.0f MHz', Temp-273, B0);
                    axis_label   = sprintf('\\Phi_{ex}^{X} (Hz^2)');

                case 'KA'
                    DATA_MIN    = 0;
                    title    = sprintf('A\\rightarrowB exchange rate at %d^oC', Temp-273);
                    axis_label   = sprintf('k_A (/sec)');

                case 'KB'
                    DATA_MIN    = 0;
                    title    = sprintf('B\\rightarrowA exchange rate at %d^oC', Temp-273);
                    axis_label   = sprintf('k_B (/sec)');

                case 'K'
                    DATA_MIN    = 0;
                    title    = sprintf('Equilibrium constant [A]/[B] at %d^oC', Temp-273);
                    axis_label   = sprintf('K (-)');

                case 'DH'
                    title    = sprintf('A\\rightarrowB enthalpy change');
                    axis_label   = sprintf('\\DeltaH (kcal/mol)');

                case 'DS'
                    title    = sprintf('A\\rightarrowB entropy change');
                    axis_label   = sprintf('\\DeltaS (cal/mol/K)');

                case {'EAB', 'EA(A->B)'}
                    %DATA_MIN    = 0;
                    title    = sprintf('A\\rightarrowB activation energy');
                    axis_label   = sprintf('E_A(A\\rightarrowB) (kcal/mol)');

                case {'EBA', 'EA(B->A)'}
                    %DATA_MIN    = 0;
                    title    = sprintf('B\\rightarrowA activation energy');
                    axis_label   = sprintf('E_A(B\\rightarrowA) (kcal/mol)');
                    
                otherwise
                    error('Invalid paramter name specified. Check code for valid options');
            end
        end
        
        function [symbolChar, colorRGB] = getPlotSymbolAndColor(obj, plotNum)
            %% Return plot symbol character ('o', 's', etc.) and colorRGB vector

            NPlotColors = size(obj.plotColorOrder, 1);

            % Get symbol and color, First go through each color, then next symbol
            % Color  1 2 3 4 ... Ncolors     1 2 3 4 ... Ncolors    1 2  ...    
            color_index     = mod(plotNum-1, NPlotColors)+1;

            % Symbol 1 1 1 1 ... (xNcolors)  2 2 2 2 ... (xNcolors) 3 3 3 3 ..
            symbol_index    = floor((plotNum-1)/NPlotColors)+1;
            
            % Access the values specified above
            symbolChar = obj.plotSymbolOrder{symbol_index};            
            colorRGB    = obj.plotColorOrder(color_index,:);
        end
        
        function Temp_array = getTemperatureArray( obj )
            %% Get array of temperatures from all groups
            
            Temp_array = [];
            t = 1;
            
            for g = 1:obj.Ng
                group = obj.groups{g};
                for cs = 1:group.Ncs
                    curveset = group.curvesets{cs};
                    for c = 1:curveset.Nc
                        curve = curveset.curves{c};
                        
                        % Check / replace limits
                        Temp_array(t) = curve.Temp;
                        t = t+1;
                    end
                end
            end
        end
        
        function loadDatasets( obj, script_filename )
            %% Load 1+ datasets using script file
            % 2011/04/27 Code script file
            % 2011/05/13 Code individual RD curves file
            %
            %
            % EXAMPLE SCRIPT FILE BELOW
            % # Ian Kleckner
            % # GUARDD Load data script
            % # Apo TRAP data
            % # 
            % # 2011/04/27
            % # 2011/05/20 Updated for individual DATA in table format
            % # 
            % #
            % # Syntax
            % #  #		-> Comment (ignore rest of line)
            % #
            % #  DATASET	-> Start a new dataset
            % #	Executes command: dataset = Dataset()
            % #
            % #	NAME	-> Specify name (or [] for auto-name)
            % #	AX	-> Specify AX nucleus (13C, 15N, or 19F)
            % #	B0	-> Magnetic field strength (MHz)
            % #	TEMPC	-> Celcius temperature (can also use TEMP for Kelvin)
            % #	TCPMG	-> CPMG time (sec)
            % #	SQX	-> (Boolean) single quantum mode (true/false)#	
            % #	
            % #  SETSPECS	-> Set prior specifications for dataset
            % #    		-> Must set these arguments first:
            % #	Executes command: dataset.setSpecs( NAME, AX, B0, TEMP, TCPMG, SQX )
            % #
            % #  INDEX 	-> Specify NMR signal index (AA#, peak#) for subsequent data
            % #  ATOM		-> Specify atom string (NH, C, CO, \delta_1)
            % #  OBS		-> Turns observation mode ON
            % #		-> Next column is VCPMG (the X-MODE)
            % #		-> Next column is R2 or INTENSITY (the Y-MODE)
            % #		-> Next column is ERROR (this is optional)
            % #		-> Observation mode will store data points until ADDDATA
            % #		-> Errors can be specified explicitly, calculated via repeat VCPMG values or not specified at all
            % #		-> All curves in the same dataset must have the same VCPMG values
            % #
            % #
            % #  ADDDATA 	-> Add all the data from last OBS mode
            % #  	Executes command: dataset.addData(INDEX, ATOM, X, Y, Y_E, Y_MODE);
            % #
            % #  NLINFILE 	-> Specify nlin.tab file (NMRPipe)
            % #  VCPMGFILE 	-> Specify vcpmg.txt file (VCPMG values for to nlin.tab)
            % #
            % #  READNLIN	-> Read the previously specified NLIN file
            % #	Executes command: dataset.readNlin(NLINFILE, VCPMGFILE);
            % # 		
            % #  SEQUENCEFILE	-> Specify and read sequence file
            % #  	This can be used anywhere in the script (only needed once)

            OUTPUT_DEBUG = obj.OUTPUT_DEBUG_LOAD_SCRIPT;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end
            
            % Read lines from script file
            try
                FILE    = fopen(script_filename, 'r');
            catch
                error('Could not open script file %s', script_filename);
            end
            
            % Start timing execution of function
            fprintf('Loading file...');
            t0 = tic;            
            
            lines   = textscan(FILE, '%s', 'delimiter', '\n');
            fclose(FILE);
            
            OBS_MODE = false;
            
            Nlines  = length(lines{1});
            
            for l = 1:Nlines
                % Get the line string
                line = lines{1}(l);
                line = line{1};
                
                if( OUTPUT_DEBUG )
                    fprintf('\n\nLine %d: %s', l, line);
                end
                
                % If there is something to read
                if( ~isempty(line) )
                    
                    tokens  = textscan(char(line), '%s');
                    Ntokens = size(tokens{1},1);
                    
                    %firstToken = char(tokens{1}(1));
                    firstToken = tokens{1}(1);
                    firstToken = upper( firstToken{1} );
                    
                    firstChar = firstToken(1);                    
                    
                    % # -> Comment, ignore
                    if( ~strcmp(firstChar, '#') )
                        if( Ntokens > 1 )
                            secondToken = tokens{1}(2);
                            secondToken = secondToken{1};
                        else
                            secondToken = '';
                        end
                        
                        % Check that second token is supplied for those
                        % that require it
                        switch upper(firstToken)
                            case {'NAME', 'AX', 'B0', 'TCPMG', 'SQX', ...
                                  'TEMP', 'TEMPK', 'TEMPC', 'MODE', ...
                                  'NLINEFILE', 'VCPMGFILE', ...
                                  ...
                                  'INDEX', 'ATOM', 'RESIDUE', 'OBS', ...
                                  'SEQUENCEFILE', ...
                                  'READSCRIPT'}
                                if( isempty(secondToken) )                                    
                                    error('Must supply second token for %s', firstToken);
                                end
                        end
                        
                        % Process the input
                        switch upper(firstToken)
                            case 'DATASET'
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tdataset = Dataset( obj );');
                                    fprintf('\n\tobj.addDataset(dataset, obj);');
                                end
                                
                                dataset = Dataset(obj);
                                obj.addDataset(dataset);
                                                                
                            case 'SETSPECS'
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tdataset.setSpecs(NAME, AX, B0, TEMP, TCPMG, SQX);');                                
                                end
                                dataset.setSpecs(NAME, AX, B0, TEMP, TCPMG, SQX);                               
                                                    
                            case {'NAME'}
                                % Accept second token as [] or as string
                                if( strcmp(secondToken,'[]') )
                                    if( OUTPUT_DEBUG )
                                        fprintf('\n\t%s = [];', firstToken);                                    
                                    end
                                    eval(sprintf('%s = [];', firstToken));
                                else
                                    if( OUTPUT_DEBUG )
                                        fprintf('\n\t%s = \''%s\'';', firstToken, secondToken);                                    
                                    end
                                    eval(sprintf('%s = \''%s\'';', firstToken, secondToken));
                                end
                                
                            case {'AX', 'ATOM', 'RESIDUE'}
                                % Second token is a STRING
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\t%s = \''%s\'';', firstToken, secondToken);                                
                                end
                                eval(sprintf('%s = \''%s\'';', firstToken, secondToken));
                                
                            case {'B0', 'TCPMG', ...
                                    'INDEX'}
                                % Second token is a NUMBER
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\t%s = %s;', firstToken, secondToken);                                
                                end
                                eval(sprintf('%s = %s;', firstToken, secondToken));
                                
                            case 'SQX'
                                % Second token is a BOOLEAN
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\t%s = %s;', firstToken, lower(secondToken));
                                end
                                eval(sprintf('%s = %s;', firstToken,  lower(secondToken)));
                                
                            case {'TEMP', 'TEMPK'}
                                % Second token is a number converted C <->K
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tTEMP = %s;', secondToken);                                
                                    fprintf('\n\tTEMPC = %s-273;', secondToken);
                                end
                                
                                eval(sprintf('TEMP = %s;', secondToken));
                                eval(sprintf('TEMPC = %s-273;', secondToken));
                                
                            case 'TEMPC'
                                % Second token is a number converted C <->K
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tTEMPC = %s;', secondToken);
                                    fprintf('\n\tTEMP = %s+273;', secondToken);
                                end
                                
                                eval(sprintf('TEMPC = %s;', secondToken));
                                eval(sprintf('TEMP = %s+273;', secondToken));
                                
                            case 'READNLIN'
                                % Command to read the NLIN file
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tdataset.readNlin(NLINFILE, VCPMGFILE);');                                
                                end
                                dataset.readNlin(NLINFILE, VCPMGFILE);
                                
                            case {'MODE', 'NLINFILE', 'VCPMGFILE'}
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\t%s = \''%s\'';', firstToken, secondToken);                                
                                end
                                eval(sprintf('%s = \''%s\'';', firstToken, secondToken));
                                
                            case 'OBS'
                                % Observation mode, to read series of data
                                % points in a series
                                OBS_MODE = true;
                                
                                % Get X_MODE
                                if( Ntokens < 3 )
                                    error('Must specify OBS VCPMG and Y_MODE (with optional ERROR column)');
                                end
                                
                                % Get X_MODE (column 2)
                                X_MODE = tokens{1}(2);
                                X_MODE = X_MODE{1};
                                    
                                if( ~strcmpi(X_MODE, 'VCPMG') )
                                    error('X_MODE must be "VCPMG"');
                                end
                                
                                % Get Y_MODE (column 3)
                                Y_MODE = tokens{1}(3);
                                Y_MODE = Y_MODE{1};
                                
                                if( strcmpi(Y_MODE, 'INTENSITY') || ...
                                    strcmpi(Y_MODE, 'I') )
                                
                                    Y_MODE = 'INTENSITY';
                                    
                                elseif( strcmpi(Y_MODE, 'R2') || ...
                                        strcmpi(Y_MODE, 'R2EFF') )
                                    
                                    Y_MODE = 'R2';
                                    
                                else                                    
                                    error('Invalid Y_MODE specified, see code for options');
                                end
                                
                                % Are there errors?
                                if( Ntokens == 4 )
                                    SUPPLY_ERRORS = true;                                    
                                else
                                    SUPPLY_ERRORS = false;                                    
                                end
                                
                                % Reset data arrays
                                X       = [];
                                Y       = [];
                                Y_E     = [];
                                Nobs    = 0;
                                
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tObservation mode, X_MODE=%s, Y_MODE=%s, ERRORS_SUPPLIED=%d', ...
                                    X_MODE, Y_MODE, SUPPLY_ERRORS);
                                end
                                
                            case {'ADDPEAK', 'ADDDATA'}
                                OBS_MODE = false;
                                                               
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tdataset.addData(INDEX, ATOM, RESIDUE, X, Y, Y_E, Y_MODE);');
                                end
                                
                                % Call the command to add the curve
                                dataset.addData(INDEX, ATOM, RESIDUE, X, Y, Y_E, Y_MODE);
                                
                                % Reset specifications for a potential next
                                % curve
                                INDEX   = NaN;
                                ATOM    = '';
                                RESIDUE = '';
                                
                            case 'SEQUENCEFILE'
                                % Load the sequence file
                                % Extract the sequence_array, if possible
                                try
                                    FILE_SEQ = fopen(secondToken);
                                catch
                                    error('Could not open sequence file "%s"', secondToken);
                                end

                                % Read the line of text from the file
                                inputText_Seq   = textscan( FILE_SEQ, '%s' );
                                sequence_array  = inputText_Seq{1};

                                fclose(FILE_SEQ);

                                % Commit these changes to the session
                                obj.setSequence( sequence_array );    
    
                                if( OUTPUT_DEBUG )
                                    fprintf('\n\tLoaded sequence file "%s"', secondToken);
                                end
                                
                            case 'SCRIPTFILE'
                                % Make sure you don't try to read the
                                % current scriptfile (infinite loop)
                                if( strcmp(secondToken,script_filename) )
                                    error('Cannot read the current script file as a script (infinite loop)');
                                else
                                    % Load the scriptfile
                                    if( OUTPUT_DEBUG )
                                        fprintf('\n\tExecuting script file "%s"', secondToken);
                                    end
                                    obj.loadDatasets(secondToken);                                    
                                end
                                
                            otherwise
                                if( OBS_MODE )
                                    % Get the observation
                                    Nobs = Nobs+1;
                                    
                                    DATUM = tokens{1}(2);
                                    DATUM = DATUM{1};                                                                        
                                    eval(sprintf('X(Nobs)   = %s;', DATUM));
                                    
                                    DATUM = tokens{1}(3);
                                    DATUM = DATUM{1};
                                    eval(sprintf('Y(Nobs)   = %s;', DATUM));
                                    
                                    if( SUPPLY_ERRORS )
                                        DATUM = tokens{1}(4);
                                        DATUM = DATUM{1};
                                        eval(sprintf('Y_E(Nobs)   = %s;', DATUM));
                                        
                                        if( OUTPUT_DEBUG )
                                            fprintf('\n\t(X,Y,Y_E) = (%f, %f, %f)', X(Nobs), Y(Nobs), Y_E(Nobs));
                                        end
                                    else
                                        if( OUTPUT_DEBUG )
                                            fprintf('\n\t(X,Y,Y_E) = (%f, %f, %s)', X(Nobs), Y(Nobs), 'N/A');
                                        end
                                    end
                                    
                                    
                                    
                                else
                                    error('Invalid command %s', firstToken);
                                end
                        end
                        
                    else
                        if( OUTPUT_DEBUG )
                            fprintf('\n\tComment');
                        end
                    end
                end
            end            
            fprintf('Done (%0.1f sec)\n', toc(t0));
        end
        
        
        function outputGroups(obj, FILE)
            %% Output groups to a readable file
            % FILE is an already opened handle, or =1 for command window
            
            for g = 1:obj.Ng
                fprintf(FILE, '\n\nGroup %d/%d', g, obj.Ng);
                
                group = obj.groups{g};
                group.outputSpecs(FILE, '  ');
            end
        end
        
        function removeDataset( obj, dataset )            
            %% Find the dataset, then remove it            
            for ds = 1:length(obj.datasets)
                if( dataset == obj.datasets{ds} )                    
                    % Here, "ds" is the element number to be eliminated
                    %  E.g., Remove eement #3 from a set of 5 elements
                    %  Set #3=#4, #4=#5, and #5=[], then Ntot=Ntot-1
                    for ds1 = ds:obj.Nds-1
                        obj.datasets{ds1} = obj.datasets{ds1+1};
                    end
                    
                    % Set final element to null
                    obj.datasets(obj.Nds) = [];
                    
                    % Update the number of elements
                    obj.Nds = obj.Nds-1;
                    return
                end
            end            
        end
        
        function removeDisplayCluster( obj, displayCluster )            
            %% Find the group, then remove it            
            for dc = 1:obj.Ndc
                if( displayCluster == obj.displayClusters{dc} )
                    % Here, "dc" is the element number to be eliminated
                    %  E.g., Remove eement #3 from a set of 5 elements
                    %  Set #3=#4, #4=#5, and #5=[], then Ntot=Ntot-1
                    for dc1 = dc:obj.Ndc-1
                        obj.displayClusters{dc1} = obj.displayClusters{dc1+1};
                    end

                    % Set final element to null
                    obj.displayClusters(obj.Ndc) = [];

                    % Update the number of elements
                    obj.Ndc = obj.Ndc-1;
                    return
                end
            end
        end
        
        function removeGroup( obj, group )            
            %% Find the group, then remove it            
            for g = 1:obj.Ng
                if( group == obj.groups{g} )
                    % Here, "g" is the element number to be eliminated
                    %  E.g., Remove eement #3 from a set of 5 elements
                    %  Set #3=#4, #4=#5, and #5=[], then Ntot=Ntot-1
                    for g1 = g:obj.Ng-1
                        obj.groups{g1} = obj.groups{g1+1};
                    end

                    % Set final element to null
                    obj.groups(obj.Ng) = [];

                    % Update the number of elements
                    obj.Ng = obj.Ng-1;
                    return
                end
            end            
        end
        
        function removeNote( obj, n )
            %% Remove the note from the list
            % Here, "n" is the element number to be eliminated
            %  E.g., Remove eement #3 from a set of 5 elements
            %  Set #3=#4, #4=#5, and #5=[], then Ntot=Ntot-1
            for n1 = n:obj.Nnotes-1
                obj.notes{n1}       = obj.notes{n1+1};               
            end
            
            % Set final elements to null
            obj.notes{obj.Nnotes}       = [];

            % Update the number of elements
            obj.Nnotes = obj.Nnotes-1;
            return            
        end
        
        function resetParamDisplay(obj)
            %% Debug: reset the param display instance (2011/06/04)
            obj.paramDisplay = ParamDisplay(obj);
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
            
            % Check to see if certain variables have been defined
            % appropriately
            if( ~strcmpi(obj.MC_RANDOMIZE_MODE,'A') && ...
                ~strcmpi(obj.MC_RANDOMIZE_MODE,'B') )
                errordlg('MC_RANDOMIZE_MODE must be "A" or "B" only. It has been reset to default (A).');
                
                obj.MC_RANDOMIZE_MODE = 'A';
            end
        end
        
        function setChi2TopFraction( obj, chi2TopFraction )
            %% Set the parameter
            obj.chi2TopFraction = chi2TopFraction;
        end
        
        function setChi2MapParamters( obj, pChi2Map )
            %% Set the parameter
            obj.pChi2Map = pChi2Map;
        end
        
        function setDefault( obj, paramName )
            %% Set the paramter to default value
            % TODO % Finish setDefault() function
            
        end
        
        function setName( obj, name )
            %% Set the name of the group
            obj.name    = name;
        end
        
        function setNote(obj, n, text, name )
            %% Update note number n with the new note
            % Note or name can be empty [] => do not set it
            if( n > obj.Nnotes || n < 1 )
               error('Invalid note number specified');
            end
            
            % Set the note
            if( ~isempty(text) )
                %obj.notes{n}.setText( multilineStringToCellArray(text) );
                obj.notes{n}.setText( text );
            end
            
            % Set the name
            if( ~isempty(name) )
                obj.notes{n}.setName(name);
            end
        end
                
        function setSequence( obj, sequence_array )
            %% Set residue names
            obj.sequence_array = sequence_array;
            
            % Propagate this sequence to each dataset, which propagates to each curve           
            for ds = 1:obj.Nds
                obj.datasets{ds}.setSequence( sequence_array );
            end
        end
        
        function setXNucleus( obj, AX_String )
            %% Set the X nucleus for the dataset
            
            if( strcmp(AX_String,'13C') )
                obj.gamma_X                = 10.705;     % 13C MHz/T
                obj.param_name          = { 'dwH', 'dwC', 'Pa', 'kex', 'R20' };
                obj.param_name_plot     = { '\Delta\omega_H', '\Delta\omega_C', 'P_A', 'k_e_x', 'R_2^0' };
                obj.AX_name_plot       = '^{13}C';

            elseif( strcmp(AX_String,'15N') )    
                obj.gamma_X                = 4.3156;       % 15N MHz/T
                obj.param_name          = { 'dwH', 'dwN', 'Pa', 'kex', 'R20' };
                obj.param_name_plot     = { '\Delta\omega_H', '\Delta\omega_N', 'P_A', 'k_e_x', 'R_2^0' };
                obj.AX_name_plot       = '^{15}N';
                
            elseif( strcmp(AX_String,'19F') )    
                obj.gamma_X                = 40.052;       % 19F MHz/T
                obj.param_name          = { 'dwH', 'dwF', 'Pa', 'kex', 'R20' };
                obj.param_name_plot     = { '\Delta\omega_H', '\Delta\omega_F', 'P_A', 'k_e_x', 'R_2^0' };
                obj.AX_name_plot       = '^{15}N';

            else    
                obj.gamma_X                = 999;       % ?
                obj.param_name          = { 'dwH', 'dw?', 'Pa', 'kex', 'R20' };
                obj.param_name_plot     = { '\Delta\omega_H', '\Delta\omega_?', 'P_A', 'k_e_x', 'R_2^0' };
                obj.AX_name_plot       = '^{A}X';

                error('AX_String must be either 13C, 15N, or 19F');
            end
        end
        
        function sortGroups( obj )
            %% Sort the groups{} by index (?)
            % TODO % Must also sort other _arrays that index by "c"
            % TODO % Perhaps this should go to higher-order class / function?
            % TODO % Also may want to sort curves in curveset, curveset in group, etc.
            
            % Sort the indices in the groups
            %gIndex  = zeros(1,obj.Ng);
            
            gIndexName  = cell(obj.Ng,2);
            for g = 1:obj.Ng
               %gIndex(g)    = obj.groups{g}.index;                
               
               % Build a cell with COL1 = Index(String), COL2 = Name
               gIndexName{g,1}  = sprintf('%10d', obj.groups{g}.index);
               gIndexName{g,2}  = obj.groups{g}.name;
            end            
            %[VOID, gsort] = sort(gIndex);
            
            % First sort by index, then by name
            [VOID, gsort2] = sortrows(gIndexName, [1,2]);
                        
            % Use that sorted array to rearrange order of groups
            new_groups = cell(1,obj.Ng);
            for gs = 1:obj.Ng
                % The index of gIndex which should go next
                %new_groups{gs} = obj.groups{ gsort(gs) };
                new_groups{gs} = obj.groups{ gsort2(gs) };
            end
            
            obj.groups = new_groups;
        end
        
        function updateGroups( obj )
            %% Call updateFitParams() on each group
            CONSTRAIN_RATE_ANALYSIS = false;
            for g = 1:obj.Ng
               obj.groups{g}.updateFitParams(CONSTRAIN_RATE_ANALYSIS); 
            end
        end
    end
    
    methods( Static = true )        
        function valueConverted = convertUnits( value, param_name, NATURAL_OR_DISPLAY )
            %% Convert the parameter units for arbitrary parameter for natural or display units
            % Value can be a vector
            % Natural: rad/sec  (faction)   cal/mol     K
            % Display: Hz       (percent)   kcal/mol    C
            
            multiply_to_get_display_units   = NaN;
            add_to_get_display_units        = NaN;
            
            switch( upper(param_name) )   
                case 'RESIDUE'
                    multiply_to_get_display_units   = 1;
                    add_to_get_display_units        = 0;
                    
                case {'TEMP', 'TEMPERATURE'};
                    multiply_to_get_display_units   = 1;
                    add_to_get_display_units        = -273;
                    
                case {'DWH', 'DWX'}
                    multiply_to_get_display_units   = 1/(2*3.1416);
                    add_to_get_display_units        = 0;                
                    
                case 'PA'
                    multiply_to_get_display_units   = 100;
                    add_to_get_display_units        = 0;
                    
                case {'DWH_PPM', 'DWHPPM', 'DWX_PPM', 'DWXPPM', ...
                      'KEX', 'R20', 'REX', 'ALPHA', 'PHIEXX', 'KA', 'KB', 'K', 'DS', ...
                      'PAB', 'P(A->B)', 'PBA', 'P(B->A)'}
                    multiply_to_get_display_units   = 1;
                    add_to_get_display_units        = 0;

                case {'DH', 'EAB', 'EA(A->B)', 'EBA', 'EA(B->A)'}
                    multiply_to_get_display_units   = 1/1000;
                    add_to_get_display_units        = 0;
                    
                otherwise
                    error('Invalid paramter name specified, "%s". Check code for valid options', param_name);
            end
            
            % Get size of input vector
            [Nrow, Ncol] = size(value);            
            
            % Apply conversion
            switch( upper(NATURAL_OR_DISPLAY) )
                case 'NATURAL'
                    valueConverted = value ./ multiply_to_get_display_units ...
                                    - add_to_get_display_units * ones(Nrow,Ncol); 
                    
                case 'DISPLAY'
                    valueConverted = value .* multiply_to_get_display_units ...
                                    + add_to_get_display_units * ones(Nrow,Ncol);
                    
                otherwise
                    error('Invalid conversion specified, "%s". Use either NATURAL or DISPLAY', NATURAL_OR_DISPLAY);
            end            
        end
        
        function gamma_relative = getGammaRelative( AX_String )
            %% Return the relative gamma GAMMA/GAMMA_H for the nucleus
            
            switch( upper(AX_String) )
                case '1H'
                    gamma_relative = 1;
                    
                case '13C'
                    gamma_relative = 0.25143;
                    
                case '15N'
                    gamma_relative = 0.10136;
                    
                case '19F'
                    gamma_relative = 0.9407;
                    
                otherwise
                    error('Specify AX nucleus 1H, 13C, 15N, or 19F');
            end
        end        
       
        function [needsTemp, needsB0] = getParamRequirements( param_name )
            %% Does the parameter need a particular Temp and/or B0?
            switch( upper(param_name) )
                case { 'PA', 'KEX', 'R20', 'REX', 'ALPHA', 'PHIEXX', 'KA', 'KB', 'K' }
                    needsTemp = true;
                otherwise
                    needsTemp = false;
            end
            
            switch( upper(param_name) )
                case { 'R20', 'REX', 'PHIEXX' };
                    needsB0 = true;
                otherwise
                    needsB0 = false;
            end 
        end
    end
end
classdef Group < handle
    % Group class: for a group of 1+ RD curvesets, each of which includes 1+ RD curve
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/01/17 Start coding
    %
    % Properties and methods for a group of RD curvesets (1 or more)
    
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        
        name        = 'GroupName';      % Name of group
        index       = 0;        % Index of the group (for sorting purposes)
        Ncs         = 0;        % Number of curve sets in group
        curvesets   = {};       % Array of curvesets in group
        
        notes       = '';       % Deprecated
        note        = {};       % Notes about the group (e.g., from fitting)
        
        parentSession = [];     % Parent session that owns this group
        
        % Cell matrix for parameters to describe this class via table
        %  Used in generateParamTable()        
        param_info = {...
            'name:(String) Name of group (e.g., Leu 22\delta_1, or NT Domain)' ...
            'index:(Number) Number of group (for sorting purposes)', ...
                      };
        
        % Fitting parameters
        %  There are an array of Exchange fit results, another array for
        %  excahnge fit results via grid search, and a single fit result
        %  defined as best for exchange, and a fit result for no exchange
        fitResults_Grid     = {};   % Array of FitResults for the grid search
        fitResults          = {};   % Array of FitResults for arbitrary fits (ex and no-ex)
        fitResult_NoEx      = [];   % Fit result to no exchange model
        fitResult_Best      = [];   % Best fit result out of all fits (ex and no-ex)
        Nf                  = 0;    % Number of fit results
        %f_Best              = NaN; % Don't hold best fit result because it
                                    %  is easy to lose account of
                                    %  Instead, use getFitResultString()
                                    %  or getBestFitIndex
        
        bestFitIsOK         = false;
        exhibitsExchange    = false;
        
        % Limits for the grid search
        % dwH(ppm), dwX(ppm), PA0(%), kex0(/s), dH(kcal/mol), Eab(kcal/mol)
        grid_p_min      = [ 0.01 0.1 60.0 10.0 -10 -10];
        grid_p_max      = [ 0.1 5.0 99.9 6000  10  10];
        grid_p_steps    = [ 3   3   3     4    3   3];
        
        % Temperature at which PA0 and kex0 are defined
        T0_Grid         = 298;
                  
        % Matrices have rows itemizing curves in the group 
        %  ctot=1,2,3,... by starting with cs=1, then c=1, c=2, ..., then cs=2, c=1, c=2, ...        
        % cscMatrix (ctot) -> [cs, c] which corresponds to that ctot#
        cscMatrix               = [];
                
        % indepParamIndexMatrix (ctot,p) -> Independent parameter number
        %  Ranges from zero to the total number of indepenent parameters
        %  There are often many duplicate entries
        indepParamIndexMatrix   = [];
        
        % paramScalingMatrix (ctot,p) -> Factor to multiply indep. param
        %   number by to get desired parameter
        paramScalingMatrix      = [];
        
        % Nonlinear constraint function
        %CONSTRAIN_RATE_ANALYSIS     = true;
        indepParam_dH               = NaN;
        indepParam_Eab              = NaN;
        indepParam_PA0              = NaN;
        indepParam_kex0             = NaN;
        T0                          = NaN;        
    end
    
    properties (SetAccess = private, GetAccess = private)
        % Cannot read nor write via ".", except within class functions               
    end
    
    properties(SetAccess = private, Dependent = true)
        % Calculated each time the property is accessed
    end
    
    methods
        function obj = Group( parentSession )
            %% Constructor function
            obj.parentSession = parentSession;
            
            obj.note = Note('Text', 'Name');
        end
        
        function addCurveset( obj, curveset )
            %% Add the new curveset to the end of the list
            % Since Group is a subclass of handle, the added curveset
            % points to its place in memory, and subsequent changes to
            % either source curveset or obj.curvesets{obj.Ncs+1} are linked
            obj.curvesets{obj.Ncs+1} = curveset;
            obj.Ncs = obj.Ncs+1;
            
            % This will invalidate all the fits because it changes the group
            obj.clearFits();
            CONSTRAIN_RATE_ANALYSIS = false;
            obj.updateFitParams(CONSTRAIN_RATE_ANALYSIS);
        end
        
        function addFitResult( obj, fitResult )
            %% Add fit result to the end of list (ensures NoEx fit is done first)
            
            % Must have a no-exchange fit first
            %  If there is NOT a no-exchange fit yet...
            if( isempty(obj.fitResult_NoEx) )
                % ...but it is supplied now
                if( strcmpi(fitResult.modelName, 'NOEXCHANGE') )
                    % ...then set it properly
                    obj.fitResult_NoEx = fitResult;                    
                else
                    % ...otherwise create a new one
                    fitResult_NoEx = FitResult(obj, 'NOEXCHANGE', 'NoEx');
                    updateGroupParamsFlag=false;
                    fitResult_NoEx.fitMe(updateGroupParamsFlag);
                    fitResult_NoEx.analyzeMe();
                    obj.addFitResult(fitResult_NoEx);
                    fprintf('\nAdded no-exchange fit to group');
                end
            end
            
            % Now that the no-exchange fit is available, add this one
            obj.fitResults{obj.Nf+1} = fitResult;           
            obj.Nf = obj.Nf+1;
            
            % If there is onlyi one fit result, it is the best
            if( obj.Nf == 1 )
                obj.fitResult_Best  = fitResult;
                %obj.f_Best          = 1;
            end
            
            % Update statistics
            obj.updateStatistics();
        end
        
        function clearFits( obj )
            %% Clear all the fits (in case the group has changed)
            % TODO % Call clearFits() when curveset in this group has been changed
            % Or when curve in curveset in this group has been changed
            obj.fitResults_Grid     = {};
            obj.fitResults          = {};
            obj.fitResult_NoEx      = [];
            obj.fitResult_Best      = [];
            obj.Nf                  = 0;
            %obj.f_Best              = NaN;

            obj.bestFitIsOK         = false;
            obj.exhibitsExchange    = false;
        end
        
        function curve_in_group = containsCurve( obj, curve )
            %% (Boolean) The supplied curve is in the group                        
            % Check each curve in each curveset in this group for a match
            curve_in_group = false;
            for c = 1:obj.Ncs
                if( obj.curvesets{c}.containsCurve(curve) )
                    curve_in_group = true;
                    return
                end
            end
        end
        
        function curveset_in_group = containsCurveset( obj, curveset )
            %% (Boolean) The supplied curveset is in the group                        
            % Check each curveset in this group for a match
            curveset_in_group = false;
            for cs = 1:obj.Ncs
                % If curveset is found, then return early
                if( curveset == obj.curvesets{cs} )
                    curveset_in_group = true;
                    return
                end
            end
        end
        
        
        function group_copy = copy( obj )
            %% Create a new group with all properties of the current curve
            % Except as if the group has not been fit
            % Copy each curveset too
            % DO NOT COPY the curves (these are basic data, and linked to dataset)
            %  If original curves are changed, so are those in the copy group
            group_copy              = Group(obj.parentSession);
            group_copy.name         = sprintf('Cp(%s)',obj.name);
            group_copy.index        = obj.index;
            
            group_copy.Ncs          = obj.Ncs;      
            
            % Copy each curveset, and maintain reference to original curves
            for cs = 1:obj.Ncs
                group_copy.curvesets{cs} = obj.curvesets{cs}.copy();
            end
            
            %group_copy.notes        = obj.notes;
                        
            group_copy.parentSession    = obj.parentSession;
            group_copy.param_info       = obj.param_info;
            
            group_copy.clearFits();
            CONSTRAIN_RATE_ANALYSIS = false;
            group_copy.updateFitParams(CONSTRAIN_RATE_ANALYSIS);
        end
        
        function pMatrix = delinearizePFmincon( obj, p_fmincon, fitResult )
            %% De-linearize parameter array to matrix form
            % Use group-specific infomation
            
            % Calculate temperature dependence of rates, if desired
            OUTPUT_DEBUG    = false;            
            NTemps          = obj.getNumTemps();
            
            Nctot           = obj.getNumCurves();
            Np              = obj.parentSession.Np;
                        
            if( fitResult.CONSTRAIN_RATE_ANALYSIS && NTemps > 1 )
                % Get parameters required for calculations
                R       = obj.parentSession.R;    
                T0      = obj.T0;
                PA0     = p_fmincon(obj.indepParam_PA0);
                kex0    = p_fmincon(obj.indepParam_kex0);
                dH      = p_fmincon(obj.indepParam_dH);
                Eab     = p_fmincon(obj.indepParam_Eab);

                % Calculate kinetic parameters of interest
                % kA and kB at temperature T0
                kA0 = (1-PA0) * kex0;
                kB0 = PA0 * kex0;

                % van't Hoff yields dS using PA=PA0 at temperature T=T0
                dS = R * log( (1-PA0)/PA0 ) + dH/T0;

                % van't Hoff also determines PA(T) for all T
                PA = @(T) 1 ./ ( 1+exp( dS/R - dH/R ./ T ) );

                % Arrhenius yields Pab using kA=kA0 at temperature T0
                Pab = kA0 * exp( Eab / (R*T0) );

                % Arrhenius also determines kA(T) for all T
                kA = @(T) Pab * exp( -1*Eab/R ./ T );

                % Kinetic parameters determine remaining exchange rates for all T
                kex = @(T) kA(T) ./ (1-PA(T));
                kB  = @(T) kex(T) - kA(T);
            end

            % De-linearize parameter array to matrix form
            pMatrix = zeros(Nctot, Np);
            for ctot = 1:Nctot
                for p = 1:Np

                    % If there are multiple temperatures, then reconstruct PA and kex
                    if( fitResult.CONSTRAIN_RATE_ANALYSIS && NTemps > 1 && (p==3 || p==4) )

                        [cs,c] = obj.getCurvesetCurve(ctot);

                        if( p == 3 )
                            % Calculate PA using van't Hoff analysis
                            pMatrix(ctot,p) = PA( obj.curvesets{cs}.curves{c}.Temp );
                            if( OUTPUT_DEBUG )
                                fprintf('\n\tCalculating PA(%dK) = %f%%', ...
                                    obj.curvesets{cs}.curves{c}.Temp, 100*PA( obj.curvesets{cs}.curves{c}.Temp ) );
                            end

                        elseif( p == 4 )
                            pMatrix(ctot,p) = kex( obj.curvesets{cs}.curves{c}.Temp );

                            if( OUTPUT_DEBUG )
                                fprintf('\n\tCalculating kex(%dK) = %f /s', ...
                                    obj.curvesets{cs}.curves{c}.Temp, kex( obj.curvesets{cs}.curves{c}.Temp ) );
                            end
                        end        

                    else        
                        % If parameter points to a fmincon index (not zero) then use it
                        %  Otherwise, the pMatrix value remains zero for simulation (e.g., dwH)
                        if( obj.indepParamIndexMatrix(ctot,p) ~= 0 )
                            pMatrix(ctot,p) = obj.paramScalingMatrix(ctot,p) * ...
                                p_fmincon( obj.indepParamIndexMatrix(ctot,p) );
                        end
                    end
                end
            end            
        end
        
        %{
        function multiline_note = displayNote(obj)
            %% Return multiline string for showing note in GUI edit box
            multiline_note = cellArrayToMultilineString( obj.notes );
        end
        %}
        
        
        function isequal = eq(g1, g2)
            %% (Boolean) The two groups are equal in all properties (Used as "g1==g2")            
            isequal = strcmp(g1.name, g2.name) && ...
                      g1.Ncs == g2.Ncs && ...
                      g1.Nf == g2.Nf;
                  
            % Compare each curve in the set too, if necessary
            if( isequal )
                for cs = 1:length(g1.curvesets)
                    % If curvesets are not equal, then return early
                    if( ~(g1.curvesets{cs} == g2.curvesets{cs}) )
                        isequal = false;
                        return
                    end
                end
            end
        end
        
        function [fitResult_Best, f_Best] = getBestFitResult(obj)
            %% Return the best fit result and its index in the array            
            for f = 1:obj.Nf                
                if( obj.fitResults{f} == obj.fitResult_Best )
                    f_Best          = f;
                    fitResult_Best  = obj.fitResults{f};
                    return
                end
            end
        end
    
        function curve_array = getCurves( obj )
            %% Return array of pointers to each curve in session
            % TODO % Do NOT create copy of curve using "Curve"            
            curves = {};
        end
        
        function curve_string_array = getCurveStringArray(obj)
            %% Return cell array of curve strings
            % List each curve in each curve set in this group
            curve_string_array = {''};
            ctot = 1;
            for cs = 1:obj.Ncs
                curveset = obj.curvesets{cs};
                for c = 1:curveset.Nc
                    curve = curveset.curves{c};
                    curve_string_array{ctot} = sprintf('%s[%s]', curveset.name, curve.getSpecsString());
                    ctot = ctot+1;
                end
            end
        end
                
        function curveset_array = getCurvesets(obj)
            %% Return array of pointers to each curveset in session
            % TODO % Do NOT create copy of curve using "Curveset"
            curvesets = {};
        end
        
        function [cs,c] = getCurvesetCurve( obj, ctot )
            %% Return curveset and curve number given the total curve num.
            if( isempty(obj.cscMatrix) )
                obj.cscMatrix = zeros(obj.getNumCurves(),2);
                ctot_counter = 0;
                for cs = 1:obj.Ncs
                    curveset = obj.curvesets{cs};
                    for c = 1:curveset.Nc
                        % Mark which cs and c is in the current row
                        ctot_counter = ctot_counter+1;
                        obj.cscMatrix(ctot_counter,1) = cs;
                        obj.cscMatrix(ctot_counter,2) = c;
                    end
                end
            end
            
            cs = obj.cscMatrix(ctot, 1);
            c  = obj.cscMatrix(ctot, 2);
        end
        
        function ctot = getCurveTot(obj, cs, c)
            %% Find the total curve number associated with the curveset and curve
            % Access 2D matrix: cscMatrix, with headers [ctot, cs, c]
            % Find the row numbers with column 2 equal to "cs"                        
            ctot_in_cs  = obj.cscMatrix(:,1) == cs;
            ctot_in_c   = obj.cscMatrix(:,2) == c;
                        
            % Find the ctot number at the curveset AND curve number
            ctot = find( and(ctot_in_cs, ctot_in_c) );
        end
        
        function ctot_array = getCurveTotArray(obj, cs)
            %% Find the curve total number associated with the curveset
            % Access 2D matrix: cscMatrix, with headers [ctot, cs, c]
            % Find the row numbers with column 2 equal to "cs"                        
            ctot_array = find( obj.cscMatrix(:,1) == cs );            
        end
        
        function curveset_for_probe = getCurvesetForProbe( obj, index, atom )
            %% Return the curveset instance with specified index and atom

            % Check each curveset in this group for a match
            curveset_for_probe = [];
            for cs = 1:obj.Ncs
                curveset = obj.curvesets{cs};
                % If a match is found, return early
                if( curveset.index == index && strcmp(curveset.atom, atom) )                    
                    curveset_for_probe = curveset;
                    return
                end
            end            
        end
        
        function curveset_string_array = getCurvesetStringArray(obj)
            %% Return cell array of curveset strings
            curveset_string_array = {''};            
            for cs = 1:obj.Ncs
                curveset_string_array{cs} = sprintf('%s %s',  obj.curvesets{cs}.name);                
            end
        end
        
        function [DATA, DATA_E, IS_OK] = getData(obj, param_name, Temp, B0, QC_String, CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY )
            %% Return data point (NATURAL UNITS) for the desired parameter, temperature, B0, and Quantum Coherence
            %  Only specify the Temp, B0 or QC_String if desired (otherwise [])
            % CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY  'ANY'       use cs = 1
            %                                       'ALL'       use cs = 1:Ncs
            %                                       [numeric]  use cs = numeric value specified
            
            [needsTemp, needsB0] = Session.getParamRequirements( param_name );
            
            Nctot = obj.getNumCurves();
            
            if( needsTemp && isempty(Temp) )
                error('Must specify temperature');
            end
            
            if( needsB0 && isempty(B0) )
                error('Must specify B0 field strength');
            end
            
            % Assume the value cannot be found (NaN)
            %  If it is found, then these numbers are replaced
            DATA    = NaN;
            DATA_E  = NaN;
            IS_OK   = false;
            
            % No fits available
            if( obj.Nf == 0 )
                return
            end
            
            % Check to see how the curvesets are specified
            if( isnumeric(CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY) )
                % Its just a number or array
                cs_array = CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY;
                
                % Ensure it is a row vector
                if( iscolumn(cs_array) )
                    cs_array = cs_array';
                end
                
            elseif( ischar(CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY) )
                % A string specifies a particular array of curvesets
                switch upper(CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY)
                    case 'ANY'
                        cs_array = 1;
                        
                    case 'ALL'
                        cs_array = 1:obj.Ncs;
                        
                    otherwise
                        error('Please specify CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY as ANY or ALL or a numeric array');
                end
                
            else
                error('Please specify CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY as ANY or ALL or a numeric array');
            end
            
            % Get data for each curveset
            if( length(cs_array) > 1 )
                for cs = cs_array
                    [DATA(cs), DATA_E(cs), IS_OK(cs)] = ...
                        obj.getData(param_name, Temp, B0, QC_String, cs);
                end
                return
                
            else
                % Only one element in the cs_array
                cs = cs_array;
                curveset = obj.curvesets{cs};
            end
                
            
            %{
            % Check to see if curveset is specified
            switch( upper(param_name) )
                case {'RESIDUE', 'INDEX', 'AA', 'DWH', 'DWX'}                    
                    if( obj.Ncs > 1 )
                        if( nargin == 6 )
                            cs          = varargin{1};
                            curveset	= obj.curvesets{cs};
                        else
                            % Get data for each curveset
                            for cs = 1:obj.Ncs
                                [DATA(cs), DATA_E(cs), IS_OK(cs)] = ...
                                    obj.getData(param_name, Temp, B0, QC_String, cs);
                            end                                
                            return
                            %error('Must specify curve set for %s', param_name);
                        end
                    else
                        cs          = 1;
                        curveset	= obj.curvesets{cs};
                    end
            end
            %}
            
            switch( upper(param_name) )
                case {'RESIDUE', 'INDEX', 'AA'}                    
                    DATA    = curveset.curves{1}.index;
                    DATA_E  = NaN;
                    IS_OK   = true;
                    
                case 'DWH'                        
                    % Check each curve for a MQ one
                    for c = 1:curveset.Nc
                        curve = curveset.curves{c};
                        
                        if( ~curve.SQX )
                            ctot = obj.getCurveTot(cs,c);
                            p = 1;
                            DATA   = obj.fitResult_Best.resultsMatrix(ctot,p);
                            DATA_E = obj.fitResult_Best.resultsMatrixErrors(ctot,p);
                            IS_OK  = obj.fitResult_Best.resultsMatrixIsOK(ctot,p);
                        end
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);
                    
                case {'DWHPPM','DWH_PPM','DWH PPM'}
                    % Check each curve for a MQ one
                    for c = 1:curveset.Nc
                        curve = curveset.curves{c};
                        
                        if( ~curve.SQX )
                            [DATA, DATA_E, IS_OK] = obj.fitResult_Best.getdwHppm( cs );
                            break;
                        end
                    end
                    
                case 'DWX'
                    ctot = obj.getCurveTot(cs,1);
                    p = 2;
                    DATA   = obj.fitResult_Best.resultsMatrix(ctot,p);
                    DATA_E = obj.fitResult_Best.resultsMatrixErrors(ctot,p);
                    IS_OK  = obj.fitResult_Best.resultsMatrixIsOK(ctot,p);
                    
                case {'DWXPPM','DWX_PPM','DWX PPM'}
                    [DATA, DATA_E, IS_OK] = obj.fitResult_Best.getdwXppm( cs );
                    
                case 'PA'
                    for ctot = 1:Nctot
                        [cs,c] = obj.getCurvesetCurve(ctot);
                        curve = obj.curvesets{cs}.curves{c};
                       
                        if( curve.Temp == Temp )
                            p = 3;
                            DATA   = obj.fitResult_Best.resultsMatrix(ctot,p);
                            DATA_E = obj.fitResult_Best.resultsMatrixErrors(ctot,p);
                            IS_OK  = obj.fitResult_Best.resultsMatrixIsOK(ctot,p);
                        end
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);
                    
                case 'KEX'
                    for ctot = 1:Nctot
                        [cs,c] = obj.getCurvesetCurve(ctot);
                        curve = obj.curvesets{cs}.curves{c};
                       
                        if( curve.Temp == Temp )
                            p = 4;
                            DATA   = obj.fitResult_Best.resultsMatrix(ctot,p);
                            DATA_E = obj.fitResult_Best.resultsMatrixErrors(ctot,p);
                            IS_OK  = obj.fitResult_Best.resultsMatrixIsOK(ctot,p);
                        end
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);
                    
                case 'R20'
                    for ctot = 1:Nctot
                        [cs,c] = obj.getCurvesetCurve(ctot);
                        curve = obj.curvesets{cs}.curves{c};
                       
                        if( curve.Temp == Temp && curve.B0 == B0 )
                            p = 5;
                            DATA   = obj.fitResult_Best.resultsMatrix(ctot,p);
                            DATA_E = obj.fitResult_Best.resultsMatrixErrors(ctot,p);
                            IS_OK  = obj.fitResult_Best.resultsMatrixIsOK(ctot,p);
                        end
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);

                case 'REX'
                    for ctot = 1:Nctot
                        [cs,c] = obj.getCurvesetCurve(ctot);
                        curve = obj.curvesets{cs}.curves{c};
                       
                        if( curve.Temp == Temp && curve.B0 == B0 )
                            p = 6;
                            DATA   = obj.fitResult_Best.resultsMatrix(ctot,p);
                            DATA_E = obj.fitResult_Best.resultsMatrixErrors(ctot,p);
                            IS_OK  = obj.fitResult_Best.resultsMatrixIsOK(ctot,p);
                        end
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);
                    
                case 'ALPHA'
                    for ctot = 1:Nctot
                        [cs,c] = obj.getCurvesetCurve(ctot);
                        curve = obj.curvesets{cs}.curves{c};
                       
                        if( curve.Temp == Temp && curve.B0 == B0 )
                            p = 7;
                            DATA   = obj.fitResult_Best.resultsMatrix(ctot,p);
                            DATA_E = obj.fitResult_Best.resultsMatrixErrors(ctot,p);
                            IS_OK  = obj.fitResult_Best.resultsMatrixIsOK(ctot,p);                            
                        end
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);
                    
                case 'PHIEXX'
                    for ctot = 1:Nctot
                        [cs,c] = obj.getCurvesetCurve(ctot);
                        curve = obj.curvesets{cs}.curves{c};
                       
                        if( curve.Temp == Temp && curve.B0 == B0 )
                            p = 8;
                            DATA   = obj.fitResult_Best.resultsMatrix(ctot,p);
                            DATA_E = obj.fitResult_Best.resultsMatrixErrors(ctot,p);
                            
                            p_DWX   = 2;
                            p_PA    = 3;
                            IS_OK  = and( obj.fitResult_Best.resultsMatrixIsOK(ctot,p_PA), ...
                                          obj.fitResult_Best.resultsMatrixIsOK(ctot,p_DWX) );
                                      
                            IS_OK = true;
                        end
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);
                    
                case 'KA'
                    % Find first temperature which is a match
                    t = find( Temp == obj.fitResult_Best.rateAnalysis.Temps, 1 );
                    if( ~isempty(t) )
                        DATA   = obj.fitResult_Best.rateAnalysis.kA(t);
                        DATA_E = obj.fitResult_Best.rateAnalysis.kA_E(t);
                        IS_OK  = obj.fitResult_Best.rateAnalysis.kA_isOK(t);
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);

                case 'KB'
                    % Find first temperature which is a match
                    t = find( Temp == obj.fitResult_Best.rateAnalysis.Temps, 1 );
                    if( ~isempty(t) )
                        DATA   = obj.fitResult_Best.rateAnalysis.kB(t);
                        DATA_E = obj.fitResult_Best.rateAnalysis.kB_E(t);
                        IS_OK  = obj.fitResult_Best.rateAnalysis.kB_isOK(t);
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);

                case 'K'
                    % Find first temperature which is a match
                    t = find( Temp == obj.fitResult_Best.rateAnalysis.Temps, 1 );
                    if( ~isempty(t) )
                        DATA   = obj.fitResult_Best.rateAnalysis.K(t);
                        DATA_E = obj.fitResult_Best.rateAnalysis.K_E(t);
                        IS_OK  = obj.fitResult_Best.rateAnalysis.K_isOK(t);
                    end
                    %error('Could not find %s at %0.1fK at %0.1fMHz in %s', param_name, Temp, B0, QC_String);

                case 'DH'
                    if( obj.fitResult_Best.CONSTRAIN_RATE_ANALYSIS )
                        % Global parameter for the group, just get it
                        DATA   = obj.fitResult_Best.dH;
                        DATA_E = obj.fitResult_Best.dH_E;
                        IS_OK  = obj.fitResult_Best.dH_isOK;
                    else
                        % Otherwise, obtain it from fit post-hoc fit kex/Pa
                        DATA   = obj.fitResult_Best.rateAnalysis.dH;
                        DATA_E = obj.fitResult_Best.rateAnalysis.dH_E;
                        IS_OK  = obj.fitResult_Best.rateAnalysis.vantHoff_isOK;
                    end
                    
                case 'DS'                    
                    if( obj.fitResult_Best.CONSTRAIN_RATE_ANALYSIS )
                        % Global parameter for the group, just get it
                        % dS is only OK if dH is
                        DATA   = obj.fitResult_Best.rateAnalysis.dS;
                        DATA_E = obj.fitResult_Best.rateAnalysis.dS_E;
                        IS_OK  = obj.fitResult_Best.dH_isOK;
                    else
                        % Otherwise, obtain it from fit post-hoc fit kex/Pa
                        DATA   = obj.fitResult_Best.rateAnalysis.dS;
                        DATA_E = obj.fitResult_Best.rateAnalysis.dS_E;
                        IS_OK  = obj.fitResult_Best.rateAnalysis.vantHoff_isOK;
                    end

                case {'EAB', 'EA(A->B)'}                    
                    if( obj.fitResult_Best.CONSTRAIN_RATE_ANALYSIS )
                        % Global parameter for the group, just get it
                        DATA   = obj.fitResult_Best.Eab;
                        DATA_E = obj.fitResult_Best.Eab_E;
                        IS_OK  = obj.fitResult_Best.Eab_isOK;
                    else
                        % ...otherwise need Arrhenius analysis
                        DATA   = obj.fitResult_Best.rateAnalysis.Eab;
                        DATA_E = obj.fitResult_Best.rateAnalysis.Eab_E;
                        IS_OK  = obj.fitResult_Best.rateAnalysis.arrhenius_isOK;
                    end                    

                case {'EBA', 'EA(B->A)'}
                    if( obj.fitResult_Best.CONSTRAIN_RATE_ANALYSIS )
                        % Global parameter for the group, just get it
                        % Only OK if Eab is OK
                        DATA   = obj.fitResult_Best.rateAnalysis.Eba;
                        DATA_E = obj.fitResult_Best.rateAnalysis.Eba_E;
                        IS_OK  = obj.fitResult_Best.Eab_isOK;
                    else
                        % ...otherwise need Arrhenius analysis
                        DATA   = obj.fitResult_Best.rateAnalysis.Eba;
                        DATA_E = obj.fitResult_Best.rateAnalysis.Eba_E;
                        IS_OK  = obj.fitResult_Best.rateAnalysis.arrhenius_isOK;
                    end
                    
                case {'PAB', 'P(A->B)'}
                    if( obj.fitResult_Best.CONSTRAIN_RATE_ANALYSIS )
                        % Global parameter for the group, just get it
                        % Only OK if Eab is OK
                        DATA   = obj.fitResult_Best.rateAnalysis.Pab;
                        DATA_E = obj.fitResult_Best.rateAnalysis.Pab_E;
                        IS_OK  = obj.fitResult_Best.Eab_isOK;
                    else                 
                        % ...otherwise need Arrhenius analysis
                        DATA   = obj.fitResult_Best.rateAnalysis.Pab;
                        DATA_E = obj.fitResult_Best.rateAnalysis.Pab_E;
                        IS_OK  = obj.fitResult_Best.rateAnalysis.arrhenius_isOK;
                    end
                    
                case {'PBA', 'P(B->A)'}
                    if( obj.fitResult_Best.CONSTRAIN_RATE_ANALYSIS )
                        % Global parameter
                        % Only OK if Eab is OK
                        DATA   = obj.fitResult_Best.rateAnalysis.Pba;
                        DATA_E = obj.fitResult_Best.rateAnalysis.Pba_E; 
                        IS_OK  = obj.fitResult_Best.Eab_isOK;                        
                    else
                        % ...otherwise need Arrhenius analysis
                        DATA   = obj.fitResult_Best.rateAnalysis.Pba;
                        DATA_E = obj.fitResult_Best.rateAnalysis.Pba_E; 
                        IS_OK  = obj.fitResult_Best.rateAnalysis.arrhenius_isOK;
                    end
                    
                otherwise
                    error('Invalid paramter name specified "%s". Check code for valid options', upper(param_name));
            end
            
            % Make sure the best fit is OK too
            IS_OK = and(IS_OK, obj.bestFitIsOK);
        end
        
        function [vcpmg_min, vcpmg_max, R2eff_min, R2eff_max] = getDataLimits( obj )
            %% Return vcpmg and R2eff limits on data by parsing each curve            
            vcpmg_min   = Inf;        % Smallest vcpmg value in all curves
            vcpmg_max   = -Inf;        % Largest vcpmg value in all curves
            R2eff_min   = Inf;        % Smallest R2eff value in all curves
            R2eff_max   = -Inf;        % Largest R2eff value in all curves
            
            % Check each curveset in the group
            for cs = 1:obj.Ncs
                curveset = obj.curvesets{cs};
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
        
        function [pMatrixLower, pMatrixUpper] = getDefaultFitLimits(obj)
            %% Return default fit limits for a group fit
                                   
            % Initialize upper and lower bound matrices
            Np = obj.parentSession.Np;
            pMatrixLower = zeros(obj.getNumCurves(), obj.parentSession.Np);
            pMatrixUpper = zeros(obj.getNumCurves(), obj.parentSession.Np);
                        
            % Iterate through each curve, and calculate limits
            ctot   = 0;
            for cs = 1:obj.Ncs
                curveset = obj.curvesets{cs};                
                for c = 1:curveset.Nc
                    curve = curveset.curves{c};
                    
                    % Mark which cs and c is in the current row
                    ctot = ctot+1;
            
                    pMatrixLower(ctot,1:Np) = ...
                        [ obj.parentSession.min_dwH_ppm * curve.B0 * 2*pi, ...
                          obj.parentSession.min_dwX_ppm * curve.B0 * 2*pi * curve.gammaX_relative, ...                            
                          obj.parentSession.min_pa, ...
                          obj.parentSession.min_kex, ...
                          obj.parentSession.min_R20 ];
                      
                    pMatrixUpper(ctot,1:Np) = ...
                        [ obj.parentSession.max_dwH_ppm * curve.B0 * 2*pi, ...
                          obj.parentSession.max_dwX_ppm * curve.B0 * 2*pi * curve.gammaX_relative, ...                            
                          obj.parentSession.max_pa, ...
                          obj.parentSession.max_kex, ...
                          obj.parentSession.max_R20 ];
                end
            end
        end
        
        function fileNameSuffix = getFileNameSuffix(obj)
            %% Return formatted string for the end of a filename
            
            % Remove any slash from the assignment (e.g., \delta_1)
            name = strrep(char(obj.name),'\','');
            % Replace space with underscore
            name = strrep(char(name),' ','_');
            
            fileNameSuffix = sprintf('%04d-%s', obj.index, name);            
        end
        
        function fitIndex = getFitIndex(obj, fitResult)
            %% Returns the index in the fit results array corresponding to
            % the argument
            
            % Check each element in the array and return the first hit
            for f = 1:obj.Nf
                if( obj.fitResults{f} == fitResult )
                    fitIndex = f;
                    return;
                end
            end
            % If it was not found, then return NaN
            fitIndex = NaN;
        end
        
        function [fitResult, isBestFit] = getFitResult( obj, f )
            %% Get the fit result
            fitResult   = obj.fitResults{f};
            isBestFit   = fitResult == obj.fitResult_Best;
        end
        
        function [fitStringArray, f_Best] = getFitNameArray(obj)
            %% Return cell array of curve strings
            % List each curve in each curve set in this group
            fitStringArray = {''};
            for f = 1:obj.Nf
                
                if( obj.fitResults{f} == obj.fitResult_Best )
                    f_Best      = f;
                    %obj.f_Best  = f;
                    marker      = '* ';
                else
                    marker      = '';
                end
                
                fitStringArray{f} = sprintf('%s%s', marker, obj.fitResults{f}.name);
            end
        end
        
        function [chi2_Grid, resultsMatrix_Grid, dH_Grid, Eab_Grid] = getGridResults( obj, INITIAL_OR_FINAL )
            %% Return the results from the grid search
            %  INITIAL_OR_FINAL = 'INITIAL' or 'FINAL' for parameter values
            %  chi2_Grid(igrid) is a 1D vector length Ngrid
            %  resultsMatrix_Grid(ctot,p,igrid) is a 3D matrix
            Ngrid               = length(obj.fitResults_Grid);
            Nctot               = obj.getNumCurves();
            Np                  = obj.parentSession.Np;
            Nr                  = obj.parentSession.Nr;
            
            chi2_Grid           = NaN*ones(Ngrid,1);
            resultsMatrix_Grid  = NaN*ones(Nctot,Nr,Ngrid);
            dH_Grid             = NaN*ones(Ngrid,1);
            Eab_Grid            = NaN*ones(Ngrid,1);
            
            switch upper(INITIAL_OR_FINAL)
                case 'INITIAL'
                    for igrid = 1:Ngrid
                        chi2_Grid(igrid)                        = obj.fitResults_Grid{igrid}.chi2;                
                        resultsMatrix_Grid(1:Nctot,1:Np,igrid)  = obj.fitResults_Grid{igrid}.pMatrixInit(1:Nctot,1:Np);
                        dH_Grid(igrid)                          = obj.fitResults_Grid{igrid}.dH_Init;
                        Eab_Grid(igrid)                         = obj.fitResults_Grid{igrid}.Eab_Init;
                    end
                    
                case 'FINAL'
                    for igrid = 1:Ngrid
                        chi2_Grid(igrid)                        = obj.fitResults_Grid{igrid}.chi2;                
                        resultsMatrix_Grid(1:Nctot,1:Nr,igrid)  = obj.fitResults_Grid{igrid}.resultsMatrix(1:Nctot,1:Nr);
                        dH_Grid(igrid)                          = obj.fitResults_Grid{igrid}.dH;
                        Eab_Grid(igrid)                         = obj.fitResults_Grid{igrid}.Eab;
                    end
                    
                otherwise
                    error('Specify either INITIAL or FINAL only');
            end
        end
        
        function Nc = getNumCurves( obj )
            %% Return the number of curves in the group (add from each curveset)            
            Nc = 0;
            for cs = 1:obj.Ncs
                Nc = Nc + obj.curvesets{cs}.Nc;
            end
        end
                        
        function Nip = getNumIndepParams( obj )
            %% Return the number of independent parameters in the group
            Nip = max(max( obj.indepParamIndexMatrix ));
        end
        
        function NTemps = getNumTemps( obj )
            %% Return number of temperatures
            NTemps = length(unique(obj.getTemperatureArray()));
        end
                
        function [Temp_array, cs_array, c_array] = getTemperatureArray( obj )
            %% Get array of temperatures from all curvesets in this group            
            Temp_array  = [];
            cs_array    = [];
            c_array     = [];
            ctot = 1;

            for cs = 1:obj.Ncs
                curveset = obj.curvesets{cs};
                for c = 1:curveset.Nc
                    curve = curveset.curves{c};

                    % Check / replace limits
                    Temp_array(ctot)   = curve.Temp;
                    cs_array(ctot)     = cs;
                    c_array(ctot)      = c;
                    ctot = ctot+1;
                end
            end            
        end
        
        function gridIsDone = gridSearch( obj, taskString, CONSTRAIN_RATE_ANALYSIS, handle_progress, handle_abort )
            %% Perform grid search to fit RD data with variety of initial conditions, and return updated fit_results
            % Ian Kleckner
            % Foster Lab
            % 2010/01/08 Create code
            % 2011/04/20 Update for GUARDD, classes
            % 2011/05/06 Constrain rates using dH and Eab
            % 
            % FUNCTION
            %  Perform grid search and fitting on RD data
            % 
            %  This will replace the prior grid search with the updated grid search
            %  Its a bit messy
            %  sse_local and p_matrix_fit are matrices which are as large as the largest
            %  grid search in all the curve sets. This is a waste of some space.
            % 
            % INPUT VARS
            %  taskString       = Either 'FIT' or 'SIM'
            %  handle_progress  = handle to progress bar from GUARDD GUI (allows real-time update)
            %  handle_abort     = handle to abort button from GUARDD GUI (allows real-time abort)
            %
            % OUTUPT VARS
            %  fit_results      = holds program fitting results (updated post grid-search)
            %  gridIsDone    = (boolean) is the grid search completed?
            %  
            % TO DO
            %  Permit grid to be run mutltiple times (i.e., do not overwrite prior grid, but simply add to it)
            %  Persistent variables?
            %  Inititalization element?
            % 
            %% Set up
            % Enable/disable output useful for debugging
            %{
            OUTPUT_DEBUG = obj.parentSession.OUTPUT_DEBUG_GRID_SEARCH;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end
            %}
            
            fprintf('\nStarting grid search...');
            
            % Update the fit parameters (matrices for indep. params)
            obj.updateFitParams(CONSTRAIN_RATE_ANALYSIS);

            % Copy some frequently used variables in the interest of code brevity
            Nctot   = obj.getNumCurves();            
            Np      = obj.parentSession.Np;
            
            F_RAND_MC_GUESS = obj.parentSession.F_RAND_MC_GUESS;
            F_RAND_MC_GUESS = 0
            
            % TODO % Add to grid, do not clear it first
            obj.fitResults_Grid = {};

            %% Initialize grid search
            gridIsDone = false;
            
            % Set up grid arrays for the four parameters in the grid
            dwHppm_values   = linspace( obj.grid_p_min(1), ...
                                        obj.grid_p_max(1), ...
                                        obj.grid_p_steps(1) );
                                    
            dwXppm_values   = linspace( obj.grid_p_min(2), ...
                                        obj.grid_p_max(2), ...
                                        obj.grid_p_steps(2) );
            % Convert % -> fraction
            Pa_values       = (1/100) .* linspace( obj.grid_p_min(3), ...
                                        obj.grid_p_max(3), ...
                                        obj.grid_p_steps(3) );
                                    
            kex_values      = linspace( obj.grid_p_min(4), ...
                                        obj.grid_p_max(4), ...
                                        obj.grid_p_steps(4) );
                                
            if( CONSTRAIN_RATE_ANALYSIS )
                % Only use the dH x Eab grid if necessary
                
                % Convert kcal/mol -> cal/mol
                dH_values       = 1e3 .* linspace( obj.grid_p_min(5), ...
                                            obj.grid_p_max(5), ...
                                            obj.grid_p_steps(5) );

                % Convert kcal/mol -> cal/mol
                Eab_values      = 1e3 .* linspace( obj.grid_p_min(6), ...
                                            obj.grid_p_max(6), ...
                                            obj.grid_p_steps(6) );
                                        
            else
                % Only use the dH x Eab grid if necessary
                dH_values   = 0;
                Eab_values  = 0;
            end

            % Number of steps in entire grid
            Ngrid   = length(dwHppm_values) * ...
                      length(dwXppm_values) * ...
                      length(Pa_values) * ...
                      length(kex_values) * ...
                      length(dH_values) * ...
                      length(Eab_values);                  
            
            % Update grid progress bar each Ngrid iterations
            Nupdates        = 20;
            igrid_Update    = floor(Ngrid/Nupdates);
            
            %TempArray   = obj.getTemperatureArray();
            %TempMin     = min(TempArray);
            %TempMax     = max(TempArray);
            
            % For finding the best fit from the grid
            chi2_Best   = Inf;
            igrid_Best  = NaN;            
            
            %% Set up progress bar (either internal or external)
            if( ishandle( handle_progress) )
                value = 0;
                h=handle_progress;        
                cla(h);        
                patch([0,value,value,0],[0,0,1,1],'k', 'Parent', h);
                axis(h,[0,1,0,1]); axis(h,'off'); drawnow;

            else
                % Open a new popup progress bar
                h_wb = waitbar(0, 'Performing grid search', ...
                    'CreateCancelBtn', 'setappdata(gcbf, ''cancel_button'', 1)', ...
                    'Name', 'Grid search progress');
                setappdata(h_wb, 'cancel_button', 0);    
            end
            
            %% Begin grid search            
            igrid   = 0;
            t0      = tic;
            for idwH = 1:length(dwHppm_values)
                for idwX = 1:length(dwXppm_values)
                    for iPa = 1:length(Pa_values)
                        for ikex = 1:length(kex_values)
                            for idH = 1:length(dH_values)
                                for iEab = 1:length(Eab_values)
                                    %% Update progress bar (only every igrid_Update iterations)
                                    igrid = igrid+1;
                                    % Calculate i_grid_local from other grid iterators and use parfor
                                    %fprintf('\b\b\b\b\b\b\b\b\b\b\b%05d/%05d', igrid, Ngrid);

                                    if( mod(igrid,igrid_Update)==0 )
                                        if( ishandle(handle_progress) )
                                            value = igrid/Ngrid;
                                            h=handle_progress;        
                                            cla(h);        
                                            patch([0,value,value,0],[0,0,1,1],'k', 'Parent', h);
                                            axis(h,[0,1,0,1]); axis(h,'off'); drawnow;

                                            % If user wants to abort
                                            if( ishandle(handle_abort) && get(handle_abort, 'Value') )
                                                fprintf('\nCanceled grid search\n');                                 
                                                gridIsDone = false;                        
                                                return
                                            end

                                        else
                                            % Update wait bar and provide opportunity for cancel
                                            waitbar(igrid/Ngrid,h_wb);
                                            if ( getappdata(h_wb,'cancel_button') )
                                                fprintf('\nCanceled grid search\n');
                                                gridIsDone = false;
                                                delete(h_wb);           
                                                return
                                            end
                                        end
                                    end                            

                                    %% Set up initial parameter matrix for each curve                                    
                                    switch upper(taskString)
                                        case 'FIT'
                                            fitResult = FitResult(obj, 'EXCHANGE', 'GRIDFIT', CONSTRAIN_RATE_ANALYSIS);
                                        case 'SIM'
                                            fitResult = FitResult(obj, 'EXCHANGE', 'GRIDSIM', CONSTRAIN_RATE_ANALYSIS);
                                        otherwise
                                            error('Bad taskString supplied, use either \''FIT\'' or \''SIM\''');
                                    end

                                    % Set bounds on fitting
                                    [pMatrixLower, pMatrixUpper] = obj.getDefaultFitLimits();
                                    fitResult.setLimits_pMatrix( pMatrixLower, pMatrixUpper );
                                                                        
                                    % CHEMICAL SHIFTS are same for each
                                    % curveset in the group
                                    % TODO % Alter this randomly for each cs ?
                                    dwHppm_csArray  = ones(obj.Ncs,1) .* dwHppm_values(idwH);
                                    dwXppm_csArray  = ones(obj.Ncs,1) .* dwXppm_values(idwX);
                                    fitResult.setInitial_Shifts(dwHppm_csArray, dwXppm_csArray)
                                    
                                    % Set the initial kinetics (kex and PA)
                                    if( CONSTRAIN_RATE_ANALYSIS )
                                        % Get initial rate constraints from
                                        % grid
                                        T0      = obj.T0_Grid;
                                        PA0     = Pa_values(iPa);
                                        kex0    = kex_values(ikex);
                                        dH      = dH_values(idH);
                                        Eab     = Eab_values(iEab);

                                        % Set initial conditions, and temperature-dependent parameters
                                        fitResult.setInitial_Kinetics_ConstrainedRates(T0, PA0, kex0, dH, Eab);

                                    else
                                        % Set kex and PA at each temperature
                                        Temp_Array  = unique( obj.getTemperatureArray() );
                                        TempMin     = min(Temp_Array);
                                        NTemp       = length(Temp_Array);                                        
                                        
                                        % kex doubles for each 10C (Arrhenius)
                                        N10C_Array  = (Temp_Array - TempMin)/10;
                                        kex_Array   = zeros(NTemp, 1);                                        
                                        kex_Array(:)= kex_values(ikex) .* 2.^N10C_Array;
                                        
                                        % Keep PA the same at each temperature
                                        PA_Array    = zeros(NTemp, 1);
                                        PA_Array(:) = Pa_values(iPa);

                                        % Set the initial conditions
                                        fitResult.setInitial_Kinetics_UnconstrainedRates(Temp_Array, PA_Array, kex_Array);
                                    end
                                    
                                    % Note: R20 is updated in setInitial_Shifts()

                                    %{
                                    % DEPRECATED CODE -- This was to
                                    % randomly alter the initial parameters
                                    % somewhat, and to reconstruct the
                                    % pMatrixInit
                                    % (I can't remember now why I thought I
                                    % needed to do the latter task)
                                    
                                    pMatrixInit = zeros(Nctot,Np);
                                    Nip = obj.getNumIndepParams();
                                    for ip = 0:Nip
                                        % Find each independent element in the matrix
                                        [ctotArray,pArray] = find(obj.indepParamIndexMatrix==ip);

                                        % Now set values in those rows and columns
                                        for iip = 1:length(ctotArray)
                                            ctot1   = ctotArray(1);
                                            ctot    = ctotArray(iip);
                                            p       = pArray(iip);

                                            [cs,c]  = obj.getCurvesetCurve(ctot);
                                            curve   = obj.curvesets{cs}.curves{c};

                                            % The first iteration is defined
                                            if( iip == 1 )
                                                % dwH ppm -> rad/sec
                                                if( p==1 )                                            
                                                    if( curve.SQX )
                                                        pMatrixInit(ctot,p) = 0;
                                                    else
                                                        pMatrixInit(ctot,p) = 2*pi*dwHppm_values(idwH) * curve.B0;
                                                    end

                                                % dwX ppm -> rad/sec
                                                elseif( p==2 )                                            
                                                    pMatrixInit(ctot,p) = 2*pi*dwXppm_values(idwX) * ...
                                                        curve.B0 * curve.gammaX_relative;

                                                % PA (%) -> (fraction)
                                                elseif( p==3 )
                                                    % PA randomly adjusted
                                                    PA = (PaPercent_values(iPa)/100) * ...
                                                        (1+F_RAND_MC_GUESS*rand()*(-1)^mod(ceil(10*rand()),2));

                                                    if( PA < 0.5 )
                                                        PA = 0.5;
                                                    elseif( PA > 1 )
                                                        PA = 0.999;
                                                    end
                                                    pMatrixInit(ctot,p) = PA;                                            

                                                % kex increases with temperature
                                                elseif( p==4 )
                                                    %{
                                                    % Randomly alter lowest temperature curve's kex
                                                    %  by as much as F_RAND_MC_GUESS
                                                    % Higher temperature curves
                                                    % will follow
                                                    % TODO % Fix initial calculation of kex
                                                    if( 0 && curve.Temp == TempMin )
                                                       pMatrixInit(ctot,p) = kex_values(ikex) * ...
                                                           (1+F_RAND_MC_GUESS*rand()*(-1)^mod(ceil(10*rand()),2));

                                                       % Store the kex at min temp
                                                       kex_TempMin = pMatrixInit(ctot,p);

                                                    else
                                                    %}
                                                        % kex doubles for each 10C (Arrhenius)
                                                        N10C = (curve.Temp-TempMin)/10;
                                                        pMatrixInit(ctot,p) = kex_values(ikex) * 2^N10C;
                                                    %end                                                

                                                % R20 is always independent, and easily estimated
                                                elseif( p==5 )
                                                    pMatrixInit(ctot,p) = min(curve.R2eff);
                                                end

                                             %  Remaining values are simply scaled
                                             %  from the first curve "ctot1" and
                                             %  their own scaling factor
                                            else
                                                pMatrixInit(ctot,p) = pMatrixInit(ctot1,p) * obj.paramScalingMatrix(ctot,p);
                                            end                                    

                                            %fprintf('\nCTOT=%d gX=%0.3f dwH=%0.1f dwX=%0.1f', ...
                                            %    ctot, curve.gammaX_relative, dwH, dwX);                                   
                                        end
                                    end

                                    fitResult.setFitInitialConditions( pMatrixInit );
                                    %}                                   
                                    

                                    %% Fit the group
                                    % Start the fit optimization
                                    updateGroupParamsFlag = false;

                                    switch upper(taskString)
                                        case 'FIT'
                                            fitResult.fitMe( updateGroupParamsFlag );
                                        case 'SIM'
                                            fitResult.simMe( updateGroupParamsFlag );
                                        otherwise
                                            error('Bad taskString supplied, use either \''FIT\'' or \''SIM\''');
                                    end                                    
                                    fitResult.analyzeMe();

                                    % Is this fit better than the current best?
                                    if( fitResult.chi2 < chi2_Best )
                                        chi2_Best       = fitResult.chi2;
                                        igrid_Best      = igrid;
                                    end

                                    %% Store it in the grid cell array
                                    obj.fitResults_Grid{igrid} = fitResult;
                                end
                            end
                        end
                    end
                end
            end
            
            %% Find the best fit from grid search and add it to the list
            fitResult_GridBest = obj.fitResults_Grid{igrid_Best};
            obj.addFitResult( fitResult_GridBest );
                        
            fprintf('\nGrid search done (%0.1f sec)\n', toc(t0));
            if( exist('h_wb','var') && ishandle(h_wb) )
                delete(h_wb);
            end
        end
        
        function gridIsDone = gridSearchIsDone( obj )
            %% (Boolean) the grid search is done
            gridIsDone = ~isempty(obj.fitResults_Grid);
        end
        
        function [isMQ, csFirstMQ] = isMQ( obj )
            %% (Boolean) The group contains a MQ curve
            
            isMQ        = false;
            csFirstMQ   = NaN;
            for cs = 1:obj.Ncs
                % If there is one MQ curveset, then the whole group is MQ, and return early
                if( obj.curvesets{cs}.isMQ() )
                    isMQ = true;
                    csFirstMQ = cs;
                    return
                end
            end
        end
        
        function noExFitIsDone = isNoExFitDone( obj )
            %% (Boolean) the no-exchange fit is done
            noExFitIsDone = ~isempty(obj.fitResult_NoEx);            
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
            fprintf(FILE, '\n%s%s\t%d', TAB,     'NumCurvesets', obj.Ncs);
            
            for cs = 1:obj.Ncs
                fprintf(FILE, '\n\n%sCurveset %d/%d', TAB, cs, obj.Ncs);
                
                curveset = obj.curvesets{cs};
                curveset.outputSpecs(FILE, '    ');
            end
        end
        
        function plotCurves( obj, handle, varargin )
            %% Plot each curve from the group on figure handle
            
            %LineSpec = getproperties(varargin);
            
            % Plot each curve from each curveset
            for cs = 1:obj.Ncs
                curveset = obj.curvesets{cs};
                for c = 1:curveset.Nc
                    curve = curveset.curves{c};
                    
                    errorbar( handle, curve.vcpmg, curve.R2eff, curve.eR2eff, varargin{:} );

                    % Enable curve overlay if there is more than one
                    if( cs == 1 && c == 1 )
                        hold(handle, 'on');
                    end
                end
            end            
        end
        
        function removeCurveset( obj, curveset )            
            %% Find the curve, then remove it
            for cs = 1:length(obj.curvesets)
                if( curveset == obj.curvesets{cs} )
                    % Here, "cs" is the element number to be eliminated
                    %  E.g., Remove eement #3 from a set of 5 elements
                    %  Set #3=#4, #4=#5, and #5=[], then Ntot=Ntot-1
                    for cs1 = cs:obj.Ncs-1
                        obj.curvesets{cs1} = obj.curvesets{cs1+1};
                    end

                    % Set final element to null
                    obj.curvesets(obj.Ncs) = [];

                    % Update the number of elements
                    obj.Ncs = obj.Ncs-1;
                    
                    % This will invalidate all the fits because it changes the group
                    obj.clearFits();
                    CONSTRAIN_RATE_ANALYSIS = false;
                    obj.updateFitParams(CONSTRAIN_RATE_ANALYSIS);
                    return
                end
            end
        end
        
        function removeFitResult( obj, fitResult )            
            %% Find the fitResult, then remove it
            for f = 1:obj.Nf
                if( fitResult == obj.fitResults{f} )
                    % Here, "f" is the element number to be eliminated
                    %  E.g., Remove eement #3 from a set of 5 elements
                    %  Set #3=#4, #4=#5, and #5=[], then Ntot=Ntot-1
                    for f1 = f:obj.Nf-1
                        obj.fitResults{f1} = obj.fitResults{f1+1};
                    end

                    % Set final element to null
                    obj.fitResults(obj.Nf) = [];

                    % Update the number of elements
                    obj.Nf = obj.Nf-1;
                    
                    % If the best fit was removed, set it to first fit
                    if( fitResult == obj.fitResult_Best )
                        obj.fitResult_Best  = obj.fitResults{1};
                        %obj.f_Best          = 1;
                    end

                    % Update statistics
                    obj.updateStatistics();                    
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
         
         function setBestFitIsOK( obj, bestFitIsOK )
             %% (Boolean) The best fit is OK
             obj.bestFitIsOK         = bestFitIsOK;            
         end         
         
         function setBestFitResult( obj, fitResult )
            %% Set the fit result as best (either numeric index OR class instance)
            if( isobject(fitResult) )
                % Check if it is listed already
                fIndex = obj.getFitIndex(fitResult);
                
                % If it is not yet listed...
                if( isnan(fIndex) )
                    % ...then add it
                    obj.addFitResult( fitResult );
                    obj.fitResult_Best = fitResult;
                    
                    % ...and mark it as best
                    obj.setBestFitResult(fitResult);
                else
                    % ...Otherwise just mark it as best
                    obj.fitResult_Best  = fitResult;
                    %obj.f_Best          = obj.Nf;
                end
                
            elseif( isnumeric(fitResult) )
                f = fitResult;
                
                if( f <= obj.Nf )
                    obj.fitResult_Best  = obj.fitResults{f};
                    %obj.f_Best          = f;
                end
            else
                error('Specify best fit result as FitResult object or integer index');
            end
         end
         
         function setConstrainRateAnalysis( obj, CONSTRAIN_RATE_ANALYSIS )
             %% Constrain the rate analysis (true/false)
             %obj.CONSTRAIN_RATE_ANALYSIS = CONSTRAIN_RATE_ANALYSIS;
             error('CONSTRAIN_RATE_ANALYSIS is now only a property of fitResult, not group');
         end
        
         function setExhibitsExchange( obj, exhibitsExchange )
             %% (Boolean) The group exhibits exchange             
             obj.exhibitsExchange    = exhibitsExchange;             
         end
        
         function obj = setIndex( obj, index )
            %% Set the index of the group
            obj.index = index;
         end
        
         function setGridLimits( obj, grid_p_min, grid_p_max, grid_p_steps )
             %% Set the limits on the grid search
             obj.grid_p_min     = grid_p_min;
             obj.grid_p_max     = grid_p_max;
             obj.grid_p_steps   = grid_p_steps;
         end
        
        function obj = setName( obj, name )
            %% Set the name of the group
            obj.name = name;
        end        
        
        function obj = setNote( obj, text )
            %% Set the notes for the group
            % Change newlines to new rows in the cell array
            if( isempty(obj.note) )
                obj.note = Note( multilineStringToCellArray(text) );
                %obj.note = Note( text );
            
            else
                obj.note.setText( multilineStringToCellArray(text) );
                %obj.note.setText( text );
            end
        end        
        
        function setT0_Grid(obj, T0_Grid )
            %% Set the temperature at which the grid search PA0 and kex0 are defined
            obj.T0_Grid = T0_Grid;
        end
        
        function sortCurvesets( obj )
            %% Sort the curvesets{} by index and atom
                        
            csIndexName  = cell(obj.Ncs,2);
            for cs = 1:obj.Ncs                              
               % Build a cell with COL1 = Index(String), COL2 = Name
               csIndexName{cs,1}  = sprintf('%10d', obj.curvesets{cs}.index);
               csIndexName{cs,2}  = obj.curvesets{cs}.name;
            end
                        
            % First sort by index, then by name
            [VOID, cssort2] = sortrows(csIndexName, [1,2]);
                        
            SORT_BACKWARDS = false;
            if( SORT_BACKWARDS )
               fprintf('\n\nSORTING BACKWARDS, See Group.sortCurvesets()');
               cssort2 = flipud(cssort2);               
               cssort2
            end
            
            % Use that sorted array to rearrange order of items
            new_curvesets = cell(1,obj.Ncs);
            for cs = 1:obj.Ncs
                % The index of csIndex which should go next                
                new_curvesets{cs} = obj.curvesets{ cssort2(cs) };
            end
            
            obj.curvesets = new_curvesets;
        end
        
        function updateFitParams( obj, CONSTRAIN_RATE_ANALYSIS )
            %% Identify the independent parameters and dependent scaling facotrs for the group fit            
            OUTPUT_DEBUG = obj.parentSession.OUTPUT_DEBUG_UPDATE_FIT_PARAMS;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end
            
            % Parameter types
            %  dwH -> 1H shift
            %  dwX -> AX shift
            %  PA  -> Population of A state
            %  kex -> Exchange rate A<->B
            %  R20 -> Intrinsic R2 relaxation rate
            
            % Group         -> define PA and kex for each temperature
            %  Curveset     -> define dwX
            %   Curve       -> define R20
            %    Curve(SQ)  -> define dwH = 0
            %    Curve(MQ)  -> define dwH =/= 0
            
            Np = 5;
            
            % Holds the independent parameter number for the cs,c,p in the group
            %  Some have the same number because they are shared
            obj.indepParamIndexMatrix   = zeros(obj.getNumCurves(), Np);
            
            % This holds the scaling factor for independent parameters
            obj.paramScalingMatrix      = zeros(obj.getNumCurves(), Np);
            ctot = 0;
            
            % Include each row with cs, and c
            % (or a link to curveset and curve)
            obj.cscMatrix = zeros(obj.getNumCurves(),2);
            
            % Start with matrix of curvesets and curves and parameters
            % Goal: count each independent parameter
            indepParam = 1;
            
            % List the temperatures so far encountered
            TempsEncountered    = [];
            NTemps              = length(unique(obj.getTemperatureArray()));
            
            %CONSTRAIN_RATE_ANALYSIS = obj.CONSTRAIN_RATE_ANALYSIS;
            
            %% dH and Eab - Temperature-dependence of thermodynamics and kintics
            if( CONSTRAIN_RATE_ANALYSIS && NTemps > 1 )
                % ...then itemize dΗ                
                if( OUTPUT_DEBUG )
                    fprintf('\n\t\tNumber of temperatures %d > 1', NTemps);
                    fprintf('\n\t\tItemizing parameter %d (dH)', indepParam);
                end
                obj.indepParam_dH = indepParam;
                indepParam = indepParam+1;
                
                % ...then itemize Eab          
                if( OUTPUT_DEBUG )
                    fprintf('\n\t\tItemizing parameter %d (Eab)', indepParam);
                end
                obj.indepParam_Eab = indepParam;
                indepParam = indepParam+1;
            end
                        
            for cs = 1:obj.Ncs
                curveset = obj.curvesets{cs};
                
                if( OUTPUT_DEBUG )
                    fprintf('\nWorking on CS=%d, %s', cs, curveset.name);
                end
                
                for c = 1:curveset.Nc
                    curve = curveset.curves{c};
                    
                    % Mark which cs and c is in the current row
                    ctot = ctot+1;
                    obj.cscMatrix(ctot,1) = cs;
                    obj.cscMatrix(ctot,2) = c;
                    
                    if( OUTPUT_DEBUG )
                        fprintf('\n\tWorking on C=%d (CTOT=%d), %s', c, ctot, curve.getSpecsString());                    
                    end
                        
                    % Check if the parameter has been itemized yet
                    % If independent
                    %  => add it to the linearized list
                    % Else
                    %  => mark which element number it is in the independent list
                    %  => mark the scaling factor that multiples the list value
                    
                    %%  (1) dwH -> 1H shift difference between A and B
                    % Unique dwH and dwX for each curveset 
                    
                    % Check if the curveset is MQ, and what the first
                    % curve is that makes it such, from that curveset
                    [csIsMQ, cFirstMQ ] = curveset.isMQ();
                    
                    % dwH only for MQ curve (i.e., NOT SQ)
                    if( csIsMQ && ~curve.SQX )
                        
                        % If this is the first MQ curve in the set...
                        if( c == cFirstMQ )
                            % ...then itemize dwH
                            p                               = 1;
                            obj.indepParamIndexMatrix(ctot, p)  = indepParam;

                            scale_factor                    = 1;
                            obj.paramScalingMatrix(ctot,p)      = scale_factor;             
                            if( OUTPUT_DEBUG )
                                fprintf('\n\t\tItemizing parameter %d (dwH @ CS%d, C%d)', indepParam, cs, c);
                            end
                            indepParam = indepParam+1;

                        else
                            % ...otherwise refer to the proper dwH from cFirstMQ
                            
                            % Convert the cFirstMQ to a ctot number in the group
                            ctotFirstMQ = obj.getCurveTot(cs, cFirstMQ);
                            
                            p                               = 1;
                            linkParam                       = obj.indepParamIndexMatrix(ctotFirstMQ, p);
                            obj.indepParamIndexMatrix(ctot, p)  = linkParam;

                            scale_factor                    = curve.B0 / curveset.curves{cFirstMQ}.B0;
                            obj.paramScalingMatrix(ctot,p)      = scale_factor;
                            if( OUTPUT_DEBUG )
                                fprintf('\n\t\tLinking dwH to parameter %d, scaled by %fx', ...
                                    linkParam, scale_factor);
                            end
                        end
                        
                    else
                        p                               = 1;
                        linkParam                       = 0;
                        obj.indepParamIndexMatrix(ctot, p)  = linkParam;

                        scale_factor                    = 0;
                        obj.paramScalingMatrix(ctot,p)  = scale_factor;
                        if( OUTPUT_DEBUG )
                            fprintf('\n\t\tFixing dwH to zero (SQ curve)');
                        end
                        
                    end
                                   
                    %% (2) dwX -> AX shift difference between A and B
                    % Unique dwH and dwX for each curveset 
                    
                    % If this is a new curveset...
                    if( c == 1 )
                        % ...then itemize dwX
                        p                               = 2;
                        obj.indepParamIndexMatrix(ctot, p)  = indepParam;
                        
                        scale_factor                    = 1;
                        obj.paramScalingMatrix(ctot,p)      = scale_factor;    
                        if( OUTPUT_DEBUG )
                            fprintf('\n\t\tItemizing parameter %d (dwX @ CS%d, C%d)', indepParam, cs, c);
                        end
                        indepParam = indepParam+1;
                        
                    else
                        % ...otherwise refer to the proper dwX from c=1
                        
                        % Convert the c=1 curve in this curvset to the ctot
                        %  curve number in the group
                        ctot_c1 = obj.getCurveTot(cs, 1);
                            
                        p                               = 2;
                        linkParam                       = obj.indepParamIndexMatrix(ctot_c1, p);
                        obj.indepParamIndexMatrix(ctot, p)  = linkParam;
                        
                        scale_factor                    = curve.B0 / curveset.curves{1}.B0;
                        obj.paramScalingMatrix(ctot,p)      = scale_factor;
                        if( OUTPUT_DEBUG )
                            fprintf('\n\t\tLinking dwX to parameter %d, scaled by %fx', ...
                                linkParam, scale_factor);
                        end
                    end
                    
                    %%  (3) PA  -> Population of A state
                    %   (4) kex -> Exchange rate A<->B
                    % Unique PA and kex for each temperature                    
                    
                    % Check if this temperature has been encountered
                    isSameTemp = (curve.Temp == TempsEncountered);
                    
                    % What "ctot" index was the first appearance of this
                    % temperature?
                    ctotSameTemp = find(isSameTemp);
                    if( isempty(ctotSameTemp) )
                        ctotSameTemp = 0;
                    else
                        ctotSameTemp = ctotSameTemp(1);
                    end
                    
                    % Itemize that this temperature has been encountered
                    TempsEncountered(ctot) = curve.Temp;
                    
                    % If this is a new temperature...
                    if( ctotSameTemp == 0 )
                        
                        % If there is only one temperature
                        % OR if this is the first temperature
                        if( ~CONSTRAIN_RATE_ANALYSIS || ...
                            ( NTemps == 1 || length(TempsEncountered)==1 ) )
                            % ...then itemize PA
                            p                                   = 3;
                            obj.indepParamIndexMatrix(ctot, p)  = indepParam;

                            scale_factor                        = 1;
                            obj.paramScalingMatrix(ctot,p)      = scale_factor;
                            if( OUTPUT_DEBUG )
                                fprintf('\n\t\tItemizing parameter %d (PA @ %dK) *PA0*', indepParam, curve.Temp);
                            end
                            obj.indepParam_PA0 = indepParam;
                            indepParam = indepParam+1;

                            % ... and kex
                            p                                   = 4;
                            obj.indepParamIndexMatrix(ctot, p)  = indepParam;

                            scale_factor                        = 1;
                            obj.paramScalingMatrix(ctot,p)      = scale_factor;
                            if( OUTPUT_DEBUG )
                                fprintf('\n\t\tItemizing parameter %d (kex @ %dK) *kex0*', indepParam, curve.Temp);
                            end
                            obj.indepParam_kex0 = indepParam;
                            indepParam = indepParam+1;
                            
                            % Mark the temperature for PA0 and kex0
                            obj.T0 = curve.Temp;
                        else
                            % Multiple temperatures require a parameters dH
                            % and Ea to determine PA and kex, given the Temp
                            
                            % For PA
                            p                                   = 3;
                            linkParam                           = obj.indepParam_PA0;
                            obj.indepParamIndexMatrix(ctot, p)  = linkParam;

                            scale_factor                        = NaN;
                            obj.paramScalingMatrix(ctot,p)      = scale_factor;
                            
                            if( OUTPUT_DEBUG )
                                fprintf('\n\t\tLinking PA @ %dK to PA0 @ T0=%dK (param %d) via Temp (%dK), dH (param %d), and Eab (param %d)', ...
                                    curve.Temp, obj.T0, linkParam, curve.Temp, obj.indepParam_dH, obj.indepParam_Eab);
                            end

                            % and kex
                            p                                   = 4;
                            linkParam                           = obj.indepParam_kex0;
                            obj.indepParamIndexMatrix(ctot, p)  = linkParam;

                            scale_factor                        = NaN;
                            obj.paramScalingMatrix(ctot,p)      = scale_factor;
                            if( OUTPUT_DEBUG )
                                fprintf('\n\t\tLinking kex @ %dK to kex0 @ T0=%dK (param %d) via Temp (%dK), dH (param %d), and Eab (param %d)', ...
                                    curve.Temp, obj.T0, linkParam, curve.Temp, obj.indepParam_dH, obj.indepParam_Eab);
                            end
                        end
                        
                    else
                        % ...otherwise refer to the proper PA and kex
                        % Get the first curve at which this temperature was
                        % found, then link to that parameter
                        
                        % For PA
                        p                                   = 3;
                        linkParam                           = obj.indepParamIndexMatrix(ctotSameTemp, p);
                        obj.indepParamIndexMatrix(ctot, p)  = linkParam;
                        
                        scale_factor                        = 1;
                        obj.paramScalingMatrix(ctot,p)      = scale_factor;
                        if( OUTPUT_DEBUG )
                            fprintf('\n\t\tLinking PA to parameter %d, scaled by %fx', ...
                                linkParam, scale_factor);
                        end
                        
                        % and kex
                        p                                   = 4;
                        linkParam                           = obj.indepParamIndexMatrix(ctotSameTemp, p);
                        obj.indepParamIndexMatrix(ctot, p)  = linkParam;
                        
                        scale_factor                        = 1;
                        obj.paramScalingMatrix(ctot,p)      = scale_factor;
                        if( OUTPUT_DEBUG )
                            fprintf('\n\t\tLinking kex to parameter %d, scaled by %fx', ...
                                linkParam, scale_factor);
                        end
                    end
                    
                    %%  (5) R20 -> Intrinsic R2 relaxation rate 
                    % Unique R20 for each curve
                    p                                   = 5;
                    obj.indepParamIndexMatrix(ctot, p)  = indepParam;

                    scale_factor                        = 1;
                    obj.paramScalingMatrix(ctot,p)      = scale_factor;         
                    if( OUTPUT_DEBUG )
                        fprintf('\n\t\tItemizing parameter %d (R20 @ CS%d, C%d)', indepParam, cs, c);
                    end
                    indepParam = indepParam+1;
                end
            end            
        end
        
        function updateStatistics( obj )
            %% Calculate statistics compared to no-exchange fit
            
        end        
    end    
end

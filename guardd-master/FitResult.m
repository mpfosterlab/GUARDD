classdef FitResult < handle
    % FitResult class: for storing results of fit session
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/04/12 Start coding
    % 2011/09/02 Update resultsMatrix to include PhiexX
    % 2012/01/11 Quick and dirty output the chi2 for each curve in fitMe()
            
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set        
        name            = ''    % Name of the fit (E.g., IC Chi2=37.69 [Apr 18,2011-17:35:34])
        modelName       = '';   % Name of model used to fit (E.g., Exch, or NoEx)
        fitModeString   = '';   % Name of fitting mode (IC, Sim, Grid)
        
        % Pre-fitting
        pMatrixInit    = [];   % Initial fitting parameters
        pMatrixLower   = [];   % Paramter matrix (cstot,p) for lower bound
        pMatrixUpper   = [];   % Paramter matrix (cstot,p) for upper bound
        
        % (TDR) Temperature-dependence of rates should be calculated from dH and Eab
        % TODO % Set FitResult CONSTRAIN_RATE_ANALYSIS somehow
        CONSTRAIN_RATE_ANALYSIS = false;

        % (TDR) Initial guess for enthalpy and activation energy
        dH_Init         = NaN;
        Eab_Init        = NaN;
        
        % Post-fitting
        timeRequired        = 0;    % Time required for fitting
        pMatrixBest         = [];   % Parameters to fit all curves in all curvesets in the group        
        chi2                = NaN;  % The chi2 goodness of fit metric for the fit
        df                  = NaN;  % Degrees of freedom
        chi2red             = NaN;  % chi2/df
        
        pMatrixBestErrors   = [];   % Parameters to fit all curves in all curvesets in the group
        
        % (TDR) Final result for enthalpy and activation energy
        kex0    = NaN;
        PA0     = NaN;
        dH      = NaN;
        Eab     = NaN;

        % Statistics
        Fstatistic      = NaN;
        Pvalue          = NaN;
        
        % Post-fitting analysis
        Np                  = NaN;
        Nr                  = NaN;
        resultsMatrix       = [];   % (ctot, p)
                                    % dwH, dwX, PA, kex, R20, Rex, alpha
                                    % Rex = Exchange broadening
                                    % alpha = B0-dependence of Rex                                    
                                    % PhiexX = PA(1-PA)dwX^2
        rateAnalysis        = []    % RateAnalysis instance for temperature-dependent analysis

        % Error analysis from Monte Carlo
        fitResults_MC       = {};   % Array of FitResults for storing Monte Carlo error results
        resultsMatrixErrors = [];   % Errors in fitted parameters        
        dH_E                = NaN;
        Eab_E               = NaN;
        
        % (Boolean) the parameter is OK
        resultsMatrixIsOK   = [];
        kex0_isOK           = false;
        PA0_isOK            = false;      
        dH_isOK             = false;
        Eab_isOK            = false;
        
        % Generic
        parentGroup         = [];  % Link to the group that is fit
        pNames              = {'dwH', 'dwX', 'PA', 'kex', 'R20'};
        
        
    end
    
    properties (SetAccess = private, GetAccess = private)
        % Cannot read nor write via ".", except within class functions
    end

    methods
        function obj = FitResult( group, modelName, fitModeString, CONSTRAIN_RATE_ANALYSIS )
            %% Constructor function (set group, name, create rate analysis)
            % CONSTRAIN_RATE_ANALYSIS => (Bool) use constrain on dH and Eab
            % instead of fitting (kex,PA) at each temperature
            obj.parentGroup     = group;
            
            obj.Np = group.parentSession.Np;
            obj.Nr = group.parentSession.Nr;
            
            % Format the model name for consistency            
            switch upper(modelName)
                case { 'NOEXCHANGE', 'NOEX' }
                    obj.modelName = 'NOEXCHANGE';
                    
                case {'EXCHANGE', 'TWOSTATE', 'A<->B'}
                    obj.modelName = 'EXCHANGE';
                    
                otherwise
                    error('Please specify proper modelName EXCHANGE or NOEXCHANGE');
            end
            
            % Format the model name for consistency
            switch upper(fitModeString)
                case {'NOEX'}
                    obj.fitModeString = 'NoEx';
                    
                case {'FIT'}
                    obj.fitModeString = 'FIT-1';
                    
                case {'SIM'}
                    obj.fitModeString = 'SIM-1';
                    
                case {'GRIDFIT'}
                    obj.fitModeString = 'FIT-G';
                    
                case {'GRIDSIM'}
                    obj.fitModeString = 'SIM-G';
                    
                otherwise
                    error('Please specify proper fitModeString NoEx, FIT, SIM, FITGRID, SIMGRID');
            end
            
            if(CONSTRAIN_RATE_ANALYSIS)
                obj.fitModeString = sprintf('%s[CR]', obj.fitModeString);
            else
                obj.fitModeString = sprintf('%s[--]', obj.fitModeString);
            end
            
            
            % Constrain the rate analysis?
            obj.CONSTRAIN_RATE_ANALYSIS = CONSTRAIN_RATE_ANALYSIS;
            
            % Create the rate analysis
            obj.rateAnalysis = RateAnalysis( obj );
        end
        
        function analyzeMe( obj )
            %% Analyze the fitResult (usually called after fitMe()
            % 2011/06/30 Update for calculating alpha properly (using dwX
            % not B0)
            
            % Copy some frequently used variables in the interest of code brevity            
            Np      = obj.Np;
            Nr      = obj.Nr;
            Nctot   = obj.parentGroup.getNumCurves();
            
            % Start timing
            tic0    = tic();
            fprintf('\nStarting analysis...');
            
            switch upper(obj.modelName)                
                case 'NOEXCHANGE'
                    % Statistics (compare to no-exchange fit, itself)
                    obj.Fstatistic  = 1;
                    obj.Pvalue      = fpdf(obj.Fstatistic, obj.df, obj.df);

                    % Build the results matrices (some already set)
                    % dwH, dwX, PA, kex, R20, Rex, alpha
                    obj.resultsMatrixErrors         = zeros(Nctot,Nr);
                    obj.resultsMatrixIsOK           = zeros(Nctot,Nr);

                    obj.resultsMatrix               = zeros(Nctot,Nr);
                    obj.resultsMatrix(1:Nctot,1:Np) = obj.pMatrixBest(1:Nctot,1:Np);
                    
                    % Rex is zero for no-exchang emodel
                    obj.resultsMatrix(1:Nctot,6)        = 0;
                    obj.resultsMatrixErrors(1:Nctot,6)  = 0;
                    
                    % alpha is undefined for no-exchange model
                    obj.resultsMatrix(1:Nctot,7)        = NaN;
                    obj.resultsMatrixErrors(1:Nctot,7)  = NaN;
                    
                    % PhiexX is zero
                    obj.resultsMatrix(1:Nctot,8)        = NaN;
                    obj.resultsMatrixErrors(1:Nctot,8)  = NaN;

                    % Perform rate analysis
                    obj.rateAnalysis.analyzeMe();                    
                    
                case 'EXCHANGE'                    
                    %% Calculate post-fitting quantities (ANOVA,Rex,alpha,Temp-analysis)

                    % Statistics (compare to no-exchange fit)
                    df_NoEx         = obj.parentGroup.fitResult_NoEx.df;
                    chi2red_NoEx    = obj.parentGroup.fitResult_NoEx.chi2red;            
                    obj.Fstatistic  = chi2red_NoEx / obj.chi2red;
                    obj.Pvalue      = fpdf(obj.Fstatistic, df_NoEx, obj.df);

                    % Build the results matrices
                    % dwH, dwX, PA, kex, R20, Rex, alpha
                    obj.resultsMatrixErrors         = zeros(Nctot,Nr);
                    obj.resultsMatrixIsOK           = zeros(Nctot,Nr);

                    obj.resultsMatrix               = zeros(Nctot,Nr);
                    obj.resultsMatrix(1:Nctot,1:Np) = obj.pMatrixBest(1:Nctot,1:Np);
                    
                    % Iterate through each curve for Rex and alpha
                    ctot = 0;
                    for cs = 1:obj.parentGroup.Ncs
                        curveset = obj.parentGroup.curvesets{cs};   
                                                
                        % Initialize arrays for calculating alpha
                        %  alpha calculated for each curveset
                        RexArray                        = zeros(curveset.Nc,1);
                        RexArray_E                      = zeros(curveset.Nc,1);
                        dwXArray                        = zeros(curveset.Nc,1);
                        B0Array                         = zeros(curveset.Nc,1);
                        TempArray                       = zeros(curveset.Nc,1);
                        QuantumCoherenceStringArray     = cell(curveset.Nc,1);

                        % Rex is calculated for each curve
                        for c = 1:curveset.Nc
                            ctot = ctot+1;
                            curve = curveset.curves{c};

                            % Calculate Rex
                            obj.resultsMatrix(ctot,6) = ...
                                calculate_Rex( curve.TCPMG, obj.pMatrixBest(ctot,1:Np) );                    

                            % Build arrays to calculate alpha
                            RexArray(c)                          = obj.resultsMatrix(ctot,6);
                            RexArray_E(c)                        = NaN;
                            B0Array(c)                           = curve.B0;
                            %dwXArray(c)                          = obj.resultsMatrix(ctot,2) / 2 / 3.14;   
                            % Same results obtained using B0 or rad/s for dwX
                            %dwXArray(c)                          = obj.resultsMatrix(ctot,2);                            
                            dwXArray(c)                          = curve.B0;                     
                            TempArray(c)                         = curve.Temp;
                            if( curve.SQX )
                                QuantumCoherenceStringArray{c}   = 'SQX';
                            else
                                QuantumCoherenceStringArray{c}   = 'MQX';
                            end
                        end
                        
                        % Alpha is calculated via Rex as function of B0 at each temperature and coherence
                        % Needs an array of Rex, B0, Temp, and quantum coherence            
                        [alphaArray, alphaArrayErrors] = calculate_alpha( RexArray, ...
                            RexArray_E, dwXArray, B0Array, TempArray, QuantumCoherenceStringArray, ...
                            obj.parentGroup.parentSession.OUTPUT_DEBUG_CALCULATE_ALPHA, curveset.name );

                        % Get all the ctot's associated with the cs
                        ctot_array = obj.parentGroup.getCurveTotArray(cs);
                        obj.resultsMatrix(ctot_array,7)        = alphaArray;
                        obj.resultsMatrixErrors(ctot_array,7)  = alphaArrayErrors;
                        
                        
                        % Calculate PhiexX [2011/09/02]
                        dwX     = Session.convertUnits( obj.resultsMatrix(ctot_array,2), 'DWX', 'DISPLAY');
                        dwX_E   = Session.convertUnits( obj.resultsMatrixErrors(ctot_array,2), 'DWX', 'DISPLAY');
                        
                        Pa      = obj.resultsMatrix(ctot_array,3);
                        Pa_E    = obj.resultsMatrixErrors(ctot_array,3);
                        
                        obj.resultsMatrix(ctot_array,8)        = Pa.*(1-Pa).*dwX.^2;
                        obj.resultsMatrixErrors(ctot_array,8)  = Pa.*(1-Pa).*dwX.^2 .* ...
                            sqrt( (Pa_E./Pa).^2 + (Pa_E./Pa).^2 + (dwX_E./dwX).^2 + (dwX_E./dwX).^2 );
                    end
                    
                    % Perform rate analysis
                    obj.rateAnalysis.analyzeMe();
                    
                otherwise
                    error('Bad modelName specified, found during analyzeMe()');            
            end
            fprintf('done (%0.1f sec)\n', toc(tic0));
        end
        
        function calculateErrors( obj, handle_progress, handle_abort )
            %% Estimate error in dispersion fit using Monte Carlo bootstrapping
            % Ian Kleckner
            % Foster Lab
            % 2010/01/27 Create code
            % 2011/04/25 Upgrade for GUARDD with classes
            % 2011/05/06 Update for errors on dH and Ea
            % 
            % FUNCTION
            %  Estimate error in dispersion fit using bootstrapping            
            % 
            % INPUT VARS
            %  handle_progress = handle to progress bar from GUARDD Fit RD GUI (allows real-time update)
            %  handle_abort = handle to abort button from GUARDD Fit RD GUI (allows real-time abort)
            %
            % OUTUPT VARS
            %  p_matrix_error(c,p)  = parameter errors from standard deviation of array of fit results
            %  Rex_error(c)         = error in Rex from stdev of array of fit results
            %  chi2_error            = stdev of chi^2 values from array fit residuals
            %  errors_finished      = (boolean) did the process get completed?
            % 
            % TO DO
            %  Nothing!
            
            % Debugging
            OUTPUT_DEBUG = obj.parentGroup.parentSession.OUTPUT_DEBUG_ERRORS;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end
            
            %% Set up            
            %  Note: this can only be called once the fit has been made
            if( isempty(obj.pMatrixBest) )
                error('Cannot compute errors without a fit. Call fitMe() first');
            end
            
            % Copy some frequently used variables for code brevity            
            %errors_finished     = false;
            Nctot                   = obj.parentGroup.getNumCurves();
            Np                      = obj.parentGroup.parentSession.Np;            
            
            F_RAND_MC_GUESS         = obj.parentGroup.parentSession.F_RAND_MC_GUESS;            
            Nmc                     = obj.parentGroup.parentSession.Nmc;              
            MC_RANDOMIZE_MODE  = obj.parentGroup.parentSession.MC_RANDOMIZE_MODE;
            
            
            fprintf('\nStarting MC error analysis (%d iterations)', Nmc);
            
            if( strcmpi(MC_RANDOMIZE_MODE,'A') )
                fprintf(' using method A: Fitted residuals...');
            else
                fprintf(' using method B: Experimental errors...');
            end
            
            % Cell array for the Monte Carlo fit results
            obj.fitResults_MC = cell(Nmc,1);
            
            % Update grid progress bar each Nmc iterations
            Nupdate = floor(Nmc/20);            
                       
            if( OUTPUT_DEBUG )
                hf = figure;
                hax = axes;                
                hold(hax, 'on');
                
                set( hax, 'FontName', obj.parentGroup.parentSession.FONTNAME, ...
                    'FontSize', obj.parentGroup.parentSession.FONTSIZE_MEDIUM, ...
                    'FontWeight', 'bold', 'LineWidth', obj.parentGroup.parentSession.LINEWIDTH)

                box( hax, 'on');
                xlabel( hax,'\nu_{CPMG} (Hz)')
                ylabel( hax,'R_2^{Eff} (Hz)');
                title(  hax, sprintf('MC Errors\nGroup: %s', obj.parentGroup.name) );                
                F_RAND_MC_GUESS
            end
                        
            %% Get errors for Monte Carlo data randomization (method A or B)
            % Two methods permitted as of 2011/09/11
            residuals_Mean   = zeros(Nctot,1);
            residuals_Std    = zeros(Nctot,1);
            R2eff_Model      = cell(Nctot,1);
            
            for ctot = 1:Nctot
                % Get curveset and curve number
                [cs,c] = obj.parentGroup.getCurvesetCurve( ctot );
                
                % Get model curve fit
                vcpmg   = obj.parentGroup.curvesets{cs}.curves{c}.vcpmg;
                R2eff   = obj.parentGroup.curvesets{cs}.curves{c}.R2eff;                
                TCPMG   = obj.parentGroup.curvesets{cs}.curves{c}.TCPMG;
                params  = obj.pMatrixBest(ctot,1:Np);

                % Store the modeled data
                R2eff_Model{ctot} = model_MQRD_CRJ( vcpmg, TCPMG, params );
                
                % Method A: Residuals from fit
                if( strcmpi(MC_RANDOMIZE_MODE,'A') )                    
                    % Store the residuals for the fit (DATA - MODEL)
                    residuals = R2eff - R2eff_Model{ctot};
                    
                    %[void, void, residuals] = chi2_MQRD_CRJ( vcpmg, R2eff,
                    %eR2eff, TCPMG, params, Nfixed );                    
                    
                % Method B: Get experimental error on each curve
                else
                    % ! These are called "residuals" but only to maintain
                    % nomenclature from old code
                    %  They are just experimental errors                    
                    residuals = obj.parentGroup.curvesets{cs}.curves{c}.eR2eff;
                end
                
                % Calculate statistics on the residuals                
                residuals_Mean(ctot) = mean(residuals);
                residuals_Std(ctot)  = std(residuals);
                
                if( OUTPUT_DEBUG )
                    plot(hax, vcpmg, R2eff, 'ok');
                    plot(hax, vcpmg, residuals, 'or');                    
                end
            end
                                    
            %% Set up progress bar (either internal or external)
            if( ishandle( handle_progress) )
                value = 0;
                h=handle_progress;        
                cla(h);        
                patch([0,value,value,0],[0,0,1,1],'k', 'Parent', h);
                axis(h,[0,1,0,1]); axis(h,'off'); drawnow;

            else
                % Open a new popup progress bar
                h_wb = waitbar(0, 'Performing Monte Carlo error estimation', ...
                    'CreateCancelBtn', 'setappdata(gcbf, ''cancel_button'', 1)', ...
                    'Name', 'Error estimation progress');
                setappdata(h_wb, 'cancel_button', 0);    
            end
            
            %% Begin iterative Monte Carlo fitting process           
            
            % Start timing the process
            t0 = tic();
            for imc = 1:Nmc                
                %% Update progress bar (only every Nupdate iterations)
                % Calculate i_grid_local from other grid iterators and use parfor
                %fprintf('\b\b\b\b\b\b\b\b\b\b\b%05d/%05d', imc, Nmc);

                if( mod(imc,Nupdate)==0 )
                    if( ishandle(handle_progress) )
                        value = imc/Nmc;
                        h=handle_progress;        
                        cla(h);        
                        patch([0,value,value,0],[0,0,1,1],'k', 'Parent', h);
                        axis(h,[0,1,0,1]); axis(h,'off'); drawnow;

                        % If user wants to abort
                        if( ishandle(handle_abort) && get(handle_abort, 'Value') )
                            fprintf('\nCanceled error estimation\n');
                            errors_finished = false;                        
                            return
                        end

                    else
                        % Update wait bar and provide opportunity for cancel
                        waitbar(imc/Nmc,h_wb);
                        if ( getappdata(h_wb,'cancel_button') )
                            fprintf('\nCanceled grid search\n');
                            errors_finished = false;
                            delete(h_wb);           
                            return
                        end
                    end
                end                
                
                %% Generate synthetic data from Monte Carlo randomization
                R2eff_MC = cell(Nctot,1);
                for ctot = 1:Nctot
                    % Get curveset and curve number
                    [cs,c] = obj.parentGroup.getCurvesetCurve( ctot );
                    
                    % Resample a residual vector from that Gaussian distribution
                    Nobs                = obj.parentGroup.curvesets{cs}.curves{c}.Nobs;
                    residuals_Synthetic = random('norm', residuals_Mean(c), residuals_Std(c), 1, Nobs);

                    % Add this to the model data to yield synthetic observation set
                    R2eff_MC{ctot} = R2eff_Model{ctot} + residuals_Synthetic;
                    
                    % Debugging output to overlay the simulated curves
                    if( OUTPUT_DEBUG )
                        vcpmg   = obj.parentGroup.curvesets{cs}.curves{c}.vcpmg;
                        plot(hax, vcpmg, R2eff_MC{ctot}, '.b');                        
                    end
                end
                
                %% Set up fit for this Monte Carlo iteration
                fitResult = FitResult(obj.parentGroup, 'EXCHANGE', 'GRIDFIT', obj.CONSTRAIN_RATE_ANALYSIS);
                
                % Set current fit limits
                fitResult.setLimits_pMatrix( obj.pMatrixLower, obj.pMatrixUpper );
                
                %{
                % 12/10/09 Randomize initial guess for fitting (more accurate errors)    
                pMatrixInit_MC = zeros(Nctot,Np);
                for ctot = 1:Nctot
                    for p = 1:Np
                        % Multiply current guess by {1 +/- (Random_fraction<=F_RAND_MC_GUESS)}
                        pMatrixInit_MC(ctot,p) = obj.pMatrixInit(ctot,p) * ...
                            ( 1 + F_RAND_MC_GUESS*rand()*(-1)^mod(ceil(10*rand()),2) );

                        % Make sure that guess is within allowed limits
                        if( pMatrixInit_MC(ctot,p) > obj.pMatrixUpper(ctot,p) )
                            pMatrixInit_MC(ctot,p) = obj.pMatrixUpper(ctot,p);

                        elseif( pMatrixInit_MC(ctot,p) < obj.pMatrixLower(ctot,p) )
                            pMatrixInit_MC(ctot,p) = obj.pMatrixLower(ctot,p);
                        end
                    end
                end
                fitResult.setFitInitialConditions( pMatrixInit_MC, ...
                    obj.rateAnalysis.dH, obj.rateAnalysis.Eab );
                %}
                
                % Get the initial parameters from the best fit
                % (don't use function to calculate them, just copy them)
                fitResult.pMatrixInit   = obj.pMatrixBest;
                
                if( fitResult.CONSTRAIN_RATE_ANALYSIS )
                    fitResult.dH_Init       = obj.dH;
                    fitResult.Eab_Init      = obj.Eab;
                end
                
                %% Fit the group
                % Start the fit optimization
                updateGroupParamsFlag = false;
                fitResult.fitMe( updateGroupParamsFlag, R2eff_MC );
                fitResult.analyzeMe();
                
                if( OUTPUT_DEBUG )
                    for ctot = 1:Nctot
                        [cs,c]  = obj.parentGroup.getCurvesetCurve( ctot );
                        
                        % Plot each fit to the simulated data
                        vcpmg   = obj.parentGroup.curvesets{cs}.curves{c}.vcpmg;
                        TCPMG   = obj.parentGroup.curvesets{cs}.curves{c}.TCPMG;
                        params  = fitResult.pMatrixBest(ctot,1:Np);
                        
                        VCPMG_Sim = linspace(0, 1.1*max(vcpmg), 20);

                        plot(hax, VCPMG_Sim, model_MQRD_CRJ( VCPMG_Sim, TCPMG, params ), '-b');    
                    end
                end
                
                %% Store the result in the cell array
                obj.fitResults_MC{imc} = fitResult;                
            end
            
            %% Calculate deviations in fitted parameters across all MC iterations            
            [chi2_MC, resultsMatrix_MC, dH_MC, Eab_MC] = obj.getMCResults('FINAL');
            for ctot = 1:Nctot
                for r = 1:obj.Nr
                    obj.resultsMatrixErrors(ctot, r) = std( resultsMatrix_MC(ctot, r, 1:Nmc) );
                end                
            end
            obj.dH_E = std( dH_MC );
            obj.Eab_E = std( Eab_MC );
            
            % Propagate errors to rate analysis
            obj.rateAnalysis.analyzeMe();
            
            %% Wrap up
            errors_finished = true;
            
            fprintf('\nMC Errors done (%0.1f sec)\n', toc(t0));
            if( exist('h_wb','var') && ishandle(h_wb) )
                delete(h_wb);
            end            
        end
        
        function fitMe( obj, updateGroupParamsFlag, varargin )
            %% Fit the group either to NOEXCHANGE or EXCHANGE model
            % updateGroupParamsFlag designates if the group parameters
            % should be updated (re-calculate group.indepParamIndexMatrix)
            
            % Update the fit parameters for the group
            %  This should be done each time, except for multiple fits in a
            %  row (e.g., grid search) where one KNOWS it is not needed                    
            if( updateGroupParamsFlag )     
                % This will re-calculate group.indepParamIndexMatrix
                obj.parentGroup.updateFitParams(obj.CONSTRAIN_RATE_ANALYSIS);
            end
            
            % Processing optional argument for R2eff_MC
            %  This replaces the Y-axis data stored in the curve
            %  This is useful for error estimation wherein Y-axis data are
            %  generated and supplied to this function (Monte Carlo style)
            if( nargin == 3 )
                R2eff_MC = varargin{1};                
            else
                R2eff_MC = [];
            end            
            
            % Frequently used variables in the interest of code brevity
            Np      = 5;
            Nctot   = obj.parentGroup.getNumCurves();                 
            
            switch upper(obj.modelName)                
                case 'NOEXCHANGE'
                    %% Fit a single curve set to model of NO EXCHANGE
                    % Ian Kleckner
                    % Foster Lab
                    % 2010/02/15 Create code
                    % 2010/02/16 Update to include adding no-exchange model to list in fit_results
                    % 2011/04/19 Update for GUARDD and classes
                    % 
                    % FUNCTION
                    %  Fit all dispersion curves in a single curve set to model of NO EXCHANGE
                                                           
                    % Total chi2 and observations for the group
                    chi2_tot    = 0;
                    Nobs_tot    = 0;
                    
                    % Start the fitting timer
                    t0 = tic;
                    fprintf('\nStarting fit...');

                    ctot = 0;
                    % Fit each curve independently to R20 horizontal line
                    for cs = 1:obj.parentGroup.Ncs
                        curveset = obj.parentGroup.curvesets{cs};
                        for c = 1:curveset.Nc                            
                            ctot = ctot+1;
                            curve = curveset.curves{c};
                            
                            % Get the data for the curve
                            X           = curve.vcpmg';
                            if( isempty(R2eff_MC) )
                                Y       = curve.R2eff';
                            else
                                Y       = R2eff_MC{ctot}';
                            end
                            Y_E         = curve.eR2eff';        
                            Nobs_tot    = Nobs_tot + curve.Nobs;
                            
                            % Weight each observation by inverse of experimental error
                            if( sum( Y_E ) > 0 )
                                weights = 1 ./ Y_E;
                            else
                                weights = ones( 1, length(Y) );
                            end

                            s = fitoptions( 'Method','NonlinearLeastSquares',...
                                            'Weights', weights, ...
                                            'Lower', 0, ...
                                            'Upper', Inf, ...
                                            'Startpoint', 1 );
                            f = fittype('0*x + R20','options',s);

                            
                            
                            [params, void] = fit( X, Y, f );
                            
                            % Mark the time required for fitting
                            obj.timeRequired = toc(t0);

                            % Get error in R20 parameter
                            outputfit_conf  = confint(params, 0.683); % Get limits at 1 sigma
                            R201            = outputfit_conf(1);    % Lower limit at 1 sigma
                            R202            = outputfit_conf(2);    % Upper limit at 1 sigma
                            R20             = 0.5*(R201+R202);      % Fitted value
                            R20_E           = abs(R202-R20);        % Error in fit (1 sigma)

                            % Store the results
                            % dwH dwX PA kex R20
                            %  0   0   1  0   ?
                            obj.pMatrixBest(ctot,1:Np)      = [0, 0, 1, 0, R20];
                            obj.pMatrixBestErrors(ctot,1:Np)= [0, 0, 0, 0, R20_E];
                            obj.dH                          = 0;
                            obj.dH_E                        = 0;
                            obj.Eab                         = 0;
                            obj.Eab_E                       = 0;
                            
                            % Calculate chi2
                            chi2_curve = sum( ((Y-R20)./Y_E).^2 );
                            chi2_tot = chi2_tot + chi2_curve;
                        end
                    end
                    
                    % Mark the time required for fitting
                    obj.timeRequired = toc(t0);  
                    fprintf('done (%0.1f sec)\n', obj.timeRequired);
                    
                    obj.chi2    = chi2_tot;
                    % Degrees of freedom = Nobservations - Nparameters (one per curve)
                    obj.df      = Nobs_tot - Nctot;
                    obj.chi2red = obj.chi2 / obj.df;
                    
                    % Mark the name of the fit now that it is done
                    obj.name = sprintf('%s Chi2=%0.2f [%s]', ...
                        obj.fitModeString, obj.chi2, datestr(now, 'mmm dd,yyyy-HH:MM:SS'));
                    
                case 'EXCHANGE'
                    %% Fit all dispersion curves in the group
                    % Ian Kleckner
                    % Foster Lab
                    % 2010/01/11 Create code
                    % 2011/04/17 Modify for GUARDD and groups
                    % 
                    % FUNCTION
                    %  Fit all dispersion curves in a single group (multiple curve sets)
                    % 
                    % INPUT VARS
                    %  Dispersion data and errors (via obj.parentGroup)
                    %  pMatrixInit = Initial conditions used to start calculation of best fit
                    %   
                    % OUTUPT VARS
                    %  pMatrixFit   = fitted parameter values post-optimization
                    %  chi2         = chi^2 total for entire group post-optimization
                    %  df           = Total df for the fit to the group
                    %
                    % TO DO

                    % Must have a no-exchange fit first
                    %  If there is NOT a no-exchange fit yet...
                    if( isempty(obj.parentGroup.fitResult_NoEx) )
                        % ...otherwise create a new one
                        fitResult_NoEx = FitResult(obj.parentGroup, 'NOEXCHANGE', 'NoEx', obj.CONSTRAIN_RATE_ANALYSIS);
                        
                        % Set fit limits (this is not used though)
                        pMatrixLower        = zeros( obj.parentGroup.getNumCurves(), 5 );
                        pMatrixUpper        = zeros( obj.parentGroup.getNumCurves(), 5 );
                        pMatrixUpper(:,3)   = 1;
                        pMatrixUpper(:,5)   = Inf;
                        fitResult_NoEx.setLimits_pMatrix( pMatrixLower, pMatrixUpper );                        
                        
                        updateGroupParamsFlag = false;
                        fitResult_NoEx.fitMe(updateGroupParamsFlag);
                        
                        fitResult_NoEx.analyzeMe();
                        obj.parentGroup.addFitResult(fitResult_NoEx);
                    end

                    if( isempty( obj.pMatrixInit ) || ...
                        isempty( obj.pMatrixUpper ) || ...
                        isempty( obj.pMatrixLower ) )
                        error('Must define Initial, Lower and Upper parameter matrices');
                    end
                    
                    fprintf('\nStarting fit...');

                    % Get the number of parameters for fmincon
                    Np_fmincon = obj.parentGroup.getNumIndepParams();

                    % Preallocate arrays for speed
                    p_fmincon_init      = zeros(1,Np_fmincon);
                    p_fmincon_name      = cell(1,Np_fmincon);
                    lb_fmincon          = zeros(1,Np_fmincon);
                    ub_fmincon          = zeros(1,Np_fmincon);

                    % Transform from p_matrix_init -> p_fmincon_init (03/09/2009)
                    ip_fmin = 1;
                    
                    if( obj.CONSTRAIN_RATE_ANALYSIS && obj.parentGroup.getNumTemps() > 1 )                        
                        % TODO % Get dH and Ea limits from session class                                                
                        % Enthalpy dH
                        indepParam_dH = obj.parentGroup.indepParam_dH;
                        p_fmincon_init( indepParam_dH )   = obj.dH_Init;
                        lb_fmincon( indepParam_dH )       = -100e3;
                        ub_fmincon( indepParam_dH )       = 100e3;
                        p_fmincon_name{ indepParam_dH }   = sprintf('dH');
                                                
                        % Activation energy Eab
                        indepParam_Eab = obj.parentGroup.indepParam_Eab;
                        p_fmincon_init( indepParam_Eab )   = obj.Eab_Init;
                        lb_fmincon( indepParam_Eab )       = -100e3;
                        ub_fmincon( indepParam_Eab )       = 100e3;
                        p_fmincon_name{ indepParam_Eab }   = sprintf('Eab');
                        
                        % Designate the two parameters have been completed
                        ip_fmin = 3;
                    end
                    
                    for ctot = 1:Nctot
                        for p = 1:Np
                            % Only log the first instance of each independent parameter                                              
                            if( obj.parentGroup.indepParamIndexMatrix(ctot,p) >= ip_fmin )                                
                                % The parameter should be passed to fmincon
                                p_fmincon_init( ip_fmin ) = ...
                                    obj.pMatrixInit(ctot,p) / obj.parentGroup.paramScalingMatrix(ctot,p);

                                lb_fmincon( ip_fmin ) = ...
                                    obj.pMatrixLower(ctot,p) / obj.parentGroup.paramScalingMatrix(ctot,p);

                                ub_fmincon( ip_fmin ) = ...
                                    obj.pMatrixUpper(ctot,p) / obj.parentGroup.paramScalingMatrix(ctot,p);                            

                                p_fmincon_name{ ip_fmin } = ...
                                    sprintf('%s-CTot=%02d', obj.pNames{p}, ctot);

                                
                                FIX_PARAMETERS = false;
                                if( FIX_PARAMETERS )
                                    % Fixing parameters - quick and dirty (2011/06/10)
                                    % Probe     dwH|(ppm)	dwX|(ppm)	Pa|(Percent)	kex|(/s)
                                    % Leu 22δ1	 0.1±0.0	 1.0±0.2	 95.5±3.4	 1161.3±217.6
                                    % Leu 22δ2	 0.0±0.0	 1.0±0.1	 87.3±2.5	 1046.8±103.6
                                    % Ile 26δ1	 0.0±0.0	 1.4±0.2	 96.3±2.3	 2391.8±237.6      
                                    
                                    % Leu 22d1 (2011/06/10)
                                    PA  = 0.955;
                                    kex = 1161.3;
                                    
                                    % Ile 26d1 (2011/06/10)
                                    PA  = 0.963;
                                    kex = 2391.8;
                                    
                                    % Global average (2011/06/10)
                                    % Parameter	Min     Max     Avg     Dev
                                    % PA        87.29	98.6	94.82	4.36
                                    % kex       1047	2392	1629	366
                                    PA  = 94.82/100;
                                    kex = 1629;
                                    
                                    if( p == 3 )
                                        % The parameter should be passed to fmincon
                                        fprintf('\nFitResult.fitMe: Param number %d (PA) fixed to %f', ip_fmin, PA);
                                        p_fmincon_init( ip_fmin )   = PA;
                                        lb_fmincon( ip_fmin )       = PA;
                                        ub_fmincon( ip_fmin )       = PA;

                                    elseif( p == 4 )
                                        % The parameter should be passed to fmincon
                                        fprintf('\nFitResult.fitMe: Param number %d (kex) fixed to %f', ip_fmin, kex);
                                        p_fmincon_init( ip_fmin )   = kex;
                                        lb_fmincon( ip_fmin )       = kex;
                                        ub_fmincon( ip_fmin )       = kex;
                                    end
                                end
                                
                                % Count up for the next parameter
                                ip_fmin = ip_fmin+1;
                            end
                        end
                    end
                    
                    % TODO % Consider using Simplex algorithm for optimization (recommended by Evgenii Kovrigin)
                    options = optimset('fmincon');
                    options = optimset( ...
                        options, ...
                        ...'Algorithm', 'active-set',               ... Use LAST (large steps, to add speed)
                        ...'Algorithm', 'trust-region-reflective',  ... DO NOT USE (requires providing a gradient)
                        'Algorithm', 'interior-point',           ... Use FIRST
                        ...'Algorithm', 'sqp',                      ... Use SECOND
                        ...
                        ...'Display', 'final-detailed',           ...
                        'Display', 'off',           ...
                        ...'Display', 'iter',              ... 'iter' displays output at each iteration, and gives the default exit message.
                        ...'Display', 'iter-detailed',     ... 'iter-detailed' displays output at each iteration, and gives the technical exit message.
                        ...'Diagnostics', 'on',       ...
                        ...
                        ...'GradientType', 'basic', ... Default method for computing the gradients
                        ...'GradientType', 'refined',  ... Offers a more robust and less noisy gradient calculation method than 'basic'
                        ...
                        ...'MaximallyFeasible', 1,     ... Option to specify that the optimization continue after an initial, feasible solution has been found:
                        ...
                        ...'MaxFunEvals',1e5,       ... Maximum number of function evaluations allowed, a positive integer. The default value for all algorithms except interior-point is 100*numberOfVariables; for the interior-point algorithm the default is 3000.    
                        ...'MaxIter', 350,          ... Maximum number of iterations allowed, a positive integer. The default value for all algorithms except interior-point is 400; for the interior-point algorithm the default is 1000.
                        ...
                        ...'TolCon', ?,             ... Tolerance on the constraint violation, a positive scalar. The default is 1e-6.
                        ...'TolFun', 0.1,              ... Termination tolerance on the function value, a positive scalar. The default is 1e-6. (Set to 1-5% target function)
                        ...'TolX', ?,                  ... Termination tolerance on x, a positive scalar. The default value for all algorithms except interior-point is 1e-6; for the interior-point algorithm the default is 1e-10.
                        ...
                        ... ALGORITHM OPTIONS: trust-region-reflective
                        ...
                        ... ALGORITHM OPTIONS: active-set
                        ... 'MaxSQPIter', 1e4,       ... Maximum number of SQP iterations allowed, a positive integer. The default is 10*max(numberOfVariables, numberOfInequalities + numberOfBounds).
                        ... 'RelLineSrchBnd', 10,    ... Relative bound (a real nonnegative scalar value) on the line search step length such that the total displacement in x satisfies |Δx(i)| ≤ relLineSrchBnd· max(|x(i)|,|typicalx(i)|). This option provides control over the magnitude of the displacements in x for cases in which the solver takes steps that are considered too large. The default is no bounds ([]).
                        ... ALGORITHM OPTIONS: interior-point
                        ... ALGORITHM OPTIONS: sqp
                        ...
                        ... DEPRECATED options
                        ... 'LargeScale', 'off',        ... DEPRECATED
                        ...'Simplex', 'on',                ... When LargeScale is 'off', choose between the two remaining algorithms by setting the Simplex option to:
                        ...'PlotFcns', @optimplotfval, ... Useful for optimizing fitting options
                        ...
                        'TypicalX', p_fmincon_init ... Typical x values. The number of elements in TypicalX is equal to the number of elements in x0, the starting point. The default value is ones(numberofvariables,1). lsqcurvefit uses TypicalX for scaling finite differences for gradient estimation.                        
                         );
                     
                    % Start the fitting timer
                    t0 = tic;

                    % Now send the linearized parameters to fmincon
                    % to minimize the SSE and get the best fit                 
                    try 
                        [p_fmincon_best] = ..., chi2_tot, exitflag, output, lambda, grad, hessian] = ...
                        fmincon( @(p)chi2_MQRD_CRJ_group( obj, p, R2eff_MC ), ...
                            p_fmincon_init, [], [], [], [], lb_fmincon, ub_fmincon, [], options);
                    catch err
                        fprintf('\n!! Error in fmincon -- make sure error in data is not zero! (=> chi2=Inf)\n');
                        rethrow(err);
                    end

                    % Mark the time required for fitting
                    obj.timeRequired = toc(t0);  
                    fprintf('done (%0.1f sec)\n', obj.timeRequired);
                    
                    % De-linearize parameter array to matrix form            
                    obj.pMatrixBest = obj.delinearizePFmincon(p_fmincon_best);
                    
                    % Retrieve dH and Eab after the fit
                    if( obj.CONSTRAIN_RATE_ANALYSIS && obj.parentGroup.getNumTemps() > 1 )                        
                        obj.dH  = p_fmincon_best(indepParam_dH);
                        obj.Eab = p_fmincon_best(indepParam_Eab);
                    else
                        obj.dH  = NaN;
                        obj.Eab = NaN;
                    end
                    
                    % No errors to report yet
                    obj.dH_E                = NaN;
                    obj.Eab_E               = NaN;
                    obj.pMatrixBestErrors   = NaN*ones(Nctot, Np);
                    
                    % Store the other fit results
                    [obj.chi2, obj.df, chi2_curveArray] = chi2_MQRD_CRJ_group( obj, p_fmincon_best, R2eff_MC );
                    obj.chi2red = obj.chi2/obj.df;
                    
                    %------------------------------------------------------
                    % For Sandra by IK [2012/01/11]
                    %  Goal: output the chi2 for each curve
                    %  Method: quick and dirty workaround -- show result
                    %  from chi2_MQRD_CRJ_group
                    %
                    %  using the desired fit parameters and observe the result in
                    %  command window
                    OUTPUT_CHI2_PER_CURVE = true;
                    
                    if( OUTPUT_CHI2_PER_CURVE )            
                        fprintf('\n\nCurveset\tCurve\tCurveTot\tChi2');

                        % Output the chi2 for each curve
                        for ctot = 1:Nctot
                            % Get curveset and curve number
                            [cs,c] = obj.parentGroup.getCurvesetCurve( ctot );

                            % Output the chi2 for that curve
                            if( OUTPUT_CHI2_PER_CURVE )            
                                fprintf('\n%d\t%d\t%d\t%f', cs, c, ctot, chi2_curveArray(ctot))
                            end                
                        end
                        
                        fprintf('\nTotal Chi2\t\t%f', sum(chi2_curveArray));
                        fprintf('\n');
                    end
                    %------------------------------------------------------

                    % Mark the name of the fit now that it is done
                    obj.name = sprintf('%s Chi2=%0.2f [%s]', ...
                        obj.fitModeString, obj.chi2, datestr(now, 'mmm dd,yyyy-HH:MM:SS'));
                    
                    if( FIX_PARAMETERS )
                       obj.name = sprintf('FIX PA=%0.1f kex=%0.0f, %s', 100*PA, kex, obj.name); 
                    end
                otherwise
                    error('Bad modelName specified, found during fitMe()');            
            end
            % May want to call analyzeMe() now
            
        end
        
        function pMatrix = delinearizePFmincon(obj, p_fmincon)
            %% Call delineraize p_fmincon for the parent group
            % This requires the current fitResult
            pMatrix = obj.parentGroup.delinearizePFmincon(p_fmincon, obj);
        end
        
        function [dwHppm, dwHppm_E, IS_OK] = getdwHppm( obj, cs )
            %% Get the dwHppm result
                                    
            % Get the ctot numbers which correspond to the curveset "cs"
            ctot_array = obj.parentGroup.getCurveTotArray(cs);
                        
            % Find a curve out of these which contains dwH AND is in the curveset
            %  Is there a curve available that has dwH independent?
            p           = 1;
            ctot_dwH    = find( obj.parentGroup.indepParamIndexMatrix(ctot_array,p)~=0 );
            
            % If there is a MQ curve to probe dwH
            if( ~isempty(ctot_dwH) )
                ctot    = ctot_dwH(1);                
                [cs,c]  = obj.parentGroup.getCurvesetCurve(ctot);                                
                curve   = obj.parentGroup.curvesets{cs}.curves{c};
                
                % Get the value and error
                dwHppm      = obj.resultsMatrix(ctot,p) / (2*pi*curve.B0);
                dwHppm_E    = obj.resultsMatrixErrors(ctot,p) / (2*pi*curve.B0);
                IS_OK       = obj.resultsMatrixIsOK(ctot,p);
                
            else
                % dwH is unknown since there are no MQ curves to probe it
                dwHppm      = NaN;
                dwHppm_E    = NaN;                
            end
        end
        
        function [dwXppm, dwXppm_E, IS_OK] = getdwXppm( obj, cs )
            %% Get the dwXppm result
            
            % Get the ctot numbers which correspond to the curveset "cs"
            ctot_array = obj.parentGroup.getCurveTotArray(cs);
            
            p       = 2;
            ctot    = ctot_array(1);
            [cs,c]  = obj.parentGroup.getCurvesetCurve(ctot);                                
            curve   = obj.parentGroup.curvesets{cs}.curves{c};
                
            % Get the value and error
            dwXppm      = obj.resultsMatrix(ctot,p) / (2*pi*curve.B0*curve.gammaX_relative);
            dwXppm_E    = obj.resultsMatrixErrors(ctot,p) / (2*pi*curve.B0*curve.gammaX_relative);
            IS_OK       = obj.resultsMatrixIsOK(ctot,p);
        end
        
        function [Temp_unique_Array, PA_isOK_Array, kex_isOK_Array] = getKineticsIsOK(obj)
            %% Return arrays for kinetics IsOK at each temperature
            
            % Temperature for each "ctot"
            Temp_Array  = obj.parentGroup.getTemperatureArray();                        
            PA_isOK     = obj.resultsMatrixIsOK(:,3)';
            kex_isOK    = obj.resultsMatrixIsOK(:,4)';
            
            % Make each array unique
            [Temp_unique_Array, i, j]   = unique(Temp_Array);
            PA_isOK_Array               = PA_isOK(i)==1;
            kex_isOK_Array              = kex_isOK(i)==1;
        end
        
        function [Temp_unique_Array, PA_Array, PA_E_Array, kex_Array, kex_E_Array] = getKineticsValues(obj)
            %% Return arrays for kinetics at each temperature
            
            % Temperature for each "ctot"
            Temp_Array  = obj.parentGroup.getTemperatureArray();                        
            PA          = obj.resultsMatrix(:,3)';
            PA_E        = obj.resultsMatrixErrors(:,3)';
            kex         = obj.resultsMatrix(:,4)';
            kex_E       = obj.resultsMatrixErrors(:,4)';
            
            % Make each array unique
            [Temp_unique_Array, i, j]   = unique(Temp_Array);
            PA_Array                    = PA(i);
            PA_E_Array                  = PA_E(i);
            kex_Array                   = kex(i);
            kex_E_Array                 = kex_E(i);
        end
        
        function [chi2_MC, resultsMatrix_MC, dH_MC, Eab_MC] = getMCResults( obj, INITIAL_OR_FINAL )
            %% Return the results from the Monte Carlo bootstrapping
            %  INITIAL_OR_FINAL = 'INITIAL' or 'FINAL' for parameter values
            %  chi2_MC(imc) is a 1D vector length Nmc
            %  resultsMatrix_MC(ctot,p,imc) is a 3D matrix
            Nmc     = length(obj.fitResults_MC);
            Nctot   = obj.parentGroup.getNumCurves();
            Np      = obj.Np;
            Nr      = obj.Nr;
            
            chi2_MC             = NaN*ones(Nmc,1);
            resultsMatrix_MC    = NaN*ones(Nctot,Nr,Nmc);
            dH_MC               = NaN*ones(Nmc,1);
            Eab_MC              = NaN*ones(Nmc,1);
            
            switch upper(INITIAL_OR_FINAL)
                case 'INITIAL'
                    for imc = 1:Nmc
                        chi2_MC(imc)                        = obj.fitResults_MC{imc}.chi2;                
                        resultsMatrix_MC(1:Nctot,1:Np,imc)  = obj.fitResults_MC{imc}.pMatrixInit(1:Nctot,1:Np);
                        dH_MC(imc)                          = obj.fitResults_MC{imc}.dH_Init;
                        Eab_MC(imc)                         = obj.fitResults_MC{imc}.Eab_Init;
                    end
                    
                case 'FINAL'
                    for imc = 1:Nmc
                        chi2_MC(imc)                        = obj.fitResults_MC{imc}.chi2;                
                        resultsMatrix_MC(1:Nctot,1:Nr,imc)  = obj.fitResults_MC{imc}.resultsMatrix(1:Nctot,1:Nr);
                        dH_MC(imc)                          = obj.fitResults_MC{imc}.dH;
                        Eab_MC(imc)                         = obj.fitResults_MC{imc}.Eab;
                    end
                    
                otherwise
                    error('Specify either INITIAL or FINAL only');
            end
        end
        
        function isEqual = eq(fitResult1, fitResult2)
            %% (Boolean) The two curves equal in all properties (Used as "curve1==curve2")            
            isEqual = isobject(fitResult1) && isobject(fitResult2) && ...
                      strcmp(fitResult1.name, fitResult2.name) && ...
                      strcmp(fitResult1.modelName, fitResult2.modelName) && ...
                      strcmp(fitResult1.fitModeString, fitResult2.fitModeString) && ...
                      ...
                      isequal( fitResult1.pMatrixInit, fitResult2.pMatrixInit ) && ...
                      isequal( fitResult1.pMatrixLower, fitResult2.pMatrixLower ) && ...
                      isequal( fitResult1.pMatrixUpper, fitResult2.pMatrixUpper ) && ...
                      ...
                      fitResult1.timeRequired == fitResult2.timeRequired && ...                  
                      isequal( fitResult1.pMatrixBest, fitResult2.pMatrixBest ) && ...                      
                      fitResult1.chi2 == fitResult2.chi2 && ...
                      fitResult1.df == fitResult2.df;        
        end
        
        %{
        function [PhiexX, PhiexX_Error] = getPhiexX( obj, ctot )
            %% Return PhiexX using the results matrix
            % 2011/09/02
            if( isempty(obj.resultsMatrix) )
                error('Must perform fit before accessing resultsMatrix');                
            end
            
            % Pa in fraction (-)
            p       = 3;
            Pa      = obj.resultsMatrix(ctot, p);
            
            % dwX in Hz
            p       = 2;
            dwX     = Session.convertParameterToDisplayUnits(obj.resultsMatrix(ctot, p),p);
            
            % PhiexX = Pa(1-Pa)dwX^2
            PhiexX = Pa*(1-Pa)*dwX^2;
        end
        %}
        
        function resultsMatrixString = getResultsMatrixString( obj )
            %% Return results matrix in string format (DISPLAY UNITS Hz, rad/sec)
            %  For table display in Fit RD GUI (value[error])
            
            Nctot   = obj.parentGroup.getNumCurves();
            Nr      = obj.Nr;
            resultsMatrixString = cell(Nctot,Nr+1);
                        
            resultsMatrixDisplay        = zeros(Nctot,Nr);
            resultsMatrixErrorsDisplay  = zeros(Nctot,Nr);
                        
            resultsMatrixDisplay(1:Nctot, 1:Nr) = ...
                obj.parentGroup.parentSession.convertMatrixToDisplayUnits(obj.resultsMatrix(1:Nctot, 1:Nr));
            
            resultsMatrixErrorsDisplay(1:Nctot,1:Nr) = ...
                obj.parentGroup.parentSession.convertMatrixToDisplayUnits(obj.resultsMatrixErrors(1:Nctot, 1:Nr));
            
            % Use the names of the curvesets in the first column
            curveStringArray = obj.parentGroup.getCurveStringArray();
            
            for ctot = 1:Nctot
                % First column is the CurvesetCurveName                
                resultsMatrixString{ctot,1} = curveStringArray{ctot};
                
                for r = 1:Nr                    
                    resultsMatrixString{ctot,r+1} = ...
                        sprintf('  %s[%s]', displayNumber(resultsMatrixDisplay(ctot,r), '%0.1f'), ...
                                            displayNumber(resultsMatrixErrorsDisplay(ctot,r), '%0.1f'));                    
                end
            end
        end
        
        function setParamIsOK(obj, isOK, param_name, varargin)
            %% Set the param_isOK for the parameter name
            % Also may require cs or ctot number for some paramters
            % Also takes ALL or R*
            
            ctot    = NaN;
            cs      = NaN;
            switch upper(param_name)
                % Some paramters require extra input
                case {'DWH', 'DWX'}
                    if( nargin == 4 )
                        cs = varargin{1};
                    else
                        error('Must supply curveset number');
                    end
                case {'R20', 'REX'}
                    if( nargin == 4 )
                        ctot = varargin{1};
                    else
                        error('Must supply curveTot number');
                    end
                    
                case {'PA', 'KEX'}
                    if( nargin == 4 )                        
                        % Get the "ctot" numbers that correspond to this temp
                        Temp = varargin{1};
                        [Temp_array, VOID, VOID] = obj.parentGroup.getTemperatureArray();                        
                        ctot_array = Temp==Temp_array;
                    else
                        error('Must supply temperature (K)');
                    end
            end
                        
            switch upper(param_name)
                case 'ALL'
                    obj.kex0_isOK               = isOK;
                    obj.PA0_isOK                = isOK;
                    obj.dH_isOK                 = isOK;
                    obj.Eab_isOK                = isOK;
                    obj.resultsMatrixIsOK(:,:)  = isOK;                    
                    
                case 'R*'
                    obj.resultsMatrixIsOK(:,5:6)= isOK;
                    
                case 'KEX0'
                    obj.kex0_isOK               = isOK;
                    
                    % Change all instances of kex
                    r = 4;
                    obj.resultsMatrixIsOK(:,r)  = isOK;
                    
                case 'PA0'
                    obj.PA0_isOK                = isOK;
                    
                    % Change all instances of PA0
                    r = 3;
                    obj.resultsMatrixIsOK(:,r)  = isOK;
                    
                case 'PA'
                    % Change all instances of PA at the temperature
                    r = 3;
                    obj.resultsMatrixIsOK(ctot_array,r)  = isOK;
                    
                case 'KEX'
                    % Change all instances of kex at the temperature
                    r = 4;
                    obj.resultsMatrixIsOK(ctot_array,r)  = isOK;
                    
                case 'DH'
                    obj.dH_isOK = isOK;
                    
                case 'EAB'
                    obj.Eab_isOK = isOK;
                    
                case 'DWH'
                    % Change instances for that curveset
                    r = 1;
                    ctot_array = obj.parentGroup.getCurveTotArray(cs);                    
                    obj.resultsMatrixIsOK(ctot_array,r) = isOK;
                    
                case 'DWX'
                    % Change instances for that curveset
                    r = 2;
                    ctot_array = obj.parentGroup.getCurveTotArray(cs);
                    obj.resultsMatrixIsOK(ctot_array,r) = isOK;
                    
                case 'R20'
                    % Change instances for that curve
                    r = 5;
                    obj.resultsMatrixIsOK(ctot,r) = isOK;
                    
                case 'REX'
                    % Change instances for that curve
                    r = 6;
                    obj.resultsMatrixIsOK(ctot,r) = isOK;                    
                    
                case 'PHIEXX'
                    % Change instances for that curve
                    r = 8;
                    obj.resultsMatrixIsOK(ctot,r) = isOK;                    
                    
                otherwise
                    error('Invalid paramter specified');
            end
            
            % Perform rate analysis with updated paramter OK'ness
            obj.rateAnalysis.analyzeMe();
        end
        
        
        
        function setResultsMatrixIsOK( obj, resultsMatrixIsOK )
            %% Set the results matrix is OK
            
            % If the matrices are different number of columns, copy from
            % beginning columns
            Nrows = size(resultsMatrixIsOK,1);
            Ncols = size(resultsMatrixIsOK,2);
            
            obj.resultsMatrixIsOK(1:Nrows,1:Ncols) = resultsMatrixIsOK(1:Nrows,1:Ncols);
            
            % Update the rate analysis with new parameter OK'ness
            obj.rateAnalysis.analyzeMe();
        end
        
        function simMe( obj, updateGroupParamsFlag )
            %% Simulate the fit
            if( updateGroupParamsFlag )     
                % This will re-calculate group.indepParamIndexMatrix
                obj.parentGroup.updateFitParams(obj.CONSTRAIN_RATE_ANALYSIS);
            end
            
            % Must have a no-exchange fit first
            %  If there is NOT a no-exchange fit yet...
            if( isempty(obj.parentGroup.fitResult_NoEx) )
                % ...otherwise create a new one
                fitResult_NoEx = FitResult(obj.parentGroup, 'NOEXCHANGE', 'NoEx', obj.CONSTRAIN_RATE_ANALYSIS);
                updateGroupParamsFlag=false;
                fitResult_NoEx.fitMe(updateGroupParamsFlag);
                fitResult_NoEx.analyzeMe();
                obj.parentGroup.addFitResult(fitResult_NoEx);
            end
            
            % In the simulation, the best parameters are the initial ones
            obj.pMatrixBest = obj.pMatrixInit;
            obj.dH          = obj.dH_Init;
            obj.Eab         = obj.Eab_Init;
                        
            % Start the fitting timer
            t0 = tic;
            
            % Iterate through each curve in the group
            Nctot       = obj.parentGroup.getNumCurves();
            Np          = obj.parentGroup.parentSession.Np;
            chi2Array   = zeros(1,Nctot);
            NobsTotal   = 0;
            for ctot = 1:Nctot
                % Get curveset and curve number
                [cs,c] = obj.parentGroup.getCurvesetCurve( ctot );

                % Fit each curve
                vcpmg   = obj.parentGroup.curvesets{cs}.curves{c}.vcpmg;
                R2eff   = obj.parentGroup.curvesets{cs}.curves{c}.R2eff;
                eR2eff  = obj.parentGroup.curvesets{cs}.curves{c}.eR2eff;
                TCPMG   = obj.parentGroup.curvesets{cs}.curves{c}.TCPMG;

                params  = obj.pMatrixInit(ctot,1:Np);
                Nfixed  = 0;

                [chi2Array(ctot), void, void] = chi2_MQRD_CRJ( vcpmg, R2eff, eR2eff, TCPMG, params, Nfixed );
                NobsTotal               = NobsTotal + obj.parentGroup.curvesets{cs}.curves{c}.Nobs;
            end
            
            % Mark the time required for fitting
            obj.timeRequired = toc(t0);
            
            % Total chi2 for all the curves
            obj.chi2 = sum(chi2Array);

            % Degrees of freedom in global fit = Nobservations - Nparams
            obj.df = NobsTotal - obj.parentGroup.getNumIndepParams();
            obj.chi2red = obj.chi2/obj.df;

            % Mark the name of the fit now that it is done
            obj.name = sprintf('%s Chi2=%0.2f [%s]', ...
                obj.fitModeString, obj.chi2, datestr(now, 'mmm dd,yyyy-HH:MM:SS'));
        end
        
        function setInitial_Kinetics_UnconstrainedRates( obj, Temp_Array, PA_Array, kex_Array )
            %% Set initial fitting conditions        
            
            if( obj.CONSTRAIN_RATE_ANALYSIS )
                error('Cannot use CONSTRAIN_RATE_ANALYSIS with this function. Use setInitial_Kinetics_ConstrainedRates() instead.');
            end
            
            Nctot   = obj.parentGroup.getNumCurves();
            Np      = obj.parentGroup.parentSession.Np;
            
            % Set the parameter matrix
            if( isempty(obj.pMatrixInit) )
                obj.pMatrixInit = zeros(Nctot, Np);            
            end
            
            for ctot = 1:Nctot
                [cs, c] = obj.parentGroup.getCurvesetCurve(ctot);
                curve = obj.parentGroup.curvesets{cs}.curves{c};
                
                % Find the temperature for which PA and kex are specified
                ctot_temp = find(curve.Temp == Temp_Array);
                
                if( ~isempty(ctot_temp) )
                    % Set PA and kex manually via ctot_temp
                    obj.pMatrixInit(ctot,3)     = PA_Array(ctot_temp);
                    obj.pMatrixInit(ctot,4)     = kex_Array(ctot_temp);
                
                else
                    Temp_Array
                    error('Could not find temperature %d in Temp_Array (above). Please specify PA and kex at this temperature', curve.Temp);
                end             
            end
        end
        
        function setInitial_Kinetics_ConstrainedRates( obj, T0, PA0, kex0, dH, Eab )
            %% Set initial fitting conditions
            
            if( ~obj.CONSTRAIN_RATE_ANALYSIS )
                error('Must set CONSTRAIN_RATE_ANALYSIS to specify dH and Eab, or use setInitial_Kinetics_UnconstrainedRates() instead');
            end
            
            Nctot   = obj.parentGroup.getNumCurves();
            Np      = obj.parentGroup.parentSession.Np;
            
            % Set kinetic parameters
            obj.dH_Init     = dH;
            obj.Eab_Init    = Eab;
            
            % Use kinetic parameters to get PA and kex at other Temps
            simCurveset = SimulationCurveset();
            simCurveset.setSpecsKinetics( T0, PA0, kex0, Eab, dH )
            
            % Set the parameter matrix
            if( isempty(obj.pMatrixInit) )
                obj.pMatrixInit = zeros(Nctot, Np);            
            end
            
            for ctot = 1:Nctot
                [cs, c] = obj.parentGroup.getCurvesetCurve(ctot);
                curve = obj.parentGroup.curvesets{cs}.curves{c};
                
                % Set PA and kex via temperature-depenednece of rates
                obj.pMatrixInit(ctot,3)     = simCurveset.calc_PA(curve.Temp);
                obj.pMatrixInit(ctot,4)     = simCurveset.calc_kex(curve.Temp);                
            end
        end
        
        function setInitial_Shifts( obj, dwHppm_csArray, dwXppm_csArray, varargin )
            %% Set initial fitting conditions for dwH, dwX, and R20
            % Optional argument to set R20 to a fixed value
            Nctot   = obj.parentGroup.getNumCurves();
            Np      = obj.parentGroup.parentSession.Np;
                        
            % Set the parameter matrix
            if( isempty(obj.pMatrixInit) )
                obj.pMatrixInit = zeros(Nctot, Np);            
            end
            
            for ctot = 1:Nctot
                [cs, c] = obj.parentGroup.getCurvesetCurve(ctot);
                curve = obj.parentGroup.curvesets{cs}.curves{c};
                
                % Set dwH and dwX in pMatrixInit via ppm, B0, and AX
                if( curve.SQX )
                    obj.pMatrixInit(ctot,1) = 0;
                else
                    obj.pMatrixInit(ctot,1) = 2*pi*dwHppm_csArray(cs)*curve.B0;
                end                
                obj.pMatrixInit(ctot,2)     = 2*pi*dwXppm_csArray(cs)*curve.B0*curve.gammaX_relative;
            
                % R20
                if( nargin == 4 )
                    % Option to set R20 to a custom value
                    obj.pMatrixInit(ctot,5)     = varargin{1};
                else
                    % Otherwise, set R20 as minimum R2Eff
                    obj.pMatrixInit(ctot,5)     = min( curve.R2eff );                
                end
                    
            end
        end
        
        
        function setInitial_pMatrix( obj, pMatrixInit )
            %% Set initial fitting conditions via pMatrixInit
            % TODO % Finish setFitInitialConditions to check limits
            %{
            % Check that each initial condition is within the allowed fit limits
            for ctot = 1:parentGroup.getNumCurves()
                if( p_matrix_init_display(c,1) < data.B0_all(cs,c)*session.min_dwH_ppm )
                    fprintf('\ndwH too small, should be greater than %f\n', data.B0_all(cs,c)*session.min_dwH_ppm); 
                    ICs_OK = false;
                elseif( p_matrix_init_display(c,1) > data.B0_all(cs,c)*session.max_dwH_ppm )
                    fprintf('\ndwH too large, should be less than %f\n', data.B0_all(cs,c)*session.max_dwH_ppm);
                    ICs_OK = false;
                end

                if( p_matrix_init_display(c,2) < session.gamma_X / session.gamma_H * data.B0_all(cs,c)*session.min_dwX_ppm )
                    fprintf('\ndwX too small, should be greater than %f\n', session.gamma_X / session.gamma_H * data.B0_all(cs,c)*session.min_dwX_ppm); 
                    ICs_OK = false;        
                elseif( p_matrix_init_display(c,2) > session.gamma_X / session.gamma_H * data.B0_all(cs,c)*session.max_dwX_ppm )
                    fprintf('\ndwX too large, should be less than %f\n', session.gamma_X / session.gamma_H * data.B0_all(cs,c)*session.max_dwX_ppm);
                    ICs_OK = false;
                end

                if( p_matrix_init_display(c,3) < session.min_pa* session.convert_units_to_display(3) )
                    fprintf('\nPa too small, should be greater than %f\n', session.min_pa * session.convert_units_to_display(3)); 
                    ICs_OK = false;
                elseif( p_matrix_init_display(c,3) > session.max_pa* session.convert_units_to_display(3) )
                    fprintf('\nPa too large, should be less than %f\n', session.max_pa * session.convert_units_to_display(3));
                    ICs_OK = false;
                end

                if( p_matrix_init_display(c,4) < session.min_kex )
                    fprintf('\nkex too small, should be greater than %f\n', session.min_kex); 
                    ICs_OK = false;        
                elseif( p_matrix_init_display(c,4) > session.max_kex )
                    fprintf('\nkex too large, should be less than %f\n', session.max_kex);
                    ICs_OK = false;
                end    
            end
            %}            
            obj.pMatrixInit = pMatrixInit;
        end
        
        function setLimits_pMatrix( obj, pMatrixLower, pMatrixUpper )
            %% Set limits for fit
            obj.pMatrixLower = pMatrixLower;
            obj.pMatrixUpper = pMatrixUpper;
        end
        
        function setName(obj, name)
            %% Set the name of the fitResult
            obj.name = name;
        end
    end
end
classdef RateAnalysis < handle
    % RateAnalysis class: for storing results of temperature-dependent analysis
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/04/23 Start coding (based on code from 2010/01/28)
    % 2011/06/07 Use cal instead of kcal 
    %            Any RateAnalysis dH, dH_E, Eab, Eab_E, Eba, Eba_E must get multiplied by 1000
    
    properties (SetAccess = private)
        % Temperature-independent quantities (from analysis)
        % Arrhenius analysis k(T)
        arrhenius_isOK  = false;   % (Boolean) the analysis is OK
        Eab             = NaN;  % Activation energy A->B (cal/mol)
        Eab_E           = NaN;  % Error
        Pab             = NaN;  % Pre-exponential rate constant (/sec)
        Pab_E           = NaN;  % Error
        Eba             = NaN;  % Activation energy B->A (cal/mol)
        Eba_E           = NaN;  % Error
        Pba             = NaN;  % Pre-exponential rate constant (/sec)
        Pba_E           = NaN;  % Error
        
        % van't Hoff analysis K(T)
        vantHoff_isOK   = false;    % (Boolean) the analysis is OK
        dH              = NaN;  % Enthalpy (cal/mol)
        dH_E            = NaN;  % Error
        dS              = NaN;  % Entropy (cal/mol/K)
        dS_E            = NaN;  % Error
                
        R               = [];   % Gas constant (set via parentSession)
        T0              = 298;  % Reference temperature
        
        % The class instances related to this RateAnalysis
        parentFitResult = [];
        parentGroup     = [];
    %end
    
    %properties (SetAccess = private, GetAccess = private)
        % Cannot read nor write via ".", except within class functions
        
        % Data from the parent FitResult and Group
        Temps       = [];   % Tempertaures in the rate analysis
        Temps_Inv   = [];   % 1/Temps in the rate analysis        
        NTemps      = 0;    % Number of temperatures
        
        PA          = [];   % PA (-) at each temperature
        PA_E        = [];   % Error
        PA_isOK     = [];   % (Boolean) the parameter is OK
        kex         = [];   % kex (Hz) at each temperature
        kex_E       = [];   % Error
        kex_isOK    = [];   % (Boolean) the parameter is OK
        
        % Calculated from PA
        K           = [];   % Equibrium constant at each temperature
        K_E         = [];   % Error in equilibrium constant
        
        K_isOK      = [];   % (Boolean) the parameter is OK
        
        % Calculated from PA and kex
        kA          = [];   % Rate A->B (Hz) at each temperature
        kA_E        = [];   % Error in rate
        kB          = [];   % Rate B->A (Hz) at each temperature
        kB_E        = [];   % Error in rate
        
        kA_isOK     = [];   % (Boolean) the parameter is OK
        kB_isOK     = [];   % (Boolean) the parameter is OK
    end
    
    properties (Dependent = true, SetAccess = private, GetAccess = private)
        % Cannot write via ".", except within class functions
        % Value is calculated on an as-needed basis (via function get.VARIABLE(obj) below)        
    end
        
    
    methods
        function obj = RateAnalysis( parentFitResult )
            %% Default constructor
            obj.parentFitResult = parentFitResult;
            obj.parentGroup     = parentFitResult.parentGroup;
            
            obj.R               = obj.parentGroup.parentSession.R;
        end        
        
        function analyzeMe(obj)
            %% Update kinetic quantities using PA(T) and kex(T)
            %  Get K, kA, and kB, and their _isOK vectors
                        
            OUTPUT_DEBUG = obj.parentGroup.parentSession.OUTPUT_DEBUG_RATE_ANALYSIS;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end
            
            % Itemize temperatures from the group
            Temp_array      = obj.parentGroup.getTemperatureArray();
            obj.Temps       = sort(unique(Temp_array));
            obj.Temps_Inv   = 1 ./ obj.Temps;
            obj.NTemps      = length(obj.Temps);
            
            % Itemize the parameters
            resultsMatrix       = obj.parentFitResult.resultsMatrix;
            resultsMatrixErrors = obj.parentFitResult.resultsMatrixErrors;
            resultsMatrixIsOK   = obj.parentFitResult.resultsMatrixIsOK;            
            
            for t = 1:obj.NTemps
                Temp = obj.Temps(t);
                
                % Find a "total curve" index (ctot) for this temperature
                ctot = find(Temp==Temp_array);
                ctot = ctot(1);                
                
                % Obtain PA
                p               = 3;
                obj.PA(t)       = resultsMatrix(ctot,p);
                obj.PA_E(t)     = resultsMatrixErrors(ctot,p);
                obj.PA_isOK(t)  = resultsMatrixIsOK(ctot,p);
                
                % K = Pb / Pa = (1-Pa) / Pa
                obj.K(t)        = (1-obj.PA(t)) / obj.PA(t);
                
                % Error in K = K*sqrt( (E{Pa}/Pa)^2 + (E{Pa}/Pa)^2 )
                %            = K*sqrt( 2*(E{Pa}/Pa)^2 )        
                obj.K_E(t)      = obj.K(t) * sqrt( 2*( obj.PA_E(t)/obj.PA(t) )^2 );
                
                % Calculating K requires use of PA
                obj.K_isOK(t)   = obj.PA_isOK(t);
                
                % Obtain kex
                p               = 4;
                obj.kex(t)      = resultsMatrix(ctot,p);
                obj.kex_E(t)    = resultsMatrixErrors(ctot,p);
                obj.kex_isOK(t) = resultsMatrixIsOK(ctot,p);
                
                % kA = (1-Pa)*kex
                obj.kA(t)       = (1-obj.PA(t)) * obj.kex(t);

                % Error in kA = kA*sqrt( (E{Pa}/Pa)^2 + (E{kex}/kex)^2 )
                obj.kA_E(t)     = obj.kA(t) * ...
                    sqrt( ( obj.PA_E(t)/obj.PA(t) )^2 + ( obj.kex_E(t)/obj.kex(t) )^2 );
                
                % Calculating kA requires PA and kex
                obj.kA_isOK(t)  = and( obj.PA_isOK(t), obj.kex_isOK(t) );                
                
                % kB = Pa*kex
                obj.kB(t)       = obj.PA(t) * obj.kex(t);

                % Error in kB = kB*sqrt( (E{Pa}/Pa)^2 + (E{kex}/kex)^2 )
                obj.kB_E(t)     = obj.kB(t) * ...
                    sqrt( ( obj.PA_E(t)/obj.PA(t) )^2 + ( obj.kex_E(t)/obj.kex(t) )^2 );
                
                % Calculating kB requires PA and kex
                obj.kB_isOK(t)  = and( obj.PA_isOK(t), obj.kex_isOK(t) );
                
                if( OUTPUT_DEBUG )
                    fprintf('\n\nTemp=%0.1f\tPA=%0.1f(%0.1f)\tkex=%0.1f(%0.1f)', ...
                        Temp, obj.PA(t), obj.PA_E(t), obj.kex(t), obj.kex_E(t));
                    fprintf('\n\tK=%0.1f(%0.1f)\tkA=%0.1f(%0.1f)\tkB=%0.1f(%0.1f)', ...
                        obj.K(t), obj.K_E(t), obj.kA(t), obj.kA_E(t), obj.kB(t), obj.kB_E(t));                
                end
            end
            
            %% Analysis of populations: Fit van't Hoff using K
            % van't Hoff: ln(K) = dS/R + (-dH/R)*(1/T)
            
            % Get only the observations which are OK
            [X, Y, Y_E, void, void, void] = obj.getVantHoffPlot();
            %o_isOK  = find(obj.K_isOK);
            %X       = obj.Temps_Inv(o_isOK);
            %Y       = log( obj.K(o_isOK) );
            %Y_E     = obj.K_E(o_isOK) ./ obj.K(o_isOK);
                        
            % Fit the data to a line y = m*x + b
            [ m, m_E, b, b_E ] = fit_line( X, Y, Y_E );
            
            % Check if the fit is OK
            obj.vantHoff_isOK   = ~isnan(m);

            % Assign the properties
            obj.dS   = b*obj.R;                % Entropy from intercept
            obj.dS_E = abs(b_E*obj.R);         % Error propagation for scaling

            obj.dH  = -1*obj.R*m/1;         % Enthalpy from slope
            obj.dH_E= abs(obj.R*m_E/1);     % Error propagation for scaling
            
            %% Analysis of rates: Fit Arrhenius
            % Arrhenius: ln(k) = ln(P) + (-Ea/R)*(1/T)
            
            % A->B using kA
            % Get only the observations which are OK
            [X, Y, Y_E, void, void, void] = obj.getArrheniusPlotA();
            %o_isOK  = find(obj.kA_isOK);
            %X       = obj.Temps_Inv(o_isOK);
            %Y       = log( obj.kA(o_isOK) );
            %Y_E     = obj.kA_E(o_isOK) ./ obj.kA(o_isOK);
            
            % Fit the data to a line y = m*x + b
            [ m, m_E, b, b_E ] = fit_line( X, Y, Y_E );
            
            % Check if the fit is OK
            obj.arrhenius_isOK   = ~isnan(m);

            obj.Pab   = exp(b);
            obj.Pab_E = abs(exp(b) * b_E);     % Error propagation for exp

            obj.Eab   = -1*obj.R*m/1;
            obj.Eab_E = abs(obj.R*m_E/1);   % Error propagation for scaling
            
            % B->A using kB            
            % Get only the observations which are OK
            [X, Y, Y_E, void, void, void] = obj.getArrheniusPlotB();
            %o_isOK  = find(obj.kB_isOK);
            %X       = obj.Temps_Inv(o_isOK);
            %Y       = log( obj.kB(o_isOK) );
            %Y_E     = obj.kB_E(o_isOK) ./ obj.kB(o_isOK);
            
            % Fit the data to a line y = m*x + b
            [ m, m_E, b, b_E ] = fit_line( X, Y, Y_E );
            
            % Check if the fit is OK
            obj.arrhenius_isOK   = ~isnan(m);

            obj.Pba   = exp(b);
            obj.Pba_E = abs(exp(b) * b_E);     % Error propagation for exp

            obj.Eba   = -1*obj.R*m/1;
            obj.Eba_E = abs(obj.R*m_E/1);   % Error propagation for scaling            
            
        end
        
        function [X_ok, Y_ok, Y_ok_E, X_all, Y_all, Y_all_E] = getArrheniusPlotA(obj)
            %% Return X and Y vectors for the Arrhenius plots A
            % Arrnehius -> ln(kA) vs. 1/T (or ln(kB) vs 1/T)
                       
            % Get all the data (OK or otherwise)
            X_all   = obj.Temps_Inv;
            Y_all   = log( obj.kA );
            Y_all_E = obj.kA_E ./ obj.kA;
            
            % Get only the observations which are OK
            o_isOK  = find(obj.kA_isOK);
            X_ok    = X_all(o_isOK);
            Y_ok    = Y_all(o_isOK);
            Y_ok_E  = Y_all_E(o_isOK);                
        end
        
        function [X_ok, Y_ok, Y_ok_E, X_all, Y_all, Y_all_E] = getArrheniusPlotB(obj)
            %% Return X and Y vectors for the Arrhenius plots A
            % Arrnehius -> ln(kB) vs. 1/T (or ln(kB) vs 1/T)
            
            % Get all the data (OK or otherwise)
            X_all   = obj.Temps_Inv;
            Y_all   = log( obj.kB );
            Y_all_E = obj.kB_E ./ obj.kB;
            
            % Get only the observations which are OK
            o_isOK  = find(obj.kB_isOK);
            X_ok    = X_all(o_isOK);
            Y_ok    = Y_all(o_isOK);
            Y_ok_E  = Y_all_E(o_isOK);                                     
        end
        
        function rateAnalysisTableString = getRateAnalysisTableString( obj )
            %% Return rate analysis results in string format
            %  For table display in Fit RD GUI (value[error])
            %
            % dH dS Eab Pab Eba Pba
            % Convert to kcals (except for entropy cal/mol/K)
            
            property_isOK   = [obj.vantHoff_isOK,  obj.vantHoff_isOK, ...
                               obj.arrhenius_isOK, obj.arrhenius_isOK, ...
                               obj.arrhenius_isOK, obj.arrhenius_isOK ];
            propertyNames   = {'dH', 'dS', 'Eab', 'Pab', 'Eba', 'Pba'};
            propertyScale   = [1e-3, 1, 1e-3, 1, 1e-3, 1];
            %propertyScale   = [1, 1, 1, 1, 1, 1];
            rateAnalysisTableString = cell(1,6);
            
            for pr = 1:length(propertyNames)
                % Get the property value and its error
                eval(sprintf('propertyValue     = obj.%s * %d;',    propertyNames{pr}, propertyScale(pr)));                
                eval(sprintf('propertyValue_E   = obj.%s_E * %d;',  propertyNames{pr}, propertyScale(pr)));                
                
                if( property_isOK(pr) )
                    switch propertyNames{pr}
                        case {'Pab', 'Pba'}
                            % Check if there is an error value available
                            if( ~isnan(propertyValue_E) )
                                % To get base-10 scientific notation power of x (for formatting)
                                %  format 1.1e+10 +- 3.0e+9 -> (1.1+-0.3)e10
                                SN_power = 10^floor(log10( propertyValue ));
                                rateAnalysisTableString{1,pr} = sprintf(' %0.1f[%0.1f] E%0.0f', ...
                                    propertyValue/SN_power, propertyValue_E/SN_power, log10(SN_power));
                            else
                                rateAnalysisTableString{1,pr} = sprintf(' %+0.1e[?]', propertyValue);
                            end
                        otherwise
                        
                            % Check if there is an error value available
                            if( ~isnan(propertyValue_E) )
                                rateAnalysisTableString{1,pr} = sprintf(' %+0.1f[%0.1f]', propertyValue, propertyValue_E);
                            else
                                rateAnalysisTableString{1,pr} = sprintf(' %+0.1f[?]', propertyValue);
                            end                
                    end
                else
                    rateAnalysisTableString{1,pr} = ' Not OK';
                end
            end
        end
        
        function [X_ok, Y_ok, Y_ok_E, X_all, Y_all, Y_all_E] = getVantHoffPlot(obj)
            %% Return X and Y vectors for van't Hoff plot
            % van't Hoff -> ln(K) vs 1/T
            
            % Get all the data (OK or otherwise)
            X_all   = obj.Temps_Inv;
            Y_all   = log( obj.K );
            Y_all_E = obj.K_E ./ obj.K;
            
            % Get only the observations which are OK
            o_isOK  = find(obj.K_isOK);
            X_ok    = X_all(o_isOK);
            Y_ok    = Y_all(o_isOK);
            Y_ok_E  = Y_all_E(o_isOK);           
        end
    end
end
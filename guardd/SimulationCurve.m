classdef SimulationCurve < handle & SimulationSession
    % SimulationCurve class: holds a single simulation CURVE for CPMG RD two-state exchange
    %  Analagous to the Curve class
    % NOTE: Parameters cannot be set arbitrarily, PA and kex must be
    % calculated from the parentCurveset dH and Ea, and some temperature
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/04/15 Start coding
    % 2011/06/16 Fixed error calculating eIntensity from eR2eff in addData()
    
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        
        name        = 'Sim name';      % Name of simulation
        
        % Params in natural units for math functions (rad/sec, fraction)        
        params          = [];
        Temp            = 298;
        B0              = 800;
        TCPMG           = 0.02;
        SQX             = true;
        parentCurveset  = [];   % Parent curveset (for AX, dwHppm, dwXppm)
    end
    
    properties (SetAccess = private, GetAccess = private)
        % Cannot read nor write via ".", except within class functions
    end
        
    methods
        function obj = SimulationCurve( parentCurveset, Temp, B0, R20, TCPMG, SQX )
            %% Add a new curve
            % Calulate PA and kex based on temperature and parent curveset            
            obj.parentCurveset = parentCurveset;            
            obj.setSpecs( Temp, B0, R20, TCPMG, SQX )
        end
        
        function ccopy = copy( obj )
            %% Create a new curve with all properties of the current curve            
            ccopy           = SimulationCurve( obj.parentCurveset, ...
                obj.Temp, obj.B0, obj.params(5), obj.TCPMG, obj.SQX );
            ccopy.name      = sprintf('Cp(%s)',obj.name);            
        end
        
        function dwHppm = getdwHppm( obj )
            %% Calculate ppm values: rad/sec / (2*pi*MHz) = ppm
            %dwHppm = obj.params(1) / (2*pi*obj.B0);            
            dwHppm = obj.parentCurveset.dwHppm;
        end
        
        function dwXppm = getdwXppm( obj )
            %% Calculate ppm values: rad/sec / (2*pi*MHz) = ppm
            %dwXppm = obj.params(2) / (2*pi*obj.B0 * obj.AX_gamma_relative(obj.AX_index));            
            dwXppm = obj.parentCurveset.dwXppm;
        end
        
        function valueInDisplayUnits = getParamValueForDisplay( obj, p )
            %% Return parameter value in display units (Hz, %, /sec)
            valueInDisplayUnits = obj.params(p) * obj.convert_units_to_display(p);
        end
        
        function TempC = getTempC( obj )
            %% Return temperature in Celcius
            TempC = obj.Temp-273;
        end
        
        function setName( obj, varargin )
            %% Set the name of the simulation curve ([] -> auto)
            if( nargin == 1 )
                % Automatically set name (e.g., 800-MQ-25C)
                if( obj.SQX )
                    QC_String = 'SQ';
                else
                    QC_String = 'MQ';
                end
                obj.name = sprintf('%d-%s-%dC', round(obj.B0), QC_String, obj.Temp-273);
                
            elseif( nargin == 2 )
                % Manually set name from input argument
                obj.name = varargin{1}; 
            else
                error('setName() takes either ZERO or ONE input argument');
            end
        end
        
        function [vcpmgSim, R2effSim, eR2effSim, ISim, eISim] = simulateData(obj, vcpmgMin, vcpmgMax, Nobs, fNoise)
            %% Simulate dispersion curve with given noise
            % TODO % Calculate noise as would be propagated by NMR signal
            % TODO % R2eff = -ln(I/I0) / TCPMG
            
            % Calculate the vcpmg points to simulate
            vcpmgSim    = linspace(vcpmgMin, vcpmgMax, Nobs);
            
            % Simulate data with no noise
            mR2eff      = model_MQRD_CRJ(vcpmgSim, obj.TCPMG, obj.params);
            
            % Add noise to data
            noise_mean  = 0;
            noise_stdev = fNoise .* mR2eff;
            noise_boost = random('norm', noise_mean, noise_stdev, 1, Nobs);
            R2effSim    = abs(mR2eff + noise_boost);
            
            % Determine error bar
            eR2effSim   = fNoise .* mR2eff;
            
            % Determine NMR signal intensity from which R2eff would have
            %  been extracted
            % I = I0 * exp( -TCPMG * R2eff )
            % eI = TCPMG*eR2eff * exp( TCPMG * R2eff )
            %  fNoise error in signal intensities
            I0          = 1;
            ISim        = I0 * exp( -obj.TCPMG .* R2effSim );
            %eISim       = I0 * obj.TCPMG*eR2effSim .* exp(obj.TCPMG .* R2effSim);              
            eISim       = I0 * exp(-obj.TCPMG .* R2effSim) * obj.TCPMG .* eR2effSim;
        end
        
        function setParamsNatural( obj, paramsNaturalUnits )
            %% Set parameter values (rad/sec, fraction, /s, Hz)
            obj.params = paramsNaturalUnits;            
        end
        
        function setParamsFromDisplay( obj, paramsFromDisplay )
            %% Set parameter values (convert to natural units)
            % paramsFromDisplay: dwH(Hz) dwX(Hz) Pa(%) kex(/sec) R20(Hz)
            obj.params = paramsFromDisplay ./ obj.convert_units_to_display;            
        end
        
        function setParentCurveset( obj, parentCurveset )
            %% Set parent curveset for this curve
            obj.parentCurveset = parentCurveset;          
        end
        
        function setSpecs( obj, Temp, B0, R20, TCPMG, SQX )
            %% Set specs on curve (calculate kinetics from parentCurvest)           
            
            % Only set the parameters if they are specified
            if( ~isempty(Temp) )
                obj.Temp        = Temp;
            end
            if( ~isempty(R20) )
                obj.params(5)   = R20;
            end
            if( ~isempty(B0) )
                obj.B0          = B0;
            end            
            if( ~isempty(TCPMG) )
                obj.TCPMG       = TCPMG;
            end
            if( ~isempty(SQX) )
                obj.SQX         = SQX;
            end
            
            % These must go last because they use prior parameters
            % Calulate PA and kex based on temperature and parent curveset
            obj.updateParams();
        end
        
        function updateParams( obj )                                              
            %% Update kinetic parameters PA and kex from parent curveset
            if( obj.SQX )
                obj.params(1)   = 0;
            else
                obj.params(1)   = 2*pi*obj.parentCurveset.dwHppm*obj.B0;
            end
            obj.params(2)   = 2*pi*obj.parentCurveset.dwXppm*obj.B0*obj.parentCurveset.gammaX_relative;
            obj.params(3)   = obj.parentCurveset.calc_PA(obj.Temp);
            obj.params(4)   = obj.parentCurveset.calc_kex(obj.Temp);
        end        
    end
end
        
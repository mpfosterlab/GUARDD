classdef SimulationCurveset < handle & SimulationSession
    % SimulationCurveset class: holds a single CURVESET of simulations for CPMG RD two-state exchange
    %  Analagous to the Curveset class
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/04/15 Start coding
    %
    
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        
        name        = 'Set name';      % Name of set
        
        Nc          = 0;        % Number of simulations (single RD cuves)
        curves      = {};       % Holds each simulation (single RD curve)        
        c_selected  = 0;         % The selected simulation curve
        
        dwHppm              = 0;
        dwXppm              = 1;
        AX_String           = '13C';
        gammaX_relative     = 0.25143;
        
        Eab  = 2e3;
        Pab  = 2.9e3;
        Eba  = 12e3;
        Pba  = 5.8e11;
        dH   = -10e3;
        dS   = -37.9;
        T0   = 25+273;
        PA0  = 0.9;
        kex0 = 1000;
        
        % Units for display
        dwHppm_Units    = 'ppm';
        dwXppm_Units    = 'ppm';
        AX_String_Units = '-';
        
        Eab_Units       = 'cal/mol';
        Pab_Units       = '/sec';
        Eba_Units       = 'cal/mol';
        Pba_Units       = '/sec';
        dH_Units        = 'cal/mol';
        dS_Units        = 'cal/mol/K';
        T0_Units        = 'K';
        PA0_Units       = '-';
        kex0_Units      = '/sec';   
    end
    
    properties (SetAccess = private, GetAccess = private)
        % Cannot read nor write via ".", except within class functions

    end
    
    methods
        function obj = SimulationCurveset()
            %% Constructor function
            
        end
        
        function addCurve( obj, curve )
            %% Add the new curve to the end of the list
            % Since Curveset is a subclass of handle, the added curve
            % points to its place in memory, and subsequent changes to
            % either source curve or obj.curves{obj.Nc+1} are linked            
            obj.curves{obj.Nc+1} = curve;
            obj.Nc = obj.Nc+1;
        end
        
        function PA = calc_PA( obj, Temp )
            PA = 1 ./ ( 1+exp( obj.dS / obj.R - obj.dH / obj.R ./ Temp ) );
        end

        function kA = calc_kA( obj, Temp )
            kA = obj.Pab * exp( -1*obj.Eab/obj.R ./ Temp );
        end

        function kex = calc_kex( obj, Temp )
            kex = obj.calc_kA(Temp) ./ (1-obj.calc_PA(Temp));
        end

        function kB = calc_kB( obj, Temp )
            kB = obj.calc_kex(Temp) - obj.calc_kA(Temp);
        end
        
        function cscopy = copy( obj )
            %% Create a new curve with all properties of the current curve            
            cscopy           = SimulationCurveset;
            cscopy.name    = sprintf('Cp(%s)',obj.name);
            cscopy.Eab    = obj.Eab;
            cscopy.Pab    = obj.Pab;
            cscopy.Eba    = obj.Eba;
            cscopy.Pba    = obj.Pba;
            cscopy.dH    = obj.dH;
            cscopy.dS    = obj.dS;
            cscopy.T0    = obj.T0;
            cscopy.PA0    = obj.PA0;
            cscopy.kex0    = obj.kex0;
            
            % Copy each curve
            for c = 1:obj.Nc
                curve = obj.curves{c};
                ccopy = curve.copy();
                ccopy.setParentCurveset(obj);
                cscopy.addCurve( ccopy );
            end
            cscopy.Nc = obj.Nc;
        end
        
        function curve = newCurve( obj )
            %% Add a curve with default specs
            Temp        = 298;            
            R20         = 10;
            B0          = 800;            
            TCPMG       = 0.02;
            SQX         = true;
            
            obj.curves{obj.Nc+1} = SimulationCurve( obj, Temp, B0, R20, TCPMG, SQX );
            obj.curves{obj.Nc+1}.setName(sprintf('SimCurve %d', obj.Nc+1));
            curve = obj.curves{obj.Nc+1};
            
            obj.Nc = obj.Nc+1;
        end
        
        function removeCurve( obj, curve )
            %% Find the curve, then remove it
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
            %}
        end
        
        function setName( obj, name )
            %% Set the name of the simulation curveset
            obj.name = name; 
        end
        
        function setcSelected( obj, c_selected )
            %% Set the selected curveset
            if( c_selected <= obj.Nc )
                obj.c_selected = c_selected;
            else
                error('Invalid selection of curve number');
            end
        end
        
        
        %function setSpecsKinetics( obj, Eab, Pab, Eba, Pba, dH, dS, T0, PA0, kex0 )
        function setSpecsKinetics( obj, T0, PA0, kex0, Eab, dH )
            %% Set specifications for curveset
            
            % First, set the most basic properties
            obj.T0      = T0;
            obj.PA0     = PA0;
            obj.kex0    = kex0;            
            obj.Eab     = Eab;
            obj.dH      = dH;
            
            % Then calculate the remaining properties from this information
            
            % van't Hoff yields dS using PA=PA0 at temperature T=T0
            obj.dS = obj.R * log( (1-PA0)/PA0 ) + dH/T0;
            
            % Calculate kinetic parameters of interest
            % kA and kB at temperature T0
            kA0 = (1-PA0) * kex0;
            kB0 = PA0 * kex0;
            
            % Arrhenius yields Pab using kA=kA0 at temperature T0
            obj.Pab = kA0 * exp( Eab / (obj.R*T0) );

            % Arrhenius: Use kB at two temperatures to get Eba
            T1 = 280; T2 = 320;
            obj.Eba = obj.R * log( obj.calc_kB(T1) / obj.calc_kB(T2) ) / (1/T2 - 1/T1);

            % Arrhenius: use Eba and point kB at T0 to get Pba
            obj.Pba = kB0 * exp( obj.Eba / (obj.R*T0) );
                                    
            obj.updateCurves();
        end
        
        function setSpecsShifts( obj, dwHppm, dwXppm, AX_String )
            %% Set specifications for curveset
            obj.dwHppm      = dwHppm;
            obj.dwXppm      = dwXppm;
            obj.AX_String   = AX_String;
            obj.gammaX_relative = Session.getGammaRelative(AX_String);
            
            obj.updateCurves();
        end
        
        function updateCurves( obj )
            %% Update all curve parameters given the curveset parameters
            for c = 1:obj.Nc
               %fprintf('\n\tUpdating kinetics for %s', obj.curves{c}.name);
               %obj.curves{c}.updateKinetics();
               obj.curves{c}.updateParams();
            end            
        end            
    end
end
        
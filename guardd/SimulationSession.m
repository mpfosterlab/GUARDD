classdef SimulationSession < handle
    % SimulationSession class: holds a simulator session with global specs for all
    % sets and simulations
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/04/15 Start coding
    % 2017/04/25 Added 19F (IK)
    
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        
        Ncs         = 0;      % Number of simulation sets
        curvesets   = {};      % Holds each of the simulation sets
        cs_selected = 0;    % The selected simulation set
        
        % Settings for entire session
        AX_name_array               = {'1H', '13C', '15N', '19F'};
        AX_gamma_relative           = [1.000, 0.25143, 0.10136, 0.9407];
        convert_units_to_display    = [ 1/(2*3.1416) 1/(2*3.1416) 100 1 1];
        
        % TODO % Copy R gast constant from settings (parent class?)
        R = 1.9858775;  % Gas constant cal /mol/K

    end
    
    properties (SetAccess = private, GetAccess = private)
        % Cannot read nor write via ".", except within class functions
        
        

    end
    
    methods
        function addCurveset( obj, curveset )
            %% Add the new curveset to the end of the list                
            obj.curvesets{obj.Ncs+1} = curveset;
            obj.Ncs = obj.Ncs+1;                        
        end
        
        function Nc = getNumCurves( obj )
            %% Return the number of curves in the session (add from each curveset)            
            Nc = 0;
            for cs = 1:obj.Ncs
                Nc = Nc + obj.curvesets{cs}.Nc;
            end
        end
        
        function newCurveset( obj )
            %% Add a curve with default specs
            cs = obj.Ncs+1;
            obj.curvesets{cs} = SimulationCurveset;
            obj.curvesets{cs}.setName(sprintf('SimSet %d', cs));
            obj.Ncs = obj.Ncs+1;            
        end
        
        function setcsSelected( obj, cs_selected )
            %% Set the selected curveset
            if( cs_selected <= obj.Ncs )
                obj.cs_selected = cs_selected;
            else
                error('Invalid selection of curveset number');
            end
        end
        
        function removeCurveset( obj, curveset )            
            %% Find the curveset, then remove it
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
                    return                    
                end
            end            
        end
    end
end
        
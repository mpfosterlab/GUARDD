classdef DisplayCluster < handle
    % DisplayCluster class: for a group of 1+ RD curvesets, each of which includes 1+ RD curve
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/05/10 Start coding
    %
    
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        
        name            = 'Name';      % Name of group
        description     = 'Description';
        
        color_RGB       = [0 0 1];  % Display color
        display_Data    = true;     % (Bool) display this cluster
        display_Labels  = true;     % (Bool) show labels for this cluster
        
        Ng              = 0;    % Number of groups in cluster
        groups          = {};   % Link to each group in cluster
        %display_Group   = [];   % (Bool) display this group
              
        parentSession = [];     % Parent session that owns this group
    end
    
    methods
        function obj = DisplayCluster( parentSession )
            %% Default constructor requires parent session
            if( isa(parentSession, 'Session') )
                obj.parentSession = parentSession;
            else
                error('Must specify valid Session');
            end
        end
        
        function addGroup( obj, group )
            %% Add the new group to the end of the list
            if( obj.containsGroup(group) )
                %error('Group already contained in displaycluster');
                
            else
            	obj.groups{obj.Ng+1} = group;
                obj.Ng = obj.Ng+1;
            end
        end
        
        function containsGroup = containsGroup( obj, group )
           %% (Bool) the displayCluster contains this group
           
           % Check each group
           for g = 1:obj.Ng
               if( obj.groups{g} == group )
                   containsGroup = true;
                   return
               end
           end
           containsGroup = false;
        end
        
        function contains_ok_fit = containsOKFit(obj)
            %% (Bool) the displayCluster contains at least one fit for display
            contains_ok_fit = false;
            
            % Check each group...
            for g = 1:obj.Ng
                group = obj.groups{g};                
                % Check each fit...
                if( group.Nf > 0 && group.bestFitIsOK )
                    contains_ok_fit = true;
                    return
                end                
            end            
        end
        
        function [DATA, DATA_E, IS_OK, f_abshh] = getData(obj, param_name, Temp, B0, CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY)
            %% Return the parameter data for each group in this display cluster
            % CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY  'ANY'       use cs = 1
            %                                       'ALL'       use cs = 1:Ncs
            %                                       [numeric]  use cs = numeric value specified
            
            %DATA            = NaN*zeros(obj.Ng,1);
            %DATA_E          = NaN*zeros(obj.Ng,1);
            %IS_OK           = NaN*zeros(obj.Ng,1);
            DATA            = [];
            DATA_E          = [];
            IS_OK           = [];
            
            f_abshh         = NaN;
            
            % TODO % Take QC as argument when selecting data from DisplayCluster
            QC_String       = 'SQ';
            
            for g = 1:obj.Ng
                group = obj.groups{g};
                
                % Add data from group to the array
                [GROUP_DATA, GROUP_DATA_E, GROUP_IS_OK] = ...
                    group.getData(param_name, Temp, B0, QC_String, CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY);
                
                % Copy the data into the entire array
                N                     = length(GROUP_DATA);
                DATA(end+1:end+N)     = GROUP_DATA;
                DATA_E(end+1:end+N)   = GROUP_DATA_E;
                IS_OK(end+1:end+N)    = GROUP_IS_OK;
            end
        end
        
        function removeGroup( obj, group )            
            %% Find the group, then remove it            
            for g = 1:obj.Ng
                if( group == obj.groups{g} )
                    % Here, "g" is the element number to be eliminated
                    %  E.g., Remove element #3 from a set of 5 elements
                    %  Set #3=#4, #4=#5, and #5=[], then Ntot=Ntot-1
                    for g1 = g:obj.Ng-1
                        obj.groups{g1} = obj.groups{g1+1};
                    end

                    % Set final element to null
                    obj.groups{obj.Ng} = [];

                    % Update the number of elements
                    obj.Ng = obj.Ng-1;
                    return
                end
            end            
        end
        
        function setColor( obj, color_RGB )
            %% Set color of the group
            obj.color_RGB = color_RGB;
        end
        
        function setDescription(obj, description)
            %% Set parameter
            obj.description = description;
        end
        
        function setDisplayData(obj, display_Data)
            %% Set parameter
            obj.display_Data = display_Data;
        end
        
        function setDisplayLabels(obj, display_Labels)
            %% Set parameter
            obj.display_Labels = display_Labels;
        end
        
        function setName(obj, name)
            %% Set parameter
            obj.name = name;
        end
    end
end
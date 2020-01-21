function table_data = generateParamTable( class_instance )
% Return cell matrix of parameters HTML-formatted for table display
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/01/28 Start coding
% 2017/04/25 Added 19F (IK)
% 
% FUNCTION
%  Create cell matrix of parameters HTML-formatted for table display
%  This is frequently used to show properties of dataset, curve, curveset, and group
% 
% INPUT VARS
%  A class instance for which there is a property called "param_info"
%  param_info has the format: "PARAM_NAME:PARAM_DESCRIPTION"
%    Example:
%    param_info =    {...
%        'DESCRIPTION:Description of ALL data and/or project', ...
%        'seqFileName:Name of ASCII file containing three letter residues for molecular sequence (three letters per line)', ...
%        'X_nuc:Name of X nucleus (13C, 15N, or 19F)' ...
%                    };
%
%   class_instance is the instance of the class, required to access data value
%
% OUTUPT VARS
%  table_data can be directly set into the table
%
% TO USE
%  set(handles.table_dataset, 'Data', generateParamTable(class_instance) )
%
%  Also include a function in the class to save the data from a table
%  This must belong to the class, in order to write directly to its properties
%
% Working example:
%     function saveDataFromTable( obj, table_data )
%         %% Set parameters using values from table
% 
%         % Save each property, row by row
%         %  Column 1 contains property name
%         %  Column 2 contains value entered into the table
%         for row = 1:size(table_data,1)
%             propety_name = table_data{row,1};
%             eval( sprintf('obj.%s = table_data{row,2};', propety_name) );                
%         end
%         % Special treatment of certain parameters
%         % Convert Celcius to Kelvin
%         obj.Temp = obj.Temp_C + 273;
% 
%         obj.updateCurves()
%     end

param_info  = class_instance.param_info;
Nparams     = length(param_info);
table_data  = cell(Nparams,3);

% Format: "PARAM_NAME:PARAM_DESCRIPTION"
for p = 1:Nparams
    
    % Extract the name from the information string
    curr_param_info = param_info{p};

    % The ":" character separates the NAME from DESCRIPTION
    k = strfind( curr_param_info, ':');
    if( k>1 )
        curr_param_name    = curr_param_info(1:k-1);
    else
        curr_param_name    = '';
    end

    curr_param_value        = eval( sprintf('class_instance.%s', curr_param_name) );
    curr_param_description  = sprintf('  %s', curr_param_info(k+1:end) );

    table_data{p,1}     = curr_param_name;
    table_data{p,2}     = curr_param_value;
    table_data{p,3}     = sprintf('<html><i>%s</i></html>', curr_param_description);
end


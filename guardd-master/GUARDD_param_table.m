% View fitting results in formatted table (with export to ASCII) (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/??/?? Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/05/20 Use classes
% 2011/09/05 Add PhiexX
% 
% TO DO
%  Nothing?
% 

function varargout = GUARDD_param_table(varargin)
% GUARDD_PARAM_TABLE M-file for GUARDD_param_table.fig
%      GUARDD_PARAM_TABLE, by itself, creates a new GUARDD_PARAM_TABLE or raises the existing
%      singleton*.
%
%      H = GUARDD_PARAM_TABLE returns the handle to a new GUARDD_PARAM_TABLE or the handle to
%      the existing singleton*.
%
%      GUARDD_PARAM_TABLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_PARAM_TABLE.M with the given input arguments.
%
%      GUARDD_PARAM_TABLE('Property','Value',...) creates a new GUARDD_PARAM_TABLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_param_table_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_param_table_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_param_table

% Last Modified by GUIDE v2.5 04-Jun-2011 12:26:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_param_table_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_param_table_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUARDD_param_table is made visible.
function GUARDD_param_table_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_param_table (see VARARGIN)

% Choose default command line output for GUARDD_param_table
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_param_table wait for user response (see UIRESUME)
% uiwait(handles.param_table_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
% AND same for fit_r2eff window which called this GUI
setappdata(handles.param_table_gui, 'handles_main', handles_main);

refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_param_table_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Called to update display elements
function refresh_display(handles)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.param_table_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

if( session.Ng == 0 )
    return
end

%% Fill up listbox_table_params
param_name = {...
'AA|Num', ...
'Residue', ...
'Atom', ...
'Curve|Name', ...
'Curveset|Name', ...
'Group|Name', ...
'Group|Num', ...
'Curveset|Num', ...
'Curve|Num', ...
'B0|(MHz)', ...
'Temp|(C)', ...
'Num|Curvesets', ...
'Num|Curves', ...
'Num|Temps', ...
'Dataset|Name', ...
'Notes', ...
'Exchange?', ...
'Best fit|OK?', ...
'Constrain|Rates?', ...
'DoF|(NoEx)', ...
'DoF|(Ex)', ...
'Chi2|(NoEx)', ...
'Chi2|(Ex)', ...
'Chi2Red|(NoEx)', ...
'Chi2Red|(Ex)', ...
'F(Ex)', ...
'P(Ex)|(For F>1)', ...
'dwH|(ppm)', ...
'dwX|(ppm)', ...
'Pa|(Percent)', ...
'kex|(/s)', ...
'R20|(Hz)', ...
'Rex|(Hz)', ...
'alpha|(-)', ...
'PhiexX|(Hz^2)', ...
'kA|(/s)', ...
'kB|(/s)', ...
'dH|(kcal/mol)', ...
'dS|(cal/mol/K)', ...
'Ea(A->B)|(kcal/mol)', ...
'P(A->B)|(/s)', ...
'Ea(B->A)|(kcal/mol)', ...
'P(B->A)|(/s)'...
};

param_width = {...
40, ... AA num
60, ... Residue
40, ... Atom
90, ... Curve name
90, ... Curveset name
90, ... Group name
40, ... Group|Num
40, ... Curveset|Num
40, ... Curve|Num
60, ... B0
40, ... Temp
60, ... Num curvests
60, ... Num curves
60, ... Num temps
90, ... Dataset Name
100, ... Notes
70, ...
50, ... Best fit OK?
70, ... Constrain Rates?
50, ...
50, ...
60, ... Chi2(NoEx)
60, ...
60, ...
60, ... Chi2Red(Ex)
60, ... F(Ex)
100, ... P(Ex)
80, ... dwH
80, ...
80, ...
80, ... kex
80, ...
80, ...
80, ... alpha
100, ... PhiexX
80, ... kA
80, ... kB
80, ...
100, ... dS
100, ...
100, ...
100, ... Ea(B->A)
100 ...
};

set(handles.listbox_table_params, 'String', param_name);

%% Listboxes
paramDisplay    = session.paramDisplay;

% B0 listbox
uniqueB0_String = paramDisplay.getUniqueB0String();
set(handles.listbox_B0, 'String', uniqueB0_String);

% Use any value of B0?
any_B0 = get(handles.checkbox_any_B0, 'Value');

% If so, then disable the listbox
if( any_B0 )
    set(handles.listbox_B0, 'Enable', 'off');
else
    set(handles.listbox_B0, 'Enable', 'on');
end

% Quantum coherence listbox
set(handles.listbox_QC, 'String', {'SQ', 'MQ'});

% Use any value of B0?
any_QC = get(handles.checkbox_any_QC, 'Value');

% If so, then disable the listbox
if( any_QC )
    set(handles.listbox_QC, 'Enable', 'off');
else
    set(handles.listbox_QC, 'Enable', 'on');
end   

% Temperature listbox
uniqueTemp_String = paramDisplay.getUniqueTempertureString();
set(handles.listbox_temp, 'String', uniqueTemp_String);

% Use any temperature?
any_temp = get(handles.checkbox_any_temp, 'Value');

% If so, then disable the listbox
if( any_temp )
    set(handles.listbox_temp, 'Enable', 'off');
else
    set(handles.listbox_temp, 'Enable', 'on');
end

%% Set up display into table

% Get user input from listboxes
param_index_selected    = get(handles.listbox_table_params, 'Value');
Nparam_selected         = length(param_index_selected);

B0_index_selected       = get(handles.listbox_B0, 'Value');
B0_value_selected       = paramDisplay.B0_values(B0_index_selected);

QC_index_selected       = get(handles.listbox_QC, 'Value');
QC_values               = {'SQ', 'MQ'};
QC_value_selected       = QC_values(QC_index_selected);

Temp_index_selected     = get(handles.listbox_temp, 'Value');
Temp_value_selected     = paramDisplay.Temp_values(Temp_index_selected);

% F array stores the fit number for best fit (or zero if fit is a bad fit)
%Fbest(1:data.Ncs) = fit_results.list_iFit_best(1:data.Ncs);

Nctot   = session.getNumCurves();

if( Nparam_selected > 0 )
    
    % Clear the table_data_string
    table_data_string = cell(Nctot, Nparam_selected);
    
    % Clear table data for export (extra columns for fitting errors)
    table_data_export = cell(Nctot, Nparam_selected*2);
    
    % Iterate through each column in the table
    for column = 1:Nparam_selected
        
        % Select the index for this column's parameter parameter
        param_index = param_index_selected(column);
        
        % Set title of column
        table_data_export{1, 2*column-1} = param_name{param_index};
        
        % For parameters after dwH, there should be an error column too
        USE_TWO_COLUMNS = param_index >= find(ismember(param_name, 'dwH|(ppm)')==1);
        if( USE_TWO_COLUMNS )
            table_data_export{1, 2*column} = sprintf('E{ %s }', param_name{param_index});            
        end
        
        % Iterate through each group
        row = 0;
        for g = 1:session.Ng
            group = session.groups{g};
            
            % Reset display flags to ensure single usage
            B0_displayed    = false;
            QC_displayed    = false;
            Temp_displayed  = false;
            
            % Iterate through each curveset
            for cs = 1:group.Ncs
                curveset = group.curvesets{cs};
                
                % Iterate through each curve
                for c = 1:curveset.Nc
                    curve = curveset.curves{c};
                    
                    if( curve.SQX )
                        QC_String = 'SQ';
                    else
                        QC_String = 'MQ';
                    end
                    
                    % GOAL: Only use curve if B0 and QC and Temp are as user selected                    
                    % Somewhat complex if statement set to multi-line for readability
                    B0_is_valid     = ~any_B0 && sum(curve.B0 == B0_value_selected);
                    QC_is_valid     = ~any_QC && sum( strcmpi(QC_String, QC_value_selected));
                    Temp_is_valid   = ~any_temp && sum(curve.Temp == Temp_value_selected);

                    if(     (any_B0&&~B0_displayed)    && (any_QC&&~QC_displayed)  && (any_temp&&~Temp_displayed) ...
                        || ...
                            ((any_B0&&~B0_displayed)    && (any_QC&&~QC_displayed)  && Temp_is_valid) ...
                        || ...
                            ((any_B0&&~B0_displayed)    && QC_is_valid              && (any_temp&&~Temp_displayed)) ...
                        || ...
                            ((any_B0&&~B0_displayed)    && QC_is_valid              && Temp_is_valid) ...
                        || ...
                            (B0_is_valid                && (any_QC&&~QC_displayed)  && (any_temp&&~Temp_displayed)) ...
                        || ...
                            (B0_is_valid                && (any_QC&&~QC_displayed)  && Temp_is_valid) ...
                        || ...
                            (B0_is_valid                && QC_is_valid              && (any_temp&&~Temp_displayed)) ...
                        || ...
                            (B0_is_valid                && QC_is_valid              && Temp_is_valid) ...
                        )
                    
                        % Update row number
                        row = row+1;

                        % Limit any B0 or any Temp to a single occurrance
                        if( any_B0 )
                            B0_displayed = true;
                        end
                        if( any_QC )
                            QC_displayed = true;
                        end
                        if( any_temp )
                            Temp_displayed = true;
                        end

                        % Display the appropriate paramter
                        switch param_name{param_index}
                            case 'AA|Num'
                                table_data_string{row, column} = sprintf('  %d', curve.index);
                                table_data_export{row+1, 2*column-1} = curve.index;

                            case 'Residue'
                                table_data_string{row, column} = sprintf('  %s', curve.residue);
                                table_data_export{row+1, 2*column-1} = curve.residue;

                            case 'Atom'
                                table_data_string{row, column} = sprintf('  %s', curve.atom);
                                table_data_export{row+1, 2*column-1} = curve.atom;

                            case 'Curve|Name'
                                table_data_string{row, column} = curve.name;
                                table_data_export{row+1, 2*column-1} = curve.name;
                                
                            case 'Curveset|Name'
                                table_data_string{row, column} = curveset.name;
                                table_data_export{row+1, 2*column-1} = curveset.name;

                            case 'Group|Name'
                                table_data_string{row, column} = group.name;
                                table_data_export{row+1, 2*column-1} = group.name;

                            case 'B0|(MHz)'
                                table_data_string{row, column} = sprintf('  %0.2f', curve.B0);
                                table_data_export{row+1, 2*column-1} = curve.B0;

                            case 'Temp|(C)'
                                table_data_string{row, column} = sprintf('  %0.1f', curve.Temp_C);
                                table_data_export{row+1, 2*column-1} = curve.Temp_C;
                            
                            case 'Num|Curvesets'
                                table_data_string{row, column} = sprintf('       %d', group.Ncs);
                                table_data_export{row+1, 2*column-1} = group.Ncs;                                

                            case 'Num|Curves'
                                table_data_string{row, column} = sprintf('       %d', group.getNumCurves());
                                table_data_export{row+1, 2*column-1} = group.getNumCurves();

                            case 'Num|Temps'
                                table_data_string{row, column} = sprintf('       %d', group.getNumTemps());
                                table_data_export{row+1, 2*column-1} = group.getNumTemps();

                            case 'Dataset|Name'
                                table_data_string{row, column} = curve.dataset.name;
                                table_data_export{row+1, 2*column-1} = curve.dataset.name;
                                
                            case 'Notes'
                                % TODO % Replace with getTextOneLiner();
                                text_oneliner = '';
                                Nlines = max(size(group.note.text));
                                for l = 1:Nlines
                                    if( l == 1 )
                                        text_oneliner = sprintf('%s', group.note.text{l});
                                    else
                                        text_oneliner = sprintf('%s [CR] %s', text_oneliner, group.note.text{l});
                                    end
                                end                                
                                %group.note.getTextOneLiner();
                                table_data_string{row, column} = text_oneliner;
                                table_data_export{row+1, 2*column-1} = text_oneliner;
                                
                            case 'Group|Num'
                                table_data_string{row, column} = sprintf('    %d', g);
                                table_data_export{row+1, 2*column-1} = g;
                                
                            case 'Curveset|Num'
                                table_data_string{row, column} = sprintf('    %d', cs);
                                table_data_export{row+1, 2*column-1} = cs;

                            case 'Curve|Num'
                                table_data_string{row, column} = sprintf('    %d', c);
                                table_data_export{row+1, 2*column-1} = c;

                            case 'Exchange?'
                                table_data_string(row, column) = num2cell(group.exhibitsExchange==1);
                                table_data_export{row+1, 2*column-1} = group.exhibitsExchange==1;

                            case 'Best fit|OK?'
                                table_data_string(row, column) = num2cell(group.bestFitIsOK==1);
                                table_data_export{row+1, 2*column-1} = group.bestFitIsOK==1;
                                
                            otherwise                                
                                % Make sure there is a fit AND at least one of the parameters is OK                            
                                if( group.Nf == 0 )
                                    % No fit available
                                    table_data_string{row, column} = sprintf('Not Fit');
                                    table_data_export{row+1, 2*column-1} = sprintf('Not Fit');
                                    if( USE_TWO_COLUMNS )
                                        table_data_export{row+1, 2*column} = sprintf('Not Fit');
                                    end
                                    
                                elseif( ~group.bestFitIsOK )
                                    % The best fit is not OK
                                    table_data_string{row, column} = sprintf('Fit Not OK');
                                    table_data_export{row+1, 2*column-1} = sprintf('Fit Not OK');
                                    if( USE_TWO_COLUMNS )
                                        table_data_export{row+1, 2*column} = sprintf('Fit Not OK');
                                    end
                                    
                                else
                                    switch param_name{param_index}
                                       
                                        case 'Constrain|Rates?'
                                            table_data_string(row, column) = num2cell(group.fitResult_Best.CONSTRAIN_RATE_ANALYSIS==1);
                                            table_data_export{row+1, 2*column-1} = group.fitResult_Best.CONSTRAIN_RATE_ANALYSIS==1;
                                        
                                       case 'DoF|(NoEx)'
                                            table_data_string{row, column} = sprintf('    %d', group.fitResult_NoEx.df );
                                            table_data_export{row+1, 2*column-1} = group.fitResult_NoEx.df;
                                        
                                        % Special treatment if exchange is favored
                                        case {'DoF|(Ex)', 'Chi2|(Ex)', 'Chi2Red|(Ex)', 'F(Ex)', 'P(Ex)|(For F>1)'}
                                            if( group.exhibitsExchange )
                                                switch param_name{param_index}
                                                    case 'DoF|(Ex)'
                                                        table_data_string{row, column} = sprintf('    %d', group.fitResult_Best.df);
                                                        table_data_export{row+1, 2*column-1} = group.fitResult_Best.df;

                                                    case 'Chi2|(Ex)'
                                                        table_data_string{row, column} = sprintf('  %0.2f', group.fitResult_Best.chi2 );
                                                        table_data_export{row+1, 2*column-1} = group.fitResult_Best.chi2;

                                                    case 'Chi2Red|(Ex)'
                                                        table_data_string{row, column} = sprintf('  %0.2f', group.fitResult_Best.chi2red );
                                                        table_data_export{row+1, 2*column-1} = group.fitResult_Best.chi2red;

                                                    case 'F(Ex)'
                                                        table_data_string{row, column} = sprintf('  %0.2f', group.fitResult_Best.Fstatistic );
                                                        table_data_export{row+1, 2*column-1} = group.fitResult_Best.Fstatistic;

                                                    case 'P(Ex)|(For F>1)'
                                                        table_data_string{row, column} = sprintf(' %0.2e', group.fitResult_Best.Pvalue );
                                                        table_data_export{row+1, 2*column-1} = group.fitResult_Best.Pvalue;
                                                end

                                            % No exchange
                                            else
                                                table_data_string{row, column} = '    No Exch';
                                                table_data_export{row+1, 2*column-1} = 'No Exch';
                                                if( USE_TWO_COLUMNS )
                                                    table_data_export{row+1, 2*column} = 'No Exch';
                                                end
                                            end

                                        case 'Chi2|(NoEx)'
                                            table_data_string{row, column} = sprintf('    %0.2f', group.fitResult_NoEx.chi2);
                                            table_data_export{row+1, 2*column-1} = group.fitResult_NoEx.chi2;

                                        case 'Chi2Red|(NoEx)'
                                            table_data_string{row, column} = sprintf('    %0.2f', group.fitResult_NoEx.chi2red);
                                            table_data_export{row+1, 2*column-1} = group.fitResult_NoEx.chi2red;
                                            
                                        % Special case for these parameters
                                        % which can be obtained via
                                        % group.getData() function
                                        case {'dwH|(ppm)', 'dwX|(ppm)', 'Pa|(Percent)', 'kex|(/s)', 'R20|(Hz)', ...
                                              'Rex|(Hz)', 'alpha|(-)', 'PhiexX|(Hz^2)', 'kA|(/s)', 'kB|(/s)', 'dH|(kcal/mol)', 'dS|(cal/mol/K)', ...
                                              'Ea(A->B)|(kcal/mol)', 'P(A->B)|(/s)', 'Ea(B->A)|(kcal/mol)', 'P(B->A)|(/s)' }
                                          
                                          switch param_name{param_index}                                          
                                            case 'dwH|(ppm)'
                                                PARAM_NAME = 'DWH_PPM';                                                

                                            case 'dwX|(ppm)'
                                                PARAM_NAME = 'DWX_PPM';

                                            case 'Pa|(Percent)'
                                                PARAM_NAME = 'PA';

                                            case 'kex|(/s)'
                                                PARAM_NAME = 'KEX';

                                            case 'R20|(Hz)'
                                                PARAM_NAME = 'R20';

                                            case 'Rex|(Hz)'
                                                PARAM_NAME = 'REX';

                                            case 'alpha|(-)'
                                                PARAM_NAME = 'ALPHA';
                                                
                                            case 'PhiexX|(Hz^2)'
                                                PARAM_NAME = 'PHIEXX';

                                            case 'kA|(/s)'
                                                PARAM_NAME = 'KA';

                                            case 'kB|(/s)'
                                                PARAM_NAME = 'KB';

                                            case 'dH|(kcal/mol)'
                                                PARAM_NAME = 'DH';                                                

                                            case 'dS|(cal/mol/K)'
                                                PARAM_NAME = 'DS';

                                            case 'Ea(A->B)|(kcal/mol)'
                                                PARAM_NAME = 'EAB';

                                            case 'P(A->B)|(/s)'
                                                PARAM_NAME = 'PAB';

                                            case 'Ea(B->A)|(kcal/mol)'
                                                PARAM_NAME = 'EBA';

                                            case 'P(B->A)|(/s)'
                                                PARAM_NAME = 'PBA';
                                          end
                                          
                                          % Obtain data via PARAM_NAME
                                          [DATUM, DATUM_E, IS_OK] = group.getData(PARAM_NAME, curve.Temp, curve.B0, QC_String, cs);
                                          DATUM   = Session.convertUnits(DATUM, PARAM_NAME, 'DISPLAY');
                                          DATUM_E = Session.convertUnits(DATUM_E, PARAM_NAME, 'DISPLAY');
                                          
                                          %{
                                          % Hz -> ppm
                                          if( strcmpi(PARAM_NAME,'DWH') )
                                              DATUM     = DATUM / (curve.B0);
                                              DATUM_E   = DATUM_E / (curve.B0);
                                              
                                          elseif( strcmpi(PARAM_NAME,'DWX') )
                                              DATUM     = DATUM / (curve.B0*curve.gammaX_relative);
                                              DATUM_E   = DATUM_E / (curve.B0*curve.gammaX_relative);                                              
                                          end
                                          %}
                                          
                                          % Commit these data to the table
                                          if( IS_OK )
                                              table_data_string{row, column}        = sprintf(' %s±%s', ...
                                                  displayNumber(DATUM, '%0.1f'), displayNumber(DATUM_E, '%0.1f') );
                                              table_data_export{row+1, 2*column-1}  = displayNumber(DATUM, '%0.5f');
                                              
                                              if( USE_TWO_COLUMNS )
                                                table_data_export{row+1, 2*column}    = displayNumber(DATUM_E, '%0.5f');
                                              end
                                              
                                          else
                                              table_data_string{row, column}        = sprintf('Param Not OK');
                                              table_data_export{row+1, 2*column-1}  = sprintf('Param Not OK');
                                              if( USE_TWO_COLUMNS )
                                                table_data_export{row+1, 2*column}    = sprintf('Param Not OK');
                                              end
                                          end
                                    end
                                end
                        end
                    end
                end
            end
        end
    end
       
    % Trim off any excess rows
    table_data_string = table_data_string(1:row, 1:column);
    table_data_export = table_data_export(1:row+1, 1:column*2);
       
    %fprintf('\n\nIn body\tRows = %d, Cols = %d', size(table_data_export,1), size(table_data_export,2) )
    
    %% Populate table with data
    param_name_selected = param_name(param_index_selected);

    set(handles.table_results, 'ColumnName', param_name_selected);
    set(handles.table_results, 'ColumnWidth', param_width(param_index_selected));
    set(handles.table_results, 'Data', table_data_string);
    
    %% Save the formatted data in case user wants to export to data file
    
    % This contains the formatted values with "±"
    % Must add titles to each column though
    table_data_string_head                  = cell(row+1, column);
    table_data_string_head(1,1:column)      = param_name_selected(1:column);
    table_data_string_head(2:row+1,1:column)= table_data_string(1:row,1:column);
    
    setappdata(handles.param_table_gui, 'table_data_string_head', table_data_string_head);
    setappdata(handles.param_table_gui, 'table_data_export', table_data_export);
    
% No data to show in the table
else    
    set(handles.table_results, 'ColumnName', []);
    set(handles.table_results, 'Data', []);
end


% --- Executes on button press in button_export_ASCII.
function button_export_ASCII_Callback(hObject, eventdata, handles)
% hObject    handle to button_export_ASCII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save results of table to ASCII file

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.param_table_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');
VERSION         = getappdata(handles_main.main_gui, 'VERSION');

% Ask user for file to save as
default_filename    = sprintf('GUARDD-Paramter_Table-%s.csv', datestr(now, 'yyyy.mm.dd-HH-MM'));

if( ~exist( session.outputDir, 'dir' ) )
    mkdir(session.outputDir);
end

[filename, filepath] = uiputfile('*.csv', 'Save table results to file', ...
    sprintf('./%s/%s', session.outputDir, default_filename) );

% Save it
if( ~isequal(filename,0) )
    
    % Read data from memory
    table_data_string_head  = getappdata(handles.param_table_gui, 'table_data_string_head');
    table_data_export       = getappdata(handles.param_table_gui, 'table_data_export');
        
    % Strip the end off the filename
    filename_base = strrep(filename, '.csv', '');
    
    % File 1: Formatted
    fprintf('\nWriting %s-formatted.csv...', filename_base);
    
    FILE = fopen(sprintf('%s/%s-formatted.csv', filepath, filename_base), 'w');
    fprintf(FILE, 'GUARDD %s', VERSION);
    fprintf(FILE, '\n(C) Ian Kleckner 2010-2011');
    fprintf(FILE, '\nExported results table');
    fprintf(FILE, '\nCreated on %s', datestr(now, 'yyyy.mm.dd-HH-MM'));
    fprintf(FILE, '\n\n');
    fclose(FILE);
    
    % Write table
    cell2csv_ik( sprintf('%s/%s-formatted.csv', filepath, filename_base), table_data_string_head, ',');
    fprintf('Done!\n');
    
    
    % File 2: Unformatted        
    fprintf('Writing %s.csv...', filename_base);
    
    FILE = fopen(sprintf('%s/%s.csv', filepath, filename_base), 'w');
    fprintf(FILE, 'GUARDD %s', VERSION);
    fprintf(FILE, '\n(C) Ian Kleckner 2010-2011');
    fprintf(FILE, '\nExported results table');
    fprintf(FILE, '\nCreated on %s', datestr(now, 'yyyy.mm.dd-HH-MM'));
    fprintf(FILE, '\n\n');
    fclose(FILE);
    
    % Write table
    cell2csv_ik( sprintf('%s/%s.csv', filepath, filename_base), table_data_export, ',');        
    fprintf('Done!\n');
    
    %fprintf('\n\nIn button\tRows = %d, Cols = %d', size(table_data_export,1), size(table_data_export,2) )
end

% --- Executes on selection change in listbox_table_params.
function listbox_table_params_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_table_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_table_params contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_table_params

%% Update display when value in listbox is selected
refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function listbox_table_params_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_table_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when param_table_gui is resized.
function param_table_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to param_table_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Resize each element of the GUI

% Get size of the GUI window which is just resized
fP = get(handles.param_table_gui, 'Position');


% Main table (scale size, fix relative position)
left    = 0.5;
width   = fP(3)-2*left;
height  = 8.5;
bottom  = fP(4)-height-0.1;
set(handles.panel_manage, 'Position', [left, bottom, width, height]);

% Main table (scale size, fix relative position)
left    = 0.5;
bottom  = 0.1;
width   = fP(3)-2*left;
height  = fP(4)-8.5-0.1-0.1;
set(handles.table_results, 'Position', [left, bottom, width, height]);

%{

%% Listbox
% Adjust listbox (scale position, fix size)
width   = 25;
height  = 8.5;
left    = 0.5;
bottom  = fP(4)-height-0.5;
set(handles.listbox_table_params, 'Position', [left, bottom, width, height]);

% Adjust listbox (scale position, fix size)
width   = 20;
height  = 5;
left    = 25+1;
bottom  = fP(4)-height-0.5;
set(handles.listbox_B0, 'Position', [left, bottom, width, height]);

% Adjust listbox (scale position, fix size)
width   = 20;
height  = 5;
left    = 25+1+20+1;
bottom  = fP(4)-height-0.5;
set(handles.listbox_temp, 'Position', [left, bottom, width, height]);

%% Checkbox
% Adjust button_export_ASCII (scale position, fix size)
width   = 20;
height  = 3;
left    = 25+1;
bottom  = fP(4)-5-3.5;
set(handles.checkbox_any_B0, 'Position', [left, bottom, width, height]);

% Adjust button_export_ASCII (scale position, fix size)
width   = 20;
height  = 3;
left    = 25+1+20+1;
bottom  = fP(4)-5-3.5;
set(handles.checkbox_any_temp, 'Position', [left, bottom, width, height]);

%% Button
% Adjust button_export_ASCII (scale position, fix size)
width   = 20;
height  = 2.5;
left    = 25+1+20+1+20+1;
bottom  = fP(4)-height-0.5;
set(handles.button_export_ASCII, 'Position', [left, bottom, width, height]);
%}


% --- Executes on selection change in popup_B0.
function popup_B0_Callback(hObject, eventdata, handles)
% hObject    handle to popup_B0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_B0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_B0


% --- Executes during object creation, after setting all properties.
function popup_B0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_B0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_B0.
function listbox_B0_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_B0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_B0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_B0
%% Update display when value in listbox is selected
refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function listbox_B0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_B0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_temp.
function listbox_temp_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_temp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_temp
%% Update display when value in listbox is selected
refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function listbox_temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_any_B0.
function checkbox_any_B0_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_any_B0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_any_B0

%% Disable/enable any magnetic field
refresh_display(handles);


% --- Executes on button press in checkbox_any_temp.
function checkbox_any_temp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_any_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_any_temp

%% Disable/enable any temperature
refresh_display(handles);


% --- Executes on selection change in listbox_QC.
function listbox_QC_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_QC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_QC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_QC
refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function listbox_QC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_QC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_any_QC.
function checkbox_any_QC_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_any_QC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_any_QC
refresh_display(handles);

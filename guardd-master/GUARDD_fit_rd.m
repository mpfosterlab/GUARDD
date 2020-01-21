% User-directed fit of RD curve set, set parameter validity (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/01/01
% 2011/01/11 Convert to GUARDD program
% 2011/04/18 Implement classes for data structures
% 
% TO DO
%  User defines error
%  Set custom limits on fitted parameter values pre-optimization

function varargout = GUARDD_fit_rd(varargin)
% GUARDD_FIT_RD M-file for GUARDD_fit_rd.fig
%      GUARDD_FIT_RD, by itself, creates a new GUARDD_FIT_RD or raises the existing
%      singleton*.
%
%      H = GUARDD_FIT_RD returns the handle to a new GUARDD_FIT_RD or the handle to
%      the existing singleton*.
%
%      GUARDD_FIT_RD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_FIT_RD.M with the given input arguments.
%
%      GUARDD_FIT_RD('Property','Value',...) creates a new GUARDD_FIT_RD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_fit_rd_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_fit_rd_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_fit_rd

% Last Modified by GUIDE v2.5 10-Jun-2011 11:51:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_fit_rd_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_fit_rd_OutputFcn, ...
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


% --- Executes just before GUARDD_fit_rd is made visible.
function GUARDD_fit_rd_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_fit_rd (see VARARGIN)

% Choose default command line output for GUARDD_fit_rd
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_fit_rd wait for user response (see UIRESUME)
% uiwait(handles.fit_rd_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.fit_rd_gui, 'handles_main', handles_main);

%% Read data from main GUI
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};
Nctot   = group.getNumCurves();
Np      = session.Np;

%% Initialize progress bar axis
% Update total progress wait bar (Johannes Korsawe, 2007)
value = 0;
h=handles.axes_current_progress;
cla(h);
patch([0,value,value,0],[0,0,1,1],'r', 'Parent', h);
axis(h,[0,1,0,1]); axis(h,'off'); drawnow;

%% Table: grid search
column_string = {'dwH|(ppm)', 'dwX|(ppm)', 'PA0|(%)', 'kex0|(/s)', 'dH|(kcal/mol)', 'Eab|(kcal/mol)'};
set(handles.table_grid_limits, 'ColumnName', column_string );
set(handles.table_grid_limits, 'ColumnWidth', {70});

% List each curve in each curve set in this group
row_string = {'Min', 'Max', 'Steps'};
set(handles.table_grid_limits, 'RowName', row_string);
table_data = zeros(length(row_string), length(column_string));
set(handles.table_grid_limits, 'Data', table_data);

% Temperature
column_string = [];
set(handles.table_Temp0_Grid, 'ColumnName', column_string );
set(handles.table_Temp0_Grid, 'ColumnWidth', {70});
row_string = 'Temp0(C)';
set(handles.table_Temp0_Grid, 'RowName', row_string);
set(handles.table_Temp0_Grid, 'Data', group.T0_Grid-273);

%% Table: Initial fit conditions
%{
% Set up column titles
column_string = {''};
for p = 1:session.Np
    column_string{p} = sprintf('%s|(%s)', session.param_name{p}, session.param_units_display{p});
end
set(handles.table_ics, 'ColumnName', column_string );
set(handles.table_ics, 'ColumnWidth', {60,60,50,60,50});

% List each curve in each curve set in this group
rownames = group.getCurveStringArray();        
set(handles.table_ics, 'RowName', rownames);

set(handles.table_ics, 'Data', zeros(Nctot, Np));
%}

%% Table: Initial temperature-depdenent conditions
% Set up column titles
column_string = {'T0|(C)', 'PA0|(%)', 'kex0|(/s)', 'dH|(kcal/mol)', 'Eab|(kcal/mol)'};
set(handles.table_ics_group_kinetics, 'ColumnName', column_string );
set(handles.table_ics_group_kinetics, 'ColumnWidth', {50,50,60,70,70});

% List each curve in each curve set in this group
set(handles.table_ics_group_kinetics, 'RowName', []);
table_data = [25, 90, 1000, -10, 2];
set(handles.table_ics_group_kinetics, 'Data', table_data);

%% Table: Initial conditions kinetics at each temperature
column_string = {'Temp(C)', 'PA(%)','kex(/s)'};
set(handles.table_ics_temp_kinetics, 'ColumnName', column_string );
set(handles.table_ics_temp_kinetics, 'ColumnWidth', {80,80,80});
set(handles.table_ics_temp_kinetics, 'RowName', []);

Temp_array  = unique( group.getTemperatureArray() );
Ntemp       = length(Temp_array);
table_data  = zeros(Ntemp, 3);

% Temp kex(/s) PA(%)
for t = 1:Ntemp
    table_data(t,1) = Temp_array(t)-273;
    table_data(t,2) = 90;
    table_data(t,3) = 1500;
end
set(handles.table_ics_temp_kinetics, 'Data', table_data);


%% Table: Initial conditions for group ppm shifts
column_string = {'Curveset', 'dwH(ppm)', 'dwX(ppm)'};
set(handles.table_ics_curveset_ppm, 'ColumnName', column_string );
set(handles.table_ics_curveset_ppm, 'ColumnWidth', {120,70,70});

% List each curve in each curve set in this group
%row_string = group.getCurvesetStringArray();
set(handles.table_ics_curveset_ppm, 'RowName', []);

curveset_names = group.getCurvesetStringArray();
table_data       = cell(group.Ncs,3);
for cs = 1:group.Ncs
    table_data{cs,1} = curveset_names{cs};       
    table_data{cs,2}  = 0.01;    
    table_data{cs,3}  = 1;
end
set(handles.table_ics_curveset_ppm, 'Data', table_data );

   
%% Table: View initial conditions (no editing)
column_string = {'Curveset[Curve]', ...
                'Temp(C)', 'B0(MHz)', 'SQX', ...
                'dwH(Hz)', 'dwX(Hz)', ...
                 'PA(%)', 'kex(/s)', 'R20(Hz)'};
set(handles.table_ics_view_all, 'ColumnName', column_string );
set(handles.table_ics_view_all, 'ColumnWidth', {200,   60,60,40,   60,60,   50,60,60});
set(handles.table_ics_view_all, 'ColumnFormat', {'char' 'numeric', 'numeric', 'logical', ...
    'numeric', 'numeric', 'numeric', 'numeric'});

% List each curve in each curve set in this group
row_string = [];
set(handles.table_ics_view_all, 'RowName', row_string);
table_data = [];
%zeros(Nctot, length(column_string))
set(handles.table_ics_view_all, 'Data', table_data);

% Set the vales in the table
update_table_ics(handles);


%% Table: rate analysi
%{
% Set up column titles
row_string = {'T0'; 'PA0'; 'kex0'; 'dH'; 'Eab'};
set(handles.table_constrain_rates, 'RowName', row_string );

set(handles.table_constrain_rates, 'ColumnName', [] );
set(handles.table_constrain_rates, 'ColumnWidth', {60,60});

table_data = { 25, 'C'; ...
               90, '%'; ...
               1000, '/s'; ...
               -10, 'kcal/m'; ...
               2, 'kcal/m' };
set(handles.table_constrain_rates, 'Data', table_data);
set(handles.table_constrain_rates, 'Enable', 'on');
%}

%% Table: Group Fitting results (dwH, dwX, dH, Eab)
column_string = {'T0|(C)', 'PA0|(%)', 'kex0|(/s)', 'dH|(kcal/mol)', 'Eab|(kcal/mol)'};
set(handles.table_results_group_kinetics, 'ColumnName', column_string );
set(handles.table_results_group_kinetics, 'ColumnWidth', {50,70,80,70,70});
set(handles.table_results_group_kinetics, 'Data', []);
set(handles.table_results_group_kinetics, 'RowName', []);

%% Table: Fit results for kinetics at each temperature
column_string = {'Temp(C)', 'PA(%)','kex(/s)'};
set(handles.table_results_temp_kinetics, 'ColumnName', column_string );
set(handles.table_results_temp_kinetics, 'ColumnWidth', {80,80,100});
set(handles.table_results_temp_kinetics, 'RowName', []);

table_data = [];
set(handles.table_results_temp_kinetics, 'Data', table_data);


%% Table: Results for curveset ppm
column_string = {'Curveset', 'dwH(ppm)', 'dwX(ppm)'};
set(handles.table_results_curveset_ppm, 'ColumnName', column_string );
set(handles.table_results_curveset_ppm, 'ColumnWidth', {120,70,70});

% List each curve in each curve set in this group
%row_string = group.getCurvesetStringArray();
set(handles.table_results_curveset_ppm, 'RowName', []);
%set(handles.table_results_curveset_ppm, 'Data', zeros(length(row_string), length(column_string)));

curveset_names = group.getCurvesetStringArray();
table_data       = cell(group.Ncs,3);
for cs = 1:group.Ncs
    table_data{cs,1} = curveset_names{cs};       
    table_data{cs,2}  = 0;    
    table_data{cs,3}  = 0;
end
set(handles.table_results_curveset_ppm, 'Data', table_data );


%% Fit results table
%{
column_string = {''};
for p = 1:session.Np
    column_string{p} = sprintf('%s (%s)', session.param_name{p}, session.param_units_display{p});
end
column_string{6} = sprintf('Rex (Hz)');
column_string{7} = sprintf('alpha (-)');
%}
column_string = {'Curveset[Curve]', 'dwH(Hz)', 'dwX(Hz)', 'Pa(%)', 'kex(/s)', 'R20(Hz)', 'Rex(Hz)', 'alpha(-)', 'PhiExX(Hz^2)'};
set(handles.table_fit_results, 'ColumnName', column_string );
set(handles.table_fit_results, 'ColumnWidth', {200,70,70,70,70,60,60,60});

% List each curve in each curve set in this group
%rownames = group.getCurveStringArray();
rownames = [];
set(handles.table_fit_results, 'RowName', rownames);

%% ANOVA table
set(handles.table_ANOVA, 'ColumnWidth', {45,40,80,90,60,70});
set(handles.table_ANOVA, 'ColumnName', {'Favor?','DoF','Chi2','Chi2Red','F','P(For F>1)'});
set(handles.table_ANOVA, 'RowName', {'Current fit','No Exchange'} );

%% Table: rate analysis
% Set up column titles
column_string = {'dH|(kcal/mol)', 'dS|(cal/mol/K)', 'Ea(A->B)|(kcal/mol)', ...
                 'P(A->B)|(/s)', 'Ea(B->A)|(kcal/mol)', 'P(B->A)|(/s)'};

% Put this all into the table
set(handles.table_rate_analysis, 'ColumnWidth', {90,90,90,100,90,100});
set(handles.table_rate_analysis, 'ColumnName', column_string );

%% Table: fit parameter is OK?
%{
column_string = {''};
for p = 1:session.Np
    column_string{p} = sprintf('%s', session.param_name{p});
end
column_string{6} = 'Rex';
set(handles.table_param_is_ok, 'ColumnName', column_string);
set(handles.table_param_is_ok, 'ColumnWidth', {30});

rownames = group.getCurveStringArray();
set(handles.table_param_is_ok, 'RowName', rownames);
%}

%% Tables: fit paramters OK?
% (1/3) Table for group kinetics
columnNames = [];
set(handles.table_isok_group_kinetics, 'ColumnName', columnNames);
set(handles.table_isok_group_kinetics, 'ColumnWidth', {30});
rowNames = {'kex0', 'PA0', 'dH', 'Eab'};
set(handles.table_isok_group_kinetics, 'RowName', rowNames);

columnNames = {'Temp(C)', 'PA', 'kex'};
set(handles.table_isok_temp_kinetics, 'ColumnName', columnNames);
set(handles.table_isok_temp_kinetics, 'ColumnWidth', {60,30,30});
rowNames = [];
set(handles.table_isok_temp_kinetics, 'RowName', rowNames);
set(handles.table_isok_temp_kinetics, 'ColumnFormat', {'Char', 'Logical', 'Logical'});

% (2/3) Table for curveset dw (ppm)
columnNames = {'Curveset', 'dwH', 'dwX'};
set(handles.table_isok_curveset_ppm, 'ColumnName', columnNames);
set(handles.table_isok_curveset_ppm, 'ColumnWidth', {150, 30, 30});
set(handles.table_isok_curveset_ppm, 'ColumnFormat', {'Char', 'Logical', 'Logical'});
rowNames = [];
set(handles.table_isok_curveset_ppm, 'RowName', rowNames);

% (3/3) Table for curve R20 and Rex
columnNames = {'Curve', 'R20', 'Rex'};
set(handles.table_isok_curve_r, 'ColumnName', columnNames);
set(handles.table_isok_curve_r, 'ColumnWidth', {150, 30, 30});
set(handles.table_isok_curve_r, 'ColumnFormat', {'Char', 'Logical', 'Logical'});
rowNames = [];
set(handles.table_isok_curve_r, 'RowName', rowNames);


%% Select the best fit, by default
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

if( group.Nf > 0 )
    [fitStringArray, f_Best] = group.getFitNameArray();
    set(handles.popup_fit_results, 'Value', f_Best);
end

%% Update
refresh_display(handles);


% This is called to update all display elements
function refresh_display(handles)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

%% Basic setup

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

Nctot   = group.getNumCurves();
Nf      = group.Nf;

% Curve set name
set(handles.text_title, 'String', sprintf('Group: %s',group.name));

%% Fit command options and INITIAL FIT CONDITIONS

% Cannot constrain rates with only one temperature
if( group.getNumTemps() == 1 )
    set(handles.checkbox_constrain_rates, 'Value', 0);
    set(handles.checkbox_constrain_rates, 'Enable', 'off');
    
else
    set(handles.checkbox_constrain_rates, 'Enable', 'on');    
end
CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;

% SINGLE_FIT or GRID_SEARCH
GRID_SEARCH             = get(handles.radiobutton_grid_search, 'Value') == 1;

% Show the appropriate GUI elements for the fitting task
if( GRID_SEARCH )
    % Set up for a grid search
    set(handles.panel_gridsearch, 'Visible', 'on');
    set(handles.panel_ICs, 'Visible', 'off');
    
    if( CONSTRAIN_RATE_ANALYSIS )
        column_string = {'dwH|(ppm)', 'dwX|(ppm)', 'PA0|(%)', 'kex0|(/s)', 'dH|(kcal/mol)', 'Eab|(kcal/mol)'};
        set(handles.table_grid_limits, 'ColumnName', column_string );
        set(handles.table_grid_limits, 'ColumnWidth', {70});

        % List each curve in each curve set in this group
        row_string = {'Min', 'Max', 'Steps'};
        set(handles.table_grid_limits, 'RowName', row_string);
        %set(handles.table_grid_limits, 'Data', zeros(length(row_string), length(column_string)));
        
        % Table: Grid search
        Np_Grid = 6;
        table_data = [  group.grid_p_min(1:Np_Grid)  ; ...
                        group.grid_p_max(1:Np_Grid)  ; ...
                        group.grid_p_steps(1:Np_Grid) ];
        set(handles.table_grid_limits, 'Data', table_data);
        
    else
        column_string = {'dwH|(ppm)', 'dwX|(ppm)', 'PA0|(%)', 'kex0|(/s)'};
        set(handles.table_grid_limits, 'ColumnName', column_string );
        set(handles.table_grid_limits, 'ColumnWidth', {70});

        % List each curve in each curve set in this group
        row_string = {'Min', 'Max', 'Steps'};
        set(handles.table_grid_limits, 'RowName', row_string);
        %set(handles.table_grid_limits, 'Data', zeros(length(row_string), length(column_string)));
        
        % Table: Grid search
        Np_Grid = 4;
        table_data = [  group.grid_p_min(1:Np_Grid)  ; ...
                        group.grid_p_max(1:Np_Grid)  ; ...
                        group.grid_p_steps(1:Np_Grid) ];
        set(handles.table_grid_limits, 'Data', table_data);
    end
    
    % Calibration temperature table
    table_data = group.T0_Grid-273;
    set(handles.table_Temp0_Grid, 'Data', table_data);
    
else
    % Invidiual fit, not a grid search
    set(handles.panel_gridsearch, 'Visible', 'off');
    set(handles.panel_ICs, 'Visible', 'on');

    if( CONSTRAIN_RATE_ANALYSIS )
        
        set(handles.table_ics_group_kinetics, 'Visible', 'on');        
        set(handles.table_ics_temp_kinetics, 'Visible', 'off');
    else
        
        set(handles.table_ics_group_kinetics, 'Visible', 'off');
        set(handles.table_ics_temp_kinetics, 'Visible', 'on');
    end
end

%% Panel: SELECT and view results
if( Nf > 0 )
    set(handles.panel_results, 'Visible', 'on');
    
    %% Fit result list
    % Set the list of fit results
    [fitStringArray, f_Best] = group.getFitNameArray();
    set(handles.popup_fit_results, 'String', fitStringArray );    
    
    % Get the fit number from the list
    f = get(handles.popup_fit_results, 'Value');
    
    % If the value is invalid, set it to the best one on the list
    if( f > Nf )
        set(handles.popup_fit_results, 'Value', f_Best);
        f = f_Best;
    end
    
    % Get the fit result
    [fitResult, isBestFit]  = group.getFitResult(f);    
    
    % Fit manager button
    if( isBestFit )
        set(handles.button_set_best, 'Enable', 'off');
        set(handles.button_set_best_02, 'Enable', 'off');
    else
        set(handles.button_set_best, 'Enable', 'on');
        set(handles.button_set_best_02, 'Enable', 'on');
    end
    
    % Table for rate analysis    
    if( fitResult.CONSTRAIN_RATE_ANALYSIS )
        % Display/Hide GUI elements
        set(handles.text_rates_constrained_unconstrained, 'String', 'Rates CONSTRAINED');
        set(handles.table_results_group_kinetics, 'Visible', 'on');
        set(handles.table_results_temp_kinetics, 'Visible', 'off');
        set(handles.table_rate_analysis, 'Visible', 'off');
        
        % Display table results
        % T0(C)    PA0(%)  kex0(/s)    dH(kcal/mol)    Eab(kcal/mol)        
        T0      = session.convertUnits(group.curvesets{1}.curves{1}.Temp, 'TEMP', 'DISPLAY');
        PA0     = session.convertUnits( fitResult.resultsMatrix(1,3), 'PA', 'DISPLAY');
        PA0_E   = session.convertUnits( fitResult.resultsMatrixErrors(1,3), 'PA', 'DISPLAY');
        
        kex0    = session.convertUnits( fitResult.resultsMatrix(1,4), 'KEX', 'DISPLAY');
        kex0_E  = session.convertUnits( fitResult.resultsMatrixErrors(1,4), 'KEX', 'DISPLAY');
        
        dH      = session.convertUnits( fitResult.dH, 'DH', 'DISPLAY');
        dH_E    = session.convertUnits( fitResult.dH_E, 'DH', 'DISPLAY');
        
        Eab     = session.convertUnits( fitResult.Eab, 'EAB', 'DISPLAY');
        Eab_E   = session.convertUnits( fitResult.Eab_E, 'EAB', 'DISPLAY');
        
        table_data_string       = cell(1,5);
        table_data_string{1,1}  = T0;
        table_data_string{1,2}  = sprintf('  %s[%s]', displayNumber(PA0), displayNumber(PA0_E));
        table_data_string{1,3}  = sprintf('  %s[%s]', displayNumber(kex0), displayNumber(kex0_E));
        table_data_string{1,4}  = sprintf('  %s[%s]', displayNumber(dH), displayNumber(dH_E));
        table_data_string{1,5}  = sprintf('  %s[%s]', displayNumber(Eab), displayNumber(Eab_E));
        set(handles.table_results_group_kinetics, 'Data', table_data_string );

    else
        % Display/Hide GUI elements
        set(handles.text_rates_constrained_unconstrained, 'String', 'Rates NOT Constrained');
        set(handles.table_results_group_kinetics, 'Visible', 'off');
        set(handles.table_results_temp_kinetics, 'Visible', 'on');
        set(handles.table_rate_analysis, 'Visible', 'on');
                
        % Table: PA and kex at each temperature
        [Temp_unique_Array, PA_Array, PA_E_Array, kex_Array, kex_E_Array] = fitResult.getKineticsValues();
        NTemps = length(Temp_unique_Array);
        
        table_data_string       = cell(NTemps,3);
        for t = 1:NTemps
            PA      = session.convertUnits(PA_Array(t), 'PA', 'DISPLAY');
            PA_E    = session.convertUnits(PA_E_Array(t), 'PA', 'DISPLAY');        
            kex      = session.convertUnits(kex_Array(t), 'KEX', 'DISPLAY');
            kex_E    = session.convertUnits(kex_E_Array(t), 'KEX', 'DISPLAY');     
            
            table_data_string{t,1}  = session.convertUnits(Temp_unique_Array(t), 'TEMP', 'DISPLAY');
            table_data_string{t,2}  = sprintf('  %s[%s]', displayNumber(PA), displayNumber(PA_E));
            table_data_string{t,3}  = sprintf('  %s[%s]', displayNumber(kex), displayNumber(kex_E));            
        end
        set(handles.table_results_temp_kinetics, 'Data', table_data_string );
        
        
        % Table: Rate analysis (vant Hoff and Arrhenius)
        % Get results table from rate analysis data structure
        table_data = fitResult.rateAnalysis.getRateAnalysisTableString();        
        set(handles.table_rate_analysis, 'Data', table_data );        
    end
    
    % Table: summary of results for each curve
    % _NAME_        dwH(ppm)  dwX(ppm)  Pa  kex  R20  Rex  alpha
    % 600MHz, 25C      0        1.5     0.8 1000  20   10    1.1
    %rownames = group.getCurveStringArray();    
    %set(handles.table_fit_results, 'RowName', rownames);
    table_data_string = fitResult.getResultsMatrixString();
    set(handles.table_fit_results, 'Data', table_data_string );
    
    % Table: PPM Shifts 
    % CurvesetName  dwH(ppm)  dwX(ppm)    
    table_data_string = cell(1,3);    
    for cs = 1:group.Ncs
        curveset_names = group.getCurvesetStringArray();
        table_data_string{cs,1} = curveset_names{cs};

        [dwHppm, dwHppm_E] = fitResult.getdwHppm(cs);
        table_data_string{cs,2}  = sprintf('  %s[%s]', displayNumber(dwHppm), displayNumber(dwHppm_E));    

        [dwXppm, dwXppm_E] = fitResult.getdwXppm(cs);    
        table_data_string{cs,3}  = sprintf('  %s[%s]', displayNumber(dwXppm), displayNumber(dwXppm_E));    
    end    
    set(handles.table_results_curveset_ppm, 'Data', table_data_string );
    
    
else
    % No fits for this curve set
    f = 0;
    set(handles.panel_results, 'Visible', 'off');
end



%% Table: ANOVA model selection
%{
% If the current fit is available
if( f > 0 )
    % Enable the table
    set(handles.table_ANOVA, 'Enable', 'on');
    
    % Initialize table_data_string    
    table_data_string = cell( 2,6 );
    
    % Check if the exchange model MSE(chi2Red) is less than no-exchange
    %  => Exchange model is favored
    if( fitResult.chi2red <= group.fitResult_NoEx.chi2red )
        table_data_string{1,1}  = '    X';
        table_data_string{2,1}  = '';        
    % No-exchange model is favored
    else
        table_data_string{1,1}  = '';
        table_data_string{2,1}  = '    X';        
    end

    table_data_string{1,2}  = sprintf('  %d',     fitResult.df );
    table_data_string{1,3}  = sprintf(' %0.2f',  fitResult.chi2 );
    table_data_string{1,4}  = sprintf('   %0.2f',  fitResult.chi2red );
    table_data_string{1,5}  = sprintf(' %0.2f',  fitResult.Fstatistic );
    table_data_string{1,6}  = sprintf(' %0.1e',  fitResult.Pvalue );
    
    table_data_string{2,2}  = sprintf('  %d',     group.fitResult_NoEx.df );
    table_data_string{2,3}  = sprintf(' %0.2f',  group.fitResult_NoEx.chi2 );
    table_data_string{2,4}  = sprintf('   %0.2f',  group.fitResult_NoEx.chi2red );
    table_data_string{2,5}  = sprintf('  -');
    table_data_string{2,6}  = sprintf('  -');
    set(handles.table_ANOVA, 'Data', table_data_string );
   
% No exchange fit to compare to
else
    % Disable the table
    set(handles.table_ANOVA, 'Enable', 'off');    
    set(handles.table_ANOVA, 'Data', []);
end
%}




%% Plot: fit results on inset axes
figure_handle = handles.panel_display_dispersion;

% Clear prior panel contents
children = get(figure_handle, 'Children');
for ch = 1:length(children)    
    % Delete it if its not the checkbox (auto-Y)
    if( children(ch) ~= handles.checkbox_y_scale )
        delete( children(ch) );
    end
end

% Make new axes in panel
x = 0.05;   % No room for label
y = 0.1;    % No room for label
%x = 0.11;  % Fits label
%y = 0.21;  % Fits label
h = axes('position', [x y 0.99-x 0.98-y],'Parent', figure_handle);
cla(h);
hold(h, 'off');

% Legend
lstring     = group.getCurveStringArray();

ctot = 0;
for cs = 1:group.Ncs
    curveset = group.curvesets{cs};
    for c = 1:curveset.Nc
        ctot = ctot+1;
        curve = curveset.curves{c};
        
        X   = curve.vcpmg;
        Y   = curve.R2eff;
        YE  = curve.eR2eff;
    
        errorbar( h, X, Y, YE, ...
            'o','MarkerSize',session.MARKERSIZE, ...
            'LineWidth', session.LINEWIDTH_SMALL );

        % Overlay subsequent plots    
        if( ctot == 1 );
            hold(h, 'all');
        end
    end
end

% Now plot the fits to the curves if available
% Done separately from data so the legend shows one entry for line
% and one entry for dispersion model curve
[void, vcpmg_max, void, void] = group.getDataLimits();

if( f > 0 )    
    Nx = 50;
    X = linspace(0, 1.1*vcpmg_max, Nx);
    
    for ctot = 1:Nctot
        curve = group.curvesets{cs}.curves{c};
        
        % Check if the fit corresponds to exchange-model (Pa ~= 1)
        if( fitResult.pMatrixBest(ctot,3) ~= 1 )
            mR2eff = model_MQRD_CRJ(X, curve.TCPMG, fitResult.pMatrixBest(ctot,:));
            
        % No exchange looks like constant R2Eff(vCPMG)=R20
        else            
            mR2eff = fitResult.pMatrixBest(ctot,5) .* ones(1,length(X));
        end

        % Plot the MQ R2eff fit with the rest of the data        
        plot(h, X, mR2eff, 'lineStyle','-', ...
            'LineWidth', session.LINEWIDTH_SMALL, 'Color', 'black' );       
    end
    
    % Only show legend entry once for line and one for dispersion       
    lstring{end+1} = sprintf('Fit');
end


% Set title and axes for CPMG data
set( h,'FontName', session.FONTNAME, 'FontSize', session.FONTSIZE_SMALL, ...
    'FontWeight', 'Normal', 'LineWidth', session.LINEWIDTH)
%xlabel(h, '\nu_C_P_M_G (Hz)');
%ylabel(h, 'R_2^e^f^f (Hz)')

if( get(handles.checkbox_y_scale, 'Value') )
    [void, void, R2eff_min, R2eff_max] = group.getDataLimits();
    ylim(  h, [0.9*R2eff_min, 1.1*R2eff_max] );
else
    [void, void, R2eff_min, R2eff_max] = session.getDataLimits();
    ylim(  h, [0.9*R2eff_min, 1.1*R2eff_max] );
end

xlim(  h, [ 0 1.1*vcpmg_max ] );
%set(   h, 'YGrid','on');

hl=legend(h, lstring{1,:}, 'Location', 'EastOutside' );
set(hl, 'FontSize', session.FONTSIZE_EXTRA_SMALL-2);


%% Notes on crosspeak
%nlines = size(group.notes,1);
%for l = 1:nlines
%    group.notes{l}
%end
set(handles.edit_notes, 'String', group.note.text );

% Disable save and notes buttons
set(handles.button_save_notes, 'Enable', 'off');
set(handles.button_restore_notes, 'Enable', 'off');

%% Fit parameter OK?
Nctot   = group.getNumCurves();
Np      = group.parentSession.Np;



% If fit exists
if( f > 0 )    
    set(handles.panel_param_isok, 'Visible', 'on');
    
    % (1/3) Table for group kinetics
    % Table type (A)
    if( fitResult.CONSTRAIN_RATE_ANALYSIS )
        table_data = [  fitResult.kex0_isOK; ...
                        fitResult.PA0_isOK; ...
                        fitResult.dH_isOK; ...
                        fitResult.Eab_isOK ];    
        set(handles.table_isok_group_kinetics, 'Data', table_data);
        
        set(handles.table_isok_group_kinetics, 'Visible', 'on');
        set(handles.table_isok_temp_kinetics, 'Visible', 'off');
    else
        % Table type (B)
        % Each temperature's PA and kex can be OK or not
        [Temp_unique_Array, PA_isOK_Array, kex_isOK_Array] = fitResult.getKineticsIsOK();
        NTemps = length(Temp_unique_Array);
        
        % Temp kex(/s) PA(%)
        table_data = cell(NTemps,3);
        for t = 1:NTemps
            table_data{t,1} = Temp_unique_Array(t)-273;
            table_data{t,2} = PA_isOK_Array(t)==1;
            table_data{t,3} = kex_isOK_Array(t)==1;
        end
        %table_data = [ (Temp_unique_Array-273)' PA_isOK_Array' kex_isOK_Array'];        
        set(handles.table_isok_temp_kinetics, 'Data', table_data);
        
        set(handles.table_isok_group_kinetics, 'Visible', 'off');
        set(handles.table_isok_temp_kinetics, 'Visible', 'on');
    end

    % (2/3) Table for curveset dw (ppm)    
    % Check if the dw is OK for each curveset
    table_data      = cell(group.Ncs,3);
    curveset_names  = group.getCurvesetStringArray();
    for cs = 1:group.Ncs
        % Find the first "cot" for the curveset
        ctot_array = group.getCurveTotArray(cs);
        ctot = ctot_array(1);
            
        % Set the curveset
        table_data{cs,1} = curveset_names{cs};
        table_data{cs,2} = fitResult.resultsMatrixIsOK(ctot, 1)==1;
        table_data{cs,3} = fitResult.resultsMatrixIsOK(ctot, 2)==1;
    end
    set(handles.table_isok_curveset_ppm, 'Data', table_data);

    % (3/3) Table for curve R20 and Rex
    table_data  = cell(Nctot,3);
    curve_names = group.getCurveStringArray();
    for ctot = 1:Nctot
        table_data{ctot,1} = curve_names{ctot};
        table_data{ctot,2} = fitResult.resultsMatrixIsOK(ctot, 5)==1;
        table_data{ctot,3} = fitResult.resultsMatrixIsOK(ctot, 6)==1;
    end
    set(handles.table_isok_curve_r, 'Data', table_data);
    
    set(handles.checkbox_exhibits_exchange,     'Value', group.exhibitsExchange);
    set(handles.checkbox_best_fit_is_ok, 'Value', group.bestFitIsOK);

% No fit
else
    set(handles.panel_param_isok, 'Visible', 'off');    
end


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_fit_rd_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function table_grid_limits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to table_grid_limits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%{
% --- Executes on button press in button_grid_search.
function button_grid_search_Callback(hObject, eventdata, handles)
% hObject    handle to button_grid_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Fit via grid search
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
data            = getappdata(handles_main.main_gui, 'data');
settings        = getappdata(handles_main.main_gui, 'settings');
fit_results     = getappdata(handles_main.main_gui, 'fit_results');
cs              = get(handles_main.listbox_cs, 'Value');

%% Check if user wants to override prior grid search
% If the fit IS available AND the user wants to cancel, then abort (return)
if( isstruct(fit_results) && isfield(fit_results, 'grid_is_available') && fit_results.grid_is_available(cs) && ...
    strcmp(questdlg('This will replace prior grid search results. Do you wish to proceed?', ...
        'Replace prior grid search?', 'Proceed', 'Cancel', 'Cancel'),'Cancel') ...
   )
    return
end

%% Execute grid search

% Perform grid search
[fit_results, grid_finished] = grid_search(data, settings, cs, fit_results, ...
                handles.axes_current_progress, handles.toggle_abort );

if( grid_finished )
    %fprintf('\n\nAfter grid search subroutine fit_results.fits_local');
    %fit_results.fits_local.sse_local    
    
    %fprintf('\n-- POST GRID SEARCH -- ');
    %fit_results.analysis_is_available

    % Store the results to the main GUI    
    setappdata(handles_main.main_gui, 'fit_results', fit_results);
else
    
    % Turn abort button OFF and to color BLACK
    set(handles.toggle_abort, 'Value', false);
    set(handles.toggle_abort, 'ForegroundColor', [0 0 0]);    
end

refresh_display(handles);
%}

% --- Executes when entered data in editable cell(s) in table_grid_limits.
function table_grid_limits_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_grid_limits (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Enable use of the grid search table buttons
set(handles.button_save_grid_search, 'Enable', 'on');
set(handles.button_restore_grid_search, 'Enable', 'on');


% --- Executes on button press in button_save_grid_search.
function button_save_grid_search_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_grid_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Read the table, set new values and disable use of buttons

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Read data from table
table_data = get(handles.table_grid_limits, 'Data');

% Constraining the rate analysis determines what of the grid is changed
CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;

if( CONSTRAIN_RATE_ANALYSIS )
    Np_Grid = 6;
    
else
    %  Without constraints, only a subset of the grid is updated
    Np_Grid = size(table_data,2);
end

% Get prior grid so it can be partially updated, if needed
grid_p_min              = group.grid_p_min;
grid_p_max              = group.grid_p_max;
grid_p_steps            = group.grid_p_steps;

% Make update
grid_p_min(1:Np_Grid)    = table_data(1,1:Np_Grid);
grid_p_max(1:Np_Grid)    = table_data(2,1:Np_Grid);
grid_p_steps(1:Np_Grid)  = table_data(3,1:Np_Grid);
group.setGridLimits( grid_p_min, grid_p_max, grid_p_steps )

% Update display and disable use of buttons
set(handles.button_save_grid_search, 'Enable', 'off');
set(handles.button_restore_grid_search, 'Enable', 'off');
refresh_display(handles);


% --- Executes on button press in button_restore_grid_search.
function button_restore_grid_search_Callback(hObject, eventdata, handles)
% hObject    handle to button_restore_grid_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Update display and disable use of buttons
refresh_display(handles);
set(handles.button_save_grid_search, 'Enable', 'off');
set(handles.button_restore_grid_search, 'Enable', 'off');


% --- Executes on selection change in popup_fit_results.
function popup_fit_results_Callback(hObject, eventdata, handles)
% hObject    handle to popup_fit_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_fit_results contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_fit_results
%% Update display for newly selected fit
refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function popup_fit_results_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_fit_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_set_best.
function button_set_best_Callback(hObject, eventdata, handles)
% hObject    handle to button_set_best (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Set current fit as best fit
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

if( group.Nf > 0 )
    % Get the fit number from the list
    f = get(handles.popup_fit_results, 'Value');
    
    % If the value is invalid, set it to the first one on the list
    if( f > group.Nf )
        set(handles.popup_fit_results, 'Value', 1);
        f = 1;
    end
        
    % Commit the selection
    group.setBestFitResult(f);
end

refresh_display(handles);


% --- Executes on button press in button_delete_fit.
function button_delete_fit_Callback(hObject, eventdata, handles)
% hObject    handle to button_delete_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Delete current selected fit

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

if( group.Nf > 0 )    
    % Get the fit number from the list
    f           = get(handles.popup_fit_results, 'Value');
    fitResult   = group.fitResults{f};    
    
    % If this is the NoEx fit then abort
    if( fitResult == group.fitResult_NoEx )
        errordlg( {'Cannot remove no-exchange fit'; 'It is required for model-comparison'}, 'No exchange fit' )
        
    % This is NOT the NoEx fit, so it can be deleted
    else
    
        % Ask user if they want to delete fit
        if( strcmp(questdlg('Are you sure you want to delete this fit?', ...
            'Delete fit?', 'Delete fit', 'Cancel','Cancel'),'Delete fit') )
        
            group.removeFitResult(fitResult);            
        end
    end
end

refresh_display(handles);

% --- Executes on button press in button_errors_in_dispersion.
function button_errors_in_dispersion_Callback(hObject, eventdata, handles)
% hObject    handle to button_errors_in_dispersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};


if( group.Nf > 0 )    
    % Get the fit number from the list
    f           = get(handles.popup_fit_results, 'Value');
    fitResult   = group.fitResults{f};   
    
    fprintf('\nMC error analysis (1-6000 sec, see progress bar)');    
    fitResult.calculateErrors( handles.axes_current_progress, handles.toggle_abort );
    
end

refresh_display(handles);

% --- Executes on button press in button_save_notes.
function button_save_notes_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Save notes

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Save the new text notes the user entered
getNotes = get(handles.edit_notes, 'String');

% % Convert each line of text into a cell array
% Nlines      = size(getNotes,1);
% notesString = cell(Nlines,1);
% for l = 1:Nlines
%    notesString{l} = getNotes{l}; 
% end

% Save it
group.setNote( getNotes );

% Disable save and notes buttons
set(handles.button_save_notes, 'Enable', 'off');
set(handles.button_restore_notes, 'Enable', 'off');

refresh_display(handles);

% --- Executes on button press in button_restore_notes.
function button_restore_notes_Callback(hObject, eventdata, handles)
% hObject    handle to button_restore_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Undo changes to notes

% Disable save and notes buttons
set(handles.button_save_notes, 'Enable', 'off');
set(handles.button_restore_notes, 'Enable', 'off');

% This will erase what the user typed before saving it
refresh_display(handles);


% --- Executes when entered data in editable cell(s) in table_ics.
function table_ics_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_ics (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_set_ics.
function button_set_ics_Callback(hObject, eventdata, handles)
% hObject    handle to button_set_ics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% To set fit results as initial conditions in table
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

if( group.Nf > 0 )    
    % Get the fit number from the list
    f = get(handles.popup_fit_results, 'Value');    
    fitResult = group.fitResults{f};
    
    % Set the GROUP temperature dependent properties
    T0      = group.curvesets{1}.curves{1}.Temp;
    PA0     = session.convertParameterToDisplayUnits( fitResult.pMatrixBest(1,3), 'PA' );
    kex0    = session.convertParameterToDisplayUnits( fitResult.pMatrixBest(1,4), 'KEX' );
    dH      = session.convertParameterToDisplayUnits( fitResult.dH, 'DH' );
    Eab     = session.convertParameterToDisplayUnits( fitResult.Eab, 'EAB' );    
    set(handles.table_ics_group_kinetics, 'Data', [T0-273, PA0, kex0, dH, Eab]);
    
    % Set the kex and PA at each temperature
    [Temp_unique_Array, PA_Array, PA_E_Array, kex_Array, kex_E_Array] = fitResult.getKineticsValues();
    set(handles.table_ics_temp_kinetics, 'Data', [(Temp_unique_Array-273)' 100*PA_Array' kex_Array']);
    
    % Set the CURVESET chemical shifts
    table_data = cell(group.Ncs,3);
    curveset_names = group.getCurvesetStringArray();
    for cs = 1:group.Ncs
        table_data{cs,1}    = curveset_names{cs};
        table_data{cs,2}    = fitResult.getdwHppm(cs); 
        table_data{cs,3}    = fitResult.getdwXppm(cs);
    end
    set(handles.table_ics_curveset_ppm, 'Data', table_data);
    
    % Set the CONSTRAIN RATES
    CONSTRAIN_RATE_ANALYSIS = fitResult.CONSTRAIN_RATE_ANALYSIS;
    set(handles.checkbox_constrain_rates, 'Value', CONSTRAIN_RATE_ANALYSIS);
    group.updateFitParams(CONSTRAIN_RATE_ANALYSIS);
end

refresh_display(handles);


% --- Executes on button press in checkbox_best_fit_is_ok.
function checkbox_best_fit_is_ok_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_best_fit_is_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_best_fit_is_ok
%% Update curve set as exhibiting exchange or not

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

group.setBestFitIsOK( get(handles.checkbox_best_fit_is_ok, 'Value') );
refresh_display(handles);

% --- Executes on button press in checkbox_exhibits_exchange.
function checkbox_exhibits_exchange_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_exhibits_exchange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_exhibits_exchange
%% Update curve set as exhibiting exchange or not

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

group.setExhibitsExchange( get(handles.checkbox_exhibits_exchange, 'Value')==1 );
refresh_display(handles);


% --- Executes on key press with focus on edit_notes and none of its controls.
function edit_notes_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%% If user pressed key in notes window
% Enable save and notes buttons
set(handles.button_save_notes, 'Enable', 'on');
set(handles.button_restore_notes, 'Enable', 'on');



% --- Executes on button press in toggle_abort.
function toggle_abort_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_abort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_abort

%% Abort current progress?

% Color black/red for unaborted/aborted
color = [0 0 0];
if( get(handles.toggle_abort, 'Value') )
    color = [1 0 0];
end

set(handles.toggle_abort, 'ForegroundColor', color);


% --- Executes on button press in button_all_param_ok.
function button_all_param_ok_Callback(hObject, eventdata, handles)
% hObject    handle to button_all_param_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Set all parameters to OK

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Get the fit number from the list
f = get(handles.popup_fit_results, 'Value');
[fitResult, isBestFit]  = group.getFitResult(f);


fitResult.setParamIsOK( true, 'ALL');
refresh_display(handles);



% --- Executes on button press in button_invert_param_ok.
function button_invert_param_ok_Callback(hObject, eventdata, handles)
% hObject    handle to button_invert_param_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Set all parameters to NOT( current value )
%{

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Get the data after the change
table_data              = cell2mat( get(handles.table_param_is_ok, 'Data') );

% Set all to opposite
table_data(1:end,1:end) = ~table_data(1:end,1:end);

% Get the fit number from the list
f           = get(handles.popup_fit_results, 'Value');
fitResult   = group.fitResults{f};

fitResult.setResultsMatrixIsOK( table_data );
refresh_display(handles);
%}

% --- Executes on button press in button_none_param_ok.
function button_none_param_ok_Callback(hObject, eventdata, handles)
% hObject    handle to button_none_param_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Set all parameters to NOT OK

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Get the fit number from the list
f = get(handles.popup_fit_results, 'Value');
[fitResult, isBestFit]  = group.getFitResult(f);

fitResult.setParamIsOK( false, 'ALL');
refresh_display(handles);



% --- Executes on button press in checkbox_y_scale.
function checkbox_y_scale_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_y_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_y_scale

%% Change scaling on Y-Axis for R2eff
refresh_display(handles);

%{
% --- Executes when entered data in editable cell(s) in table_param_is_ok.
function table_param_is_ok_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_param_is_ok (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Save changes to param_is_ok

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Check which curve and parameter has been changed
data_changed    = eventdata.Indices;
ctot_changed    = data_changed(1);
p_changed       = data_changed(2);

columnNames     = get(handles.table_param_is_ok, 'ColumnName');

% Get the data after the change
table_data      = cell2mat( get(handles.table_param_is_ok, 'Data') );

% Change all other co-dependent variables in table_data
%   If user did not pick Rex
switch upper(columnNames{p_changed})
    
    case {'DWH', 'DWX', 'PA', 'KEX', 'R20'}
        % Get the independent parameter number
        indepParam = group.indepParamIndexMatrix(ctot_changed, p_changed);

        % Find other curves and paramters which are dependent upon this value
        [C,P,V] = find( group.indepParamIndexMatrix == indepParam );

        % Change these elements to new value, too
        for cp = 1:length(C)
            table_data(C(cp), P(cp) ) = eventdata.NewData;
        end
        
    case 'REX'
        % Nothing else to do (One Rex per curve)
        
    case 'DH'
        % Change all other DH values
        
    case 'EAB'
        % Change all other EAB values
        
    otherwise
        error('Invalid parameter changed');        
end
set(handles.table_param_is_ok, 'Data', table_data);

% Commit this change to the current fitResult
f           = get(handles.popup_fit_results, 'Value');
fitResult   = group.fitResults{f};
fitResult.setResultsMatrixIsOK( table_data );

refresh_display(handles);
%}

% --- Executes on button press in button_r_star.
function button_r_star_Callback(hObject, eventdata, handles)
% hObject    handle to button_r_star (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Set all parameters to NOT OK

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Get the fit number from the list
f = get(handles.popup_fit_results, 'Value');
[fitResult, isBestFit]  = group.getFitResult(f);

fitResult.setParamIsOK( true, 'R*');
refresh_display(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_notes.
function edit_notes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% If user pressed key in notes window (2010/08/30)
% Enable save and notes buttons
set(handles.button_save_notes, 'Enable', 'on');
set(handles.button_restore_notes, 'Enable', 'on');


% --- Executes on button press in button_go.
function button_go_Callback(hObject, eventdata, handles)
% hObject    handle to button_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');    
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% SINGLE_FIT or GRID_SEARCH
GRID_SEARCH     = get(handles.radiobutton_grid_search, 'Value') == 1;

% OPTIMIZE or SIMULATE
OPTIMIZE_FIT    = get(handles.radiobutton_fit, 'Value') == 1;

% Fit or simulate the data using initial conditions from table
CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;

if( GRID_SEARCH )    
    % Fit via grid search
    
    %% Check if the grid in the table matches the grid in the group
    %  (Sometimes the user forgets to *Save* the table)
    
    % Read data from table
    table_data = get(handles.table_grid_limits, 'Data');
    grid_p_min              = table_data(1,:);
    grid_p_max              = table_data(2,:);
    grid_p_steps            = table_data(3,:);

    % Constraining the rate analysis determines what of the grid is changed
    CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;

    if( CONSTRAIN_RATE_ANALYSIS )
        Np_Grid = 6;
    else
        %  Without constraints, only a subset of the grid is updated
        Np_Grid = size(table_data,2);
    end

    % Check to see if the there is a mismatch
    if( any( grid_p_min(1:Np_Grid)      ~= group.grid_p_min(1:Np_Grid) ) || ...
        any( grid_p_max(1:Np_Grid)      ~= group.grid_p_max(1:Np_Grid) ) || ...
        any( grid_p_steps(1:Np_Grid)    ~= group.grid_p_steps(1:Np_Grid) ) )        
        msgbox('Grid search has been changed, but not saved. Please use "Save" or "Restore" before starting');
        return
    end    
    
    %% Check if user wants to override prior grid search    
    if( group. gridSearchIsDone() && ...
        strcmp(questdlg('This will replace prior grid search results. Do you wish to proceed?', ...
            'Replace prior grid search?', 'Proceed', 'Cancel', 'Cancel'),'Cancel') )
        return
    end   
    
    if( OPTIMIZE_FIT )
        taskString = 'FIT';
    else
        taskString = 'SIM';
    end
    
    fprintf('\nGrid search (1-6000 sec, see progress bar)');
    
    % Call grid search
    group.gridSearch( taskString, CONSTRAIN_RATE_ANALYSIS, handles.axes_current_progress, handles.toggle_abort );        
    
    % Turn abort button OFF and to color BLACK
    set(handles.toggle_abort, 'Value', false);
    set(handles.toggle_abort, 'ForegroundColor', [0 0 0]);

else
    %% Create new fit result to perform the desired function
    if( OPTIMIZE_FIT )
        fitResult   = FitResult(group, 'EXCHANGE', 'FIT', CONSTRAIN_RATE_ANALYSIS);        
    else
        fitResult   = FitResult(group, 'EXCHANGE', 'SIM', CONSTRAIN_RATE_ANALYSIS);
    end

    % Get initial ppm values for each curve
    table_data  = get(handles.table_ics_curveset_ppm, 'Data');
    Nrows       = size(table_data,1);

    dwHppm_csArray  = zeros(Nrows, 1);
    dwXppm_csArray  = zeros(Nrows, 1);
    for r = 1:Nrows
        % CurvesetName in column 1, dwH and dwX in 2 and 3
        dwHppm_csArray(r)  = table_data{r,2};
        dwXppm_csArray(r)  = table_data{r,3};
    end
    fitResult.setInitial_Shifts(dwHppm_csArray, dwXppm_csArray)
    
    % Debugging (2011/06/17)
    %  This was used to make SIM fits that exactly matched simulation
    %  conditions (otherwise I cannot set R20's initial value)
    SET_CONSTANT_R20 = false;
    if( SET_CONSTANT_R20 )
        R20_Fixed = 10;
        fprintf('\n\n!!GUARDD_fit_rd.m: Setting initial R20 to a fixed value of %f\n\n', R20_Fixed);
        fitResult.setInitial_Shifts(dwHppm_csArray, dwXppm_csArray, R20_Fixed);
    end
         

    % Set the initial kinetics (kex and PA)
    if( CONSTRAIN_RATE_ANALYSIS )
        % Get initial rate constraints from table
        table_data = get(handles.table_ics_group_kinetics, 'Data');                
        T0      = table_data(1,1) + 273;
        PA0     = table_data(1,2) / 100;
        kex0    = table_data(1,3);
        dH      = table_data(1,4) * 1e3;
        Eab     = table_data(1,5) * 1e3;

        % Set initial conditions, and temperature-dependent parameters
        fitResult.setInitial_Kinetics_ConstrainedRates(T0, PA0, kex0, dH, Eab);

    else
        % Get kex and PA at each temperature
        table_data = get(handles.table_ics_temp_kinetics, 'Data');
        Nrows       = size(table_data,1);

        Temp_Array  = zeros(Nrows, 1);
        PA_Array    = zeros(Nrows, 1);
        kex_Array   = zeros(Nrows, 1);

        % One temperature at each row
        for r = 1:Nrows
            Temp_Array(r)   = table_data(r,1)+273;
            PA_Array(r)     = table_data(r,2)/100;
            kex_Array(r)    = table_data(r,3);
        end
        % Set initial conditions
        fitResult.setInitial_Kinetics_UnconstrainedRates(Temp_Array, PA_Array, kex_Array);
    end
    
    [pMatrixLower, pMatrixUpper] = group.getDefaultFitLimits();
    fitResult.setLimits_pMatrix( pMatrixLower, pMatrixUpper );
        
        %{
    try
    catch e
        msgbox('Initial condition(s) out of range, check command window', ...
               'ICs out of range','error');
       return;
    end
        %}
    
    fprintf('\nFitting (1-60 sec)');
        
    % Start either optimizing or simulating
    updateGroupParamsFlag = true;
    if( OPTIMIZE_FIT )
        fitResult.fitMe( updateGroupParamsFlag );        
    else
        fitResult.simMe( updateGroupParamsFlag );
    end
    fitResult.analyzeMe();
        
    % Add this result to the group
    group.addFitResult( fitResult );
    
    %fprintf('\n\tDone (%0.1f sec)\n', fitResult.timeRequired);
end

%% Save app data
%{
variables       = getappdata(handles_main.main_gui, 'variables');
default_filename    = 'RDA-Session--Batch_Progress.mat';
default_filepath    = session.outputDir;

% Make the directory if it does not yet exist
if( ~exist(session.outputDir, 'dir') )
    mkdir(session.outputDir);
end

fprintf('\n\nFinished group, wrote to HDD');

% Check each variable and save
for v = 1:length(variables)
    fprintf('\n%s\t', variables{v} );
    if(length(variables{v})<8);   fprintf('\t');  end

    % If variable exists
    if( isappdata(handles_main.main_gui, variables{v}) )

        eval(sprintf('%s = getappdata(handles_main.main_gui, variables{v});',variables{v}));

        % If outupt file exists, append it otherwise create a new one
        if( exist(sprintf('%s/%s', default_filepath,default_filename), 'file') )
            save( sprintf('%s/%s', default_filepath, default_filename), variables{v}, '-append' );
        else
            save( sprintf('%s/%s', default_filepath, default_filename), variables{v} );
        end

        fprintf('\tFound and saved');
    else
        fprintf('\tNo info to save');
    end
end    
%}

% Select it from the list
set(handles.popup_fit_results, 'Value', group.Nf );

refresh_display(handles);

% --- Executes when selected object is changed in panel_procedure.
function panel_procedure_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_procedure 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

refresh_display(handles);


% --- Executes when selected object is changed in panel_task.
function panel_task_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_task 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_constrain_rates.
function checkbox_constrain_rates_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_constrain_rates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_constrain_rates

%% Constrain the rate analysis?
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');    
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;
group.updateFitParams(CONSTRAIN_RATE_ANALYSIS);

update_table_ics(handles);
refresh_display(handles);


% --- Executes when entered data in editable cell(s) in table_constrain_rates.
function table_constrain_rates_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_constrain_rates (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

refresh_display(handles);

% --- Executes when entered data in editable cell(s) in table_ics_group_kinetics.
function table_ics_group_kinetics_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_ics_group_kinetics (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
update_table_ics(handles);
refresh_display(handles);

% --- Executes when entered data in editable cell(s) in table_ics_curveset_ppm.
function table_ics_curveset_ppm_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_ics_curveset_ppm (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
update_table_ics(handles);
refresh_display(handles);


% --- Executes when entered data in editable cell(s) in table_Temp0_Grid.
function table_Temp0_Grid_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_Temp0_Grid (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% Change calibration temperature
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');    
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

T0_Grid = get(handles.table_Temp0_Grid, 'Data')+273;

% Must be a valid temperature
if( T0_Grid > 0 )
    group.setT0_Grid( T0_Grid );    
end

refresh_display(handles);


function update_table_ics( handles )
%% update_table_ics

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');    
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};
Nctot   = group.getNumCurves();

CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;

% Create a new FitResult so that the initial conditions can be set
fitResult   = FitResult(group, 'EXCHANGE', 'FIT', CONSTRAIN_RATE_ANALYSIS);

% Get initial ppm values for each curve
table_data  = get(handles.table_ics_curveset_ppm, 'Data');
Nrows       = size(table_data,1);

dwHppm_csArray  = zeros(Nrows, 1);
dwXppm_csArray  = zeros(Nrows, 1);
for r = 1:Nrows
    % CurvesetName in column 1, dwH and dwX in 2 and 3
    dwHppm_csArray(r)  = table_data{r,2};
    dwXppm_csArray(r)  = table_data{r,3};
end
fitResult.setInitial_Shifts(dwHppm_csArray, dwXppm_csArray)

% Set the initial kinetics (kex and PA)
if( CONSTRAIN_RATE_ANALYSIS )
    % Get initial rate constraints from table
    table_data = get(handles.table_ics_group_kinetics, 'Data');                
    T0      = table_data(1,1) + 273;
    PA0     = table_data(1,2) / 100;
    kex0    = table_data(1,3);
    dH      = table_data(1,4) * 1e3;
    Eab     = table_data(1,5) * 1e3;
    
    % Set initial conditions, and temperature-dependent parameters
    fitResult.setInitial_Kinetics_ConstrainedRates(T0, PA0, kex0, dH, Eab);
    
else
    % Get kex and PA at each temperature
    table_data = get(handles.table_ics_temp_kinetics, 'Data');
    Nrows       = size(table_data,1);

    Temp_Array  = zeros(Nrows, 1);
    PA_Array    = zeros(Nrows, 1);
    kex_Array   = zeros(Nrows, 1);
    
    % One temperature at each row
    for r = 1:Nrows
        Temp_Array(r)   = table_data(r,1)+273;
        PA_Array(r)     = table_data(r,2)/100;
        kex_Array(r)    = table_data(r,3);
    end
    
    fitResult.setInitial_Kinetics_UnconstrainedRates(Temp_Array, PA_Array, kex_Array);
end

% Build the information for the table
% 'Curveset[Curve]', 'Temp(C)', 'B0(MHz)', 'SQX', 'dwH(Hz)', 'dwX(Hz)', 'PA(%)', 'kex(/s)', 'R20(Hz)'
pMatrixInit = fitResult.pMatrixInit;

curveString_Array = group.getCurveStringArray();

table_data_string = cell(Nctot,9);
for ctot = 1:Nctot
    [cs,c] = group.getCurvesetCurve(ctot);
    curve = group.curvesets{cs}.curves{c};
    
    table_data_string{ctot, 1} = curveString_Array{ctot};
    %table_data_string{ctot, 1}  = sprintf(' %d', cs);
    %table_data_string{ctot, 2}  = sprintf(' %d', c);
    table_data_string{ctot, 2}  = sprintf(' %0.1f', curve.Temp-273);
    table_data_string{ctot, 3}  = sprintf(' %0.1f', curve.B0);
    table_data_string{ctot, 4}  = curve.SQX==1;
    
    p = 1; table_data_string{ctot, 4+p}  = sprintf(' %0.1f', ...
        session.convertParameterToDisplayUnits(pMatrixInit(ctot,p),p));
    p = 2; table_data_string{ctot, 4+p}  = sprintf(' %0.1f', ...
        session.convertParameterToDisplayUnits(pMatrixInit(ctot,p),p));
    p = 3; table_data_string{ctot, 4+p}  = sprintf(' %0.1f', ...
        session.convertParameterToDisplayUnits(pMatrixInit(ctot,p),p));
    p = 4; table_data_string{ctot, 4+p}  = sprintf(' %0.1f', ...
        session.convertParameterToDisplayUnits(pMatrixInit(ctot,p),p));
    p = 5; table_data_string{ctot, 4+p}  = sprintf(' %0.1f', ...
        session.convertParameterToDisplayUnits(pMatrixInit(ctot,p),p));    
end
set(handles.table_ics_view_all, 'Data', table_data_string);



function table_isok_group_kinetics_CellEditCallback(hObject, eventdata, handles)
% --- Executes when entered data in editable cell(s) in table_isok_group_kinetics.
% hObject    handle to table_isok_group_kinetics (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% User changed ISOK in kex,PA0,dH, or Eab

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Check which curve and parameter has been changed
data_changed    = eventdata.Indices;
row_changed     = data_changed(1);
rowNames        = get(handles.table_isok_group_kinetics, 'RowName');
pName_changed   = rowNames{row_changed};

% Commit this change to the current fitResult
f           = get(handles.popup_fit_results, 'Value');
fitResult   = group.fitResults{f};
fitResult.setParamIsOK( eventdata.NewData, pName_changed );

refresh_display(handles);


% --- Executes when entered data in editable cell(s) in table_isok_curveset_ppm.
function table_isok_curveset_ppm_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_isok_curveset_ppm (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% User changed ISOK in kex,PA0,dH, or Eab

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Check which curve and parameter has been changed
data_changed    = eventdata.Indices;
row_changed     = data_changed(1);
cs_changed      = row_changed;

col_changed     = data_changed(2);
columnNames     = get(handles.table_isok_curveset_ppm, 'ColumnName');
pName_changed   = columnNames{col_changed};

% Commit this change to the current fitResult
f           = get(handles.popup_fit_results, 'Value');
fitResult   = group.fitResults{f};
fitResult.setParamIsOK( eventdata.NewData, pName_changed, cs_changed );

refresh_display(handles);


% --- Executes when entered data in editable cell(s) in table_isok_curve_r.
function table_isok_curve_r_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_isok_curve_r (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% User changed ISOK in kex,PA0,dH, or Eab

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Check which curve and parameter has been changed
data_changed    = eventdata.Indices;
row_changed     = data_changed(1);
ctot_changed     = row_changed;

col_changed     = data_changed(2);
columnNames     = get(handles.table_isok_curve_r, 'ColumnName');
pName_changed   = columnNames{col_changed};

% Commit this change to the current fitResult
f           = get(handles.popup_fit_results, 'Value');
fitResult   = group.fitResults{f};
fitResult.setParamIsOK( eventdata.NewData, pName_changed, ctot_changed );

refresh_display(handles);

% --- Executes when entered data in editable cell(s) in table_isok_temp_kinetics.
function table_isok_temp_kinetics_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_isok_temp_kinetics (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% User changed ISOK in kex,PA0,dH, or Eab

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Check which curve and parameter has been changed
data_changed    = eventdata.Indices;
row_changed     = data_changed(1);
table_data      = get(handles.table_isok_temp_kinetics, 'Data');
Temp_changed    = table_data{row_changed,1}+273;


col_changed     = data_changed(2);
columnNames     = get(handles.table_isok_temp_kinetics, 'ColumnName');
pName_changed   = columnNames{col_changed};

eventdata.NewData
pName_changed
Temp_changed

% Commit this change to the current fitResult
f           = get(handles.popup_fit_results, 'Value');
fitResult   = group.fitResults{f};
fitResult.setParamIsOK( eventdata.NewData, pName_changed, Temp_changed );

refresh_display(handles);


% --- Executes when entered data in editable cell(s) in table_ics_temp_kinetics.
function table_ics_temp_kinetics_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_ics_temp_kinetics (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
update_table_ics(handles);
refresh_display(handles);


% --- Executes on button press in button_rename_fit.
function button_rename_fit_Callback(hObject, eventdata, handles)
% hObject    handle to button_rename_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Call the load curves GUI using variables from main handle

handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% Get the current fitResult
f           = get(handles.popup_fit_results, 'Value');
fitResult   = group.fitResults{f};

% Pause GUI execution until the data are loaded
h = GUARDD_edit_fit_name('GUARDD', handles_main, 'fitResult', fitResult);
waitfor(h);

refresh_display(handles);

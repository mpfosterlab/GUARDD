% Simulate two-state kinetic rates and populations to study nature of exchange phenomena (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/01/13 Start coding
%

function varargout = GUARDD_kinetic_simulator(varargin)
% GUARDD_KINETIC_SIMULATOR M-file for GUARDD_kinetic_simulator.fig
%      GUARDD_KINETIC_SIMULATOR, by itself, creates a new GUARDD_KINETIC_SIMULATOR or raises the existing
%      singleton*.
%
%      H = GUARDD_KINETIC_SIMULATOR returns the handle to a new GUARDD_KINETIC_SIMULATOR or the handle to
%      the existing singleton*.
%
%      GUARDD_KINETIC_SIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_KINETIC_SIMULATOR.M with the given input arguments.
%
%      GUARDD_KINETIC_SIMULATOR('Property','Value',...) creates a new GUARDD_KINETIC_SIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_kinetic_simulator_OpeningFcn gets
%      called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_kinetic_simulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_kinetic_simulator

% Last Modified by GUIDE v2.5 16-Apr-2011 11:01:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_kinetic_simulator_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_kinetic_simulator_OutputFcn, ...
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


% --- Executes just before GUARDD_kinetic_simulator is made visible.
function GUARDD_kinetic_simulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_kinetic_simulator (see VARARGIN)

% Choose default command line output for GUARDD_kinetic_simulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_kinetic_simulator wait for user response (see UIRESUME)
% uiwait(handles.kinetic_simulator_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.kinetic_simulator_gui, 'handles_main', handles_main);


%% Initialize the table for simulation parameters

column_string = {'Temp|(C)', 'PA|(%)', 'kex|(/s)', 'dH|(kcal/mol)', 'Eab|(kcal/mol)'};
set(handles.table_kinetic_specs, 'RowName', []);
set(handles.table_kinetic_specs, 'ColumnWidth', {40,40,60,70,70});
set(handles.table_kinetic_specs, 'ColumnName', column_string );
set(handles.table_kinetic_specs, 'Data', [25,90,1000,-10, 2]);

set(handles.table_input_temp, 'RowName', {'Temp (C)'});
set(handles.table_input_temp, 'ColumnWidth', {40});
set(handles.table_input_temp, 'ColumnName', []);
set(handles.table_input_temp, 'Data', 25);

%set(handles.table_kinetic_report1, 'RowName', {'PA (%)', 'kex (/sec)', 'kA (/sec)', 'kB (/sec)'});
set(handles.table_kinetic_report1, 'RowName', {'PA', 'kex', 'kA', 'kB'});
set(handles.table_kinetic_report1, 'ColumnWidth', {80});
set(handles.table_kinetic_report1, 'ColumnName', []);
set(handles.table_kinetic_report1, 'Data', []);

set(handles.table_kinetic_report2, 'RowName', {'Eab', 'Pab', 'Eba', 'Pba', 'dH', 'dS'});
set(handles.table_kinetic_report2, 'ColumnWidth', {120});
set(handles.table_kinetic_report2, 'ColumnName', []);
set(handles.table_kinetic_report2, 'Data', []);

% Select the same set as in the RD simulator GUI
simulationSession = getappdata(handles_main.main_gui, 'simulationSession');
set(handles.popup_cs, 'Value', simulationSession.cs_selected);

% TODO % Tooltips


refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_kinetic_simulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% This is called to update all display elements
function refresh_display(handles)
%%
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

if( simulationSession.Ncs > 0 )    
    % Check which set is selected from the list
    Ncs         = simulationSession.Ncs;
    cs_selected = get(handles.popup_cs, 'Value');
    
    if( cs_selected > simulationSession.Ncs )
        cs_selected = 1;
        set(handles.popup_cs, 'Value', cs_selected);
    end
    curveset    = simulationSession.curvesets{cs_selected};

    %% Populate pull-down menu with sets
    set_names = cell(Ncs,1);
    for cs = 1:Ncs
        if( cs == cs_selected )
            prefix      = '* ';
        else
            prefix      = '';
        end

        set_names{cs} = sprintf('%s[%02d] %s', prefix, cs, simulationSession.curvesets{cs}.name);
    end
    
    set(handles.popup_cs, 'String', set_names );
    set(handles.popup_cs, 'Enable', 'on');

    %% Set name in edit box
    set(handles.edit_cs_name, 'String', curveset.name);

    %% Output to tables
    
    table_data = [ curveset.T0-273, ...
                   curveset.PA0*100, ...
                   curveset.kex0, ...
                   curveset.dH/1000, ...
                   curveset.Eab/1000 ...
                 ];
    set(handles.table_kinetic_specs, 'Data', table_data);

    % Get temperature from table
    Tt = get(handles.table_input_temp, 'Data')+273;

    table_data = { sprintf('%0.1f%%',  curveset.calc_PA(Tt)*100   ); ...
                   sprintf('%0.1f /s', curveset.calc_kex(Tt)      ); ...
                   sprintf('%0.1f /s', curveset.calc_kA(Tt)       ); ...
                   sprintf('%0.1f /s', curveset.calc_kB(Tt)       ) };
    set(handles.table_kinetic_report1, 'Data', table_data);

    table_data = { sprintf('%0.1f kcal/mol',    curveset.Eab / 1000); ...
                   sprintf('%0.1e /s',          curveset.Pab); ...
                   sprintf('%0.1f kcal/mol',    curveset.Eba / 1000); ...
                   sprintf('%0.1e /s',          curveset.Pba); ...
                   sprintf('%0.1f kcal/mol',    curveset.dH / 1000); ...
                   sprintf('%0.1f cal/mol/K',   curveset.dS) };
    set(handles.table_kinetic_report2, 'Data', table_data);

    % Plot the simulated kinetic data
    plot_sim( handles.panel_plot_sim, handles);
    
    
else
    % No simulations to show
    
    
end

function plot_sim( figure_handle, handles )

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');
session        = getappdata(handles_main.main_gui, 'session');

% Get temperature from table
Tt = get(handles.table_input_temp, 'Data')+273;
%{
% TODO % Remove hard-coded session
session.FONTSIZE_SMALL         = 12;
session.FONTSIZE_MEDIUM        = 14;
session.FONTSIZE_LARGE         = 16;
session.FONTNAME               = 'Arial';
session.MARKERSIZE         = 8;
session.LINEWIDTH          = 2;
session.R = 1.9858775;            % Gas constant cal/mol/K
%}

%% Goal: plot PA(T) and kex(T) given certain temperature-dependences

% Check which set is selected from the list
cs_selected = get(handles.popup_cs, 'Value');
curveset    = simulationSession.curvesets{cs_selected};

% TODO % Get user to set TMin, TMax, TSteps
T_grid = linspace(273, 328, 12);

h = subplot(1,1,1, 'Parent', figure_handle);
cla(h);
hold(h,'off');

% Plot equilibrium data
h = subplot(2,1,1, 'Parent', figure_handle);
set( h, 'FontName', session.FONTNAME, 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'bold', 'LineWidth', session.LINEWIDTH)
hold(h, 'on');

%PA_grid = simulation.calc_PA(T_grid,cs_selected);
PA_grid =  curveset.calc_PA(T_grid);

plot(h, T_grid-273, PA_grid, '-ok', 'MarkerSize', session.MARKERSIZE, 'LineWidth', session.LINEWIDTH)
ylim([0,1]);

% Add a vertical line marking the selected temperature
plot(h, [Tt Tt]-273, [0 1], '--k', 'LineWidth', session.LINEWIDTH);
ylabel('P_A', 'FontSize', session.FONTSIZE_MEDIUM);

% TODO % Move title to somewhere useful
title(h, sprintf('Kinetic Simulator'), 'FontName', session.FONTNAME, 'FontSize', session.FONTSIZE_LARGE);

% Kinetic plot
h = subplot(2,1,2, 'Parent', figure_handle);
set( h, 'FontName', session.FONTNAME, 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'bold', 'LineWidth', session.LINEWIDTH)

%kex_grid = simulation.calc_kex(T_grid, cs_selected);
%kA_grid  = simulation.calc_kA(T_grid, cs_selected);
%kB_grid  = simulation.calc_kB(T_grid, cs_selected);

kex_grid = curveset.calc_kex(T_grid);
kA_grid  = curveset.calc_kA(T_grid);
kB_grid  = curveset.calc_kB(T_grid);

hold(h, 'on');
plot(h, T_grid-273, kex_grid,'-ok', 'MarkerFaceColor', 'k', 'MarkerSize', session.MARKERSIZE, 'LineWidth', session.LINEWIDTH);
plot(h, T_grid-273, kA_grid, '->k', 'MarkerSize', session.MARKERSIZE, 'LineWidth', session.LINEWIDTH);
plot(h, T_grid-273, kB_grid, '-<r', 'MarkerSize', session.MARKERSIZE, 'LineWidth', session.LINEWIDTH);
ylabel('Exch. rate (/s)', 'FontSize', session.FONTSIZE_MEDIUM);
xlabel('Temperature (^oC)', 'FontSize', session.FONTSIZE_MEDIUM);

% Add a vertical line marking the selected temperature (don't alter default Y-scale)
YLIM = get(h, 'YLim');
plot(h, [Tt Tt]-273, YLIM, '--k', 'LineWidth', session.LINEWIDTH);
set(h,'YLim', YLIM);

hl=legend(h, {'k_{ex}', 'k_A', 'k_B'}, 'Location', 'NorthWest' );
set(hl, 'FontName', session.FONTNAME, 'FontSize', session.FONTSIZE_SMALL, ...
    'FontWeight', 'Normal', 'LineWidth', session.LINEWIDTH);


% --- Executes when entered data in editable cell(s) in table_input_temp.
function table_input_temp_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_input_temp (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
refresh_display(handles);


% --- Executes on selection change in popup_cs.
function popup_cs_Callback(hObject, eventdata, handles)
% hObject    handle to popup_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_cs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_cs
%% Commit the selection
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

simulation.cs_selected = get(handles.popup_cs, 'Value');
setappdata(handles_main.main_gui, 'simulationSession', simulationSession);

refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function popup_cs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function ui_save_figure_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save figure to file

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
session        = getappdata(handles_main.main_gui, 'session');

if( ~isfield(session, 'outputDir') )
    session.outputDir = '/';
end

% Ask user for file to save as
default_filename     = sprintf('simulation.%s', '');
[filename, filepath] = uiputfile('*.fig; *.png; *.ps', 'Save graphic file', ...
    sprintf('./%s/%s', session.outputDir, default_filename) );

% If the user selected a file
if( ~isequal(filename,0) )
    
    % Make new figure window
    h = figure;

    % Plot the data to a new figure    
    plot_sim( h, handles )
    
    % Get file extension from file name if there is one
    if( findstr(filename, '.') < length(filename) )
        file_ext = filename(findstr(filename, '.')+1:end);
    else
        file_ext = '?';
    end
    
    if( strcmp(file_ext, 'fig') )
        hgsave( h, sprintf('%s/%s',filepath,filename) );
    elseif( strcmp(file_ext, 'png') )
        print( h, '-r300', '-dpng', sprintf('%s/%s',filepath,filename) );
    else
        print( h, '-dpsc', sprintf('%s/%s',filepath,filename) );
    end    

    close(h);
end




% --------------------------------------------------------------------
function ui_popout_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_popout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Display plot in new figure window
h = figure;
plot_sim( h, handles )


% --- Executes during object creation, after setting all properties.
function listbox_params_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_new_cs.
function button_new_cs_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Add a new simulation set with default parameters
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

% Add new set
simulationSession.newCurveset();

% Select the new set
simulationSession.setcsSelected(simulationSession.Ncs);

setappdata(handles_main.main_gui, 'simulationSession', simulationSession);
set(handles.popup_cs, 'Value', simulationSession.Ncs);
refresh_display(handles);


% --- Executes on button press in button_delete_cs.
function button_delete_cs_Callback(hObject, eventdata, handles)
% hObject    handle to button_delete_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Remove simulation set from list
% Ask user if they want to delete fit
if( strcmp(questdlg('Are you sure you want to delete this simulation set?', ...
    'Delete set?', 'Delete set', 'Cancel','Cancel'),'Delete set') )

    % Access main window handle which is stored in this window's main handle
    handles_main        = getappdata(handles.kinetic_simulator_gui, 'handles_main');
    simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');    
    cs                  = get(handles.popup_cs, 'Value');
    curveset            = simulationSession.curvesets{cs};

    simulationSession.removeCurveset( curveset );

    % Select a new simulation
    if( cs > 1 )
        set(handles.popup_cs, 'Value', cs-1);
    else
        set(handles.popup_cs, 'Value', 1);
    end

    % Save data
    setappdata(handles_main.main_gui, 'simulationSession', simulationSession);
end

refresh_display(handles);


% --- Executes on button press in button_view_2d.
function button_view_2d_Callback(hObject, eventdata, handles)
% hObject    handle to button_view_2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_view_2d

refresh_display(handles);


% --- Executes on button press in button_view_3d.
function button_view_3d_Callback(hObject, eventdata, handles)
% hObject    handle to button_view_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_view_3d

refresh_display(handles);



% --- Executes on button press in button_save_cs.
function button_save_cs_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save the name of the set
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

cs      = get(handles.popup_cs, 'Value');
name    = get(handles.edit_cs_name, 'String');

simulationSession.curvesets{cs}.setName( name );

% Save it
setappdata(handles_main.main_gui, 'simulationSession', simulationSession);
refresh_display(handles);



% --- Executes on key press with focus on edit_cs_name and none of its controls.
function edit_cs_name_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_cs_name (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%% User typed in set name box, so allow user to save changes
set(handles.button_save_cs, 'Enable', 'on');


% --- Executes on selection change in popup_y_var.
function popup_y_var_Callback(hObject, eventdata, handles)
% hObject    handle to popup_y_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_y_var contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_y_var

%% Select new variable for y-axis

y_axis_string = {'Sim num', 'dwH', 'dwX', 'Pa', 'kex', 'R20'};
y_axis_index = get(handles.popup_y_var,'Value');

if( strcmp('Sim num', y_axis_string{y_axis_index}) )
    Y_LIM = [0,0];
    set(handles.checkbox_surface, 'Value', 0);

elseif( strcmp('dwH', y_axis_string{y_axis_index}) )
    Y_LIM = [0,1000];

elseif( strcmp('dwX', y_axis_string{y_axis_index}) )
    Y_LIM = [0,500];

elseif( strcmp('Pa', y_axis_string{y_axis_index}) )
    Y_LIM = [50,100];

elseif( strcmp('kex', y_axis_string{y_axis_index}) )
    Y_LIM = [0,10000];

elseif( strcmp('R20', y_axis_string{y_axis_index}) )
    Y_LIM = [0,100];
end

% Read in prior Nx and Ny values to maintain them
table_data  = get(handles.table_kinetic_specs, 'Data');

set(handles.table_kinetic_specs, 'Data', [Y_LIM, table_data(3:4)]);


refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function popup_y_var_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_y_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in table_kinetic_specs.
function table_kinetic_specs_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_kinetic_specs (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Read data from tables and calculate kinetic parameters

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

%session        = getappdata(handles_main.main_gui, 'session');
R               = simulationSession.R;           % Gas constant cal/mol/K

% Obtain basic kinetic specifications, which determine remaining two-site exchange parameters
table_data = get(handles.table_kinetic_specs, 'Data');
T0      = table_data(1) + 273;  % Celcius -> Kelvin
PA0     = table_data(2) / 100;  % Percent -> Fraction
kex0    = table_data(3);
dH      = table_data(4) * 1000; % kcal/mol -> cal/mol
Eab     = table_data(5) * 1000; % kcal/mol -> cal/mol

% Check which set is selected from the list
cs          = get(handles.popup_cs, 'Value');
curveset    = simulationSession.curvesets{cs};

% Save to data structure
curveset.setSpecsKinetics( T0, PA0, kex0, Eab, dH );

%{
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

% Arrhenius: Use kB at two temperatures to get Eba
T1 = 280; T2 = 320;
Eba = R * log( kB(T1) / kB(T2) ) / (1/T2 - 1/T1);

% Arrhenius: use Eba and point kB at T0 to get Pba
Pba = kB0 * exp( Eba / (R*T0) );

%[(T-273)' PA' kex' kA' kB']
%[Eab Pab 0 0 dH dS]

% Save to data structure
curveset.setSpecsKinetics( Eab, Pab, Eba, Pba, dH, dS, T0, PA0, kex0 );
%}

setappdata(handles_main.main_gui, 'simulationSession', simulationSession);

% The changed set IS SAVED by default
set(handles.button_save_cs, 'Enable', 'off');

refresh_display(handles);


% --- Executes on button press in checkbox_surface.
function checkbox_surface_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_surface

refresh_display(handles);


% --- Executes when selected object is changed in panel_view_angle.
function panel_view_angle_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_view_angle 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

refresh_display(handles);


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export_menu_Callback(hObject, eventdata, handles)
% hObject    handle to export_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Call export GUI using variables from main handle
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
GUARDD_rd_simulator_export('GUARDD', handles_main);


% --- Executes on button press in button_copy_cs.
function button_copy_cs_Callback(hObject, eventdata, handles)
% hObject    handle to button_copy_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Copy simulation set 

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.kinetic_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

% Select current set
cs          = get(handles.popup_cs, 'Value');
curveset    = simulationSession.curvesets{cs};

% Copy it
curveset_cp = curveset.copy();
simulationSession.addCurveset( curveset_cp );

setappdata(handles_main.main_gui, 'simulationSession', simulationSession);

% Select the new simulation set
set(handles.popup_cs, 'Value', simulationSession.Ncs );

refresh_display(handles);


% --- Executes when kinetic_simulator_gui is resized.
function kinetic_simulator_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to kinetic_simulator_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ui_refresh_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

refresh_display(handles);



function edit_cs_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cs_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cs_name as text
%        str2double(get(hObject,'String')) returns contents of edit_cs_name as a double
set(handles.button_save_cs, 'Enable', 'on');

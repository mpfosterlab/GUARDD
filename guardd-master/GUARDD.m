% Controls essential functions (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% Started January, 2010
% 2011/01/18 Update for classes
% 2017/04/25 Added 19F (IK)
% 2017/05/26 Fix bug loading prior session

function varargout = GUARDD(varargin)
% GUARDD M-file for GUARDD.fig
%      GUARDD, by itself, creates a new GUARDD or raises the existing
%      singleton*.
%
%      H = GUARDD returns the handle to a new GUARDD or the handle to
%      the existing singleton*.
%
%      GUARDD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD.M with the given input arguments.
%
%      GUARDD('Property','Value',...) creates a new GUARDD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD

% Last Modified by GUIDE v2.5 05-Sep-2011 21:04:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_OutputFcn, ...
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

% TODO % Set all font sizes via session.FONTSIZE and session.FONTNAME (axes med, legend small, title large)
% TODO % Backup code used to run


% --- Executes just before GUARDD is made visible.
function GUARDD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD (see VARARGIN)

% Choose default command line output for GUARDD
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD wait for user response (see UIRESUME)
% uiwait(handles.main_gui);

%% Initialize program
% TODO % Render with OpenGL for increased speed ?
echo off all;

clc;

% Start timing the program's execution
tic;

%% Initialize data structures
%variables   = {'data', 'settings', 'fit_results', 'sse_fit', 'group_analysis', 'simulationSession', 'session' };
variables   = {'session', 'simulationSession' };
setappdata(handles.main_gui, 'variables', variables);

session = Session;
setappdata(handles.main_gui, 'session', session);

VERSION = 'v.2017.04.25';
set(handles.text_version, 'String', sprintf('%s',VERSION));
setappdata(handles.main_gui, 'VERSION', VERSION);

% Set contents of listbox for all curve sets
set(handles.listbox_g, 'String', '(Load data)');

update_gui_status(handles);

% Erase all of the variable values
clear all;


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% This is called to update all display elements for selecting current g
% g        The current group selected for display
function refresh_display(handles)

%% Get the data structure
%settings    = getappdata(handles.main_gui, 'settings');
session     = getappdata(handles.main_gui, 'session');


%% Nothing to do if data is empty
if( isempty(session) || session.Ng == 0 )
    set(handles.listbox_g, 'String', '(No groups)');
    set(handles.listbox_g, 'Value', 1);
    set(handles.checkbox_exhibits_exchange, 'Enable', 'off');
    
    % Clear old contents and create new axes (instead of subplot(1,1,1))
    delete( get(handles.panel_dispersion, 'Children') );
    %h = axes('Parent', handles.panel_dispersion);
    %hold(h, 'all');
    %h = subplot(1,1,1,'Parent', handles.panel_dispersion );    
    %cla(h);
    return
end

%% Fill list box with AA names formatted for FIT/EXCHANGE status
listbox_string  = cell(1, session.Ng);
for g = 1:session.Ng
    group = session.groups{g};
    
    % Format string for exchange status
    if( group.exhibitsExchange )
        marker_string = '[Ex] ';
    else
        marker_string = '[--] ';
    end

    % Format string for fit status
    %  Fit available and OK         Bold + Black
    %  Fit available and not OK     Bold + Gray
    %  No fits available            Italics

    % Fit available and OK
    if( group.Nf > 0 && group.bestFitIsOK )
        listbox_string{g} = sprintf('<html><b>%s%s</b></html>', ...
            marker_string, group.name );

    % Fit available but NOT OK
    elseif( group.Nf > 0 && ~group.bestFitIsOK )
        listbox_string{g} = sprintf('<html><font color="Gray"><b>%s%s</b></font></html>', ...
            marker_string, group.name );

    % Fit is not available
    else
        listbox_string{g} = sprintf('<html><i>%s%s</i></html>', ...
            marker_string, group.name );
    end
end

% Set contents of listbox for all curve sets
set(handles.listbox_g, 'String', listbox_string);

% Update list box
if( get(handles.listbox_g, 'Value') > session.Ng )    
    set(handles.listbox_g, 'Value', 1);
end

%% Plot the data in the panel for selected group
%h = subplot(1,1,1, 'Parent', handles.panel_dispersion);
%cla(h);
%hold(h, 'off');

% Clear old contents and create new axes (instead of subplot(1,1,1))
delete( get(handles.panel_dispersion, 'Children') );
h = axes('Parent', handles.panel_dispersion);
hold(h, 'all');

% Get the selected group
g       = get(handles.listbox_g, 'Value');
group   = session.groups{g};

group.plotCurves(h, 'ok', 'MarkerSize', session.MARKERSIZE, 'LineWidth', session.LINEWIDTH);
%group.plotCurves(h);

[vcpmg_min, vcpmg_max, R2eff_min, R2eff_max] = session.getDataLimits();
set(h, 'XLim', [0 100*ceil(vcpmg_max*1.1/100)] );
set(h, 'YLim', [0 50*ceil(R2eff_max/50)] );

% Set title and axes for CPMG data
%set(h,'FontName', session.FONTNAME, 'FontSize', session.FONTSIZE_LARGE, 'LineWidth', session.LINEWIDTH)
set(h, 'FontSize', session.FONTSIZE_MEDIUM, 'LineWidth', session.LINEWIDTH)
title(h, group.name, 'FontSize', session.FONTSIZE_LARGE, 'FontWeight', 'Bold');
xlabel(h, '\nu_{CPMG} (Hz)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
ylabel(h, 'R_2^{eff} (Hz)',  'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold')

set(h, 'Box', 'on');
set(h, 'YGrid','on');

% Exhibits exchange?
set(handles.checkbox_exhibits_exchange, 'Enable', 'on');
set(handles.checkbox_exhibits_exchange, 'Value', group.exhibitsExchange)


%% Update other displays
%{
% If the user wants crosspeak notes
if( get(handles.checkbox_notes, 'Value') )
    handles.notes_gui = GUARDD_notes('GUARDD',handles);
end
%}

% If the user wants 2D fitting plot
if( get(handles.checkbox_fit_rd, 'Value') )
    handles.fit_rd_gui = GUARDD_fit_rd('GUARDD',handles);
end


% If the user wants 3D fitting plot
if( get(handles.checkbox_display_rd, 'Value') )
    handles.display_rd_gui = GUARDD_display_rd('GUARDD',handles);
end


% SSE scatter display from grid search
if( get(handles.checkbox_display_chi2_map, 'Value') )
    handles.display_chi2_map_gui = GUARDD_display_chi2_map('GUARDD',handles);
end

% Temperature dependent analysis of current fit
if( get(handles.checkbox_display_rates, 'Value') )
    handles.display_rates_gui = GUARDD_display_rates('GUARDD',handles);                                
end

% Save this window handle for future use
% I don't think this works like this (2010/01/18)
%setappdata(handles.main_gui, 'fit_rd_gui', handles.fit_rd_gui);





% --------------------------------------------------------------------
function menu_input_Callback(hObject, eventdata, handles)
% hObject    handle to menu_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
% File...Load config file...
function menu_file_new_session_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_new_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Erase data and start over
if( strcmp(questdlg('Are you sure you want to clear this session?', ...
        'Clear?', 'Clear', 'Cancel', 'Cancel'),'Clear') )
    
    session = Session;
    setappdata(handles.main_gui, 'session', session);

    simulationSession = SimulationSession;
    setappdata(handles.main_gui, 'simulationSession', simulationSession);

    % Set contents of listbox for all curve sets
    set(handles.listbox_g, 'String', '(Load data)');

    update_gui_status(handles);
end


% --------------------------------------------------------------------
% File...Save config file...
function menu_file_save_session_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save variables to Matlab .mat file

% TODO % MATLAB .MAT file is SLOW to load (3-5 min), and v.7.3 is LARGE 170 Mb
% TODO % How can this be addressed??

variables   = getappdata(handles.main_gui, 'variables');
session    = getappdata(handles.main_gui, 'session');

if( ~exist(session.outputDir, 'dir') )
    mkdir(session.outputDir);
end

default_filepath    = session.outputDir;

% Ask user for file to save as
default_filename    = 'GUARDD-Session.mat';

[filename, filepath] = uiputfile('*.mat', 'Save session', ...
    sprintf('./%s/%s', default_filepath, default_filename) );

if( ~isequal(filename,0) )
    fprintf('\nSaving (1-30 sec)');
    
    % Start timing
    tic0 = tic();
    
    % If outupt file exists, delete it and create a new one
    if( exist(sprintf('%s/%s', filepath,filename), 'file') )
        delete(sprintf('%s/%s', filepath,filename));
    end
        
    % Check each variable and save
    for v = 1:length(variables)
        fprintf('\n\t%s\t', variables{v} );
        if(length(variables{v})<8)
            fprintf('\t');
        end

        % If variable exists
        if( isappdata(handles.main_gui, variables{v}) && ~isempty(variables{v}) )
            
            eval(sprintf('%s = getappdata(handles.main_gui, variables{v});',variables{v}));

            % First variable creates a new file
            if( ~exist(sprintf('%s/%s', filepath,filename), 'file') )
                %save( sprintf('%s/%s', filepath, filename), variables{v}, '-v7.3' );
                save( sprintf('%s/%s', filepath, filename), variables{v} );

            % Subsequent variables append the existing file
            else                
                %save( sprintf('%s/%s', filepath, filename), variables{v}, '-append', '-v7.3' );                
                save( sprintf('%s/%s', filepath, filename), variables{v}, '-append' );
            end

            fprintf('Saved');
                
        else
            fprintf('\tNo info to save');
        end
    end
    fprintf('\n\tWrote "%s" (%0.1f sec)\n', filename, toc(tic0));
end


%{
% --------------------------------------------------------------------
% File...Quit
function menu_file_quit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Exit program
if( strcmp(questdlg('Are you sure you want to exit GUARDD?', ...
        'Exit?', 'Exit', 'Cancel', 'Cancel'),'Exit') )
    
    % TO DO: Perform a test and close every associated window
    if( isfield(handles, 'fit_rd_gui') )
        close(handles.fit_rd_gui);
    end

    if( isfield(handles, 'display_chi2_map_gui') )
        close(handles.display_chi2_map_gui);
    end


    close(handles.main_gui);
    
end
%}

% --- Executes on selection change in listbox_g.
function listbox_g_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_g contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_g

%% Update all displays for the currently selected curve set
refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function listbox_g_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

%% Listbox creation
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in button_next_g.
function button_next_g_Callback(hObject, eventdata, handles)
% hObject    handle to button_next_g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Next group number
% Get the current group number
g = get(handles.listbox_g, 'Value');

% If its not at the end of this list, then increment
session = getappdata(handles.main_gui, 'session');
if( g < session.Ng )
    set(handles.listbox_g, 'Value',g+1);
    refresh_display(handles);
end



% --- Executes on button press in button_prev_g.
function button_prev_g_Callback(hObject, eventdata, handles)
% hObject    handle to button_prev_g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Next group number
% Get the current group number
g = get(handles.listbox_g, 'Value');

% If its not at the end of this list, then increment
session = getappdata(handles.main_gui, 'session');
if( g > 1 )
    set(handles.listbox_g, 'Value',g-1);
    refresh_display(handles);
end

%{
% --- Executes on button press in button_update_displays.
function button_update_displays_Callback(hObject, eventdata, handles)
% hObject    handle to button_update_displays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Update all desired plots
refresh_display(handles);
%}

% --- Executes on button press in checkbox_fit_rd.
function checkbox_fit_rd_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fit_rd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fit_rd



% --------------------------------------------------------------------
function menu_file_load_session_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_load_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Load variables

% Ask user to select file
[filename,pathname] = uigetfile('*.mat', 'Load Matlab variables file');

if( ~isequal(filename,0) )
    fprintf('\nLoading (30-600 sec)');
    % Start timing
    tic0 = tic();
    
    % Get list of variable names from main GUI
    variables           = getappdata(handles.main_gui, 'variables');
    %file_variables      = whos('-file', sprintf('%s/%s', pathname, filename));
    
    % Clear contents of current variables
    for v = 1:length( variables )
        
        fprintf('\n\t%s', variables{v});
        
        % If variable exists, clear it
        if( isappdata(handles.main_gui, variables{v}) )
            try
                eval( sprintf('rmappdata(handles.main_gui, \''%s\'');', variables{v}) );
                eval( sprintf('clear %s;', variables{v}) );
                fprintf('\tCleared');
            catch
                fprintf('\tNothing to clear');
            end
        end
        
        % Check if the variable is in the file
        %for fv = 1:length(file_variables)
            
            % Only try to load the variable if its in the file
            %if( strcmp(file_variables(fv).name, variables{v}) )
                % Load variable from file
                try
                    load(sprintf('%s/%s', pathname, filename), variables{v});
                    
                    % Set variable as application data
                    eval( sprintf('setappdata(handles.main_gui, \''%s\'', %s);', variables{v}, variables{v}) );
                    fprintf('\tLoaded');
                catch
                    fprintf('\tNothing to load');
                end
            %end            
        %end
    end
    fprintf('\n\tLoaded "%s" (%0.1f sec)\n', filename, toc(tic0));
end

update_gui_status(handles);

function update_gui_status(handles)
%% This will check whether data structures exist and then update GUI accordingly

% Get data from main handle
%data    = getappdata(handles.main_gui, 'data');
session = getappdata(handles.main_gui, 'session');

MODE        = 'off';
MENU_MODE   = 'off';

if( ~isempty(session) &&  session.Ng > 0 )    
    % Set contents of listbox for all groups
    MODE = 'on';
    MENU_MODE = 'on';
    
    g = get(handles.listbox_g, 'Value');
    if( g > session.Ng )
        set(handles.listbox_g, 'Value', 1);
    end
end 

% Set GUI elements now that data are loaded
set(handles.button_prev_g,                     'Enable', MODE);
set(handles.button_next_g,                     'Enable', MODE);
set(handles.listbox_g,                         'Enable', MODE);

%set(handles.button_update_displays,             'Enable', MODE);
%set(handles.checkbox_notes,                     'Enable', MODE);
set(handles.checkbox_fit_rd,                    'Enable', MODE);
set(handles.checkbox_display_rd,                'Enable', MODE);
set(handles.checkbox_display_chi2_map,          'Enable', MODE);
set(handles.checkbox_display_rates,             'Enable', MODE);

set(handles.menu_analysis,                      'Enable', MENU_MODE);
set(handles.menu_output,                        'Enable', MENU_MODE);

refresh_display(handles);


% --- Executes on button press in checkbox_display_chi2_map.
function checkbox_display_chi2_map_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display_chi2_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display_chi2_map


% --- Executes on button press in checkbox_display_rates.
function checkbox_display_rates_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display_rates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display_rates


% --- Executes on button press in checkbox_display_rd.
function checkbox_display_rd_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display_rd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display_rd


% --------------------------------------------------------------------
function menu_output_Callback(hObject, eventdata, handles)
% hObject    handle to menu_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tools_batch_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Batch set up window
handles.batch_gui= GUARDD_batch('GUARDD',handles);

%{
% --- Executes on button press in checkbox_notes.
function checkbox_notes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_notes
%}

% --------------------------------------------------------------------
function menu_tools_group_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Group analysis window
handles.groups_gui = GUARDD_groups('GUARDD',handles);


% --------------------------------------------------------------------
function menu_tools_seq_map_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools_seq_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Sequence map figure window
handles.seq_map_gui = GUARDD_param_exam('GUARDD',handles);



% --------------------------------------------------------------------
function menu_run_code_Callback(hObject, eventdata, handles)
% hObject    handle to menu_run_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Check each variable
variables = getappdata(handles.main_gui, 'variables');

% Initialize each variable to zero and set as app data
for v = 1:length(variables)    
    eval(sprintf('%s = getappdata(handles.main_gui, \''%s\'')', ...
        variables{v}, variables{v}) );    
end

clc

%% 2011/09/05 Create new ParamDisplay() for session
%session.paramDisplay = ParamDisplay(session);

%% 2011/09/02 Perform analyzeMe() on all fitResult instances
% Re-compute the resultsMatrix for PhiexX
%{
for g = 1:session.Ng
%for g = 1:1
    fprintf('\nWorking on group %d/%d', g, session.Ng);
    group = session.groups{g};
    for f = 1:group.Nf
        fitResult = group.fitResults{f};        
        fitResult.analyzeMe();
    end
    
    for fG = 1:length(group.fitResults_Grid)
        fprintf('\n\tWorking on grid %d/%d', fG, length(group.fitResults_Grid));
        fitResult = group.fitResults_Grid{fG};
        fitResult.analyzeMe();
    end
    
    if( ~isempty(group.fitResult_NoEx) )
        group.fitResult_NoEx.analyzeMe();
    end
    
    if( ~isempty(group.fitResult_Best) )
        group.fitResult_Best.analyzeMe();    
    end
end
%}

%% 2011 / 05 / 26 - Enforce minimum errors in data
%{
for ds = 1:session.Nds
    fprintf('\nWorking on dataset %d, %s', ds, session.datasets{ds}.name);
    session.datasets{ds}.enforceMinimumError(session.MIN_F_ERROR);    
end
%}


%% 2011 / 06 / 01 - PErform rate analysis on all fits
%{
% Convert kcal -> cal
for g = 1:session.Ng
    fprintf('\nWorking on group %d/%d', g, session.Ng);
    group = session.groups{g};
    for f = 1:group.Nf
        fitResult = group.fitResults{f};        
        fitResult.rateAnalysis.analyzeMe();
    end
    
    for fG = 1:length(group.fitResults_Grid)
        fprintf('\n\tWorking on grid %d/%d', fG, length(group.fitResults_Grid));
        fitResult = group.fitResults_Grid{fG};
        fitResult.rateAnalysis.analyzeMe();
    end
    
    if( ~isempty(group.fitResult_NoEx) )
        group.fitResult_NoEx.rateAnalysis.analyzeMe();
    end
    
    if( ~isempty(group.fitResult_Best) )
        group.fitResult_Best.rateAnalysis.analyzeMe();    
    end
end
%}


%% Correct stereo-assignments (2011 / 05 / 22)
% Put this in GUARDD.m -> check_variables()
%{
OLD_LEU_DELTA_1     = '\delta_1';
NEW_LEU_DELTA_1     = '\delta_2';

OLD_LEU_DELTA_2     = '\delta_2';
NEW_LEU_DELTA_2     = '\delta_1';

OLD_VAL_GAMMA_1     = '\gamma_1';
NEW_VAL_GAMMA_1     = '\gamma_2';

OLD_VAL_GAMMA_2     = '\gamma_2';
NEW_VAL_GAMMA_2     = '\gamma_1';

fprintf('\n\nWORKING ON DATASETS');

% Rename all curves in datasets
for ds = 1:session.Nds
    dataset = session.datasets{ds};
    fprintf('\n\n(%d / %d) Working on %s', ds, session.Nds, dataset.name);
    
    for c = 1:dataset.Nc
        curve = dataset.curves{c};
        old_atom = curve.atom;
        
        switch upper( curve.residue )
            case 'LEU'
                switch( old_atom )
                    case OLD_LEU_DELTA_1
                        new_atom = NEW_LEU_DELTA_1;
                    case OLD_LEU_DELTA_2
                        new_atom = NEW_LEU_DELTA_2;
                    otherwise
                        new_atom = '???';
                end
            case 'VAL'
                switch( old_atom )
                    case OLD_VAL_GAMMA_1
                        new_atom = NEW_VAL_GAMMA_1;
                    case OLD_VAL_GAMMA_2
                        new_atom = NEW_VAL_GAMMA_2;
                    otherwise
                        new_atom = '???';
                end
                
            case 'ILE'
                % Nothing to do
                new_atom = curve.atom;
                    
            otherwise
        end
        
        fprintf('\n\t(%d / %d) Working on %s', c, dataset.Nc, curve.name);
        curve.setAssignment( curve.index, new_atom, curve.residue )        
        fprintf(' -> %s', curve.name);
    end
end


fprintf('\n\nWORKING ON GROUPS');
for g = 1:session.Ng
    group = session.groups{g};
    
    % Find the residue via name
    group_residue = group.name(1:3);
    old_name = group.name;
    
    switch upper( group_residue )
            case 'LEU'
                if( ~isempty(strfind(old_name, OLD_LEU_DELTA_1)) )
                    old_atom = OLD_LEU_DELTA_1;
                    new_atom = NEW_LEU_DELTA_1;
                    
                elseif( ~isempty(strfind(old_name, OLD_LEU_DELTA_2)) )
                    old_atom = OLD_LEU_DELTA_2;
                    new_atom = NEW_LEU_DELTA_2;
                else
                    '??';
                end
                
            case 'VAL'
                if( ~isempty(strfind(old_name, OLD_VAL_GAMMA_1)) )
                    old_atom = OLD_VAL_GAMMA_1;
                    new_atom = NEW_VAL_GAMMA_1;
                    
                elseif( ~isempty(strfind(old_name, OLD_VAL_GAMMA_2)) )
                    old_atom = OLD_VAL_GAMMA_2;
                    new_atom = NEW_VAL_GAMMA_2;
                else
                    '??';
                end
                
            case 'ILE'
                % Nothing to do
                old_atom = 'X';
                new_atom = old_atom;
                    
            otherwise
    end
    
    fprintf('\n\n(%d / %d) Working on %s', g, session.Ng, group.name);
   
    % Sometimes the name is formatted differently than default
    new_name = strrep(old_name, old_atom, new_atom);
    group.setName(new_name);
    fprintf(' -> %s', group.name);
        
    
    for cs = 1:group.Ncs
        curveset = group.curvesets{cs};        
        old_name = curveset.name;
        old_atom = curveset.atom;
        
        switch upper( curveset.residue )
            case 'LEU'
                switch( old_atom )
                    case OLD_LEU_DELTA_1
                        new_atom = NEW_LEU_DELTA_1;
                        
                    case OLD_LEU_DELTA_2
                        new_atom = NEW_LEU_DELTA_2;
                        
                    otherwise
                        new_atom = '???';
                end
            case 'VAL'
                switch( old_atom )
                    case OLD_VAL_GAMMA_1
                        new_atom = NEW_VAL_GAMMA_1;
                        
                    case OLD_VAL_GAMMA_2
                        new_atom = NEW_VAL_GAMMA_2;
                        
                    otherwise
                        new_atom = '???';
                end
                
            case 'ILE'
                % Nothing to do
                new_atom = curveset.atom;
                    
            otherwise
        end
        
        fprintf('\n\t(%d / %d) Working on %s', cs, group.Ncs, curveset.name);
        
        % Do NOT change curves within, because those have already been
        % changed above
        curveset.setAtom_DoNotChangeCurves( new_atom );
        
        % Sometimes the name is formatted differently than default
        new_name = strrep(old_name, old_atom, new_atom);
        curveset.setName(new_name);
        fprintf(' -> %s', curveset.name);
    end
end
%}


% --------------------------------------------------------------------
function menu_tools_param_table_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools_param_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Run results_table GUI
handles.param_table_gui = GUARDD_param_table(...
    'GUARDD',handles);


% --------------------------------------------------------------------
function menu_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to menu_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_data_table_Callback(hObject, eventdata, handles)
% hObject    handle to menu_data_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_dataset_notes_Callback(hObject, eventdata, handles)
% hObject    handle to menu_dataset_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% If the user wants dataset notes
handles.notes_dataset_gui = ...
    GUARDD_notes_dataset('GUARDD',handles);


% --------------------------------------------------------------------
function menu_simulator_Callback(hObject, eventdata, handles)
% hObject    handle to menu_simulator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.rd_simulator_gui = ...
    GUARDD_rd_simulator('GUARDD',handles);


% --------------------------------------------------------------------
function menu_data_manager_Callback(hObject, eventdata, handles)
% hObject    handle to menu_data_manager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%{
handles.manager_gui = ...
    GUARDD_manager('GUARDD',handles);
%}

handles.data_manager_gui = ...
    GUARDD_data_manager('GUARDD',handles);



% --------------------------------------------------------------------
function menu_backup_data_Callback(hObject, eventdata, handles)
% hObject    handle to menu_backup_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Something
% TODO % Finish backup data


% --- Executes on button press in button_manual_refresh.
function button_manual_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to button_manual_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_gui_status(handles);


% --- Executes on button press in checkbox_display_rates.
function checkbox_display_rate_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display_rates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display_rates

% --------------------------------------------------------------------
function menu_display_rd_Callback(hObject, eventdata, handles)
% hObject    handle to menu_display_rd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.display_rd_gui = GUARDD_display_rd('GUARDD',handles);

% --------------------------------------------------------------------
function menu_display_chi2_map_Callback(hObject, eventdata, handles)
% hObject    handle to menu_display_chi2_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.display_chi2_map_gui = GUARDD_display_chi2_map('GUARDD',handles);

% --------------------------------------------------------------------
function menu_display_rates_Callback(hObject, eventdata, handles)
% hObject    handle to menu_display_rates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.display_rates_gui = GUARDD_display_rates('GUARDD',handles);

% --------------------------------------------------------------------
function menu_fitting_window_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fit_rd_gui = GUARDD_fit_rd('GUARDD',handles);


% --- Executes on button press in checkbox_exhibits_exchange.
function checkbox_exhibits_exchange_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_exhibits_exchange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_exhibits_exchange
%% Update curve set as exhibiting exchange or not

% Access main window handle which is stored in this window's main handle
%handles_main    = getappdata(handles.fit_rd_gui, 'handles_main');
session         = getappdata(handles.main_gui, 'session');

% Get the selected group
g       = get(handles.listbox_g, 'Value');
group   = session.groups{g};

group.setExhibitsExchange( get(handles.checkbox_exhibits_exchange, 'Value')==1 );
refresh_display(handles);


% --------------------------------------------------------------------
function menu_settings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data_manager_gui = GUARDD_settings('GUARDD',handles);


% --- Executes when user attempts to close main_gui.
function main_gui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( strcmp(questdlg('Are you sure you want to exit GUARDD?', ...
        'Exit?', 'Exit', 'Cancel', 'Cancel'),'Exit') )

    % Hint: delete(hObject) closes the figure
    delete(hObject);
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
VERSION = getappdata(handles.main_gui, 'VERSION');
GUARDD_about('VERSION', VERSION);


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

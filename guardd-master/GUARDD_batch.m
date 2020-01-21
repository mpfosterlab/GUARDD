% Set-up and execute batch tasks for grid search and error analysis (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/02/01 Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/04/28 Update for classes
% 2011/05/12 Update for modifying grid limits


function varargout = GUARDD_batch(varargin)
% GUARDD_BATCH M-file for GUARDD_batch.fig
%      GUARDD_BATCH, by itself, creates a new GUARDD_BATCH or raises the existing
%      singleton*.
%
%      H = GUARDD_BATCH returns the handle to a new GUARDD_BATCH or the handle to
%      the existing singleton*.
%
%      GUARDD_BATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_BATCH.M with the given input arguments.
%
%      GUARDD_BATCH('Property','Value',...) creates a new GUARDD_BATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_batch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_batch_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_batch

% Last Modified by GUIDE v2.5 07-Jun-2011 10:39:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_batch_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_batch_OutputFcn, ...
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


% --- Executes just before GUARDD_batch is made visible.
function GUARDD_batch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_batch (see VARARGIN)

% Choose default command line output for GUARDD_batch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_batch wait for user response (see UIRESUME)
% uiwait(handles.batch_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
% AND same for fit_r2eff window which called this GUI
setappdata(handles.batch_gui, 'handles_main', handles_main);

% Initialize some variables for the batch grid
g_selected = 1;
setappdata(handles.batch_gui, 'g_selected', g_selected);

%% Initialize progress bar axis
% Update total progress wait bar (Johannes Korsawe, 2007)
value = 0;
h=handles.axes_total_progress;
cla(h);        
patch([0,value,value,0],[0,0,1,1],'k', 'Parent', h);
axis(h,[0,1,0,1]); axis(h,'off'); drawnow;

h=handles.axes_current_progress;
cla(h);        
patch([0,value,value,0],[0,0,1,1],'k', 'Parent', h);
axis(h,[0,1,0,1]); axis(h,'off'); drawnow;

%% Table
set(handles.table_batch, 'ColumnFormat', {'char', 'numeric', 'logical'});
set(handles.table_batch, 'ColumnEditable', [false, false, false]);
set(handles.table_batch, 'ColumnWidth', {350, 70, 40});
set(handles.table_batch, 'ColumnName', {'Name' 'Steps' 'In?'} );

table_data{1,1} = 'Test';
table_data{1,2} = 0;
table_data{1,3} = true;
set(handles.table_batch, 'Data', table_data);

set(handles.table_Temp0_Grid, 'RowName', 'Temp0 (C)');
set(handles.table_Temp0_Grid, 'ColumnName', []);

set(handles.table_sec_per_step, 'ColumnName', []);
set(handles.table_sec_per_step, 'RowName', 'sec/step');
set(handles.table_sec_per_step, 'ColumnWidth', {30});
set(handles.table_sec_per_step, 'Data', 30);

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
set(handles.table_Temp0_Grid, 'Data', 0);


refresh_display(handles)

% Update elements in display
function refresh_display(handles)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.batch_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Functions to perform
PERFORM_GRID        = get(handles.checkbox_grid, 'Value')==1;
CALCULATE_ERORRS    = get(handles.checkbox_errors, 'Value')==1;

% Does the user want to constrain the rate analysis (use dH and Eab)
CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;

if( ~isempty(session) && session.Ng > 0 )    
    %% GUI functions
    set(handles.checkbox_grid, 'Enable', 'on');
    set(handles.checkbox_errors, 'Enable', 'on');
    
    if( PERFORM_GRID )
        set(handles.panel_grid_options, 'Visible', 'on');
        set(handles.checkbox_constrain_rates, 'Visible', 'on');        
    else
        set(handles.panel_grid_options, 'Visible', 'off');
        set(handles.checkbox_constrain_rates, 'Visible', 'off');
    end
    
    if( CALCULATE_ERORRS )
        set(handles.panel_error_options, 'Visible', 'on');
    else
        set(handles.panel_error_options, 'Visible', 'off');
    end

    %% Populate table with groups    
    
    % Enable GUI elements
    set(handles.button_all, 'Enable', 'on');
    set(handles.button_none, 'Enable', 'on');
    set(handles.button_invert, 'Enable', 'on');
    set(handles.button_go, 'Enable', 'on');    
    %set(handles.panel_progress, 'Visible', 'on');
    
    table_data_prior  = get(handles.table_batch, 'Data');    
    Nrows_prior       = size(table_data_prior,1);
    
    % Fill out data table in roundabout way to satisfy formatting of cell array :-/
    table_data = cell(session.Ng, 3);
    for g = 1:session.Ng
        group = session.groups{g};
        
        % Format name to indiate if grid search has been done yet
        if( group.gridSearchIsDone() )            
            in_cluster_string = '[G] ';
        else            
            in_cluster_string = '[-] ';            
        end       
        
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
            table_data{g,1} = sprintf('<html><b>%s%s%s</b></html>', ...
                in_cluster_string, marker_string, group.name );

        % Fit available but NOT OK
        elseif( group.Nf > 0 && ~group.bestFitIsOK )
            table_data{g,1} = sprintf('<html><font color="Gray"><b>%s%s%s</b></font></html>', ...
                in_cluster_string, marker_string, group.name );

        % Fit is not available
        else
            table_data{g,1} = sprintf('<html><i>%s%s%s</i></html>', ...
                in_cluster_string, marker_string, group.name );
        end
        
        
        % Number of grid steps
        Nsteps = 0;
        if( PERFORM_GRID )            
            if( CONSTRAIN_RATE_ANALYSIS )
                Np_Grid = 6;
            else
                Np_Grid = 4;
            end
            
            Nsteps = Nsteps + prod(group.grid_p_steps(1:Np_Grid));
        end
        
        if( CALCULATE_ERORRS )
            Nsteps = Nsteps + session.Nmc;            
        end
        
        table_data{g,2} = Nsteps;
        
        % New rows added        
        if( Nrows_prior < session.Ng )
            %table_data{g,3} = '   X';
            table_data{g,3} = true;
        else
            % Otherwise, preserve what was in the table
            table_data{g,3} = table_data_prior{g,3};
        end
    end       
    set(handles.table_batch, 'Enable', 'on');
    set(handles.table_batch, 'Data', table_data);
    
    %% Show grid search limits
    g_selected  = getappdata(handles.batch_gui, 'g_selected');
    group       = session.groups{g_selected};
    
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
    %{
    
    % Calibration temperature table
    table_data = group.T0_Grid-273;
    set(handles.table_Temp0_Grid, 'Data', table_data);
    
    table_data = [  group.grid_p_min(1:6)  ; ...
                    group.grid_p_max(1:6)  ; ...
                    group.grid_p_steps(1:6) ];                
    set(handles.table_grid_limits, 'Data', table_data);
    %}
    
    set(handles.panel_grid_search, 'Title', sprintf('Grid search: %s', group.name));

    % Calibration temperature table
    set(handles.table_Temp0_Grid, 'Data', group.T0_Grid-273);
    
    %% Estimated total time    
    table_data      = get(handles.table_batch, 'Data');
    searchGroup     = zeros(session.Ng,1);
    Npoints_Group   = zeros(session.Ng,1);
    for g = 1:session.Ng
        % (Boolean) The group should be searched
        %searchGroup(g)      = ~isempty(table_data{g,3});
        searchGroup(g)      = table_data{g,3} == true;

        % Number of grid points in each curve set
        Npoints_Group(g)   = table_data{g,2};
    end

    % Number of groups to search
    NgSearch                    = sum(searchGroup);
    NpointsInTotal              = sum( Npoints_Group .* searchGroup );
    
    set(handles.text_total_steps, 'String', sprintf('Steps: %d', NpointsInTotal));
    
    table_data      = get(handles.table_sec_per_step, 'Data');
    sec_per_step    = table_data(1);
    
    [hour, minute, second] = sec2hms( sec_per_step * NpointsInTotal );
    set(handles.text_estimated_time, 'String', sprintf('%0.0f hr, %0.0f min, %0.0f sec', ...
        hour, minute, second));
    
    

% If there is no data structure to display yet
else
    % Disable GUI elements
    set(handles.button_all, 'Enable', 'off');
    set(handles.button_none, 'Enable', 'off');
    set(handles.button_invert, 'Enable', 'off');
    set(handles.button_go, 'Enable', 'off');
    set(handles.table_batch, 'Enable', 'off');
    set(handles.table_batch, 'Data', []);    
    
    set(handles.checkbox_grid, 'Enable', 'off');
    set(handles.checkbox_errors, 'Enable', 'off');
    
    set(handles.panel_error_options, 'Visible', 'off');
    set(handles.panel_grid_options, 'Visible', 'off');
        
    %set(handles.panel_progress, 'Visible', 'off');
end


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_batch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_all.
function button_all_Callback(hObject, eventdata, handles)
% hObject    handle to button_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Set each curve set to true for batch run
table_data = get(handles.table_batch, 'Data');
Nrows = size(table_data,1);
for r = 1:Nrows
    %table_data{r,3} = '   X';
    table_data{r,3} = true;
end

set(handles.table_batch, 'Data', table_data);

%refresh_display(handles);



% --- Executes on button press in button_none.
function button_none_Callback(hObject, eventdata, handles)
% hObject    handle to button_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Set each curve set to false for batch run
table_data = get(handles.table_batch, 'Data');
Nrows = size(table_data,1);
for r = 1:Nrows
    %table_data{r,3} = '';
    table_data{r,3} = false;
end

set(handles.table_batch, 'Data', table_data);

% --- Executes on button press in button_invert.
function button_invert_Callback(hObject, eventdata, handles)
% hObject    handle to button_invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Set each curve set to false for batch run
table_data = get(handles.table_batch, 'Data');
Nrows = size(table_data,1);
for r = 1:Nrows
    %{
    prior_value = table_data{r,3};        
    if( isempty(prior_value) )
        table_data{r, 3} = '   X';
    else
        table_data{r, 3} = '';
    end
    %}
    table_data{r, 3} = ~table_data{r, 3};
end

set(handles.table_batch, 'Data', table_data);

% --- Executes on button press in button_go.
function button_go_Callback(hObject, eventdata, handles)
% hObject    handle to button_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Start the grid search on each desired curve set
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.batch_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');
variables       = getappdata(handles_main.main_gui, 'variables');


table_data      = get(handles.table_batch, 'Data');
searchGroup     = zeros(session.Ng,1);
Npoints_Group   = zeros(session.Ng,1);
for g = 1:session.Ng
    % (Boolean) The group should be searched
    %searchGroup(g)      = ~isempty(table_data{g,3});
    searchGroup(g)      = table_data{g,3} == true;
    
    % Number of grid points in each curve set
    Npoints_Group(g)   = table_data{g,2};
end

% Number of groups to search
NgSearch                    = sum(searchGroup);
NpointsInTotal              = sum( Npoints_Group .* searchGroup );

% Running number of total points completed
Npoints_Completed           = 0;
Npoints_Completed_per_sec   = 0;
Ng_Completed                = 0;

% Update total progress wait bar (Johannes Korsawe, 2007)
value = Npoints_Completed / NpointsInTotal;
h=handles.axes_total_progress;
cla(h);        
patch([0,value,value,0],[0,0,1,1],'k', 'Parent', h);
axis(h,[0,1,0,1]); axis(h,'off'); drawnow;

% Start timing execution of batch grid search
Ttot = 0;
t0  = tic;

% Does the user want to FIT or SIMULATE the data?
if( get(handles.radiobutton_fit, 'Value')==1 )
    taskString = 'FIT';
else
    taskString = 'SIM';
end

% Functions to perform
PERFORM_GRID        = get(handles.checkbox_grid, 'Value')==1;
CALCULATE_ERORRS    = get(handles.checkbox_errors, 'Value')==1;

% Does the user want to constrain the rate analysis (use dH and Eab)
CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;

for g = 1:session.Ng
    if( searchGroup(g) )
        group = session.groups{g};
        
        % Update text strings with progress
        set(handles.text_cs_name, 'String', ...
            sprintf('%s', group.name));            

        set(handles.text_total_progress, 'String', ...
            sprintf('Total progress [%d/%d]', Ng_Completed+1, NgSearch));
        
        % If grid search is desired
        if( PERFORM_GRID )            
            % Perform grid search
            group.gridSearch(taskString, CONSTRAIN_RATE_ANALYSIS, handles.axes_current_progress, handles.toggle_abort );            
        end
        
        % If error estimation is desired
        if( CALCULATE_ERORRS )
            % Does the user want to calculate errors on the BEST fit as marked
            %  or on the best result from the grid search?
            if( get(handles.radiobutton_error_on_best, 'Value')==1 )
                % ... calculate errors on best fit
                fitResult = group.fitResult_Best;
            else
                % ... calculate errors from the most recent fit (the GRID)
                fitResult = group.fitResults{group.Nf};
            end
            
            fitResult.calculateErrors( handles.axes_current_progress, handles.toggle_abort );
        end
        
        % If user wants to abort batch run
        if( get(handles.toggle_abort, 'Value') )
            set(handles.toggle_abort, 'Value', false);
            set(handles.toggle_abort, 'ForegroundColor', [0 0 0]);            
            return
        end
        
        % Update searched curve set
        %Npoints_Completed   = Npoints_Completed + Npoints_cs(cs);
        Npoints_Completed   = Npoints_Completed + Npoints_Group(g);
        Ng_Completed        = Ng_Completed + 1;

        Ttot                        = toc(t0);
        Npoints_Completed_per_sec   = Npoints_Completed / Ttot;
        
        %% Save results to HDD after each curve set
        fprintf('\n\nFinished group, writing to HDD...');
               
        setappdata(handles_main.main_gui, 'session', session);
        
        filename    = 'GUARDD-Session--Batch_Progress.mat';
        filepath    = session.outputDir;
        
        % Make the directory if it does not yet exist
        if( ~exist(filepath, 'dir') )
            mkdir(filepath);
        end
        
        % If outupt file exists, delete it and create a new one
        if( exist(sprintf('%s/%s', filepath,filename), 'file') )
            delete(sprintf('%s/%s', filepath,filename));
        end

        % Check each variable and save
        for v = 1:length(variables)
            fprintf('\n%s\t', variables{v} );
            if(length(variables{v})<8)
                fprintf('\t');
            end

            % If variable exists
            if( isappdata(handles_main.main_gui, variables{v}) )
                eval(sprintf('%s = getappdata(handles_main.main_gui, variables{v});',variables{v}));

                % First variable creates a new file
                if( ~exist(sprintf('%s/%s', filepath,filename), 'file') )
                    save( sprintf('%s/%s', filepath, filename), variables{v} );

                % Subsequent variables append the existing file
                else
                    save( sprintf('%s/%s', filepath, filename), variables{v}, '-append' );
                end

                fprintf('\tFound and saved');
            else
                fprintf('\tNo info to save');
            end
        end
        fprintf('\n Wrote file "%s"\n', filename);        
    end
    
    % Update total progress wait bar (Johannes Korsawe, 2007)
    value = Npoints_Completed / NpointsInTotal;
    h=handles.axes_total_progress;
    cla(h);        
    patch([0,value,value,0],[0,0,1,1],'k', 'Parent', h);
    axis(h,[0,1,0,1]); axis(h,'off'); drawnow;    
    
    [hour, minute, second] = sec2hms(Ttot);
    set(handles.text_time_elapsed, 'String', sprintf('%0.0f hr, %0.0f min, %0.0f sec', ...
        hour, minute, second));
    
    [hour, minute, second] = sec2hms((NpointsInTotal-Npoints_Completed)/Npoints_Completed_per_sec);
    set(handles.text_time_remaining, 'String', sprintf('%0.0f hr, %0.0f min, %0.0f sec', ...
        hour, minute, second));
end

function [hour, minute, second] = sec2hms(sec)
%% SEC2HMS  Convert seconds to hours, minutes and seconds.
%
%   [HOUR, MINUTE, SECOND] = SEC2HMS(SEC) converts the number of seconds in
%   SEC into hours, minutes and seconds.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 12:50:09 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

hour   = fix(sec/3600);      % get number of hours
sec    = sec - 3600*hour;    % remove the hours
minute = fix(sec/60);        % get number of minutes
sec    = sec - 60*minute;    % remove the minutes
second = sec;


% --- Executes on button press in toggle_abort.
function toggle_abort_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_abort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Abort current batch run?

% Color black/red for unaborted/aborted
color = [0 0 0];
if( get(handles.toggle_abort, 'Value') )
    color = [1 0 0];
end

set(handles.toggle_abort, 'ForegroundColor', color);


% --- Executes on button press in checkbox_grid.
function checkbox_grid_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_grid
refresh_display(handles);


% --- Executes on button press in checkbox_errors.
function checkbox_errors_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_errors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_errors
refresh_display(handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_save_grid_group.
function button_save_grid_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_grid_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Save grid to the selected group

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.batch_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check the selected group
g_selected  = getappdata(handles.batch_gui, 'g_selected');
group       = session.groups{g_selected};

% Does the user want to constrain the rate analysis (use dH and Eab)
CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;
if( CONSTRAIN_RATE_ANALYSIS )
    Np_Grid = 6;
else
    Np_Grid = 4;
end

% Read data from table
table_data = get(handles.table_grid_limits, 'Data');

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
set(handles.button_save_grid_group, 'Enable', 'off');
set(handles.button_restore, 'Enable', 'off');
refresh_display(handles);

% --- Executes on button press in button_restore.
function button_restore_Callback(hObject, eventdata, handles)
% hObject    handle to button_restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refresh_display(handles);


% --- Executes on button press in pushbutton_save_selected.
function pushbutton_save_selected_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Save grid to the ALL the selected groups

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.batch_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Does the user want to constrain the rate analysis (use dH and Eab)
CONSTRAIN_RATE_ANALYSIS = get(handles.checkbox_constrain_rates, 'Value')==1;
if( CONSTRAIN_RATE_ANALYSIS )
    Np_Grid = 6;
else
    Np_Grid = 4;
end

% Read data from table
table_data = get(handles.table_grid_limits, 'Data');

% Check which groups are selected
table_data_selected = get(handles.table_batch, 'Data');

for g = 1:session.Ng
    
    % Only save it to the group if it is selected
    if( table_data_selected{g,3} == true )
        group = session.groups{g};
        
        % Get prior grid so it can be partially updated, if needed
        grid_p_min              = group.grid_p_min;
        grid_p_max              = group.grid_p_max;
        grid_p_steps            = group.grid_p_steps;

        % Make update
        grid_p_min(1:Np_Grid)    = table_data(1,1:Np_Grid);
        grid_p_max(1:Np_Grid)    = table_data(2,1:Np_Grid);
        grid_p_steps(1:Np_Grid)  = table_data(3,1:Np_Grid);
        group.setGridLimits( grid_p_min, grid_p_max, grid_p_steps )      
    end    
end

refresh_display(handles);

% --- Executes when selected cell(s) is changed in table_batch.
function table_batch_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_batch (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

% Find the first group that has been selected
if( ~isempty(eventdata.Indices) )
    g_selected      = eventdata.Indices(1,1);
    g_selected_all  = eventdata.Indices(:,1);
    setappdata(handles.batch_gui, 'g_selected', g_selected);
    setappdata(handles.batch_gui, 'g_selected_all', g_selected_all);
    
    column_clicked  = eventdata.Indices(1,2);
    table_data      = get(handles.table_batch, 'Data');
    %{
    % If the user clicked the "In?" column...
    if( column_clicked == 3 )
        % ... then change its value
        prior_value = table_data{g_selected, 3};
        
        if( isempty(prior_value) )
            table_data{g_selected, 3} = '   X';
        else
            table_data{g_selected, 3} = '';
        end
        
        % Update the new value in the table
        set(handles.table_batch, 'Data', table_data);
    end
    %}
    refresh_display(handles);
end


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
%% User made change to grid search table
set(handles.button_save_grid_group, 'Enable', 'on');
set(handles.button_restore, 'Enable', 'on');

%{
% --- Executes when entered data in editable cell(s) in table_batch.
function table_batch_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_batch (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
% Find the first group that has been selected
if( ~isempty(eventdata.Indices) )
    g_selected = eventdata.Indices(1,1);
    setappdata(handles.batch_gui, 'g_selected', g_selected);
    
    
    table_data = get(handles.table_batch, 'Data');
    table_data
    table_data{eventdata.Indices(1,1), 3} = eventdata.NewData==1;
    table_data
    set(handles.table_batch, 'Data', table_data);
    
    refresh_display(handles);
end
%}


% --- Executes on button press in button_include_group.
function button_include_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_include_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Include the selected groups in the grid
table_data = get(handles.table_batch, 'Data');

% Check which groups are selected from the table
g_selected_all = getappdata(handles.batch_gui, 'g_selected_all');
for g = g_selected_all'    
    table_data{g,3} = true;
end

set(handles.table_batch, 'Data', table_data);
refresh_display(handles);

% --- Executes on button press in button_exclude_group.
function button_exclude_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_exclude_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Exclude the selected groups in the grid
table_data = get(handles.table_batch, 'Data');

% Check which groups are selected from the table
g_selected_all = getappdata(handles.batch_gui, 'g_selected_all');
for g = g_selected_all'    
    table_data{g,3} = false;
end

set(handles.table_batch, 'Data', table_data);
refresh_display(handles);


% --- Executes when entered data in editable cell(s) in table_sec_per_step.
function table_sec_per_step_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_sec_per_step (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
refresh_display(handles);


% --- Executes on button press in checkbox_constrain_rates.
function checkbox_constrain_rates_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_constrain_rates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_constrain_rates
refresh_display(handles);

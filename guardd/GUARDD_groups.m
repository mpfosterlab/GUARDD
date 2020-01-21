% Display fitted parameters for all NMR signals, address nature of group-based dynamics (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/02/01 Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/05/10 Use classes
% 2011/05/20 Plot histograms
% 
% TO DO
%  

function varargout = GUARDD_groups(varargin)
% GUARDD_GROUPS M-file for GUARDD_groups.fig
%      GUARDD_GROUPS, by itself, creates a new GUARDD_GROUPS or raises the existing
%      singleton*.
%
%      H = GUARDD_GROUPS returns the handle to a new GUARDD_GROUPS or the handle to
%      the existing singleton*.
%
%      GUARDD_GROUPS('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in GUARDD_GROUPS.M with the given input arguments.
%
%      GUARDD_GROUPS('Property','Value',...) creates a new GUARDD_GROUPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_groups_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_groups_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_groups

% Last Modified by GUIDE v2.5 08-Jun-2011 22:04:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_groups_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_groups_OutputFcn, ...
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



% --- Executes just before GUARDD_groups is made visible.
function GUARDD_groups_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_groups (see VARARGIN)

% Choose default command line output for GUARDD_groups
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_groups wait for user response (see UIRESUME)
% uiwait(handles.displaycluster_gui);

%% Update display upon creation

% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
% AND same for fit_r2eff window which called this GUI
setappdata(handles.displaycluster_gui, 'handles_main', handles_main);

% Access main window handle which is stored in this window's main handle
session = getappdata(handles_main.main_gui, 'session');

%{
% Initialize variables which hold info from GUI elements
group_analysis.plot_type_Array    = 1;
group_analysis.param_type_x = 1;
group_analysis.param_type_y = 1;
group_analysis.Temp_type_x  = 1;
group_analysis.Temp_type_y  = 1;
group_analysis.B0_type_x    = 1;
group_analysis.B0_type_y    = 1;
setappdata(handles_main.main_gui, 'group_analysis', group_analysis);
%}

% Set up default subplot settings
if( session.Ndc > 0 )
    % Default group is first one
    set(handles.popup_displaycluster, 'Value', 1);
else
    
end

% Number of subplots (rows and cols)
set(handles.table_subplot, 'ColumnWidth', {40,40});
set(handles.table_subplot, 'ColumnName', {'Rows', 'Cols'});
set(handles.table_subplot, 'Data', [1 1]);

% Set up types of plots for list
plot_type_string = session.paramDisplay.getPlotTypeString();
set(handles.popup_plot_type, 'String', plot_type_string);
set(handles.popup_plot_type, 'Value',1);

% Disable save button
set(handles.button_save_displaycluster, 'Enable', 'off');
set(handles.button_restore_displaycluster, 'Enable', 'off');

% Set default plot limits
set(handles.table_plot_limits, 'ColumnName', {'Min', 'Max', 'Log'});
set(handles.table_plot_limits, 'RowName', {'x', 'y'});
set(handles.table_plot_limits, 'Data', {NaN, NaN, false; NaN, NaN, false});

refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_groups_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% This is called to update all display elements
function refresh_display(handles)
    
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

Ndc             = session.Ndc;
if( Ndc == 0 )
    % Set title of panel
    set(handles.panel_scatter_plots, 'Title', sprintf('Add new display cluster' ));
    
    set(handles.popup_displaycluster, 'Enable', 'off');
    set(handles.button_delete_displaycluster, 'Enable', 'off');
   
    set(handles.panel_group_settings, 'Visible', 'off');
    set(handles.panel_show_all_groups, 'Visible', 'off');
    set(handles.panel_groups_in_cluster, 'Visible', 'off');
    set(handles.panel_display_settings, 'Visible', 'off');
    
    %{
    set(handles.table_subplot,      'Enable', 'off');
    set(handles.popup_subplot_num,  'Enable', 'off');
    set(handles.popup_plot_type,    'Enable', 'off');
    set(handles.popup_param_x,     'Enable', 'off');
    set(handles.popup_temp_x,       'Enable', 'off');
    set(handles.popup_B0_x,         'Enable', 'off');
    set(handles.popup_param_y,     'Enable', 'off');
    set(handles.popup_temp_y,       'Enable', 'off');
    set(handles.popup_B0_y,         'Enable', 'off');
    %}
    return
end

%% Show information on current displaycluster
set(handles.popup_displaycluster, 'Enable', 'on');
set(handles.button_delete_displaycluster, 'Enable', 'on');
set(handles.panel_group_settings, 'Visible', 'on');
set(handles.panel_show_all_groups, 'Visible', 'on');
set(handles.panel_groups_in_cluster, 'Visible', 'on');
set(handles.panel_display_settings, 'Visible', 'on');

% List all DisplayClusters in popup menu
displaycluster_names = session.getDisplayClusterStringArray();
set(handles.popup_displaycluster, 'String', displaycluster_names )


% Read current group number from popup menu
dc              = get(handles.popup_displaycluster, 'Value');
displayCluster  = session.displayClusters{dc};

% Group name and description
set(handles.edit_displaycluster_name, 'String', displayCluster.name );
set(handles.edit_displaycluster_description, 'String', displayCluster.description );

% Update group color table "table_displaycluster_color"
table_data(1:3) = displayCluster.color_RGB(1:3);
set(handles.table_displaycluster_color, 'ColumnWidth', {40,40,40});
set(handles.table_displaycluster_color, 'ColumnName', {'R', 'G', 'B'});
set(handles.table_displaycluster_color, 'Data', table_data);

% Display the current group?
set(handles.checkbox_display_data, 'Value', displayCluster.display_Data);
set(handles.checkbox_display_labels, 'Value', displayCluster.display_Labels);

%% Fill list box with ALL group names formatted for FIT/EXCHANGE status
listbox_string  = cell(1, session.Ng);
for g = 1:session.Ng
    group = session.groups{g};
    
    if( displayCluster.containsGroup(group) )
        in_cluster_string   = '- ';
    else
        in_cluster_string   = '+ ';
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
        listbox_string{g} = sprintf('<html><b>%s%s%s</b></html>', ...
            in_cluster_string, marker_string, group.name );

    % Fit available but NOT OK
    elseif( group.Nf > 0 && ~group.bestFitIsOK )
        listbox_string{g} = sprintf('<html><font color="Gray"><b>%s%s%s</b></font></html>', ...
            in_cluster_string, marker_string, group.name );

    % Fit is not available
    else
        listbox_string{g} = sprintf('<html><i>%s%s%s</i></html>', ...
            in_cluster_string, marker_string, group.name );
    end
end

% Set contents of listbox for all curve sets
set(handles.listbox_all_groups, 'String', listbox_string);

% Update list box
if( get(handles.listbox_all_groups, 'Value') > session.Ng )    
    set(handles.listbox_all_groups, 'Value', 1);
end

%% Fill list box with group names ONLY IN THE DISPLAYCLUSTER formatted for FIT/EXCHANGE status
listbox_string  = cell(1, displayCluster.Ng);
for g = 1:displayCluster.Ng
    group = displayCluster.groups{g};
    
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
set(handles.listbox_groups_in_displaycluster, 'String', listbox_string);

% Update list box
if( get(handles.listbox_groups_in_displaycluster, 'Value') > displayCluster.Ng )    
    set(handles.listbox_groups_in_displaycluster, 'Value', 1);
end

% Check if there is NO information to show
if( session.getNumGroupsInDisplayClusters() == 0 )
    set(handles.panel_scatter_plots, 'Title', sprintf('(Add a Group to a DisplayCluster)'));
    set(handles.panel_display_settings, 'Visible', 'off');
    delete( get(handles.panel_scatter_plots, 'Children') );
    return
end

if( ~session.displayClustersContainOKFit() )
    set(handles.panel_scatter_plots, 'Title', sprintf('(None of the groups contains an "OK" fit)'));
    set(handles.panel_display_settings, 'Visible', 'off');
    delete( get(handles.panel_scatter_plots, 'Children') );
    return    
end

%% Set title of panel and enable gui elements
set(handles.panel_display_settings, 'Visible', 'on');
set(handles.panel_scatter_plots, 'Title', '');
%{
set(handles.table_subplot, 'Enable', 'on');
set(handles.popup_subplot_num, 'Enable', 'on');
set(handles.popup_plot_type, 'Enable', 'on');
%}

%% Generate sub plot types
paramDisplay = session.paramDisplay;

% This should list anything the user wants to plot along the sequence
set(handles.popup_param_x, 'String', paramDisplay.param_type_string);
set(handles.popup_param_y, 'String', paramDisplay.param_type_string);

% Generate available temperatures
set(handles.popup_temp_x, 'String', paramDisplay.getUniqueTempertureString());
set(handles.popup_temp_y, 'String', paramDisplay.getUniqueTempertureString());

% Generate available B0s
set(handles.popup_B0_x, 'String', paramDisplay.getUniqueB0String());
set(handles.popup_B0_y, 'String', paramDisplay.getUniqueB0String());

%% Read in number of ROW x COL for subplotting from table_subplot
table_data = get(handles.table_subplot, 'Data');

subplot_rows    = table_data(1);
subplot_cols    = table_data(2);
Nsubplots       = subplot_rows * subplot_cols;

% List all subplots in popup_subplot_num
subplot_string{Nsubplots} = '';
for sp = 1:Nsubplots
    subplot_string{sp} = sprintf('Subplot %02d', sp);
end
set(handles.popup_subplot_num, 'String', subplot_string);

% Update the number of subplots, and the corresponding options arrays
paramDisplay.updateNumSubplots(subplot_rows, subplot_cols);

%% Show plot type for the selected number
subplot_selected = get(handles.popup_subplot_num, 'Value');

% Make sure subplot selected is not outside the limits
if( subplot_selected > Nsubplots )
    subplot_selected = 1;
    set(handles.popup_subplot_num, 'Value', 1);
end

% Set GUI elements to be show current subplot specs
set(handles.popup_plot_type,'Value', paramDisplay.plot_type_num_Array(subplot_selected));
set(handles.popup_param_x,  'Value', paramDisplay.param_num_X_Array(subplot_selected));
set(handles.popup_param_y,  'Value', paramDisplay.param_num_Y_Array(subplot_selected));
set(handles.popup_temp_x,   'Value', paramDisplay.Temp_num_X_Array(subplot_selected));
set(handles.popup_temp_y,   'Value', paramDisplay.Temp_num_Y_Array(subplot_selected));
set(handles.popup_B0_x,     'Value', paramDisplay.B0_num_X_Array(subplot_selected));
set(handles.popup_B0_y,     'Value', paramDisplay.B0_num_Y_Array(subplot_selected));
set(handles.checkbox_histogram, 'Value', paramDisplay.displayHistogram_Array(subplot_selected));

%% Get information on B0 and/or Temp if necessary

% If it is the custom field, enable additional GUI elements
plot_type           = get(handles.popup_plot_type,'Value');
plot_type_strings   = get(handles.popup_plot_type, 'String');

% Get the temperature from any of these parameters
if( strcmpi(plot_type_strings{plot_type}, 'CUSTOM') )
    set(handles.popup_param_x,     'Enable', 'on');
    set(handles.popup_param_y,     'Enable', 'on');
    set(handles.checkbox_histogram,'Enable', 'on');
    
    % Get names of selected x and y parameters
    param_string    = get(handles.popup_param_x, 'String');
    
    param_type_x            = get(handles.popup_param_x, 'Value');
    param_name_x            = param_string{param_type_x};
    [needsTemp, needsB0]    = Session.getParamRequirements( param_name_x );
        
    % CUSTOM X
    if( needsTemp )
        % Enable popup window and set value
        set(handles.popup_temp_x, 'Enable', 'on');        
        set(handles.popup_temp_x, 'Value', paramDisplay.Temp_num_X_Array(subplot_selected));

        % Get the B0 from either of these too
        if( needsB0 )
            % Enable popup window and set value
            set(handles.popup_B0_x, 'Enable', 'on');
            set(handles.popup_B0_x, 'Value', paramDisplay.B0_num_X_Array(subplot_selected));
        else
            set(handles.popup_B0_x, 'Enable', 'off');
        end

    % Extra information on T and B0 is not needed
    else
        set(handles.popup_temp_x, 'Enable', 'off');
        set(handles.popup_B0_x, 'Enable', 'off');
    end
    
    % Display of histogram only uses X-axis values
    DISPLAY_HISTOGRAM = get(handles.checkbox_histogram, 'Value') == 1;
    
    if( DISPLAY_HISTOGRAM )
        set(handles.popup_param_y, 'Enable', 'off');
        set(handles.popup_temp_y,  'Enable', 'off');    
        set(handles.popup_B0_y,    'Enable', 'off');
    
    else
        param_type_y            = get(handles.popup_param_y, 'Value');
        param_name_y            = param_string{param_type_y};
        [needsTemp, needsB0]    = Session.getParamRequirements( param_name_y );

        % CUSTOM Y
        if( needsTemp )
            % Enable popup window and set value
            set(handles.popup_temp_y, 'Enable', 'on');        
            set(handles.popup_temp_y, 'Value', paramDisplay.Temp_num_Y_Array(subplot_selected));

            % Get the B0 from either of these too
            if( needsB0 )
                % Enable popup window and set value
                set(handles.popup_B0_y, 'Enable', 'on');
                set(handles.popup_B0_y, 'Value', paramDisplay.B0_num_Y_Array(subplot_selected));
            else
                set(handles.popup_B0_y, 'Enable', 'off');
            end

        % Extra information on T and B0 is not needed
        else
            set(handles.popup_temp_y, 'Enable', 'off');
            set(handles.popup_B0_y, 'Enable', 'off');
        end
    end
    
% Not a custom parameter
else
    set(handles.popup_param_x,      'Enable', 'off');
    set(handles.popup_temp_x,       'Enable', 'off');    
    set(handles.popup_B0_x,         'Enable', 'off');
    set(handles.popup_param_y,      'Enable', 'off');
    set(handles.popup_temp_y,       'Enable', 'off');    
    set(handles.popup_B0_y,         'Enable', 'off');  
    set(handles.checkbox_histogram, 'Enable', 'off');
end

%% Plot the data in the panel window
plot_groups( handles.panel_scatter_plots, handles );


function plot_groups( figure_handle, handles )

%% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui,    'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% This hold the information for display of the parameters
%  How many subplots, and what should be displayed on X and Y
paramDisplay    = session.paramDisplay;

%% Get number of subplots from data structure
Nsubplots       = paramDisplay.Nsubplots;
subplot_rows    = paramDisplay.subplot_rows;
subplot_cols    = paramDisplay.subplot_cols;

%% Clear axes
% Clear old contents and create new axes (instead of subplot(1,1,1))
delete( get(figure_handle, 'Children') );
h = axes('Parent', figure_handle);
%hold(h, 'off');
hold(h, 'all');

%% Iterate through each subplot
for sp = 1:Nsubplots

    % Get plot type from popup menu    
    plot_type           = paramDisplay.getValue(sp, 'PLOT', 'NA');
    param_names         = get(handles.popup_param_x, 'String');    
    DISPLAY_HISTOGRAM   = paramDisplay.getValue(sp, 'HISTOGRAM', 'NA');
    
    % Build the data for all curves
    %X   = zeros(1, data.Ncs);
    %X_E = zeros(1, data.Ncs);
    %Y   = zeros(1, data.Ncs);
    %Y_E = zeros(1, data.Ncs);
    
    XMIN_DEFAULT    = NaN;    
    YMIN_DEFAULT    = NaN;
    
    % f_abshhx and f_abshhy used for setting widths on error bar handles
    % Values determined empirically for each plot type
    %f_abshhx    = 0.1;
    %f_abshhy    = 0.01;
    
    f_abshhx    = 0;
    f_abshhy    = 0;
    
    % There are a total of Ntemps+3 plot types
    
    if( plot_type == 1 )
        % van't Hoff
        param_name_x    = 'DH';             
        Temp_x          = NaN;
        B0_x            = NaN;
        
        param_name_y    = 'DS';
        Temp_y          = NaN;
        B0_y            = NaN;
        
        title_string    = 'vant Hoff scatter';
        
        %f_abshhx    = 0.1;
        %f_abshhy    = 0.01;
        
    
    elseif( plot_type == 2 )
        % Arrhenius
        param_name_x    = 'EAB';        
        Temp_x          = NaN;
        B0_x            = NaN;
        
        param_name_y    = 'EBA';
        Temp_y          = NaN;
        B0_y            = NaN;
        
        title_string    = 'Arrhenius scatter';
        
        %f_abshhx    = 0.01;
        %f_abshhy    = 0.1;        
    
    elseif( 3 <= plot_type && plot_type <= paramDisplay.NTemp+2 )
        % Plot types 3 through (Ntemps+2) are for kA and kB
        param_name_x    = 'KA';
        t               = plot_type-2;
        Temp_x          = paramDisplay.Temp_values(t);
        B0_x            = NaN;
        
        param_name_y    = 'KB';
        Temp_y          = Temp_x;
        B0_y            = NaN;
        
        title_string    = sprintf('k_A,k_B scatter at %d^oC',Temp_x-273);
    
    elseif( plot_type == paramDisplay.NTemp+3 )
        % For the custom plot type
        
        param_name_x    = paramDisplay.getValue(sp, 'PARAM', 'X');
        %param_name_x    = param_names{param_type_x};
        Temp_x          = paramDisplay.getValue(sp, 'TEMP', 'X');
        B0_x            = paramDisplay.getValue(sp, 'B0', 'X');
        
        param_name_y    = paramDisplay.getValue(sp, 'PARAM', 'Y');
        %param_name_y    = param_names{param_type_y};
        Temp_y          = paramDisplay.getValue(sp, 'TEMP', 'Y');
        B0_y            = paramDisplay.getValue(sp, 'B0', 'Y');      
    else
        error('Invalid plot type selected');
    end
    
    % Get axes labels and minimum values
    [VOID, axis_label_x, X_MIN] = ...
            session.getPlotLabels( param_name_x, Temp_x, B0_x );
        
    [VOID, axis_label_y, Y_MIN] = ...
        session.getPlotLabels( param_name_y, Temp_y, B0_y );
    
    if( plot_type == paramDisplay.NTemp+3 )
        title_string = sprintf('%s vs. %s', axis_label_y, axis_label_x);
    end
    
    % Now get the DATA for each DisplayCluster
    X       = cell(1,session.Ndc);
    X_E     = cell(1,session.Ndc);
    X_IS_OK = cell(1,session.Ndc);
    Y       = cell(1,session.Ndc);
    Y_E     = cell(1,session.Ndc);
    Y_IS_OK = cell(1,session.Ndc);
    DISPLAY = cell(1,session.Ndc);
    
    % Check to see if 1 or all curvesets should be displayed, considering
    % the paramter type    
    % If either X or Y value is unique for each curveset...
    if( strcmpi(param_name_x, 'DWH_PPM') || ...
        strcmpi(param_name_x, 'DWX_PPM') || ...
        strcmpi(param_name_y, 'DWH_PPM') || ...
        strcmpi(param_name_y, 'DWX_PPM') )    
        % ...make sure all curvesets are shown
        CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY = 'ALL';
        
    else
        % Otherwise, only show the first curveset
        CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY = 'ANY';
    end   
    
    % Build the data arrays
    for dc = 1:session.Ndc
        displayCluster = session.displayClusters{dc};
                
        [X{dc}, X_E{dc}, X_IS_OK{dc}, f_abshh_x] = ...
            displayCluster.getData(param_name_x, Temp_x, B0_x, ...
            CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY);

        [Y{dc}, Y_E{dc}, Y_IS_OK{dc}, f_abshh_y] = ...
            displayCluster.getData(param_name_y, Temp_y, B0_y, ...
            CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY);
        
        % In case user wants offset in residue number
        % IK used for TRAP BSte -> BSub
        OFFSET_RESIDUE_X    = false;        
        if( OFFSET_RESIDUE_X && strcmpi(param_name_x,'RESIDUE') )
            'Using offset in residue number'
            OFFSET_VALUE        = 2;            
            X{dc} = X{dc} + OFFSET_VALUE;
            X_E{dc} = X_E{dc} + OFFSET_VALUE;
        end
        
        % Convert units to display
        X{dc}   = Session.convertUnits(X{dc},   param_name_x, 'DISPLAY');
        X_E{dc} = Session.convertUnits(X_E{dc}, param_name_x, 'DISPLAY');
        Y{dc}   = Session.convertUnits(Y{dc},   param_name_y, 'DISPLAY');
        Y_E{dc} = Session.convertUnits(Y_E{dc}, param_name_y, 'DISPLAY');
                
        % DISPLAY holds NaN for data that are NOT OK, and "1" for OK
        XY_OK                       = and(X_IS_OK{dc}, Y_IS_OK{dc});
        DISPLAY_TEMP                = NaN*ones( 1,length(X_IS_OK{dc}) );                
        DISPLAY_TEMP( find(XY_OK) ) = 1;        
        DISPLAY{dc}                 = DISPLAY_TEMP;
        if( isempty(DISPLAY{dc}) )
             DISPLAY{dc} = [];
        end
    end
                
    %% Set log scale if desired    
    % Get limit information from table
    table_data  = get(handles.table_plot_limits, 'Data');
    log_x       = table_data{1,3};
    log_y       = table_data{2,3};
    
    if( log_x )
        set(h,'XScale', 'log');
        set(h,'XGrid', 'off');        
    else
        set(h,'XScale', 'linear');
        set(h,'XGrid','on');
    end
    
    if( log_y )
        set(h,'YScale', 'log');
        set(h,'YGrid', 'off');
    else
        set(h,'YScale', 'linear');
        set(h,'YGrid','on');
    end
    
    %% Set up size of error bars    
    % Find limits of all data to be shown
    XDATA_MIN = NaN;
    XDATA_MAX = NaN;
    YDATA_MIN = NaN;
    YDATA_MAX = NaN;
    
    for dc = 1:session.Ndc
        displayCluster = session.displayClusters{dc};
        if( displayCluster.display_Data )            
            % Calculate the new minimum among this and all prior displayClusters
            XDATA_MIN = min( [XDATA_MIN (X{dc}-X_E{dc}) .* DISPLAY{dc}] );
            XDATA_MAX = max( [XDATA_MAX (X{dc}+X_E{dc}) .* DISPLAY{dc}] );
            YDATA_MIN = min( [YDATA_MIN (Y{dc}-Y_E{dc}) .* DISPLAY{dc}] );
            YDATA_MAX = max( [YDATA_MAX (Y{dc}+Y_E{dc}) .* DISPLAY{dc}] );
        end
    end
    
    XDATA_RANGE = XDATA_MAX - XDATA_MIN;
    YDATA_RANGE = YDATA_MAX - YDATA_MIN;
        
    %% Plot data
    
    % Clear the axes in panel
    h = subplot(subplot_rows,subplot_cols,sp,'Parent', figure_handle);
    cla(h);
    hold(h, 'all');
    
    for dc = 1:session.Ndc
        displayCluster = session.displayClusters{dc};
        
        if( displayCluster.display_Data )
            
            if( DISPLAY_HISTOGRAM )
                [n,xout] = hist(h, X{dc} .* DISPLAY{dc});
                bar(h, xout, n, 'FaceColor', displayCluster.color_RGB, 'EdgeColor', 'k');
                %hh = findobj(h,'Type','patch');
                %set(hh, 'FaceColor', displayCluster.color_RGB, 'EdgeColor', 'k');
        
            else
                % May have to scale error bar widths with axis range
                if( log_x && log_y )
                    hp = ploterr_ik(h, X{dc}   .* DISPLAY{dc},   Y{dc} .* DISPLAY{dc}, ...
                                  X_E{dc} .* DISPLAY{dc}, Y_E{dc} .* DISPLAY{dc}, ...
                                'o', 'logxy', 'abshhx', 0.02, 'abshhy', 0.02 );

                elseif( log_x && ~log_y )
                    hp = ploterr_ik(h, X{dc}   .* DISPLAY{dc},   Y{dc} .* DISPLAY{dc}, ...
                                  X_E{dc} .* DISPLAY{dc}, Y_E{dc} .* DISPLAY{dc}, ...
                                'o', 'logx', 'abshhx', 0.05, 'abshhy', 0.05 );

                elseif( ~log_x && log_y )                
                    hp = ploterr_ik(h, X{dc}   .* DISPLAY{dc},   Y{dc} .* DISPLAY{dc}, ...
                                  X_E{dc} .* DISPLAY{dc}, Y_E{dc} .* DISPLAY{dc}, ...
                                'o', 'logy', 'abshhx', 0.05, 'abshhy', 0.05 );

                elseif( ~log_x && ~log_y )
                    hp = ploterr_ik(h, X{dc}   .* DISPLAY{dc},   Y{dc} .* DISPLAY{dc}, ...
                                  X_E{dc} .* DISPLAY{dc}, Y_E{dc} .* DISPLAY{dc}, ...
                                'o', 'abshhx', f_abshhx*XDATA_RANGE, 'abshhy', f_abshhy*YDATA_RANGE );
                end

                % Color the displayCluster
                set(hp, 'LineWidth', session.LINEWIDTH);
                set(hp, 'Color', displayCluster.color_RGB);
            end
        end
    end
    
    %% Update plot limits if required
    
    % Read auto-scale limits
    set(h, 'XLimMode', 'Auto');
    set(h, 'YLimMode', 'Auto');
    
    x_limits = get(h,'Xlim');
    y_limits = get(h,'Ylim');
    
    XMIN_AUTO = x_limits(1);
    XMAX_AUTO = x_limits(2);
    YMIN_AUTO = y_limits(1);
    YMAX_AUTO = y_limits(2);
        
    % Get limit information from table
    table_data = get(handles.table_plot_limits, 'Data');
    
    XMIN = table_data{1,1};
    XMAX = table_data{1,2};
    YMIN = table_data{2,1};
    YMAX = table_data{2,2};
        
    % Set only XMIN
    if( ~isnan(XMIN) && isnan(XMAX) )
        set(h, 'XLim', [XMIN XMAX_AUTO]);
    end
    % Set only XMAX
    if( isnan(XMIN) && ~isnan(XMAX) )        
        % Check if the default value should be used instead of autoscale
        if( ~isnan(XMIN_DEFAULT) )
            set(h, 'XLim', [XMIN_DEFAULT XMAX]);
        else
            set(h, 'XLim', [XMIN_AUTO XMAX]);
        end
    end
    % Set both XMIN and XMAX
    if( ~isnan(XMIN) && ~isnan(XMAX) )
        set(h, 'XLim', [XMIN XMAX]);
    end
    
    % If both are auto and default value supplied
    if( isnan(XMIN) && isnan(XMAX) && ~isnan(XMIN_DEFAULT) )
        set(h, 'XLim', [XMIN_DEFAULT XMAX_AUTO]);
    end    
    
    % Set only YMIN
    if( ~isnan(YMIN) && isnan(YMAX) )
        set(h, 'YLim', [YMIN YMAX_AUTO]);
    end
    % Set only YMAX
    if( isnan(YMIN) && ~isnan(YMAX) )
        % Check if the default value should be used instead of autoscale
        if( ~isnan(YMIN_DEFAULT) )
            set(h, 'YLim', [YMIN_DEFAULT YMAX]);
        else
            set(h, 'YLim', [YMIN_AUTO YMAX]);
        end
    end
    % Set both YMIN and YMAX
    if( ~isnan(YMIN) && ~isnan(YMAX) )
        set(h, 'YLim', [YMIN YMAX]);
    end
    % If both are auto and default value supplied
    if( isnan(YMIN) && isnan(YMAX) && ~isnan(YMIN_DEFAULT) )
        set(h, 'YLim', [YMIN_DEFAULT YMAX_AUTO]);
    end
    
    %% Add text labels
    % This has to be done after all the plotting to get the proper plot ranges
    
    % Determine appropriate fontsize
    if( Nsubplots == 1 )
        max_fontsize = session.FONTSIZE_MEDIUM;
        
    elseif( Nsubplots == 2 )
        max_fontsize = session.FONTSIZE_SMALL;
        
    elseif( Nsubplots <= 4 )
        max_fontsize = session.FONTSIZE_EXTRA_SMALL;
        
    else        
        max_fontsize = session.FONTSIZE_EXTRA_SMALL-2;
    end
    
    if( ~DISPLAY_HISTOGRAM )
    
        % Read limits again in case they changed above
        x_limits    = get(h,'Xlim');
        y_limits    = get(h,'Ylim');

        x_range     = x_limits(2)-x_limits(1);
        y_range     = y_limits(2)-y_limits(1);
        F_TEXT_LABEL = session.F_TEXT_LABEL;
        
        % Iterate through each displayCluster for plotting
        for dc = 1:session.Ndc
            displayCluster = session.displayClusters{dc};

            if( displayCluster.display_Labels )      
                % Get X and Y vectors for the displaycluster
                %  Each point in the vector is for a curveset in a group in the cluster                
                Xdc         = X{dc};
                X_IS_OKdc   = X_IS_OK{dc};
                Ydc         = Y{dc};
                Y_IS_OKdc   = Y_IS_OK{dc};
                cstot       = 0;
                
                %DISPLAYdc   = DISPLAY{dc};
                
                for g = 1:displayCluster.Ng
                    group = displayCluster.groups{g};
                    
                    % Should ALL the curvesets be shown, or not?
                    if( strcmpi(CURVESET_ANY_OR_ALL_OR_NUMERIC_ARRAY, 'ALL') )
                        cs_array = 1:group.Ncs;
                        SHOW_CS_LABEL = true;
                    else
                        cs_array = 1;
                        SHOW_CS_LABEL = false;
                    end
                    
                    % Iterate through each curveset relevant to the display
                    for cs = cs_array
                        curveset    = group.curvesets{cs};
                        cstot       = cstot+1;
                        
                        % Only show the label if it is OK
                        if( X_IS_OKdc(cstot) && Y_IS_OKdc(cstot) )

                            % Text label needs repositioning for log scale
                            if( log_x )
                                x_mod = log(x_range)*F_TEXT_LABEL;
                            else
                                x_mod =    (x_range)*F_TEXT_LABEL;
                            end

                            if( log_y )
                                y_mod = log(y_range)*F_TEXT_LABEL;
                            else
                                y_mod =    (y_range)*F_TEXT_LABEL;
                            end

                            if( SHOW_CS_LABEL )
                                TEXT_LABEL = sprintf('%s\n[%s]', group.name, curveset.name);
                            else
                                TEXT_LABEL = sprintf('%s', group.name);
                            end
                            
                            ht = text( (Xdc(cstot) + x_mod), (Ydc(cstot) + y_mod), TEXT_LABEL, ...
                                'HorizontalAlignment', 'Left', ...
                                'FontSize', max_fontsize-2, ...
                                'FontName', session.FONTNAME);
                            set(ht, 'Color', displayCluster.color_RGB);
                        end
                    end
                end
            end
        end
    end
    
    %% Set axes titles (depends on type of data plotted)    
    
    set(h,'FontName', session.FONTNAME, 'FontSize', max_fontsize-2, 'FontWeight', 'Normal', 'LineWidth', session.LINEWIDTH)
    title(sprintf('[%d] %s', sp, title_string), 'FontSize', max_fontsize, 'FontWeight', 'Bold');
    xlabel(axis_label_x,   'FontSize', max_fontsize-1, 'FontWeight', 'Bold');
    
    if( DISPLAY_HISTOGRAM )
        ylabel('Number',   'FontSize', max_fontsize-1, 'FontWeight', 'Bold');     
    else
        ylabel(axis_label_y,   'FontSize', max_fontsize-1, 'FontWeight', 'Bold');
    end
    
    % Grid turned on/off must be done last to work on final subplot
    set(h,'XGrid','on');
    if( log_x )
        set(h,'XGrid', 'off');
    end
    
    set(h,'YGrid','on');
    if( log_y )
        set(h,'YGrid', 'off');
    end
    
    set(h, 'Box', 'on');
end



% --- Executes when displaycluster_gui is resized.
function displaycluster_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to displaycluster_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Resize each element of the GUI
%{
% Get size of entire GUI window
fP = get(handles.displaycluster_gui, 'Position');

%% Main panel for figure
% Adjust (scale size, fix relative position)
left    = 0.5;
width   = 34;
height  = 4.4;
bottom  = fP(4)-height;
set(handles.panel_manage, 'Position', [left, bottom, width, height]);

% Adjust (scale size, fix relative position)
left    = 0.5;
width   = 34;
height  = 12;
bottom  = fP(4)-4.4-0.1-height;
set(handles.panel_group_settings, 'Position', [left, bottom, width, height]);

% Adjust (scale size, fix relative position)
left    = 0.5;
width   = 34;
bottom  = 0.1;
height  = fP(4)-4.4-0.1-12-0.1;

set(handles.panel_groups_in_cluster, 'Position', [left, bottom, width, height]);

% Adjust (scale size, fix relative position)
left    = 0.1;
bottom  = 5;
width   = 33;
height  = fP(4)-4.4-0.1-12-0.1-bottom-1.5;
set(handles.table_groups_in_displaycluster, 'Position', [left, bottom, width, height]);

% Adjust (scale size, fix relative position)
left    = 0.1+34+0.5;
width   = fP(3)-0.5-34;
height  = 7.7;
bottom  = fP(4)-height;
set(handles.panel_display_settings, 'Position', [left, bottom, width, height]);

% Adjust (scale size, fix relative position)
left    = 0.1+34+0.5;
width   = fP(3)-34-0.5;
bottom  = 0.1;
height  = fP(4)-7.7-0.1;
set(handles.panel_scatter_plots, 'Position', [left, bottom, width, height]);


%{
%% Buttons
% Adjust (scale position, fix size)
width   = 33;
height  = 1.5;
left    = 1;
bottom  = fP(4)-2;
set(handles.popup_displaycluster, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 8;
height  = 1.8;
left    = 0.5;
bottom  = fP(4)-2-1.5-0.5;
set(handles.button_new_displaycluster, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 8;
height  = 1.8;
left    = 0.5+8+0.5;
bottom  = fP(4)-2-1.5-0.5;
set(handles.button_save_displaycluster, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 8;
height  = 1.8;
left    = 0.5+8+0.5+8+0.5;
bottom  = fP(4)-2-1.5-0.5;
set(handles.button_restore_displaycluster, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 8;
height  = 1.8;
left    = 0.5+8+0.5+8+0.5+8+0.5;
bottom  = fP(4)-2-1.5-0.5;
set(handles.button_delete_displaycluster, 'Position', [left, bottom, width, height]);
%}

%% Tables
%{
% Adjust (scale position, fix size)
width   = 13.8;
height  = 2.8;
left    = 1+33+1+6;
bottom  = fP(4)-3.3;
set(handles.table_subplot, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 25.5;
height  = 4;
left    = 1+33+2;
bottom  = fP(4)-3.3-0.5-4;
set(handles.table_plot_limits, 'Position', [left, bottom, width, height]);

%% Popup menus

% Adjust (scale position, fix size)
width   = 36;
height  = 1.5;
left    = 1+33+1+25.5+2;
bottom  = fP(4)-2;
set(handles.popup_subplot_num, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 36;
height  = 1.5;
left    = 1+33+1+25.5+2;
bottom  = fP(4)-3.5;
set(handles.popup_plot_type, 'Position', [left, bottom, width, height]);

% --- X ---
% Adjust (scale position, fix size)
width   = 10;
height  = 1.5;
left    = 1+33+1+25.5+2;
bottom  = fP(4)-3.5-2.5;
set(handles.popup_param_x, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 8;
height  = 1.5;
left    = 1+33+1+25.5+2+10+1;
bottom  = fP(4)-3.5-2.5;
set(handles.popup_temp_x, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 15;
height  = 1.5;
left    = 1+33+1+25.5+2+10+1+8+1;
bottom  = fP(4)-3.5-2.5;
set(handles.popup_B0_x, 'Position', [left, bottom, width, height]);

% --- Y ---
% Adjust (scale position, fix size)
width   = 10;
height  = 1.5;
left    = 1+33+1+25.5+2;
bottom  = fP(4)-3.5-2.5-1.5;
set(handles.popup_param_y, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 8;
height  = 1.5;
left    = 1+33+1+25.5+2+10+1;
bottom  = fP(4)-3.5-2.5-1.5;
set(handles.popup_temp_y, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 15;
height  = 1.5;
left    = 1+33+1+25.5+2+10+1+8+1;
bottom  = fP(4)-3.5-2.5-1.5;
set(handles.popup_B0_y, 'Position', [left, bottom, width, height]);
%}

%% Edit boxes
%{

% Adjust (scale position, fix size)
width   = 33;
height  = 1.8;
left    = 0.5;
bottom  = fP(4)-2-1.5-0.5-1.8-0.5;
set(handles.edit_displaycluster_name, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 33;
height  = 5.5;
left    = 0.5;
bottom  = fP(4)-2-1.5-0.5-1.8-0.5-0.5-5.5;
set(handles.edit_displaycluster_description, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 20.45;
height  = 2.8;
left    = 0.5;
bottom  = fP(4)-2-1.5-0.5-1.8-0.5-0.5-5.5-0.5-3;
set(handles.table_displaycluster_color, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 12;
height  = 1.8;
left    = 0.5+21.6+0.5;
bottom  = fP(4)-2-1.5-0.5-1.8-0.5-0.5-5.5-0.5-1.6;
set(handles.checkbox_display_data, 'Position', [left, bottom, width, height]);

% Adjust (scale position, fix size)
width   = 12;
height  = 1.8;
left    = 0.5+21.6+0.5;
bottom  = fP(4)-2-1.5-0.5-1.8-0.5-0.5-5.5-0.5-3.1;
set(handles.checkbox_display_labels, 'Position', [left, bottom, width, height]);
%}

%% Crosspeaks
%{
% Adjust (scale size, fix relative position)
left    = 0.5;
bottom  = 6.2;
width   = 35;
height  = fP(4)-22.5;
set(handles.panel_groups_in_cluster, 'Position', [left, bottom, width, height]);

% Adjust (scale size, fix relative position)
left    = 0.6;
bottom  = 6.5;
width   = 34;
height  = fP(4)-22.5-2;
set(handles.table_groups_in_displaycluster, 'Position', [left, bottom, width, height]);
%}
%}

% --- Executes during object creation, after setting all properties.
function popup_displaycluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_displaycluster.
function popup_displaycluster_Callback(hObject, eventdata, handles)
% hObject    handle to popup_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_displaycluster contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_displaycluster

refresh_display(handles);

% --- Executes on button press in button_new_displaycluster.
function button_new_displaycluster_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

%% Add a new displaycluster
displayCluster      = DisplayCluster(session);
session.addDisplayCluster(displayCluster);

% Select this new group in the list
set(handles.popup_displaycluster, 'Value', session.Ndc);
refresh_display(handles);


% --- Executes on button press in button_save_displaycluster.
function button_save_displaycluster_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

%% Update current displaycluster
% Read current group number from popup menu
dc = get(handles.popup_displaycluster, 'Value');
displayCluster = session.displayClusters{dc};

displayCluster.setName( get(handles.edit_displaycluster_name, 'String') );
displayCluster.setDescription(  get(handles.edit_displaycluster_description, 'String') );
displayCluster.setColor( get(handles.table_displaycluster_color, 'Data') );
displayCluster.setDisplayData( get(handles.checkbox_display_data, 'Value') );
displayCluster.setDisplayLabels( get(handles.checkbox_display_labels, 'Value') );

% Disable save button since it was just saved
set(handles.button_save_displaycluster, 'Enable', 'off');
set(handles.button_restore_displaycluster, 'Enable', 'off');

refresh_display(handles);

% --- Executes on button press in button_delete_displaycluster.
function button_delete_displaycluster_Callback(hObject, eventdata, handles)
% hObject    handle to button_delete_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Delete selected group

% Ask user if they want to delete fit
if( strcmp(questdlg('Are you sure you want to delete this group?', ...
    'Delete group?', 'Delete group', 'Cancel','Cancel'),'Delete group') )

    % Access main window handle which is stored in this window's main handle
    handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
    session         = getappdata(handles_main.main_gui, 'session');

    % Read current group number from popup menu
    dc = get(handles.popup_displaycluster, 'Value');
    displayCluster = session.displayClusters{dc};

    % Remove it
    session.removeDisplayCluster( displayCluster );

    % Select the first group on the list
    set(handles.popup_displaycluster, 'Value',1);
    
    refresh_display(handles);
end


function edit_selection_code_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selection_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selection_code as text
%        str2double(get(hObject,'String')) returns contents of edit_selection_code as a double


% --- Executes during object creation, after setting all properties.
function edit_selection_code_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selection_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_make_selection.
function button_make_selection_Callback(hObject, eventdata, handles)
% hObject    handle to button_make_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Make selection based on selection code
'Code me please!'
%{
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
data            = getappdata(handles_main.main_gui, 'data');
group_analysis  = getappdata(handles_main.main_gui, 'group_analysis');

%% Update current group
% Read current group number from popup menu
g       = get(handles.popup_displaycluster, 'Value');

code    = get(handles.edit_selection_code,'String');

% Create simple syntax for selecting curve sets
% N = Residue number
N       = data.index_all(1:data.Ncs);
% ASS = Atom (NH or \delta_1, etc.)
ASS     = data.ASS_all(1:data.Ncs);
% AA = 3 letter amino acid code (e.g., Trp)
for cs = 1:data.Ncs
    AA{cs} = data.AA_all(cs, 1:3);
end

try
    cs_selection = eval(code);
catch ME
    fprintf('Impossible syntax: "%s"\n', code);
end

if( length(cs_selection) == data.Ncs )
    group_analysis.group{g}.display_cs = cs_selection;
else
    fprintf('Bad syntax: "%s"\n', code);
end

% Enable save button
set(handles.button_save_displaycluster, 'Enable', 'on');
set(handles.button_restore_displaycluster, 'Enable', 'on');

setappdata(handles_main.main_gui, 'group_analysis', group_analysis);
refresh_display(handles);
%}



% --- Executes on selection change in popup_subplot_num.
function popup_subplot_num_Callback(hObject, eventdata, handles)
% hObject    handle to popup_subplot_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_subplot_num contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_subplot_num

refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function popup_subplot_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_subplot_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_displaycluster_description_Callback(hObject, eventdata, handles)
% hObject    handle to edit_displaycluster_description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_displaycluster_description as text
%        str2double(get(hObject,'String')) returns contents of edit_displaycluster_description as a double

%% Enable save button
set(handles.button_save_displaycluster, 'Enable', 'on');
set(handles.button_restore_displaycluster, 'Enable', 'on');


% --- Executes during object creation, after setting all properties.
function edit_displaycluster_description_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_displaycluster_description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%{
function edit_displaycluster_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_displaycluster_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_displaycluster_name as text
%        str2double(get(hObject,'String')) returns contents of edit_displaycluster_name as a double
%}



% --- Executes during object creation, after setting all properties.
function edit_displaycluster_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_displaycluster_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popup_plot_type.
function popup_plot_type_Callback(hObject, eventdata, handles)
% hObject    handle to popup_plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_plot_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_plot_type
%% Change subplot type for the selected subplot number
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check which plot type and subplot are selected
subplot_selected    = get(handles.popup_subplot_num, 'Value');
feature_string      = 'PLOT';
XorY_string         = 'NA';
number              = get(handles.popup_plot_type,'Value');

% Commit these changes to the data structure
paramDisplay = session.paramDisplay;
paramDisplay.setPlottingSpecs( subplot_selected, feature_string, XorY_string, number )

refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function popup_plot_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_display_data.
function checkbox_display_data_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display_data

%% Set display data boolean
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Read displaycluster number from popup
dc = get(handles.popup_displaycluster, 'Value');
displayCluster = session.displayClusters{dc};

displayCluster.setDisplayData( get(handles.checkbox_display_data, 'Value')==1 );
refresh_display(handles);

%{
% --- Executes on button press in button_scatter_plot.
function button_scatter_plot_Callback(hObject, eventdata, handles)
% hObject    handle to button_scatter_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');

% Call GUI for fitting SSE-space and send the main window handle
GUARDD_display_scatter('GUARDD', handles_main);
%}

% --- Executes on key press with focus on edit_selection_code and none of its controls.
function edit_selection_code_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_selection_code (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%eventdata
key = eventdata.Key;
modifier = eventdata.Modifier;
%isempty(modifier)

% If the user pressed return key
if( strcmp(key, 'return') && ~isempty(modifier) && strcmp(modifier, 'control') )
    fprintf('\nRun the code!');
    % Run the selection code
    button_make_selection_Callback(handles.button_make_selection,0, handles);
end


% --- Executes when entered data in editable cell(s) in table_groups_in_displaycluster.
function table_groups_in_displaycluster_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_groups_in_displaycluster (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Update current group when crosspeak is checked/unchecked

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
data            = getappdata(handles_main.main_gui, 'data');
group_analysis  = getappdata(handles_main.main_gui, 'group_analysis');


% Read current group number from popup menu
g = get(handles.popup_displaycluster, 'Value');

table_data = get(handles.table_groups_in_displaycluster, 'Data');
for cs = 1:data.Ncs
    group_analysis.group{g}.display_cs(cs) = table_data{cs,2};
end

% Enable save button
set(handles.button_save_displaycluster, 'Enable', 'on');
set(handles.button_restore_displaycluster, 'Enable', 'on');

setappdata(handles_main.main_gui, 'group_analysis', group_analysis);
refresh_display(handles);


% --- Executes on selection change in popup_param_x.
function popup_param_x_Callback(hObject, eventdata, handles)
% hObject    handle to popup_param_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_param_x contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_param_x

%% Change subplot type for the selected subplot number
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check which plot type and subplot are selected
subplot_selected    = get(handles.popup_subplot_num, 'Value');
feature_string      = 'PARAM';
XorY_string         = 'X';
number              = get(handles.popup_param_x, 'Value');

% Commit these changes to the data structure
paramDisplay = session.paramDisplay;
paramDisplay.setPlottingSpecs( subplot_selected, feature_string, XorY_string, number )

refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function popup_param_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_param_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_param_y.
function popup_param_y_Callback(hObject, eventdata, handles)
% hObject    handle to popup_param_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_param_y contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_param_y
%% Change subplot type for the selected subplot number
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check which plot type and subplot are selected
subplot_selected    = get(handles.popup_subplot_num, 'Value');
feature_string      = 'PARAM';
XorY_string         = 'Y';
number              = get(handles.popup_param_y, 'Value');

% Commit these changes to the data structure
paramDisplay = session.paramDisplay;
paramDisplay.setPlottingSpecs( subplot_selected, feature_string, XorY_string, number )

refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function popup_param_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_param_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in table_subplot.
function table_subplot_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_subplot (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Ensure table_data does not have improper values for subplotting

table_data = get(handles.table_subplot, 'Data');

% Set this to an integer > 0 by always rounding up
table_data = ceil(table_data);
% Replace values <=0 with 1
table_data( table_data<=0 ) = 1;

set(handles.table_subplot, 'Data', ceil(table_data));

refresh_display(handles);


% --------------------------------------------------------------------
function ui_save_figure_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save figure to file
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session        = getappdata(handles_main.main_gui, 'session');

% Build an intelligent filename
% Get names of selected x and y parameters
param_string    = get(handles.popup_param_x, 'String');

param_type_x    = get(handles.popup_param_x, 'Value');
param_name_x    = param_string{param_type_x};

param_type_y    = get(handles.popup_param_y, 'Value');
param_name_y    = param_string{param_type_y};

% Ask user for file to save as
default_filename     = sprintf('%s-%s_vs_%s.%s', 'GroupScatter', param_name_y, param_name_x,'');
[filename, filepath] = uiputfile('*.fig; *.png; *.ps', 'Save graphic file', ...
    sprintf('./%s/%s', session.outputDir, default_filename) );

% If the user selected a file
if( ~isequal(filename,0) )
    
    % Open new figure window
    h = figure;
    plot_groups( h, handles )
    
    % Print using GUI proportions
    %set(handles.displaycluster_gui, 'PaperPositionMode', 'auto');
    
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

%% Open new figure window
h = figure;
plot_groups( h, handles )


% --- Executes when entered data in editable cell(s) in table_displaycluster_color.
function table_displaycluster_color_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_displaycluster_color (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Enable save button
set(handles.button_save_displaycluster, 'Enable', 'on');
set(handles.button_restore_displaycluster, 'Enable', 'on');


% --- Executes on key press with focus on edit_displaycluster_name and none of its controls.
function edit_displaycluster_name_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_displaycluster_name (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%% Enable save button if user typed a key
set(handles.button_save_displaycluster, 'Enable', 'on');
set(handles.button_restore_displaycluster, 'Enable', 'on');


% --- Executes on key press with focus on edit_displaycluster_description and none of its controls.
function edit_displaycluster_description_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_displaycluster_description (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%% Enable save button if user typed a key
set(handles.button_save_displaycluster, 'Enable', 'on');
set(handles.button_restore_displaycluster, 'Enable', 'on');


% --- Executes on button press in button_restore_displaycluster.
function button_restore_displaycluster_Callback(hObject, eventdata, handles)
% hObject    handle to button_restore_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Reload all the data
refresh_display(handles);


% --- Executes on button press in checkbox_log_x.
function checkbox_log_x_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_log_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_log_x


% --- Executes on button press in checkbox_log_y.
function checkbox_log_y_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_log_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_log_y


% --- Executes when entered data in editable cell(s) in table_plot_limits.
function table_plot_limits_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_plot_limits (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Axes limits / scale has been changed
refresh_display(handles);


% --- Executes on button press in checkbox_display_labels.
function checkbox_display_labels_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display_labels

%% Set display data boolean
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Read displaycluster number from popup
dc = get(handles.popup_displaycluster, 'Value');
displayCluster = session.displayClusters{dc};

displayCluster.setDisplayLabels( get(handles.checkbox_display_labels, 'Value')==1 );
refresh_display(handles);


% --- Executes on selection change in popup_temp_x.
function popup_temp_x_Callback(hObject, eventdata, handles)
% hObject    handle to popup_temp_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_temp_x contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_temp_x
%% Change subplot type for the selected subplot number
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check which plot type and subplot are selected
subplot_selected    = get(handles.popup_subplot_num, 'Value');
feature_string      = 'TEMP';
XorY_string         = 'X';
number              = get(handles.popup_temp_x, 'Value');

% Commit these changes to the data structure
paramDisplay = session.paramDisplay;
paramDisplay.setPlottingSpecs( subplot_selected, feature_string, XorY_string, number )

refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function popup_temp_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_temp_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_B0_x.
function popup_B0_x_Callback(hObject, eventdata, handles)
% hObject    handle to popup_B0_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_B0_x contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_B0_x
%% Change subplot type for the selected subplot number
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check which plot type and subplot are selected
subplot_selected    = get(handles.popup_subplot_num, 'Value');
feature_string      = 'B0';
XorY_string         = 'X';
number              = get(handles.popup_B0_x, 'Value');

% Commit these changes to the data structure
paramDisplay = session.paramDisplay;
paramDisplay.setPlottingSpecs( subplot_selected, feature_string, XorY_string, number )

refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function popup_B0_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_B0_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_temp_y.
function popup_temp_y_Callback(hObject, eventdata, handles)
% hObject    handle to popup_temp_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_temp_y contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_temp_y
%% Change subplot type for the selected subplot number
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check which plot type and subplot are selected
subplot_selected    = get(handles.popup_subplot_num, 'Value');
feature_string      = 'TEMP';
XorY_string         = 'Y';
number              = get(handles.popup_temp_y, 'Value');

% Commit these changes to the data structure
paramDisplay = session.paramDisplay;
paramDisplay.setPlottingSpecs( subplot_selected, feature_string, XorY_string, number )

refresh_display(handles);



% --- Executes during object creation, after setting all properties.
function popup_temp_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_temp_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_B0_y.
function popup_B0_y_Callback(hObject, eventdata, handles)
% hObject    handle to popup_B0_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_B0_y contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_B0_y
%% Change subplot type for the selected subplot number
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check which plot type and subplot are selected
subplot_selected    = get(handles.popup_subplot_num, 'Value');
feature_string      = 'B0';
XorY_string         = 'Y';
number              = get(handles.popup_B0_y, 'Value');

% Commit these changes to the data structure
paramDisplay = session.paramDisplay;
paramDisplay.setPlottingSpecs( subplot_selected, feature_string, XorY_string, number )

refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function popup_B0_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_B0_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_all_groups.
function listbox_all_groups_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_all_groups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_all_groups contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_all_groups


% --- Executes during object creation, after setting all properties.
function listbox_all_groups_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_all_groups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_add_group_to_displaycluster.
function pushbutton_add_group_to_displaycluster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_group_to_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Add group to displaycluster
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Read displaycluster number from popup
dc = get(handles.popup_displaycluster, 'Value');
displayCluster = session.displayClusters{dc};

% Read current group number from popup menu
G = get(handles.listbox_all_groups, 'Value');

% For multiple items selected
for g = G    
    % Add it to the displaycluster if it is not there yet
    group = session.groups{g};
    displayCluster.addGroup(group);
end

refresh_display(handles);



% --- Executes on selection change in listbox_groups_in_displaycluster.
function listbox_groups_in_displaycluster_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_groups_in_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_groups_in_displaycluster contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_groups_in_displaycluster


% --- Executes during object creation, after setting all properties.
function listbox_groups_in_displaycluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_groups_in_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_remove_group_from_displaycluster.
function pushbutton_remove_group_from_displaycluster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_group_from_displaycluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Remove group from displaycluster
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Read displaycluster number from popup
dc = get(handles.popup_displaycluster, 'Value');
displayCluster = session.displayClusters{dc};

% Read current group number from popup menu
G = get(handles.listbox_all_groups, 'Value');

% For multiple items selected
for g = G    
    % Remove it from the displaycluster if it is not there yet
    group = session.groups{g};
    displayCluster.removeGroup(group);
end

refresh_display(handles);

% --- Executes on button press in pushbutton_remove_group_from_displaycluster_02.
function pushbutton_remove_group_from_displaycluster_02_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_group_from_displaycluster_02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Remove group from displaycluster
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Read displaycluster number from popup
dc = get(handles.popup_displaycluster, 'Value');
displayCluster = session.displayClusters{dc};

% Read current group number from popup menu
G = get(handles.listbox_groups_in_displaycluster, 'Value');

% For multiple items selected
for g = G
    % Remove it from the displaycluster if it is not there yet
    group = session.groups{g};
    displayCluster.removeGroup(group);
end

refresh_display(handles);


% --- Executes on button press in checkbox_histogram.
function checkbox_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_histogram
%% Change display histogram checkbox
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.displaycluster_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Check which plot type and subplot are selected
subplot_selected    = get(handles.popup_subplot_num, 'Value');
feature_string      = 'HISTOGRAM';
XorY_string         = 'NA';
number              = get(handles.checkbox_histogram,'Value') == 1;

% Commit these changes to the data structure
paramDisplay = session.paramDisplay;
paramDisplay.setPlottingSpecs( subplot_selected, feature_string, XorY_string, number );

refresh_display(handles);

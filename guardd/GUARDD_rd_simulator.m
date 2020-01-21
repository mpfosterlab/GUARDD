% Simulate RD curves to study nature of RD phenomena, and export for custom use (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/??/?? Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/01/14 Add simulation sets to generate temperature-dependent RD data
% 2011/04/15 Use classes for SimulationSession, SimulationSet, SimulationCurve
% 2017/04/25 Added 19F (IK)
%
% TO DO


function varargout = GUARDD_rd_simulator(varargin)
% GUARDD_RD_SIMULATOR M-file for GUARDD_rd_simulator.fig
%      GUARDD_RD_SIMULATOR, by itself, creates a new GUARDD_RD_SIMULATOR or raises the existing
%      singleton*.
%
%      H = GUARDD_RD_SIMULATOR returns the handle to a new GUARDD_RD_SIMULATOR or the handle to
%      the existing singleton*.
%
%      GUARDD_RD_SIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_RD_SIMULATOR.M with the given
%      input arguments.
%
%      GUARDD_RD_SIMULATOR('Property','Value',...) creates a new GUARDD_RD_SIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_rd_simulator_OpeningFcn gets
%      called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_rd_simulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_rd_simulator

% Last Modified by GUIDE v2.5 17-Jun-2011 11:49:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_rd_simulator_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_rd_simulator_OutputFcn, ...
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

% TODO % Tooltips


% --- Executes just before GUARDD_rd_simulator is made visible.
function GUARDD_rd_simulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_rd_simulator (see VARARGIN)

% Choose default command line output for GUARDD_rd_simulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_rd_simulator wait for user response (see UIRESUME)
% uiwait(handles.rd_simulator_gui);

% TODO % Tooltips

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.rd_simulator_gui, 'handles_main', handles_main);


%% Initialize the table for simulation parameters
column_string = {'dwH(ppm)', 'dwX(ppm)', 'AX'};
set(handles.table_cs_ppm, 'RowName', []);
set(handles.table_cs_ppm, 'ColumnWidth', {70,70,60});
set(handles.table_cs_ppm, 'ColumnName', column_string );
set(handles.table_cs_ppm, 'Data', []);

column_string = {'T0|(C)', 'PA0|(%)', 'kex0|(/s)', 'dH|(kcal/mol)', 'Eab|(kcal/mol)'};
set(handles.table_cs_kinetic_specs, 'ColumnName', column_string );
set(handles.table_cs_kinetic_specs, 'ColumnWidth', {50,50,60,70,70});
set(handles.table_cs_kinetic_specs, 'RowName', [] );
set(handles.table_cs_kinetic_specs, 'Data', []);

column_string = {'Temp(C)', 'B0(MHz)', 'R20(Hz)', 'TCPMG(ms)', 'SQX'};
set(handles.table_c_input_specs, 'RowName', []);
set(handles.table_c_input_specs, 'ColumnWidth', {60,70,70,75, 50});
set(handles.table_c_input_specs, 'ColumnName', column_string );
set(handles.table_c_input_specs, 'Data', []);
set(handles.table_c_input_specs, 'ColumnEditable', [true]);

column_string = {'dwH'; 'dwX'; 'PΑ'; 'kex'; 'R20'};
set(handles.table_c_output_params, 'RowName', column_string);
set(handles.table_c_output_params, 'ColumnWidth', {80});
set(handles.table_c_output_params, 'ColumnName', []);
set(handles.table_c_output_params, 'Data', []);

row_string = {'YMin'; 'YMax'; 'Nx'; 'Ny'};
set(handles.table_grid_limits, 'RowName', row_string);
set(handles.table_grid_limits, 'ColumnWidth', {50});
set(handles.table_grid_limits, 'ColumnName', [] );
set(handles.table_grid_limits, 'Data', [0;0;50;20]);

% Default to auto-set the name
set(handles.checkbox_autoname, 'Value', 1);


%% Initialize data structure
simulationSession = getappdata(handles_main.main_gui, 'simulationSession');

% Only if the data structure is not yet created
if( ~strcmp(class(simulationSession),'SimulationSession') )
    simulationSession = SimulationSession();

    % Save data
    setappdata(handles_main.main_gui, 'simulationSession', simulationSession);    
end

refresh_display(handles);
% TODO % Use classes for Simulationset and Simulation 
% TODO % (lots of work, only payoff is organization)


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_rd_simulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% This is called to update all display elements
function refresh_display(handles)

%% Access main window handle which is stored in this window's main handle
handles_main        = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');

% Assorted GUI stuff
set(handles.button_save_c, 'Enable', 'off');
set(handles.button_save_cs, 'Enable', 'off');

% If there are simulation sets to display
if( simulationSession.Ncs > 0 )
    
    %% GUI
    set(handles.popup_cs,  'Enable', 'on');
    set(handles.edit_cs_name,  'Enable', 'on');
    set(handles.button_delete_cs,  'Enable', 'on');
    set(handles.button_copy_cs,     'Enable', 'on');
    set(handles.button_run_kinetics_simulator,  'Enable', 'on');
    set(handles.table_cs_kinetic_specs,  'Enable', 'on');
    set(handles.table_cs_ppm,  'Enable', 'on');
    
    set(handles.button_new_c,  'Enable', 'on');
    set(handles.checkbox_autoname,  'Enable', 'on');
    
    set(handles.button_run_kinetics_simulator,  'Enable', 'on');
    
            
    %% Populate pull-down menu with sets
    
    % Just for convenience
    Ncs         = simulationSession.Ncs;
    cs_selected = get(handles.popup_cs, 'Value');
    curveset    = simulationSession.curvesets{cs_selected};
    
    % Commit to the selected set
    %set(handles.popup_cs, 'Value', simulation.cs_selected);
        
    set_names   = cell(Ncs,1);    
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

    %% Show the selected fit paramters    
    %{
    table_data = { sprintf('%0.1f C', curveset.T0-273); ...
                   sprintf('%0.1f%%', curveset.PA0 * 100); ...
                   sprintf('%0.1f /s', curveset.kex0); ...
                   sprintf('%0.1f kcal/mol', curveset.dH/1000); ...
                   sprintf('%0.1f kcal/mol', curveset.Eab/1000) };  
    %}
	table_data = [curveset.T0-273, curveset.PA0*100, curveset.kex0, curveset.dH/1000, curveset.Eab/1000];
    set(handles.table_cs_kinetic_specs, 'Data', table_data);
    
    table_data = {curveset.dwHppm, curveset.dwXppm, curveset.AX_String};    
    set(handles.table_cs_ppm, 'Data', table_data);
   
    % If there are simulations to display for the current set
    if( curveset.Nc > 0 )
        
        %% Enable GUI elements               
        set(handles.button_delete_c,  'Enable', 'on');
        set(handles.button_copy_c,  'Enable', 'on');
        set(handles.button_view_2d,     'Enable', 'on');
        set(handles.button_view_3d,     'Enable', 'on');
        set(handles.checkbox_show_legend,     'Enable', 'on');
        set(handles.popup_c,         'Enable', 'on');
        set(handles.edit_c_name,      'Enable', 'on');
        set(handles.table_c_output_params,   'Enable', 'on');
        set(handles.table_c_input_specs,      'Enable', 'on');
        %set(handles.table_input_temp, 'Enable', 'on');

        set(handles.popup_y_var,        'Enable', 'on');

        % Some GUI elements off unless otherwise turned on
        set(handles.checkbox_surface,   'Enable', 'off');
        set(handles.table_grid_limits,  'Enable', 'off');

        if( get(handles.popup_y_var, 'Value') > 1 )
            set(handles.checkbox_surface,   'Enable', 'on');

            if( get(handles.checkbox_surface, 'Value') )
                set(handles.table_grid_limits,  'Enable', 'on');
            end
        end

        % Set title of panel if there is no data to show
        set(handles.panel_plot_sim, 'Title', '');
       
        
        %% Check which simulation is selected from the list
        c_selected = get(handles.popup_c, 'Value');
        if( c_selected > curveset.Nc )
            c_selected = curveset.Nc;
            set(handles.popup_c, 'Value', c_selected)
        end
        curveset.setcSelected( c_selected );

        % Just for convenience
        Nc          = curveset.Nc;
        c_selected  = curveset.c_selected;
        curve       = curveset.curves{c_selected};

        %% Populate pull-down menu with fits
        curve_names = cell(Nc,1);
        for c = 1:Nc
            if( c == c_selected )
                prefix      = '* ';
            else
                prefix      = '';
            end

            curve_names{c} = sprintf('%s[%02d] %s', prefix, c, curveset.curves{c}.name);
        end

        set(handles.popup_c, 'String', curve_names );
        set(handles.popup_c, 'Enable', 'on');

        % Set name in edit box
        set(handles.edit_c_name, 'String', curve.name);

        %% Show the selected fit paramters
        set(handles.table_c_output_params, 'Data', ...
            { sprintf('%0.1f Hz', curve.getParamValueForDisplay(1)); ...
              sprintf('%0.1f Hz', curve.getParamValueForDisplay(2)); ...
              sprintf('%0.1f%%', curve.getParamValueForDisplay(3)); ...
              sprintf('%0.1f /sec', curve.getParamValueForDisplay(4)); ...
              sprintf('%0.1f Hz', curve.getParamValueForDisplay(5)) });
                                           
        set(handles.table_c_input_specs, 'Data', { curve.getTempC(), ...                                           
                                                   curve.B0, ...                                           
                                                   curve.params(5), ...                                           
                                                   curve.TCPMG * 1000, ...
                                                   curve.SQX } );
                                               
        %table_data = get(handles.table_input_temp, 'Data');
        %set(handles.table_input_temp, 'Data', curve.getTempC() );

        % Plot the simulated data
        plot_sim( handles.panel_plot_sim, handles);

    else
        %% There are no simulations to show

        %% GUI elements
        set(handles.button_delete_c,  'Enable', 'off');
        set(handles.checkbox_autoname,  'Enable', 'off');
        set(handles.button_copy_c,  'Enable', 'off');
        set(handles.button_view_2d,     'Enable', 'off');
        set(handles.button_view_3d,     'Enable', 'off');
        set(handles.checkbox_show_legend,     'Enable', 'off');
        set(handles.popup_c,         'Enable', 'off');
        set(handles.edit_c_name,      'Enable', 'off');
        set(handles.table_c_output_params,   'Enable', 'off');
        set(handles.table_c_input_specs,      'Enable', 'off');
        set(handles.popup_y_var,        'Enable', 'off');
        set(handles.checkbox_surface,   'Enable', 'off');
        set(handles.table_grid_limits,  'Enable', 'off');
        %set(handles.table_input_temp, 'Enable', 'off');

        %{
        % Set title of panel if there is no data to show
        set(handles.panel_plot_sim, 'Title', sprintf('No simulations to display'));

        h = subplot(1,1,1,'Parent', handles.panel_plot_sim );
        cla(h);
        hold(h, 'off');    
        legend(h, 'off');
        title(h, []);   
        %}
    end
    
% No simulation curvesets to show
else
    
    %% GUI elements
    set(handles.button_new_c,  'Enable', 'off');
    set(handles.checkbox_autoname,  'Enable', 'off');
    set(handles.button_delete_c,  'Enable', 'off');
    set(handles.button_copy_c,  'Enable', 'off');
    set(handles.button_view_2d,     'Enable', 'off');
    set(handles.button_view_3d,     'Enable', 'off');
    set(handles.checkbox_show_legend,     'Enable', 'off');
    set(handles.popup_c,         'Enable', 'off');
    set(handles.edit_c_name,      'Enable', 'off');
    set(handles.table_c_output_params,   'Enable', 'off');
    set(handles.table_c_input_specs,      'Enable', 'off');
    set(handles.popup_y_var,        'Enable', 'off');
    set(handles.checkbox_surface,   'Enable', 'off');
    set(handles.table_grid_limits,  'Enable', 'off');
%    set(handles.table_input_temp, 'Enable', 'off');
    
    set(handles.popup_cs,  'Enable', 'off');
    set(handles.edit_cs_name,  'Enable', 'off');
    set(handles.button_delete_cs,  'Enable', 'off');
    set(handles.button_copy_cs,     'Enable', 'off');
    set(handles.button_run_kinetics_simulator,  'Enable', 'off');
    set(handles.table_cs_kinetic_specs,  'Enable', 'off');       
    set(handles.table_cs_ppm,  'Enable', 'off');       
end
    
if( simulationSession.getNumCurves() == 0 )
    % Set title of panel if there is no data to show
    set(handles.panel_plot_sim, 'Title', sprintf('No simulations to display'));

    h = subplot(1,1,1,'Parent', handles.panel_plot_sim );
    cla(h);
    hold(h, 'off');    
    legend(h, 'off');
    title(h, []);    
end


function plot_sim( figure_handle, handles )

% Access main window handle which is stored in this window's main handle
handles_main        = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');
session            = getappdata(handles_main.main_gui, 'session');

%% Display 
%figure_handle = handles.panel_plot_sim;
%h = subplot(1,1,1,'Parent', figure_handle );
delete( get(figure_handle, 'Children') );
h = axes('Parent', figure_handle);
hold(h, 'all');
%cla(h);
%hold(h, 'off');

if( simulationSession.getNumCurves() > 0 )
    numPlots    = 0;
    %lstring     = cell(1,Nsims);
    lstring     = cell(1,1);


    %% Initialize y-axis types
    y_axis_string = {'SimCurveset Num', 'SimCurve Num', 'dwH', 'dwX', 'Pa', 'kex', 'R20', 'Temp'};
    set(handles.popup_y_var, 'String', y_axis_string);

    %% Y-axis points (determine which one user selected)
    % TODO % Add tempertaure to y-axis list
    y_axis_index = get(handles.popup_y_var,'Value');

    if( strcmp('SimCurveset Num', y_axis_string{y_axis_index}) )    
        Y_label = 'Curveset Number';
        s_yaxis = -1;
    
    elseif( strcmp('SimCurve Num', y_axis_string{y_axis_index}) )    
        Y_label = 'Curve Number';
        s_yaxis = 0;

    elseif( strcmp('dwH', y_axis_string{y_axis_index}) )
        Y_label = '\Delta\omega_H (Hz)';
        s_yaxis = 1;

    elseif( strcmp('dwX', y_axis_string{y_axis_index}) )
        Y_label = '\Delta\omega_X (Hz)';
        s_yaxis = 2;

    elseif( strcmp('Pa', y_axis_string{y_axis_index}) )
        Y_label = 'P_A (%)';
        s_yaxis = 3;

    elseif( strcmp('kex', y_axis_string{y_axis_index}) )
        Y_label = 'k_{ex} (/sec)';
        s_yaxis = 4;

    elseif( strcmp('R20', y_axis_string{y_axis_index}) )
        Y_label = 'R_2^0 (Hz)';
        s_yaxis = 5;

    elseif( strcmp('Temp', y_axis_string{y_axis_index}) )
        Y_label = 'Temp (^oC)';
        s_yaxis = 6;
    end

    %% X-axis points
    % TODO % Make variable limits on x-axis (vcpmg)
    vCPMG_LIM   = [1,1200];

    %% Make a surface plot    
    SHOW_SURFACE_PLOT = get(handles.checkbox_surface, 'Value');
    if( SHOW_SURFACE_PLOT && s_yaxis ~= 0 )

        % Get the selected curveset and curve
        cs_selected = get(handles.popup_cs, 'Value');
        curveset    = simulationSession.curvesets{cs_selected};

        c_selected  = get(handles.popup_c, 'Value');
        curve       = curveset.curves{c_selected};

        % Nx x Ny grid to make 3d surface
        table_data  = get(handles.table_grid_limits, 'Data');
        Nx = table_data(3);
        Ny = table_data(4);

        % Get the grid points
        xgrid = linspace( vCPMG_LIM(1), vCPMG_LIM(2), Nx);
        ygrid = linspace( table_data(1), table_data(2), Ny);

        % Take all of the simulation parameters from selected curve    
        sim_params = curve.params(1:5);

        Z = zeros(Ny,Nx);

        % Then run the y-axis variable through the grid, simulating each point
        for iy = 1:Ny        
            % If simulating temperature (#6), must get PA and kex values
            if( s_yaxis == 6 )

                Temp = ygrid(iy)+273;            
                sim_params(3) = curveset.calc_PA(  Temp );
                sim_params(4) = curveset.calc_kex( Temp );

            else        
                % Update the y-axis point to span the grid
                % (need to convert from display units to natural units)
                sim_params(s_yaxis) = ygrid(iy) / simulationSession.convert_units_to_display(s_yaxis);
            end

            % Simulate all the x-axis points
            Z(iy,1:Nx) = model_MQRD_CRJ(xgrid, curve.TCPMG, sim_params);
        end

        surf( h, xgrid, ygrid, Z );

        colormap(h, 'Bone');
        %shading(h, 'interp');

        %hold(h,'all');
        numPlots = numPlots + 1;
        % Set the text for the figure legend
        lstring{numPlots} = sprintf('Grid from %s params', char(curve.name));
    end

    %% Plot each fit in the array

    % X-Axis points
    % TODO % User-defined limits on x-axis grid
    xgrid       = linspace(vCPMG_LIM(1), vCPMG_LIM(2), 100);

    % Iterate through each simulation set
    Ncs = simulationSession.Ncs;
    for cs = 1:Ncs
        % Iterate throuch each simulation in the current set
        curveset = simulationSession.curvesets{cs};
        Nc = curveset.Nc;
        for c = 1:Nc
            curve = curveset.curves{c};

            %% Get Y-axis data points
            % -1 => SimCurveset number
            if( s_yaxis == -1 )
                Y = ones(1,length(xgrid)) .* cs;
                
            % 0 => Simulation number
            elseif( s_yaxis == 0 )
                Y = ones(1,length(xgrid)) .* c;

            % 6 => Temperature
            elseif( s_yaxis == 6 )
                Y = ones(1,length(xgrid)) .* curve.getTempC();

            else
                % Convert from natural units to display units
                Y = ones(1,length(xgrid)) .* curve.params(s_yaxis) .* simulationSession.convert_units_to_display(s_yaxis);
            end

            %% Z-axis points (R2eff simulation)
            % Parameters are already in natural units
            sim_params = curve.params;
            mR2eff = model_MQRD_CRJ(xgrid, curve.TCPMG, sim_params);

            % Plot the MQ R2eff simulation
            plot3( h, xgrid, Y, mR2eff, 'lineStyle','-', 'LineWidth',2 );

            % Overlay subsequent plots
            if( (~SHOW_SURFACE_PLOT || s_yaxis == 0) && c == 1 );
                hold(h, 'all');
            end

            % Indicate that it has been plotted
            numPlots = numPlots + 1;
            % Set the text for the figure legend
            lstring{numPlots} = sprintf('%s, %s: %0.0fHz, %0.0fHz\n%0.1f%%, %0.0f/s, %0.0fHz, %0.0fC', ...
                char(curveset.name), ...
                char(curve.name), ...
                curve.getParamValueForDisplay(1), ...
                curve.getParamValueForDisplay(2), ...
                curve.getParamValueForDisplay(3), ...
                curve.getParamValueForDisplay(4), ...
                curve.getParamValueForDisplay(5), ...
                curve.getTempC() ...
                );
        end
    end    
else
    % Set title of panel if there is no data to show
    set(handles.panel_plot_sim, 'Title', sprintf('No simulations to display'));

    %h = subplot(1,1,1,'Parent', handles.panel_plot_sim );
    %cla(h);
    %hold(h, 'off');    
    %legend(h, 'off');
    %title(h, []);       
end



%% Set plot labels, etc.
set( h, 'FontName', session.FONTNAME, 'FontSize', session.FONTSIZE_MEDIUM, ...
    'FontWeight', 'bold', 'LineWidth', session.LINEWIDTH)

%title(h, sprintf('RD Simulator'), 'FontSize', session.FONTSIZE_LARGE);
xlabel(h, '\nu_C_P_M_G (Hz)', 'FontSize', session.FONTSIZE_MEDIUM);
ylabel(h, Y_label, 'FontSize', session.FONTSIZE_MEDIUM);
zlabel(h, 'R_2^e^f^f (Hz)', 'FontSize', session.FONTSIZE_MEDIUM)

% Set axes limits
set(h, 'XLim', vCPMG_LIM);
    
    
%set( h, 'YTick', 1:Nsims );

VIEW_3D = get(handles.button_view_3d, 'Value');

if( VIEW_3D )
    view(h, 48,20);
    set( h, 'XGrid','on', 'YGrid', 'on', 'ZGrid','on' );
else    
    view(h, 0,0);
    set(h, 'Box', 'on');    
end

% Show the legend?
SHOW_LEGEND = get(handles.checkbox_show_legend, 'Value');

if( SHOW_LEGEND )
    hl=legend(h, lstring{1,:}, 'Location', 'East' );
    set(hl, 'FontName', session.FONTNAME, 'FontSize', session.FONTSIZE_SMALL-4, ...
        'FontWeight', 'Normal', 'LineWidth', session.LINEWIDTH)
    set(hl, 'Interpreter', 'Tex');
end


% --- Executes when rd_simulator_gui is resized.
function rd_simulator_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to rd_simulator_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%{
% --- Executes when entered data in editable cell(s) in table_c_output_params.
function table_c_output_params_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_c_output_params (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Modify simulation parameters

% Access main window handle which is stored in this window's main handle
handles_main        = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');
%settings        = getappdata(handles_main.main_gui, 'settings');


% Read user settings from table
% Params: dwH(Hz) dwX(Hz) Pa(%) kex(/sec) R20(Hz)
table_data      = get(handles.table_c_output_params, 'Data');
cs              = get(handles.popup_cs, 'Value');
c               = get(handles.popup_c, 'Value');

% Commit changes
% TODO % Can't just change PA or kex (they are linked via temperature)
simulationSession.curvesets{cs}.curves{c}.setParamsFromDisplay( table_data(1:5) );
simulationSession.curvesets{cs}.curves{c}.setTCPMG( table_data(6)/1000 );

% Save data
setappdata(handles_main.main_gui, 'simulationSession', simulationSession);

refresh_display(handles);
%}


% --- Executes when entered data in editable cell(s) in table_c_input_specs.
function table_c_input_specs_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_c_input_specs (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% Modify simulation parameters in ppm mode

% Access main window handle which is stored in this window's main handle
handles_main        = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');
%settings        = getappdata(handles_main.main_gui, 'settings');

% Read user settings from table
table_data      = get(handles.table_c_input_specs, 'Data');

% Params
Temp        = table_data{1}+273;
B0          = table_data{2};
R20         = table_data{3};
TCPMG       = table_data{4}/1000;
SQX         = table_data{5};

% Update
cs              = get(handles.popup_cs, 'Value');
curveset        = simulationSession.curvesets{cs};
c               = get(handles.popup_c, 'Value');
curve           = curveset.curves{c};

% Set dw, B0, TCPMG, Temp, and calculate PA and kex
curve.setSpecs( Temp, B0, R20, TCPMG, SQX );

% Automatically set name?
if( get(handles.checkbox_autoname, 'Value') == 1 )
    curve.setName();
end

% Save data
setappdata(handles_main.main_gui, 'simulationSession', simulationSession);
refresh_display(handles);


% --- Executes on selection change in popup_c.
function popup_c_Callback(hObject, eventdata, handles)
% hObject    handle to popup_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_c contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_c

%% Commit the selection

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
simulation      = getappdata(handles_main.main_gui, 'simulation');

simulation.c_selected = get(handles.popup_c, 'Value');
refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function popup_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_c (see GCBO)
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
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
session        = getappdata(handles_main.main_gui, 'session');

% Make the output directory if it does not yet exist
if( ~exist(session.outputDir, 'dir') )
    mkdir(session.outputDir);
end

% Ask user for file to save as
default_filename     = sprintf('GUARDD-Simulation.%s', '');
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


% --- Executes on button press in button_new_c.
function button_new_c_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Add new simulation to the list

% Access main window handle which is stored in this window's main handle
handles_main        = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');

% Get currently selected simulation set
cs_selected = get(handles.popup_cs, 'Value');

curve = simulationSession.curvesets{cs_selected}.newCurve();

% Automatically set name?
if( get(handles.checkbox_autoname, 'Value') == 1 )
    curve.setName();
end

% Save data
setappdata(handles_main.main_gui, 'simulationSession', simulationSession);

% Select the new simulation
set(handles.popup_c, 'Value', simulationSession.curvesets{cs_selected}.Nc);

refresh_display(handles);

% --- Executes on button press in button_delete_c.
function button_delete_c_Callback(hObject, eventdata, handles)
% hObject    handle to button_delete_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Remove simulation from list

% Ask user if they want to delete fit
if( strcmp(questdlg('Are you sure you want to delete this simulation?', ...
    'Delete simulation?', 'Delete simulation', 'Cancel','Cancel'),'Delete simulation') )

    % Access main window handle which is stored in this window's main handle
    handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
    simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');
    
    cs              = get(handles.popup_cs, 'Value');
    curveset        = simulationSession.curvesets{cs};    
    c               = get(handles.popup_c, 'Value');
    curve           = curveset.curves{c};
    
    curveset.removeCurve( curve );

    % Select a new simulation
    if( c > 1 )
        set(handles.popup_c, 'Value', c-1);
    else
        set(handles.popup_c, 'Value', 1);
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



% --- Executes on button press in button_save_c.
function button_save_c_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save the name of the simulation
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

cs = get(handles.popup_cs, 'Value');
curveset = simulationSession.curvesets{cs};

c  = get(handles.popup_c, 'Value');
curve = curveset.curves{c};

curve.setName( get(handles.edit_c_name, 'String') );

% Save it
setappdata(handles_main.main_gui, 'simulationSession', simulationSession);

refresh_display(handles);



% --- Executes on key press with focus on edit_c_name and none of its controls.
function edit_c_name_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_c_name (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(handles.button_save_c, 'Enable', 'on');


% --- Executes on selection change in popup_y_var.
function popup_y_var_Callback(hObject, eventdata, handles)
% hObject    handle to popup_y_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_y_var contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_y_var

%% Select new variable for y-axis

y_axis_string = {'SimCurveset Num', 'SimCurve Num', 'dwH', 'dwX', 'Pa', 'kex', 'R20', 'Temp'};
y_axis_index = get(handles.popup_y_var,'Value');

if( strcmp('SimCurveset Num', y_axis_string{y_axis_index}) )
    Y_LIM = [0,0];
    set(handles.checkbox_surface, 'Value', 0);
    
elseif( strcmp('SimCurve Num', y_axis_string{y_axis_index}) )
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
    
elseif( strcmp('Temp', y_axis_string{y_axis_index}) )
    Y_LIM = [0,60];
end

% Read in prior Nx and Ny values to maintain them
table_data  = get(handles.table_grid_limits, 'Data');

set(handles.table_grid_limits, 'Data', [Y_LIM(1); Y_LIM(2); table_data(3); table_data(4)]);
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
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
GUARDD_rd_simulator_export('GUARDD', handles_main);


% --- Executes on button press in button_copy_c.
function button_copy_c_Callback(hObject, eventdata, handles)
% hObject    handle to button_copy_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Copy simulation set 

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

% Select current set
cs          = get(handles.popup_cs, 'Value');
curveset    = simulationSession.curvesets{cs};

c           = get(handles.popup_c, 'Value');
curve       = curveset.curves{c};

% Copy it
curve_cp = curve.copy();
curveset.addCurve( curve_cp );

setappdata(handles_main.main_gui, 'simulationSession', simulationSession);

% Select the new curve
set(handles.popup_c, 'Value', curveset.Nc );

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
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
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


% --- Executes on button press in button_copy_cs.
function button_copy_cs_Callback(hObject, eventdata, handles)
% hObject    handle to button_copy_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Copy simulation set 

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
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


% --- Executes on button press in button_save_cs.
function button_save_cs_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save the name of the set
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

cs      = get(handles.popup_cs, 'Value');
name    = get(handles.edit_cs_name, 'String');

simulationSession.curvesets{cs}.setName( name );

% Save it
setappdata(handles_main.main_gui, 'simulationSession', simulationSession);
refresh_display(handles);


function edit_cs_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cs_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cs_name as text
%        str2double(get(hObject,'String')) returns contents of edit_cs_name as a double

set(handles.button_save_cs, 'Enable', 'on');


% --- Executes during object creation, after setting all properties.
function edit_cs_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cs_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
    handles_main        = getappdata(handles.rd_simulator_gui, 'handles_main');
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



% --- Executes on button press in button_new_cs.
function button_new_cs_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Add a new simulation set with default parameters
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

% Add new set
simulationSession.newCurveset();

% Select the new set
simulationSession.setcsSelected(simulationSession.Ncs);

setappdata(handles_main.main_gui, 'simulationSession', simulationSession);
set(handles.popup_cs, 'Value', simulationSession.Ncs);
refresh_display(handles);


% --- Executes on button press in button_run_kinetics_simulator.
function button_run_kinetics_simulator_Callback(hObject, eventdata, handles)
% hObject    handle to button_run_kinetics_simulator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Call export GUI using variables from main handle
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
GUARDD_kinetic_simulator('GUARDD', handles_main);

% --------------------------------------------------------------------
function ui_refresh_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refresh_display(handles);

% --------------------------------------------------------------------
function ui_export_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Call export GUI using variables from main handle
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
GUARDD_rd_simulator_export('GUARDD', handles_main);

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


% --- Executes when entered data in editable cell(s) in table_cs_ppm.
function table_cs_ppm_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_cs_ppm (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% Modify simulation parameters for curve set in ppm mode

% Access main window handle which is stored in this window's main handle
handles_main        = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');
%settings        = getappdata(handles_main.main_gui, 'settings');

% Read user settings from table
table_data      = get(handles.table_cs_ppm, 'Data');

% Params
dwHppm      = table_data{1};
dwXppm      = table_data{2};
AX_String   = table_data{3};
%AX_index                    = find(strcmp( AX_name, simulationSession.AX_name_array ));

% Update
cs              = get(handles.popup_cs, 'Value');
curveset        = simulationSession.curvesets{cs};

% Set dw, B0, TCPMG, Temp, and calculate PA and kex
curveset.setSpecsShifts( dwHppm, dwXppm, AX_String );

% Save data
setappdata(handles_main.main_gui, 'simulationSession', simulationSession);
refresh_display(handles);


% --- Executes on button press in checkbox_autoname.
function checkbox_autoname_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_autoname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_autoname
%% Auto-name the curve?
% Access main window handle which is stored in this window's main handle
handles_main        = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');

% Select current set
cs          = get(handles.popup_cs, 'Value');
curveset    = simulationSession.curvesets{cs};

c           = get(handles.popup_c, 'Value');
curve       = curveset.curves{c};

% Automatically set name?
if( get(handles.checkbox_autoname, 'Value') == 1 )
    curve.setName();
end
refresh_display(handles);


% --- Executes on button press in checkbox_show_legend.
function checkbox_show_legend_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_show_legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_show_legend
refresh_display(handles);


% --- Executes when entered data in editable cell(s) in table_cs_kinetic_specs.
function table_cs_kinetic_specs_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_cs_kinetic_specs (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.rd_simulator_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');

% Obtain basic kinetic specifications, which determine remaining two-site exchange parameters
table_data = get(handles.table_cs_kinetic_specs, 'Data');
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
refresh_display(handles);

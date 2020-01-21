% Display RD curve set with many view options (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/02/01 Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/06/10 Show residuals
% 2011/06/14 Residuals in a subplot

function varargout = GUARDD_display_rd(varargin)
% GUARDD_DISPLAY_RD M-file for GUARDD_display_rd.fig
%      GUARDD_DISPLAY_RD, by itself, creates a new GUARDD_DISPLAY_RD or raises the existing
%      singleton*.
%
%      H = GUARDD_DISPLAY_RD returns the handle to a new GUARDD_DISPLAY_RD or the handle to
%      the existing singleton*.
%
%      GUARDD_DISPLAY_RD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_DISPLAY_RD.M with the given input arguments.
%
%      GUARDD_DISPLAY_RD('Property','Value',...) creates a new GUARDD_DISPLAY_RD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_display_rd_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_display_rd_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_display_rd

% Last Modified by GUIDE v2.5 10-Jun-2011 13:18:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_display_rd_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_display_rd_OutputFcn, ...
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


% --- Executes just before GUARDD_display_rd is made visible.
function GUARDD_display_rd_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_display_rd (see VARARGIN)

% Choose default command line output for GUARDD_display_rd
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_display_rd wait for user response (see UIRESUME)
% uiwait(handles.display_rd_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
% AND same for fit_r2eff window which called this GUI
setappdata(handles.display_rd_gui, 'handles_main', handles_main);


%% Select the best fit on the list (initially)
%{
fit_results     = getappdata(handles_main.main_gui, 'fit_results');
cs              = get(       handles_main.listbox_cs, 'Value');

if( isstruct(fit_results) && isfield(fit_results,'list_Nfits') && ...
        fit_results.list_Nfits(cs)>0 )
    
    % Set the fit number on the list
    set(handles.popup_fit_results, 'Value', fit_results.list_iFit_best(cs));
    set(handles.popup_fit_results, 'Enable', 'on');
end

% 3D view by default
%set(handles.checkbox_3d, 'Value', 1);
%}

refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_display_rd_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function refresh_display(handles)
%% This is called to update all display elements

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.display_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

%% Set title
set(handles.panel_display_r2eff_3d, 'Title', '');...
    %sprintf('%s: Dispersion 3D plot',data.name_string{cs}));

%% List: Fit results
if( group.Ncs > 0 )
    % List curves in the group
    rownames = group.getCurveStringArray();
    set(handles.listbox_curves, 'String', rownames);
    set(handles.listbox_curves, 'Enable', 'on');
    
    % Get the fit number from the list
    ctot_Display    = get(handles.listbox_curves, 'Value');
    
    % If the value is invalid, set it to the best fit
    if( max(ctot_Display) > group.getNumCurves() )        
        set(handles.listbox_curves, 'Value', 1);
    end    
else
    set(handles.listbox_curves, 'Enable', 'off');
end

if( group.Nf > 0 )
    [fitStringArray, f_Best] = group.getFitNameArray();
    set(handles.popup_fit_results, 'String', fitStringArray );    
    set(handles.popup_fit_results, 'Enable', 'on');
    
    % Get the fit number from the list
    f = get(handles.popup_fit_results, 'Value');
    
    % If the value is invalid, set it to the best fit
    if( f > group.Nf )
        f = f_Best;
        set(handles.popup_fit_results, 'Value', f);        
    end    
else
    f = 0;
    set(handles.popup_fit_results, 'Enable', 'off');    
end

% Plot the data to the panel
plot_r2eff_3d( handles.panel_display_r2eff_3d, handles );


function plot_r2eff_3d( figure_handle, handles, varargin )
%% Plot: fit results on inset panel "panel_display_r2eff_3d"

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.display_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Size factor for top / bottom subplot
TOP_SIZE_FACTOR = 3;

% Get info on group to display
if( nargin == 4 )
    group       = varargin{1};
    fitResult   = varargin{2};
else
    g       = get(handles_main.listbox_g, 'Value');
    group   = session.groups{g};
    
    % Get the fit number from the list
    f = get(handles.popup_fit_results, 'Value');
    
    if( f <= group.Nf )
        fitResult   = group.fitResults{f};
        
    else
        % No fits to show
        fitResult   = [];
        set(handles.checkbox_show_fit, 'Value', 0);
        set(handles.checkbox_show_residuals, 'Value', 0);
    end
end

% Scale Z-axis for the group only 
Z_SCALE_FOR_GROUP   = get(handles.checkbox_y_scale, 'Value');

% View at 3D angle
VIEW_3D         = get(handles.checkbox_3d, 'Value');
    
% Plot the residuals?
SHOW_RESIDUALS = get(handles.checkbox_show_residuals, 'Value');

% Show the fit?
SHOW_FIT = get(handles.checkbox_show_fit, 'Value');

numPlots    = 0;
%lstring     = cell(1,group.Ncs);
lstring     = {};

% Get list of curves to show
ctot_Display    = get(handles.listbox_curves, 'Value');
Ncurves_Display = length(ctot_Display);
%ctot_Names      = group.getCurveStringArray();

% Clear old contents and create new axes (instead of subplot(1,1,1))
delete( get(figure_handle, 'Children') );
h = axes('Parent', figure_handle);
hold(h, 'all');

% If there are RESIDUALS, then use a subplot (data at top, resid at bottom)
if( SHOW_RESIDUALS )
    hTop = subplot(TOP_SIZE_FACTOR,1,1:TOP_SIZE_FACTOR-1,'Parent', figure_handle);
    hold(hTop, 'all');
        
    hBottom = subplot(TOP_SIZE_FACTOR,1,TOP_SIZE_FACTOR,'Parent', figure_handle);
    hold(hBottom, 'all');
end


%% Plot data
for ic = 1:Ncurves_Display
    ctot        = ctot_Display(ic);
    [cs,c]      = group.getCurvesetCurve(ctot);
    curveset    = group.curvesets{cs};
    curve       = group.curvesets{cs}.curves{c};
    
    % Indicate that it has been plotted
    numPlots = numPlots + 1;
    % Set the text for the figure legend
    if( group.Ncs > 1 )
        lstring{numPlots} = sprintf('%s[%s]', curveset.name, curve.getSpecsString());
    else
        lstring{numPlots} = sprintf('%s', curve.getSpecsString());
    end

    % Plot the data
    X = curve.vcpmg;
    Y = ones(1,curve.Nobs).*curve.Temp_C;
    Z = curve.R2eff;
    E = curve.eR2eff;
       
    % Get the proper symbol and color for this plot
    [symbolChar, colorRGB] = session.getPlotSymbolAndColor(numPlots);
    
    % If there are RESIDUALS, then use a subplot (data at top, resid at bottom)
    if( SHOW_RESIDUALS )
        h = hTop;
    end
    
    hp_data = plot3( h, X, Y, Z, symbolChar, ...
        'Color', colorRGB, ...
        'MarkerSize', session.MARKERSIZE, 'LineWidth', session.LINEWIDTH);
    %color_data = get(hp_data, 'Color');

    % Plot error bars if they are available (2010/04/08)
    % Idea from Jeremiah Faith (Aug 19, 2007)
    %  http://code.izzid.com/2007/Aug/Matlab_3D_Plot_with_Errorbars/
    if( sum(E) > 0 )
        for i=1:length(X)
            xV = [X(i); X(i)];
            yV = [Y(i); Y(i)];
            zMin = Z(i) + E(i);
            zMax = Z(i) - E(i);

            zV = [zMin, zMax];
            % draw vertical error bar but don't show it in the legend
            hp_errorbar = plot3(h, xV, yV, zV, '-', ...
                'HandleVisibility','off', 'Color', colorRGB, ...
                'LineWidth', session.LINEWIDTH);
        end
    end
    
    %% Plot residuals    
    if( SHOW_RESIDUALS )
        if( fitResult.resultsMatrix(ctot,3) ~= 1 )
            mR2eff = model_MQRD_CRJ(X, curve.TCPMG, fitResult.resultsMatrix(ctot,:));

        % No exchange looks like constant R2Eff(vCPMG)=R20
        else            
            mR2eff = fitResult.resultsMatrix(ctot,5) .* ones(1,length(X));
        end
        
        % Bottom subplot
        h = hBottom;

        % Calculate residuals for the fit (DATA - MODEL)
        ZR = Z - mR2eff;
        hp_residuals = plot3( h, X, Y, ZR, symbolChar, ...
            'HandleVisibility','off', ...
            'Color', colorRGB, 'MarkerSize', session.MARKERSIZE, 'LineWidth', session.LINEWIDTH);        
        
        % Draw a line at Z = 0 at each Y=Temperature
        X = [0 5000];
        Y = ones(1, length(X)).*curve.Temp_C;
        Z = [0 0];
        plot3( h, X, Y, Z, '--k', ...
            'HandleVisibility','off', 'LineWidth', session.LINEWIDTH);
    
    end
end

%% Plot the fits to the curves if available
% Done separately from data so the legend shows one entry for line
% and one entry for dispersion model curve
% Plot each curve from each curveset

if( group.Nf > 0 && ~isempty(fitResult) && SHOW_FIT )
    % Already read fitResult from above
    
    [void, vcpmg_max, void, void] = group.getDataLimits();
    Nx = 50;
    X = linspace(0, 1.1*vcpmg_max, Nx);
    
    for ic = 1:Ncurves_Display
        ctot    = ctot_Display(ic);    
        [cs,c]  = group.getCurvesetCurve(ctot);    
        curve = group.curvesets{cs}.curves{c};
        
        Y       = ones(1,Nx).*curve.Temp_C;

        if( fitResult.resultsMatrix(ctot,3) ~= 1 )
            mR2eff = model_MQRD_CRJ(X, curve.TCPMG, fitResult.resultsMatrix(ctot,:));
            
        % No exchange looks like constant R2Eff(vCPMG)=R20
        else            
            mR2eff = fitResult.resultsMatrix(ctot,5) .* ones(1,length(X));
        end
        
        % If there are RESIDUALS, then use a subplot (data at top, resid at bottom)
        if( SHOW_RESIDUALS )
            h = hTop;
        end

        % Plot the MQ R2eff fit with the rest of the data
        plot3( h, X, Y, mR2eff, ...
             'lineStyle','-', 'LineWidth', session.LINEWIDTH, 'Color', 'black' );
    end

    % Only show legend entry once for line and one for dispersion       
    numPlots = numPlots + 1;    
    lstring{numPlots} = sprintf('Fit' );
end


%% Set title and axes for CPMG data

% Set up the title
if( SHOW_FIT )
    titleString = sprintf('Group: %s\nFit: %s', group.name, fitResult.name);
else
    titleString = sprintf('Group: %s', group.name);
end

% If there are RESIDUALS, then use a subplot (data at top, resid at bottom)
if( SHOW_RESIDUALS )
    % Top (data + fit)    
    set(hTop, 'FontSize', session.FONTSIZE_MEDIUM, 'LineWidth', session.LINEWIDTH)
            
    % TITLE
    % If the residuals and/or fit are shown, then move the upper subplot DOWN
    % to make room above the top subplot for a multi-line title
    p = get(hTop, 'Pos');
    if( SHOW_RESIDUALS )        
        p(4) = p(4) - 0.05;
    end
    if( SHOW_FIT )
        p(4) = p(4) - 0.05;
    end
    set(hTop, 'Pos', p);
    title(hTop, titleString, 'FontSize', session.FONTSIZE_LARGE, 'FontWeight', 'Bold');    
    
    % Axes labels
    if( VIEW_3D )
        xlabel(hTop, '\nu_{CPMG} (Hz)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
    else
        set(hTop, 'XTickLabel', []);
        xlabel(hTop, '');
    end    
    ylabel(hTop, 'Temperature (C)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
    zlabel(hTop, 'R_2^{eff} (Hz)',  'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
    
    
    % Bottom (residuals)    
    set(hBottom, 'FontSize', session.FONTSIZE_MEDIUM, 'LineWidth', session.LINEWIDTH)
    %title(hBottom, titleString, 'FontSize', session.FONTSIZE_LARGE, 'FontWeight', 'Bold');
    xlabel(hBottom, '\nu_{CPMG} (Hz)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
    ylabel(hBottom, 'Temperature (C)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
    zlabel(hBottom, 'R_2^{eff} (Hz)',  'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');   
    
else
    % Only a top plot
    set(h, 'FontSize', session.FONTSIZE_MEDIUM, 'LineWidth', session.LINEWIDTH)
    title(h, titleString, 'FontSize', session.FONTSIZE_LARGE, 'FontWeight', 'Bold');
    xlabel(h, '\nu_{CPMG} (Hz)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
    ylabel(h, 'Temperature (C)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
    zlabel(h, 'R_2^{eff} (Hz)',  'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
    
    % Even if a subplot is not used, then the code below can remain
    hTop    = h;
    hBottom = h;
end        

% X-Axis and Z-Axis
% If Z-axis scale should be for this group only
if( Z_SCALE_FOR_GROUP )
    [vcpmg_min, vcpmg_max, R2eff_min, R2eff_max] = group.getDataLimits();
   
    set(hTop, 'ZLim', [ 0.9 * R2eff_min, 1.1 * R2eff_max ] );
    set(hTop, 'XLim', [0 100*ceil(vcpmg_max*1.1/100)] );
    set(hBottom, 'XLim', [0 100*ceil(vcpmg_max*1.1/100)] );
    
% Z-scale for entire session
else
    [vcpmg_min, vcpmg_max, R2eff_min, R2eff_max] = session.getDataLimits();
        
    set(hTop, 'ZLim', [0 50*ceil(R2eff_max/50)] );
    set(hTop, 'XLim', [0 100*ceil(vcpmg_max*1.1/100)] );
    set(hBottom, 'XLim', [0 100*ceil(vcpmg_max*1.1/100)] );
end

% Y-Axis
% If there is temperature data in this group
Temp_array = group.getTemperatureArray();

% Get only the temperatures which are displayed
Temp_array  = Temp_array(ctot_Display);

if( std(Temp_array) ~= 0 )
    % Set the limits at the min and max temperatures
    set( hTop, 'YLim', [min(Temp_array)-273, max(Temp_array)-273] );
    set( hTop, 'YTick', unique(Temp_array)-273 );
    
    set( hBottom, 'YLim', [min(Temp_array)-273, max(Temp_array)-273] );
    set( hBottom, 'YTick', unique(Temp_array)-273 );
else
    ylim( hTop, 'auto' );
    set( hTop, 'YTick', unique(session.getTemperatureArray())-273 );
    
    ylim( hBottom, 'auto' );
    set( hBottom, 'YTick', unique(session.getTemperatureArray())-273 );    
end

if( VIEW_3D )
    view(hTop, 48,20);
    set(hTop, 'XGrid','on', 'YGrid', 'on', 'ZGrid','on' );
    
    view(hBottom, 48,20);
    set(hBottom, 'XGrid','on', 'YGrid', 'on', 'ZGrid','on' );
else
    view(hTop, 0, 0);
    set(hTop, 'Box', 'on');    
    
    view(hBottom, 0, 0);
    set(hBottom, 'Box', 'on');
end

show_legend = get(handles.checkbox_show_legend, 'Value');
if( show_legend )
    hl=legend(hTop, lstring{1,:}, 'Location', 'NorthEast' );
    set(hl, 'FontName', session.FONTNAME, ...
        'FontSize', session.FONTSIZE_SMALL, 'FontWeight', 'Normal', 'LineWidth', session.LINEWIDTH);
end


%{
% --- Executes when display_rd_gui is resized.
function display_rd_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to display_rd_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Resize each element of the GUI

% Get size of the GUI window which is just resized
fP = get(handles.display_rd_gui, 'Position');

% Main panel (scale size, fix relative position)
left    = 0.5;
bottom  = 0.1;
width   = fP(3)-left;
height  = fP(4)-2.5;
set(handles.panel_display_r2eff_3d, 'Position', [left, bottom, width, height]);

% Adjust position of y-scale checkbox (scale position, fix size)
width   = 13;
height  = 2;
left    = 0.5;
bottom  = fP(4)-2;
set(handles.checkbox_y_scale, 'Position', [left, bottom, width, height]);

% Adjust position of y-scale checkbox (scale position, fix size)
width   = 15;
height  = 2;
left    = 0.5+13+0.5;
bottom  = fP(4)-2;
set(handles.checkbox_3d, 'Position', [left, bottom, width, height]);

% Adjust position of popup_fit_results (scale position, fix size)
width   = 50;
height  = 1.5;
left    = 1+15+1+15+1;
bottom  = fP(4)-2;
set(handles.popup_fit_results, 'Position', [left, bottom, width, height]);
%}




% --- Executes when entered data in editable cell(s) in table_topf.
function table_topf_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_topf (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


%% Get TOP_FRACTION from the table

% Read top fraction and Nbins from "table_sse_specs"
handles_main    = getappdata(handles.display_rd_gui, 'handles_main');
settings        = getappdata(handles_main.main_gui, 'settings');

% Read user settings from table and update default settings structure
table_data                      = get(handles.table_topf, 'Data');
session.SCATTER_TOP_FRACTION   = table_data(1)/100;
setappdata(handles_main.main_gui, 'settings', settings);

refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function display_rd_gui_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_rd_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popup_fit_results.
function popup_fit_results_Callback(hObject, eventdata, handles)
% hObject    handle to popup_fit_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_fit_results contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_fit_results

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


% --------------------------------------------------------------------
function ui_save_figure_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save figure to file
% Read data from main GUI
handles_main    = getappdata(handles.display_rd_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g           = get(handles_main.listbox_g, 'Value');
group       = session.groups{g};

% Get the fit number from the list
f = get(handles.popup_fit_results, 'Value');

% Ask user for file to save as
default_filename     = sprintf('%s-%s-F%02d.%s', 'FitRD', group.getFileNameSuffix(), f, '');
[filename, filepath] = uiputfile('*.fig; *.png; *.ps', 'Save graphic file', ...
    sprintf('./%s/%s', session.outputDir, default_filename) );

% If the user selected a file
if( ~isequal(filename,0) )
    
    % Make new figure window
    h = figure;

    % Plot the data to a new figure
    plot_r2eff_3d( h, handles );
    
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

%% Create popout display on new figure

% Access main window handle which is stored in this window's main handle
%handles_main    = getappdata(handles.display_rd_gui, 'handles_main');
%data            = getappdata(handles_main.main_gui, 'data');
%fit_results     = getappdata(handles_main.main_gui, 'fit_results');
%cs              = get(       handles_main.listbox_cs, 'Value');

% Get the fit number from the list
%f = get(handles.popup_fit_results, 'Value');

% Plot the data to a new figure
h = figure;
plot_r2eff_3d( h,handles  );


% --------------------------------------------------------------------
function ui_save_all_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_save_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Save best fit for every curve set

% Ask user if they want to proceed
if( strcmp(questdlg({'This will output a figure for every curve set using default filename';'Do you wish to proceed?'}, ...
        'Output all figures?', 'Proceed', 'Cancel', 'Cancel'),'Proceed') )
    
    % Access main window handle which is stored in this window's main handle
    handles_main    = getappdata(handles.display_rd_gui, 'handles_main');
    session         = getappdata(handles_main.main_gui, 'session');   
    
    % Make a new figure window
    h = figure;
    
    for g = 1:session.Ng
        group   = session.groups{g};
        fprintf('\nWorking on group %d: %s', g, group.name);
        
        % Plot the data to a new figure
        plot_r2eff_3d( h, handles, group, group.fitResult_Best );

        % Save to file
        default_filename = sprintf('%s-%s-FBest.%s', 'fit-R2eff', group.getFileNameSuffix(),'ps');
        print( h, '-dpsc', sprintf('%s/%s',session.outputDir,default_filename ) );        

    end
    close(h);
end


% --- Executes on button press in checkbox_y_scale.
function checkbox_y_scale_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_y_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_y_scale

%% Change scaling on Y-Axis for R2eff
refresh_display(handles);


% --- Executes on button press in checkbox_3d.
function checkbox_3d_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_3d

refresh_display(handles);


% --- Executes on selection change in listbox_curves.
function listbox_curves_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_curves contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_curves
refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function listbox_curves_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_show_fit.
function checkbox_show_fit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_show_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_show_fit
refresh_display(handles);


% --- Executes on button press in checkbox_show_legend.
function checkbox_show_legend_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_show_legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_show_legend
refresh_display(handles);


% --- Executes on button press in checkbox_show_residuals.
function checkbox_show_residuals_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_show_residuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_show_residuals
refresh_display(handles);

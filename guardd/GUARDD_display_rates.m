% Display rate analyses for selected curve set (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/01/20 Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/04/24 Use classes for data
% 2011/06/07 RateAnalysis uses cal instead of kcal


function varargout = GUARDD_display_rates(varargin)
% GUARDD_DISPLAY_RATES M-file for GUARDD_display_rates.fig
%      GUARDD_DISPLAY_RATES, by itself, creates a new GUARDD_DISPLAY_RATES or raises the existing
%      singleton*.
%
%      H = GUARDD_DISPLAY_RATES returns the handle to a new GUARDD_DISPLAY_RATES or the handle to
%      the existing singleton*.
%
%      GUARDD_DISPLAY_RATES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_DISPLAY_RATES.M with the given input arguments.
%
%      GUARDD_DISPLAY_RATES('Property','Value',...) creates a new GUARDD_DISPLAY_RATES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_display_rates_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_display_rates_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_display_rates

% Last Modified by GUIDE v2.5 11-Jan-2011 15:00:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_display_rates_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_display_rates_OutputFcn, ...
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


% --- Executes just before GUARDD_display_rates is made visible.
function GUARDD_display_rates_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_display_rates (see VARARGIN)

% Choose default command line output for GUARDD_display_rates
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_display_rates wait for user response (see UIRESUME)
% uiwait(handles.display_rates_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.display_rates_gui, 'handles_main', handles_main);

% Read data from main GUI
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

% If there are some fits
if( group.Nf > 0 )
    % Set the list of fit results
    [fitStringArray, f_Best] = group.getFitNameArray();    
    set(handles.popup_fit_results, 'String', fitStringArray);
    set(handles.popup_fit_results, 'Value', f_Best);
    set(handles.popup_fit_results, 'Enable', 'on');
end

refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_display_rates_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% This is called to update all display elements
function refresh_display(handles)
    
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.display_rates_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};

Nctot   = group.getNumCurves();
Nf      = group.Nf;

%% List: Fit results
if( Nf > 0 )    
    % Set the list of fit results
    [fitStringArray, f_Best] = group.getFitNameArray();    
    set(handles.popup_fit_results, 'String', fitStringArray );
    set(handles.popup_fit_results, 'Enable', 'on');
    
    % Get the fit number from the list
    f = get(handles.popup_fit_results, 'Value');
    
    % If the value is invalid, set it to the best one on the list
    if( f > Nf )
        set(handles.popup_fit_results, 'Value', f_Best);
    end
    
    % Display data and fit results
    plot_rate_analysis( handles.panel_display_rate_analysis, handles );

else
    set(handles.popup_fit_results, 'String', {} );
    set(handles.popup_fit_results, 'Enable', 'off');
    
    % Set title of panel if there is no data to show
    set(handles.panel_display_rate_analysis, 'Title', ...
        sprintf('%s: No rate analysis to display', group.name));
end


function plot_rate_analysis( figure_handle, handles )
% Read data from main GUI
handles_main    = getappdata(handles.display_rates_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g           = get(handles_main.listbox_g, 'Value');
group       = session.groups{g};

% Get the fit number from the list
f           = get(handles.popup_fit_results, 'Value');
fitResult   = group.fitResults{f};
rateAnalysis= fitResult.rateAnalysis;
R           = session.R;

% Clear the axes in figure
h = subplot(1,1,1,'Parent', figure_handle);
cla(h);
hold(h, 'off');

% To scale the x-axis units
X_SCALE = 1000;

%% Top (1/2 subplot) is Arrhenius: ln(k) = ln(P) + (-Ea/R)*(1/T)
h = subplot(2,1,1, 'Parent', figure_handle);
set(h, 'FontSize', session.FONTSIZE_SMALL);
hold(h,'on');

% Get the data from the fitResult
[X_ok, Y_ok, Y_ok_E, X_all, Y_all, Y_all_E] = rateAnalysis.getArrheniusPlotA();

% Plot ALL data in gray with errors A->B
hp = errorbar(h, X_SCALE .* X_all, Y_all, Y_all_E, ...
    'LineWidth', session.LINEWIDTH, 'LineStyle', 'none', ...
    'Marker', '>', 'MarkerSize', session.MARKERSIZE, 'Color',[0.5 0.5 0.5] );
% Exclude line from legend
set(get(get(hp,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

% Plot OK data in black with errors A->B
errorbar(h, X_SCALE .* X_ok, Y_ok, Y_ok_E, ...
    'LineWidth', session.LINEWIDTH, 'LineStyle', 'none', ...
    'Marker', '>', 'MarkerSize', session.MARKERSIZE, 'Color',[0 0 0] );


% Get the data from the fitResult
[X_ok, Y_ok, Y_ok_E, X_all, Y_all, Y_all_E] = rateAnalysis.getArrheniusPlotB();

% Plot ALL data in red with errors B->A
hp = errorbar(h, X_SCALE .* X_all, Y_all, Y_all_E, ...
    'LineWidth', session.LINEWIDTH, 'LineStyle', 'none', ...
    'Marker', '<', 'MarkerSize', session.MARKERSIZE, 'Color',[1.0 0.5 0.5] );
% Exclude line from legend
set(get(get(hp,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

% Plot OK data in red with errors B->A
errorbar(h, X_SCALE .* X_ok, Y_ok, Y_ok_E, ...
    'LineWidth', session.LINEWIDTH, 'LineStyle', 'none', ...
    'Marker', '<', 'MarkerSize', session.MARKERSIZE, 'Color',[1 0 0] );

% Only show the legend if there is some data
if( ~all(isinf(Y_all)) )
    legend( {'A\rightarrowB', 'B\rightarrowA'}, 'FontSize', session.FONTSIZE_SMALL );
end


if( rateAnalysis.arrhenius_isOK )
    % Plot the fit on top A->B if it is OK
    xgrid = X_SCALE .* [0 1/200];
    ypoints = log(rateAnalysis.Pab) - 1*rateAnalysis.Eab/R .* [0 1/200];
    plot(h, xgrid, ypoints, 'LineStyle', '--', 'Color', 'k', 'LineWidth', session.LINEWIDTH);

    % Plot the fit on top B->A if it is OK
    xgrid = X_SCALE .* [0 1/200];
    ypoints = log(rateAnalysis.Pba) - 1*rateAnalysis.Eba/R .* [0 1/200];
    plot(h, xgrid, ypoints, 'LineStyle', '--', 'Color', 'r', 'LineWidth', session.LINEWIDTH);

    % Set sub title
    Eab     = Session.convertUnits(rateAnalysis.Eab, 'EAB', 'DISPLAY');
    Eab_E   = Session.convertUnits(rateAnalysis.Eab_E, 'EAB', 'DISPLAY');
    Pab     = Session.convertUnits(rateAnalysis.Pab, 'PAB', 'DISPLAY');
    Eba     = Session.convertUnits(rateAnalysis.Eba, 'EBA', 'DISPLAY');
    Eba_E   = Session.convertUnits(rateAnalysis.Eba_E, 'EBA', 'DISPLAY');
    Pba     = Session.convertUnits(rateAnalysis.Pba, 'PBA', 'DISPLAY');
    sub_title = sprintf('(>) E^{A\\rightarrowB}=%0.1f\\pm%0.1f kcal/mol, P=%0.1e /sec\n(<) E^{B\\rightarrowA}=%0.1f\\pm%0.1f kcal/mol, P=%0.1e /sec', ...
        Eab, Eab_E, Pab, Eba, Eba_E, Pba);

else
    sub_title = 'Fit is not OK';
end


xlim(h, X_SCALE .* [0.98 1.02] .* [ min(X_all) max(X_all) ] );
%ylim( [ rate_analysis.plot_min_ln_K rate_analysis.plot_max_ln_K ] );
ylabel(h, 'Ln(k)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
%set(h,'XGrid', 'on');
box(h, 'on');

%title(h, sprintf('\nArrhenius analysis k_A(o) and k_B([])\n%s', sub_title), ...
title(h, sprintf('Arrhenius analysis\n%s', sub_title), ...
    'FontSize',session.FONTSIZE_SMALL, 'FontWeight', 'Normal' );

%% Bottom (1/2 subplot) is van't Hoff: ln(K) = dS/R + (-dH/R)*(1/T)
h = subplot(2,1,2, 'Parent', figure_handle);
set(h, 'FontSize', session.FONTSIZE_SMALL);
hold(h,'on');

% Get the data from the fitResult
[X_ok, Y_ok, Y_ok_E, X_all, Y_all, Y_all_E] = rateAnalysis.getVantHoffPlot();

% Plot ALL data in gray with errors A->B
hp = errorbar(h, X_SCALE .* X_all, Y_all, Y_all_E, ...
    'LineWidth', session.LINEWIDTH, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', session.MARKERSIZE, 'Color',[0.5 0.5 0.5] );
% Exclude line from legend
set(get(get(hp,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

% Plot OK data in black with errors A->B
errorbar(h, X_SCALE .* X_ok, Y_ok, Y_ok_E, ...
    'LineWidth', session.LINEWIDTH, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', session.MARKERSIZE, 'Color',[0 0 0] );

% Plot the fit on top of the data already there IF its OK
if( rateAnalysis.vantHoff_isOK)
    xgrid = X_SCALE .* [0 1/200];
    ypoints = rateAnalysis.dS/R - 1*rateAnalysis.dH/R .* [0 1/200];
    plot(h, xgrid, ypoints, 'LineStyle', '--', 'Color', 'k', 'LineWidth', session.LINEWIDTH);

    
    dH     = Session.convertUnits(rateAnalysis.dH, 'DH', 'DISPLAY');
    dH_E   = Session.convertUnits(rateAnalysis.dH_E, 'DH', 'DISPLAY');
    dS     = Session.convertUnits(rateAnalysis.dS, 'DS', 'DISPLAY');
    dS_E   = Session.convertUnits(rateAnalysis.dS_E, 'DS', 'DISPLAY');
    sub_title = sprintf('\\DeltaH=%+0.1f\\pm%0.1f kcal/mol, \\DeltaS=%+0.1f\\pm%0.1f cal/mol/K', ...
        dH, dH_E, dS, dS_E );

else
    sub_title = 'Fit is not OK';
end

xlim(h, X_SCALE .* [0.98 1.02] .* [ min(X_all) max(X_all) ] );
xlabel(h, sprintf('%d / Temperature (1/K)', X_SCALE), 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold' );
ylabel(h, 'Ln(K)', 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold' );
%set(h,'XGrid', 'on');
box(h, 'on');

title(h, sprintf('van\''t Hoff analysis\n%s', sub_title), ...
    'FontSize',session.FONTSIZE_SMALL, 'FontWeight', 'Normal' );

%% Set main title
title_string = sprintf('%s', group.name);
% Set title of panel if drawing current panel
if( figure_handle == handles.panel_display_rate_analysis )
    set(handles.panel_display_rate_analysis, 'Title', title_string );

% Set super-title of subplots if using new figure window
else
    % Position [left, bottom, width, height]
    % Empirically determined for this plot
    ax_position = [0.08, 0.08, 0.84, 0.85];
    ax=axes('Units','Normal','Position',ax_position,'Visible','off', 'Parent', figure_handle);
    set(get(ax,'Title'),'Visible','on')        
    title(ax, title_string, 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
end


% --- Executes when display_rates_gui is resized.
function display_rates_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to display_rates_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Resize each element of the GUI

% Get size of the GUI window which is just resized
fP = get(handles.display_rates_gui, 'Position');

left    = 0.5;
bottom  = 0.1;
width   = fP(3)-left;
height  = fP(4)-2.5;
set(handles.panel_display_rate_analysis, 'Position', [left, bottom, width, height]);

% Adjust position of popup_fit_results (scale position, fix size)
width   = 50;
height  = 1.5;
left    = fP(3)/2 - width/2;
bottom  = fP(4)-2;
set(handles.popup_fit_results, 'Position', [left, bottom, width, height]);


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
function ui_popout_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_popout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Display plot in new figure window

% Read data from main GUI
handles_main    = getappdata(handles.display_rates_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g           = get(handles_main.listbox_g, 'Value');
group       = session.groups{g};

% Get the fit number from the list
f           = get(handles.popup_fit_results, 'Value');

% Check if there is information to get on this fit
if( f <= group.Nf )    
    % Get the fit number from the list
    h = figure;
    plot_rate_analysis( h, handles )
end

% --------------------------------------------------------------------
function ui_save_figure_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save figure to file
% Read data from main GUI
handles_main    = getappdata(handles.display_rates_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g           = get(handles_main.listbox_g, 'Value');
group       = session.groups{g};

% Get the fit number from the list
f = get(handles.popup_fit_results, 'Value');

% Check if there is information to get on this fit
if( f <= group.Nf )    
    
    

    % Ask user for file to save as
    default_filename     = sprintf('%s-%s-F%02d.%s', 'RateAnalysis', group.getFileNameSuffix(), f, '');
    [filename, filepath] = uiputfile('*.fig; *.png; *.ps', 'Save graphic file', ...
        sprintf('./%s/%s', session.outputDir, default_filename) );

    % If the user selected a file
    if( ~isequal(filename,0) )

        % Make new figure window
        h = figure;

        % Plot the data to a new figure
        plot_rate_analysis( h, handles );

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
end

% --------------------------------------------------------------------
function ui_save_all_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_save_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Save best fit for every curve set
%{
% Ask user if they want to proceed
if( strcmp(questdlg({'This will output a figure for every curve set using default filename';'Do you wish to proceed?'}, ...
        'Output all figures?', 'Proceed', 'Cancel', 'Cancel'),'Proceed') )
    
    % Access main window handle which is stored in this window's main handle
    handles_main    = getappdata(handles.display_rates_gui, 'handles_main');
    data            = getappdata(handles_main.main_gui, 'data');
    settings        = getappdata(handles_main.main_gui, 'settings');
    fit_results     = getappdata(handles_main.main_gui, 'fit_results');

    % Make a new figure window
    h = figure;
    for cs = 1:data.Ncs
        % Check if the best fit has data to show
        f = fit_results.list_iFit_best(cs);  
        
        if( f>0 && fit_results.analysis_is_available(cs,f) )
            fprintf('\nWorking on cs %d: %s', cs, data.name_string{cs});
            
            % Plot the data to a new figure
            plot_rate_analysis( h, handles, cs, f );
            
            % Save to file
            default_filename     = sprintf('%s-%s-F%02d.%s', 'rate_analysis', data.filename_suffix_cs{cs}, f, 'ps');
            print( h, '-dpsc', sprintf('%s/%s',session.outputDir,default_filename ) );     
        end
    end
    close(h);
end
%}

% Display chi2 map for selected curve set (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/02/01 Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/04/22 Convert to use of classes
% 2011/09/02 (dwX, PaPb) surface plot (overhead view)

function varargout = GUARDD_display_chi2_map(varargin)
% GUARDD_DISPLAY_CHI2_MAP M-file for GUARDD_display_chi2_map.fig
%      GUARDD_DISPLAY_CHI2_MAP, by itself, creates a new GUARDD_DISPLAY_CHI2_MAP or raises the existing
%      singleton*.
%
%      H = GUARDD_DISPLAY_CHI2_MAP returns the handle to a new GUARDD_DISPLAY_CHI2_MAP or the handle to
%      the existing singleton*.
%
%      GUARDD_DISPLAY_CHI2_MAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_DISPLAY_CHI2_MAP.M with the given input arguments.
%
%      GUARDD_DISPLAY_CHI2_MAP('Property','Value',...) creates a new GUARDD_DISPLAY_CHI2_MAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_display_chi2_map_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_display_chi2_map_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_display_chi2_map

% Last Modified by GUIDE v2.5 03-Jun-2011 18:37:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_display_chi2_map_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_display_chi2_map_OutputFcn, ...
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


% --- Executes just before GUARDD_display_chi2_map is made visible.
function GUARDD_display_chi2_map_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_display_chi2_map (see VARARGIN)

% Choose default command line output for GUARDD_display_chi2_map
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_display_chi2_map wait for user response (see UIRESUME)
% uiwait(handles.display_chi2_map_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.display_chi2_map_gui, 'handles_main', handles_main);

%% Populate listbox with candidate parameters to show
params_string = {'dwH', 'dwX', 'Pa', 'kex', 'R20', 'Rex', 'alpha', 'PhiexX', 'dH', 'Eab', 'dwX v PaPb'};
set(handles.listbox_params, 'String', params_string);
%set(handles.listbox_params, 'Value', [1 2 3 4 5 6 7 8]);

% Top fraction table
set(handles.table_topf, 'RowName', 'Top %');
set(handles.table_topf, 'ColumnName', []);
set(handles.table_topf, 'ColumnWidth', {40});

%% List of grid fits

% Do not pick value from grid
set(handles.checkbox_select_fit_from_grid, 'Value', 0);
f_Selected_Top = [];
setappdata(handles.display_chi2_map_gui, 'f_Selected_Top', f_Selected_Top);

%set(handles.table_grid_fits, 'ColumnName', 'chi2');
set(handles.table_grid_fits, 'RowName', 'chi2');
set(handles.table_grid_fits, 'ColumnWidth', {100});

%% Table for top fraction
session = getappdata(handles_main.main_gui, 'session');
set(handles.slider_topf, 'Value', session.chi2TopFraction);
set(handles.table_topf, 'Data', session.chi2TopFraction*100);

refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_display_chi2_map_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% This is called to update all display elements
function refresh_display(handles)
% TODO % Display list of fits from grid, and select each one for display
% TODO % and allow adding fit as an "official" listed fit

% Read data from main GUI
handles_main    = getappdata(handles.display_chi2_map_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};
Nf      = group.Nf;

DISPLAY_MC_ERRORS = get(handles.radiobutton_display_MC_errors, 'Value') == 1;
SELECT_GRID_FIT = get(handles.checkbox_select_fit_from_grid, 'Value') == 1;

%% Show table for grid fit results
if( SELECT_GRID_FIT && ~DISPLAY_MC_ERRORS && group.gridSearchIsDone() )

    % Display the inital or final conditions?
    if( get(handles.radiobutton_initial_conditions, 'Value') == 1 )
        INITIAL_OR_FINAL = 'Initial';
    else
        INITIAL_OR_FINAL = 'Final';
    end

    % Get info on the grid
    [chi2_Grid, VOID, VOID, VOID] = group.getGridResults(INITIAL_OR_FINAL);
    Ngrid = length(chi2_Grid);

    % Current parameters in the top fraction of fits
    Ntop   = ceil(session.chi2TopFraction * Ngrid);

    % Sort the chi2 array
    [chi2_Grid_Sorted, igrid_Sorted]  = sort( chi2_Grid );

    % Assign SSE values out of top fraction for current (cs,c,p)
    chi2_Grid_Top = chi2_Grid_Sorted(1:Ntop); 

    % Build table from list
    set(handles.table_grid_fits, 'Data', chi2_Grid_Top');
    
    set(handles.table_grid_fits, 'Visible', 'on');
    set(handles.pushbutton_add_grid_fit_to_list, 'Visible', 'on');
    
    set(handles.panel_display_map_large, 'Visible', 'off');
    set(handles.panel_display_map_small, 'Visible', 'on');
    display_panel_handle = handles.panel_display_map_small;
    
else
    
    set(handles.table_grid_fits, 'Visible', 'off');
    set(handles.pushbutton_add_grid_fit_to_list, 'Visible', 'off');
    
    set(handles.panel_display_map_large, 'Visible', 'on');
    set(handles.panel_display_map_small, 'Visible', 'off');
    display_panel_handle = handles.panel_display_map_large;
end


% Clear the axes
h = subplot(1,1,1,'Parent', display_panel_handle);
cla(h);

if( Nf > 0 )
    %% Check if there is information to show
    % Parameters selected for display
    set(handles.listbox_params, 'Value', session.pChi2Map);
    set(handles.listbox_params, 'Enable', 'on');

    % List curves in the group
    rownames = group.getCurveStringArray();
    set(handles.listbox_curves, 'String', rownames);
    set(handles.listbox_curves, 'Enable', 'on');
    
    % Make sure the selected curves are valid
    c_selected = get(handles.listbox_curves, 'Value');
    if( max(c_selected) > group.getNumCurves() )
        % ...Invalid selection => just pick the first one
        set(handles.listbox_curves, 'Value', 1);
    end
    
    % Set the list of fit results
    [fitStringArray, f_Best] = group.getFitNameArray();
    set(handles.popup_fit_results, 'String', fitStringArray );  
    set(handles.popup_fit_results, 'Enable', 'on');
    
    % Get the fit number from the list
    f = get(handles.popup_fit_results, 'Value');
    
    % If the value is invalid, set it to the best one on the list
    if( f > Nf )
        f = f_Best;
        set(handles.popup_fit_results, 'Value', f);
    end
    
    [fitResult_Selected, isBestFit]  = group.getFitResult(f);
    
    if( DISPLAY_MC_ERRORS )
        %% Show Monte Carlo error results
        % Alter GUI elements
        
        % If there are MC results to show
        if( ~isempty(fitResult_Selected.fitResults_MC) )
            set(handles.table_topf, 'Enable', 'on');
            set(handles.slider_topf, 'Enable', 'on');
            
            % Display data and fit results            
            plot_chi2_map( display_panel_handle, handles );
            
        else
            set(handles.listbox_params, 'Enable', 'off');
            set(handles.listbox_curves, 'Enable', 'off');
            set(handles.table_topf, 'Enable', 'off');
            set(handles.slider_topf, 'Enable', 'off');
            
            set(display_panel_handle, 'Title', ...
                sprintf('%s: No Monte Carlo results\n(perform via Fitting Window)', group.name));            
        end
        
    else
        %% Show grid search results
        if( group.gridSearchIsDone() )            
            % Set the table top fraction
            set(handles.table_topf, 'Data', 100*session.chi2TopFraction);
            set(handles.table_topf, 'Enable', 'on');
            set(handles.slider_topf, 'Enable', 'on');
            
            set(handles.checkbox_select_fit_from_grid, 'Enable', 'on');
        
            % Display data and fit results
            plot_chi2_map( display_panel_handle, handles );
            
        else
            set(handles.popup_fit_results, 'Enable', 'off');
            set(handles.listbox_params, 'Enable', 'off');
            set(handles.listbox_curves, 'Enable', 'off');
            set(handles.table_topf, 'Enable', 'off');
            set(handles.slider_topf, 'Enable', 'off');
            set(handles.checkbox_select_fit_from_grid, 'Enable', 'off');

            % Set title of panel if there is no data to show
            set(display_panel_handle, 'Title', ...
                sprintf('%s: No grid search\n(perform via Fitting Window)', group.name));
        end
    end
    
else
    %% There are no fit results to show, disable GUI elements
    set(handles.popup_fit_results, 'Enable', 'off');
    set(handles.listbox_params, 'Enable', 'off');
    set(handles.listbox_curves, 'Enable', 'off');
    set(handles.table_topf, 'Enable', 'off');
    set(handles.slider_topf, 'Enable', 'off');
    
    % Set title of panel if there is no data to show
    set(display_panel_handle, 'Title', ...
        sprintf('%s: No fits\n(perform via Fitting Window)', group.name));
end



function plot_chi2_map( figure_handle, handles )
% TODO % Plot the data and fit here too
% Read data from main GUI
handles_main    = getappdata(handles.display_chi2_map_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected group
g       = get(handles_main.listbox_g, 'Value');
group   = session.groups{g};
Nf      = group.Nf;

%% Get info on selected fit
f = get(handles.popup_fit_results, 'Value');
[fitResult_Selected, isBestFit]  = group.getFitResult(f);

% Display Monte Carlo errors? (Or grid search?)
DISPLAY_MC_ERRORS = get(handles.radiobutton_display_MC_errors, 'Value') == 1;

% Display the inital or final conditions?
if( get(handles.radiobutton_initial_conditions, 'Value') == 1 )
    INITIAL_OR_FINAL = 'Initial';
else
    INITIAL_OR_FINAL = 'Final';
end

% Display MC errors as histogram (or chi2 map?)
DISPLAY_HISTOGRAM = get(handles.radiobutton_histogram, 'Value') == 1;

if( DISPLAY_MC_ERRORS )
    %% Get info on Monte Carlo bootstrapping error estimation    
       
    [chi2_MC, resultsMatrix_MC, dH_MC, Eab_MC] = fitResult_Selected.getMCResults(INITIAL_OR_FINAL);
    Nmc = length(chi2_MC);
    
    % Current parameters in the top fraction of fits
    Ntop   = ceil(session.chi2TopFraction * Nmc);

    % Sort the chi2 array    
    [chi2_MC_Sorted, iMC_Sorted]  = sort( chi2_MC );

    % Assign SSE values out of top fraction for current (cs,c,p)
    chi2_MC_Top = chi2_MC_Sorted(1:Ntop);  
    
    % Set title of panel
    if( figure_handle == handles.panel_display_map_small || ...
        figure_handle == handles.panel_display_map_large )
        set(figure_handle, 'Title', ...
            sprintf('%s: Monte Carlo errors (best %d/%d)', group.name, Ntop, Nmc));
    end
    
    % Turn off GUI elements for grid results table
    
else
    %% Get info on the grid
    [chi2_Grid, resultsMatrix_Grid, dH_Grid, Eab_Grid] = group.getGridResults(INITIAL_OR_FINAL);
    Ngrid = length(chi2_Grid);

    % Current parameters in the top fraction of fits
    Ntop   = ceil(session.chi2TopFraction * Ngrid);

    % Sort the chi2 array
    [chi2_Grid_Sorted, igrid_Sorted]  = sort( chi2_Grid );

    % Assign SSE values out of top fraction for current (cs,c,p)
    chi2_Grid_Top = chi2_Grid_Sorted(1:Ntop);  
    
    % Set title of panel
    if( figure_handle == handles.panel_display_map_small || ...
        figure_handle == handles.panel_display_map_large )
        set(figure_handle, 'Title', ...
            sprintf('%s: Grid search (best %d/%d)', group.name, Ntop, Ngrid));
    end
    
    % Turn on GUI elements for grid results table
end

%% Display the fit from the grid, as selected in the list
DISPLAY_GRID_FIT = get(handles.checkbox_select_fit_from_grid, 'Value') == 1;

if( ~DISPLAY_MC_ERRORS && DISPLAY_GRID_FIT )
    f_Selected_Top      = getappdata(handles.display_chi2_map_gui, 'f_Selected_Top');
    
    % Repair an invalid selection
    if( isempty(f_Selected_Top) || f_Selected_Top > Ntop )
        f_Selected_Top = 1;
        setappdata(handles.display_chi2_map_gui, 'f_Selected_Top', f_Selected_Top);
    end
    
    %Nf_Selected_List    = length(f_Selected_Top);
    
    % Get the fit index from the unsorted list
    f_Selected_Grid         = igrid_Sorted( f_Selected_Top );
    fitResult_Selected_List = group.fitResults_Grid{f_Selected_Grid};
        
    resultsMatrix_Selected_List  = fitResult_Selected_List.resultsMatrix;
    
else
    fitResult_Selected_List     = [];
    resultsMatrix_Selected_List = [];
    f_Selected_Top             = NaN;
end


%% Get list of parameters to show
params_string = {'dwH', 'dwX', 'Pa', 'kex', 'R20', 'Rex', 'alpha', 'PhiexX', 'dH', 'Eab', 'dwX v PaPb'};
param_index_selected = get(handles.listbox_params, 'Value');

% Get list of curves to show
ctot_Display    = get(handles.listbox_curves, 'Value');
ctot_Names      = group.getCurveStringArray();

%% Display data
% Clear the axes
h = subplot(1,1,1,'Parent', figure_handle);
cla(h);
hold(h, 'off');

% For the subplots
Nrows = length(param_index_selected);
Ncols = length(ctot_Display);

% Iterate through each selected curve (COLUMN)
for ic = 1:Ncols
    ctot    = ctot_Display(ic);    
    [cs,c]  = group.getCurvesetCurve(ctot);
        
    % Iterate through each seleted paramter (ROW)
    for ip = 1:Nrows        
        p = param_index_selected(ip);
        
        % Subplot numbering (e.g., 2 rows, 3 cols)
        %
        %   1  |  2  |  3
        %  -----------------
        %   4  |  5  |  6
        %   
        subplot_number  = (ip-1)*Ncols + ic;
        h = subplot(Nrows, Ncols, subplot_number, 'Parent', figure_handle );
        set(h, 'FontSize', session.FONTSIZE_SMALL);
        hold(h,'on');
        
        %fprintf('\nSubplot=%d/%d (%d, %d)', subplot_number, Nrows*Ncols, ip, ic);

        % If it is the left-most subplot (left column)        
        if( ic == 1 )
            if( strcmpi(params_string{p}, 'dwX v PaPb') )
                ylabel('\Delta\omega_X^2', 'FontWeight', 'Bold')
            else
                ylabel('\chi^2', 'FontWeight', 'Bold')
            end
            %set(h,'YTick', []);
        else
            ylabel('');
            set(h,'YTick', []);
        end

        % Only output the title for the parameter at the top
        s_title = '';
        if( ip == 1 )
            s_title = sprintf('%s\n', ctot_Names{ctot});
        end        
        
        if( DISPLAY_MC_ERRORS )
            % Get X-values depending on what is selected
            switch( params_string{p} )
                case {'dwH', 'dwX', 'Pa', 'kex', 'R20', 'Rex', 'alpha', 'PhiexX'}
                    % Get the data from the Monte Carlo results
                    p_MC_Top        = NaN*ones(Ntop,1);
                    p_MC_Top(1:Ntop) = resultsMatrix_MC(ctot, p, iMC_Sorted(1:Ntop));

                    % Plot
                    X = session.convertParameterToDisplayUnits(p_MC_Top,p);
                    
                case 'dH'
                    p_MC_Top        = NaN*ones(Ntop,1);
                    p_MC_Top(1:Ntop) = dH_MC(iMC_Sorted(1:Ntop));
                    X               = p_MC_Top ./ 1e3;
                    
                case 'Eab'
                    p_MC_Top        = NaN*ones(Ntop,1);
                    p_MC_Top(1:Ntop) = Eab_MC(iMC_Sorted(1:Ntop));
                    X               = p_MC_Top ./ 1e3;
                    
                otherwise
                    % kex vs. pa
                    %error('Invalid selection');
            end                        
            Y = chi2_MC_Top;
            
            % This would have otherwise held the selected fit from the grid
            % search
            X_Selected_List             = NaN;
            Y_Selected_List             = NaN;
            
        else
            % Display chi2 map (not errors)
            % Get X-values depending on what is selected
            switch( params_string{p} )
                case {'dwH', 'dwX', 'Pa', 'kex', 'R20', 'Rex', 'alpha', 'PhiexX'}
                    
                    p
                    
                    % Assign parameter values out of top fraction for current (cs,c,p)        
                    p_Grid_Top          = NaN*ones(Ntop,1);
                    p_Grid_Top(1:Ntop)  = resultsMatrix_Grid(ctot, p, igrid_Sorted(1:Ntop));

                    % Plot
                    X = session.convertParameterToDisplayUnits(p_Grid_Top,p);
                    
                    % From the grid search selected list
                    if( ~isempty(fitResult_Selected_List) )
                        X_Selected_List = session.convertParameterToDisplayUnits(resultsMatrix_Selected_List(ctot,p),p);
                    else
                        X_Selected_List = NaN;
                    end
                    
                case 'dH'
                    p_Grid_Top          = NaN*ones(Ntop,1);
                    p_Grid_Top(1:Ntop)  = dH_Grid(igrid_Sorted(1:Ntop));                    
                    X                   = p_Grid_Top ./ 1e3;
                    
                    % From the grid search selected list
                    if( ~isempty(fitResult_Selected_List) )
                        X_Selected_List     = fitResult_Selected_List.dH ./ 1e3;
                    else
                        X_Selected_List = NaN;
                    end
                    
                case 'Eab'
                    p_Grid_Top          = NaN*ones(Ntop,1);
                    p_Grid_Top(1:Ntop)  = Eab_Grid(igrid_Sorted(1:Ntop));                    
                    X                   = p_Grid_Top ./ 1e3;
                    
                    % From the grid search selected list
                    if( ~isempty(fitResult_Selected_List) )
                        X_Selected_List     = fitResult_Selected_List.Eab ./ 1e3;
                    else
                        X_Selected_List = NaN;
                    end
                    
                otherwise
                    % kex vs. pa
                    %error('Invalid selection');
            end            
            Y = chi2_Grid_Top;
            
            if( ~isnan(f_Selected_Top) )
                Y_Selected_List     = chi2_Grid_Sorted(f_Selected_Top(1));
            else
                Y_Selected_List     = NaN;
            end
            
            % Moved below for proper plotting [2011/09/02]
%             if( ~DISPLAY_HISTOGRAM )
%                 % Highlight lowest SSE fit in blue square
%                 plot(h, X(1), Y(1), 'sb', 'Linewidth', 2, 'MarkerSize', 8);
%             end
        end
        
        % Get the selected value (best fit)
        switch( params_string{p} )
            case {'dwH', 'dwX', 'Pa', 'kex', 'R20', 'Rex', 'alpha', 'PhiexX'}
                % Obtain the selected fit's X and Y values
                Xs      = session.convertParameterToDisplayUnits(fitResult_Selected.resultsMatrix(ctot,p),p);
                Xs_STD  = session.convertParameterToDisplayUnits(...
                    fitResult_Selected.resultsMatrixErrors(ctot,p),p);

            case 'dH'
                % Obtain the selected fit's X and Y values
                Xs      = fitResult_Selected.dH ./ 1e3;    
                Xs_STD  = fitResult_Selected.dH_E ./ 1e3;

            case 'Eab'
                % Obtain the selected fit's X and Y values
                Xs      = fitResult_Selected.Eab ./ 1e3;  
                Xs_STD  = fitResult_Selected.Eab_E ./ 1e3;
                
            otherwise
                % kex vs. pa
                %error('Invalid selection');
        end        
        Ys = fitResult_Selected.chi2;
        
        
        % (dwX, PaPb) surface plot (overhead view)
        % 2011/09/02 Start coding
        %  This code is a bit sloppy since it is separated from the above
        %  switch(case) statements, but it works for this special case
        if( strcmpi(params_string{p}, 'dwX v PaPb') )... || strcmpi(params_string{p}, 'PhiexX') )
            if( DISPLAY_HISTOGRAM )
                X       = NaN;
                Y       = NaN;
                
                if( strcmpi(params_string{p}, 'dwX v PaPb') )
                    param_name_plot     = '\Delta\omega_X vs P_AP_B';
                    value_string        = 'Histogram N/A';
                    param_units_display = '';
                    
                elseif( strcmpi(params_string{p}, 'PhiexX') )
                    param_name_plot     = '\Phi_{ex}';
                    value_string        = 'value_string';
                    %value_string        = sprintf('(%0.1f %s vs %0.2f %s)', Ys, Y_units, Xs, X_units);
                    param_units_display = 'Hz^2';
                end
                
            else
                
                % Get values from data structures                
                pp = find( strcmp( params_string, 'dwX' ) );
                if( DISPLAY_MC_ERRORS )
                    % Get the data from the Monte Carlo results                    
                    p_MC_Top            = NaN*ones(Ntop,1);
                    p_MC_Top(1:Ntop)    = resultsMatrix_MC(ctot, pp, iMC_Sorted(1:Ntop));                    
                    dwX                 = session.convertParameterToDisplayUnits(p_MC_Top,pp);
                else
                    p_Grid_Top          = NaN*ones(Ntop,1);
                    p_Grid_Top(1:Ntop)  = resultsMatrix_Grid(ctot, pp, igrid_Sorted(1:Ntop));
                    dwX                 = session.convertParameterToDisplayUnits(p_Grid_Top,pp);
                end
                dwX_S               = session.convertParameterToDisplayUnits(fitResult_Selected.resultsMatrix(ctot,pp),pp);
                dwX_units_display   = 'Hz^2';
                
                % Next parameter
                pp = find( strcmp( params_string, 'Pa' ) );
                if( DISPLAY_MC_ERRORS )
                    % Get the data from the Monte Carlo results                    
                    p_MC_Top            = NaN*ones(Ntop,1);
                    p_MC_Top(1:Ntop)    = resultsMatrix_MC(ctot, pp, iMC_Sorted(1:Ntop));                    
                    Pa                  = p_MC_Top;
                else
                    p_Grid_Top          = NaN*ones(Ntop,1);
                    p_Grid_Top(1:Ntop)  = resultsMatrix_Grid(ctot, pp, igrid_Sorted(1:Ntop));
                    Pa                  = p_Grid_Top;
                end
                Pa_S               = fitResult_Selected.resultsMatrix(ctot,pp);
                Pa_units_display   = '-';
                

                % From the grid search selected list
                if( ~isempty(fitResult_Selected_List) )
                    dwX_Selected_List = session.convertParameterToDisplayUnits(resultsMatrix_Selected_List(ctot,pp),pp);
                else
                    dwX_Selected_List = NaN;
                end

                % From the grid search selected list
                if( ~isempty(fitResult_Selected_List) )
                    Pa_Selected_List = resultsMatrix_Selected_List(ctot,pp);
                else
                    Pa_Selected_List = NaN;
                end

                Pb      = 1 - Pa;
                Pb_S    = 1 - Pa_S;
                
                if( strcmpi(params_string{p}, 'dwX v PaPb') )
                    % Set values for plotting
                    X               = Pa .* Pb;
                    Xs              = Pa_S .* Pb_S;
                    Xs_STD          = NaN;
                    X_Selected_List = Pa_Selected_List * (1-Pa_Selected_List);
                    X_units         = Pa_units_display;

                    Y               = dwX .^ 2;
                    Ys              = dwX_S .^ 2;
                    Y_Selected_List = dwX_Selected_List;
                    Y_units         = dwX_units_display;

                    % Set specs for plotting
                    param_name_plot     = '\Delta\omega_X^2 vs P_AP_B';
                    value_string        = sprintf('(%0.1f %s vs %0.2f %s)', Ys, Y_units, Xs, X_units);
                    param_units_display = '';

                    if( DISPLAY_MC_ERRORS )
                        COLOR_ARRAY = chi2_MC_Top;                    
                    else
                        COLOR_ARRAY = chi2_Grid_Top;
                    end

                    SIZE_ARRAY = 30;
                    scatter(X, Y, SIZE_ARRAY, COLOR_ARRAY, 'o', 'LineWidth', 2, 'Parent', h)
                    colormap('Gray');
                    
                elseif( strcmpi(params_string{p}, 'PhiexX') )
                    
                    % Calculate PhiexX using pre-defined parameters
                    X   = Pa .* Pb .* dwX.^2;
                    Xs  = Pa_S .* Pb_S .* dwX_S.^2;
                    
                    if( DISPLAY_MC_ERRORS )
                        Y = chi2_MC_Top;                    
                    else
                        Y = chi2_Grid_Top;
                    end
                    
                    % Set specs for plotting
                    param_name_plot     = '\Phi_{ex}';
                    param_units_display = 'Hz^2';
                    
                    % If there are errors to show, then show them
                    if( ~isnan(Xs_STD) )
                        value_string = sprintf('%0.1f[%0.1f]', Xs, Xs_STD);
                    else
                        value_string = sprintf('%0.1f', Xs);
                    end            
                    
                end
            end

            
        else
            % Plotting specs can be specified in a general way if the
            % surface plot is not used
            param_name_plot = session.param_name_plot{p};
            param_units_display = session.param_units_display{p};
            
            % If there are errors to show, then show them
            if( ~isnan(Xs_STD) )
                value_string = sprintf('%0.1f[%0.1f]', Xs, Xs_STD);
            else
                value_string = sprintf('%0.1f', Xs);
            end            
        end
        
        
        if( DISPLAY_HISTOGRAM )            
            hist(h, X);
            hh = findobj(h,'Type','patch');
            set(hh, 'FaceColor', 'k', 'EdgeColor', 'k');
            
            if( ic == 1 )
                ylabel('Num', 'FontWeight', 'Bold')
            end
                        
            % Only set parameter limits if there is a deviation in values
            if( std(X) > 1e-10 )                
                xlim(h, [1 1] .* [ min([X' Xs]), max([X' Xs]) ] );               
            end
            
            YLIM = get(h, 'YLim');
            ylim(h, [0, YLIM(2)]);
            
            X_NORM  = linspace( min([X' Xs]), max([X' Xs]), 50 );
            Y_NORM  = pdf('Normal', X_NORM, Xs, Xs_STD);
            
            % Scale Y axis values to match current limits
            Y_NORM  = 0.9 * Y_NORM .* (YLIM(2) / max(Y_NORM));

            plot(h, X_NORM, Y_NORM, '-', 'Color', [0.5 0.5 0.5], ...
                'LineWidth', session.LINEWIDTH);
            
        else                     
            % Highlight lowest SSE fit in blue square
            plot(h, X(1), Y(1), 'sb', 'Linewidth', 2, 'MarkerSize', 8);
            
            if( ~strcmpi(params_string{p}, 'dwX v PaPb') )
                % Plot all in black
                plot(h, X, Y, '.k', 'Linewidth', 1);
            end

            % Highlight currently selected fit with red circle with dot in
            % the middle
            plot(h, Xs, Ys, 'or', 'Linewidth', 2, 'MarkerSize', 8);
            plot(h, Xs, Ys, 'or', 'Linewidth', 2, 'MarkerSize', 2);
            
            % Plot the selected fit from the grid list
            plot(h, X_Selected_List, Y_Selected_List, 'dg', 'Linewidth', 2, 'MarkerSize', 8);
            
                        
            % Only set parameter limits if there is a deviation in values
            if( std(X) > 1e-10 )
                xlim(h, [1 1] .* [ min([X' Xs]), max([X' Xs]) ] );                
            end
            ylim(h, [1 1] .* [ min([Y' Ys]), max([Y' Ys]) ] );
        end
        
        title(h, sprintf('%s%s = %s%s', s_title, param_name_plot, ...
            value_string, param_units_display), ...
            'FontSize', session.FONTSIZE_SMALL, 'FontWeight', 'Bold');

        %set(h,'XGrid', 'on');
        %set(h,'YGrid', 'on');
        box(h,'on');        
    end
    
    
    %% (Pa,kex) sufrace plot (overhead view)
    % Only display if parameter is selected
    if( 0 && sum(6==param_index_selected) )
        
        % Assign parameter values out of top fraction for current (cs,c,p)
        p = 3;
        p_top = fit_results.grid_search{cs}.p_matrix_fits(sort_igrid(1:Ntop),c,p); % (cs,igrid,c,p)
        pa_top      = p_top';
        select_pa   = select_p_matrix_fit(c,p);
        
        % Assign parameter values out of top fraction for current (cs,c,p)
        p = 4;
        p_top = fit_results.grid_search{cs}.p_matrix_fits(sort_igrid(1:Ntop),c,p); % (cs,igrid,c,p)
        kex_top     = p_top';
        select_kex = select_p_matrix_fit(c,p);
    
        %subplot_number  = (data.Np)*Ncols + ic;
        subplot_number  = subplots_used*Ncols + ic;
        subplots_used   = subplots_used+1;

        h = subplot(Nrows, Ncols, subplot_number, 'Parent', figure_handle );        
        set(h, 'FontSize', session.FONTSIZE_SMALL);
        hold(h,'on');
        
        % Only output the title for the parameter at the top
        s_title = '';
        if( subplots_used == 1 )
            s_title = sprintf('%s\n', char(data.name_all(cs,c)));
        end

        %contour3(hf, pa_top', kex_top', sse_top')
        %scatter3(pa_top, kex_top, sse_top, 20, sse_top, 'o', 'LineWidth', 2, 'Parent', h)
        scatter(pa_top*session.convert_units_to_display(3), kex_top, 20, sse_top, 'o', 'LineWidth', 2, 'Parent', h)
        colormap('Gray');

        % Highlight lowest SSE fit in blue square
        plot(h, pa_top(1)*session.convert_units_to_display(3), kex_top(1), 'sb', 'Linewidth', 2, 'MarkerSize', 8);

        % Highlight currently selected fit with red circle with dot in the middle
        %scatter3(select_pa, select_kex, select_sse, 10, [1 0 0], 'o', 'Linewidth', 2, 'Parent', h);
        %scatter3(select_pa, select_kex, select_sse, 10, [1 0 0], 'o', 'Linewidth', 2, 'Parent', h);  
        %plot3(h, select_pa, select_kex, select_sse, 'or', 'Linewidth', 2, 'MarkerSize', 8);    
        plot(h, select_pa*session.convert_units_to_display(3), select_kex, 'or', 'Linewidth', 2, 'MarkerSize', 8);
        plot(h, select_pa*session.convert_units_to_display(3), select_kex, 'or', 'Linewidth', 2, 'MarkerSize', 2);

        title(h, sprintf('%sk_{ex} vs P_A scatter', s_title), ...
                    'FontSize', session.FONTSIZE_SMALL, 'FontWeight', 'Bold');
                
        set(h,'XGrid', 'on');
        set(h,'YGrid', 'on');

        % If it IS the left-most subplot
        %if( mod(subplot_number,length(Temps)) == 1 )
        if( ic == 1 )
            ylabel(sprintf('k_{ex} (%s)',session.param_units_display(3)), 'FontWeight', 'Bold')
        end
        xlabel(sprintf('P_A (%s)', session.param_units_display(3)), 'FontWeight', 'Bold');

        xlim( [1 1] .* [ min(pa_top), max(pa_top) ] .* session.convert_units_to_display(3));
        ylim( [1 1] .* [ min(kex_top), max(kex_top) ] );
    end
end




return

% Add a line to separate each column of subplots [2011/01/07]
% Each column corresponds to a single temperature
ax=axes('Units','Normal','Position',[0, 0, 1, 1], 'Visible','off', 'Parent', figure_handle);
set( ax, 'XLim', [0 1]);
set( ax, 'YLim', [0 1]);

%line([0 1], [0 1], 'Parent', ax);

if( Ncols >= 4 )
    % Need to calibrate where the lines go
    % X_POSITION = [0.25 0.5 0.75];
elseif( Ncols == 3 )
    X_POSITION = [0.37 0.65];
elseif( Ncols == 2 )
    X_POSITION = [0.5];
elseif( Ncols == 1 )
    X_POSITION = NaN;
end

% Only output lines if I know where to put them
if( Ncols <= 4 )
    for col = 1:Ncols-1        
        line(X_POSITION(col).*[1 1], [.1 0.96], 'LineWidth', 3, 'Color', 'k', 'Parent', ax)        
        %line( 'XData', [-.25 -.25], 'YData', [.05 0.97], 'LineWidth', 3, 'Color', 'k', 'Parent', ax)
    end
end


%% Set title
title_string = sprintf('%s: Best %d/%d fits', data.name_string{cs}, Ntop, ...
                                            fit_results.Ngrid_local(cs) );

                                        % Set title of panel if drawing current panel
if( figure_handle == handles.panel_display_map_small )
    set(handles.panel_display_map_small, 'Title', title_string, 'FontSize', session.FONTSIZE_MEDIUM );
    
% Set super-title of subplots if using new figure window
else
    % Position [left, bottom, width, height]
    % Empirically determined for this plot
    if( data.Nc(cs) == 1 )
        ax_position = [0.08, 0.08, 0.84, 0.84];
    elseif( data.Nc(cs) == 2 )
        ax_position = [0.08, 0.08, 0.84, 0.84];
    elseif( data.Nc(cs) == 3 )
        ax_position = [0.08, 0.08, 0.84, 0.85];
    else
        ax_position = [0.08, 0.08, 0.84, 0.85];
    end
    ax=axes('Units','Normal','Position',ax_position,'Visible','off', 'Parent', figure_handle);
    set(get(ax,'Title'),'Visible','on')        
    title(ax, title_string, 'FontSize', session.FONTSIZE_MEDIUM, 'FontWeight', 'Bold');
end

% --- Executes when display_chi2_map_gui is resized.
function display_chi2_map_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to display_chi2_map_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%{

%% Resize each element of the GUI

% Get size of the GUI window which is just resized
fP = get(handles.display_chi2_map_gui, 'Position');

left    = 0.5;
bottom  = 0.1;
width   = fP(3)-left;
height  = fP(4)-3.5;
set(handles.panel_display_map_small, 'Position', [left, bottom, width, height]);

%% Listbox
width   = 18;
height  = 4;
left    = 1;
bottom  = fP(4)-4.2;
set(handles.listbox_params, 'Position', [left, bottom, width, height]);

%% Other elements
% Adjust position of popup_fit_results (scale position, fix size)
width   = 45;
height  = 1.5;
left    = 1+18+1;
bottom  = fP(4)-2.5;
set(handles.popup_fit_results, 'Position', [left, bottom, width, height]);

% Adjust position of top_fraction table
width   = 30;
height  = 2.5;
left    = 1+18+1+45+1;
bottom  = fP(4)-3;
set(handles.table_topf, 'Position', [left, bottom, width, height]);
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
handles_main    = getappdata(handles.display_chi2_map_gui, 'handles_main');
session        = getappdata(handles_main.main_gui, 'session');

% Read user settings from table and update default settings structure
table_data = get(handles.table_topf, 'Data');
session.setChi2TopFraction(table_data(1)/100);

% Update the slider
set(handles.slider_topf, 'Value', table_data(1)/100);

refresh_display(handles);


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
handles_main    = getappdata(handles.display_chi2_map_gui, 'handles_main');
data            = getappdata(handles_main.main_gui, 'data');
settings        = getappdata(handles_main.main_gui, 'settings');
cs              = get(       handles_main.listbox_cs, 'Value');

% Get the fit number from the list
f = get(handles.popup_fit_results, 'Value');

% Ask user for file to save as
default_filename     = sprintf('%s-%s-F%02d.%s', 'scatter-chi2', data.filename_suffix_cs{cs}, f, '');
[filename, filepath] = uiputfile('*.fig; *.png; *.ps', 'Save graphic file', ...
    sprintf('./%s/%s', session.outputDir, default_filename) );

% If the user selected a file
if( ~isequal(filename,0) )
    
    % Make new figure window
    h = figure;

    % Plot the data to a new figure    
    plot_chi2_map( h, handles )
    
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
plot_chi2_map( h, handles )


% --- Executes on selection change in listbox_params.
function listbox_params_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_params contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_params

%% Update GUI with new parameters selected
handles_main    = getappdata(handles.display_chi2_map_gui, 'handles_main');
session        = getappdata(handles_main.main_gui, 'session');

% Save the changes
session.setChi2MapParamters( get(handles.listbox_params, 'Value') );
refresh_display(handles);


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


% --- Executes when selected object is changed in panel_view_mode.
function panel_view_mode_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_view_mode 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
refresh_display(handles);


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
refresh_display(handles);


% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
refresh_display(handles);


% --- Executes on button press in pushbutton_add_grid_fit_to_list.
function pushbutton_add_grid_fit_to_list_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_grid_fit_to_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Add fit to list

handles_main    = getappdata(handles.display_chi2_map_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

g               = get(handles_main.listbox_g, 'Value');
group           = session.groups{g};

% Get info on the grid

% Display the inital or final conditions?
if( get(handles.radiobutton_initial_conditions, 'Value') == 1 )
    INITIAL_OR_FINAL = 'Initial';
else
    INITIAL_OR_FINAL = 'Final';
end
[chi2_Grid, VOID, VOID, VOID] = group.getGridResults(INITIAL_OR_FINAL);

% Sort the chi2 array
[chi2_Grid_Sorted, igrid_Sorted]  = sort( chi2_Grid );

f_Selected_Top          = getappdata(handles.display_chi2_map_gui, 'f_Selected_Top');
f_Selected_Grid         = igrid_Sorted( f_Selected_Top );
fitResult_Selected_List = group.fitResults_Grid{f_Selected_Grid};

% Add it to the group
group.addFitResult( fitResult_Selected_List );

% Select the fit from the list
set(handles.popup_fit_results, 'Value', group.Nf);
refresh_display(handles);


% --- Executes on button press in checkbox_select_fit_from_grid.
function checkbox_select_fit_from_grid_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_select_fit_from_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_select_fit_from_grid
refresh_display(handles);

% --- Executes when selected cell(s) is changed in table_grid_fits.
function table_grid_fits_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_grid_fits (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

rows_cols_selected = eventdata.Indices;

if( ~isempty(rows_cols_selected) )
    % Only get the first one selected
    f_Selected_Top = rows_cols_selected(1,2);
else
    f_Selected_Top = [];
end

setappdata(handles.display_chi2_map_gui, 'f_Selected_Top',f_Selected_Top);
refresh_display(handles);


% --- Executes on slider movement.
function slider_topf_Callback(hObject, eventdata, handles)
% hObject    handle to slider_topf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%% Get TOP_FRACTION from the slider

% Read top fraction and Nbins from "table_sse_specs"
handles_main    = getappdata(handles.display_chi2_map_gui, 'handles_main');
session        = getappdata(handles_main.main_gui, 'session');

% Read user settings from table and update default settings structure
slider_value = get(handles.slider_topf, 'Value');
session.setChi2TopFraction(slider_value);

% Update the table
set(handles.table_topf, 'Data', slider_value*100);

refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function slider_topf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_topf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

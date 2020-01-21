% Export simulated data from RD simulator to ASCII file (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/??/?? Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/01/12 Add noise to export
% 2011/04/15 Use classes for data storage
% 2011/05/06 Update CSV export
%
% TO DO
%  Check if parameters are OK pre-export?

function varargout = GUARDD_rd_simulator_export(varargin)
% GUARDD_RD_SIMULATOR_EXPORT M-file for GUARDD_rd_simulator_export.fig
%      GUARDD_RD_SIMULATOR_EXPORT, by itself, creates a new GUARDD_RD_SIMULATOR_EXPORT or raises the existing
%      singleton*.
%
%      H = GUARDD_RD_SIMULATOR_EXPORT returns the handle to a new GUARDD_RD_SIMULATOR_EXPORT or the handle to
%      the existing singleton*.
%
%      GUARDD_RD_SIMULATOR_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_RD_SIMULATOR_EXPORT.M with the given input arguments.
%
%      GUARDD_RD_SIMULATOR_EXPORT('Property','Value',...) creates a new GUARDD_RD_SIMULATOR_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_rd_simulator_export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_rd_simulator_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_rd_simulator_export

% Last Modified by GUIDE v2.5 16-Jan-2011 14:13:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_rd_simulator_export_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_rd_simulator_export_OutputFcn, ...
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


% --- Executes just before GUARDD_rd_simulator_export is made visible.
function GUARDD_rd_simulator_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_rd_simulator_export (see VARARGIN)

% Choose default command line output for GUARDD_rd_simulator_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_rd_simulator_export wait for user response (see UIRESUME)
% uiwait(handles.rd_simulator_export_gui);



%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.rd_simulator_export_gui, 'handles_main', handles_main);


%% Initialize the table for export params
row_string = {'Min( vCPMG )', 'Max( vCPMG )', 'Sim Points','Sim Noise (%)'};
set(handles.table_export, 'RowName', row_string);
%set(handles.table_export, 'ColumnWidth', {50,50,50,50});
set(handles.table_export, 'ColumnName', [] );
set(handles.table_export, 'Data', {100; 1000; 10; 5});


%refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_rd_simulator_export_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in button_export_csv.
function button_export_csv_Callback(hObject, eventdata, handles)
% hObject    handle to button_export_csv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Export simulations to text file
% Access main window handle which is stored in this window's main handle
handles_main        = getappdata(handles.rd_simulator_export_gui, 'handles_main');
simulationSession   = getappdata(handles_main.main_gui, 'simulationSession');
session             = getappdata(handles_main.main_gui, 'session');
VERSION             = getappdata(handles_main.main_gui, 'VERSION');

% Make the output directory if it does not yet exist
if( ~exist(session.outputDir, 'dir') )
    mkdir(session.outputDir);
end

% Ask user for file to save as
default_filename    = sprintf('output-GUARDD-Sim--%s.txt', datestr(now, 'yyyy.mm.dd-HH-MM'));
[filename, filepath] = uiputfile('*.txt', 'Save results to file', ...
    sprintf('./%s/%s', session.outputDir, default_filename) );

% Save it
if( isequal(filename,0) )
    return
end
fprintf('Writing %s...', filename);

FILE = fopen(sprintf('%s/%s', filepath, filename), 'w');

fprintf(FILE, 'GUARDD %s', VERSION);
fprintf(FILE, '\n(C) Ian Kleckner 2010-2011');
fprintf(FILE, '\nExported data from RD Simulator');
fprintf(FILE, '\nCreated on %s', datestr(now, 'yyyy.mm.dd-HH-MM'));

% Simulation specifications from table
table_data  = get(handles.table_export, 'Data');
vcpmgMin    = table_data{1};
vcpmgMax    = table_data{2};
Nobs        = table_data{3};            
fNoise      = table_data{4} / 100;

% Iterate through each simulation
for cs = 1:simulationSession.Ncs
    curveset = simulationSession.curvesets{cs};
    
    fprintf(FILE, '\n\nCurveset %d/%d, General Specifications', cs, simulationSession.Ncs);
    fprintf(FILE, '\n%s\t%s',     'NameCurveset', curveset.name);
    fprintf(FILE, '\n%s\t%d',     'NumCurves', curveset.Nc);
    fprintf(FILE, '\n%s\t%f\t%s', 'dwH',    curveset.dwHppm,    curveset.dwHppm_Units);
    fprintf(FILE, '\n%s\t%f\t%s', 'dwX',    curveset.dwXppm,    curveset.dwXppm_Units);
    fprintf(FILE, '\n%s\t%s\t%s', 'AX',     curveset.AX_String, curveset.AX_String_Units);
    fprintf(FILE, '\n%s\t%f\t%s', 'Eab',    curveset.Eab,       curveset.Eab_Units);
    fprintf(FILE, '\n%s\t%f\t%s', 'Pab',    curveset.Pab,       curveset.Pab_Units);
    fprintf(FILE, '\n%s\t%f\t%s', 'Eba',    curveset.Eba,       curveset.Eba_Units);
    fprintf(FILE, '\n%s\t%f\t%s', 'Pba',    curveset.Pba,       curveset.Pba_Units);
    fprintf(FILE, '\n%s\t%f\t%s', 'dH',     curveset.dH,        curveset.dH_Units);
    fprintf(FILE, '\n%s\t%f\t%s', 'dS',     curveset.dS,        curveset.dS_Units);
    fprintf(FILE, '\n%s\t%f\t%s', 'Temp0',  curveset.T0-273,    'C');
    fprintf(FILE, '\n%s\t%f\t%s', 'PA0',    curveset.PA0*100,   '%');
    fprintf(FILE, '\n%s\t%f\t%s', 'kex0',   curveset.kex0,      curveset.kex0_Units);
    
    if( curveset.Nc > 0 )

        % Iterate through each simulation in the set
        for c = 1:curveset.Nc
            curve = curveset.curves{c};
            
            fprintf(FILE, '\n\nCurveset %d/%d, Curve %d/%d', cs, simulationSession.Ncs, c, curveset.Nc);
            fprintf(FILE, '\n%s\t%s',     'NameCurve', curve.name);
            fprintf(FILE, '\n%s\t%s',     'NameParentCurveset', curveset.name);
            fprintf(FILE, '\n%s\t%f\t%s', 'Temp',   curve.getTempC,                     'C');
            fprintf(FILE, '\n%s\t%f\t%s', 'B0',     curve.B0,                           'MHz');
            fprintf(FILE, '\n%s\t%f\t%s', 'TCPMG',  curve.TCPMG,                        'sec');
            fprintf(FILE, '\n%s\t%d\t%s', 'SQX',    curve.SQX,                          '-');
            fprintf(FILE, '\n%s\t%f\t%s', 'dwH',    curve.getParamValueForDisplay(1),   'Hz');
            fprintf(FILE, '\n%s\t%f\t%s', 'dwX',    curve.getParamValueForDisplay(2),   'Hz');
            fprintf(FILE, '\n%s\t%f\t%s', 'PA',     curve.getParamValueForDisplay(3),   '%');
            fprintf(FILE, '\n%s\t%f\t%s', 'kex',    curve.getParamValueForDisplay(4),   '/sec');
            fprintf(FILE, '\n%s\t%f\t%s', 'R20',    curve.getParamValueForDisplay(5),   'Hz');
            
            fprintf(FILE, '\n\n%s\t%d',     'NObs',   Nobs);
            fprintf(FILE, '\n%s\t%f',     'Noise',  fNoise);
                        
            [vcpmgSim, R2effSim, eR2effSim, ISim, eISim] = ...
                curve.simulateData(vcpmgMin, vcpmgMax, Nobs, fNoise);
            
            
            fprintf(FILE, '\n%s\t%s\t%s\t%s\t%s\t%s', 'Obs', 'vcpmg(Hz)', ...
                'R2Eff(Hz)', 'eR2Eff(Hz)', 'Intensity', 'eIntensity');
            for o = 1:Nobs
                fprintf(FILE, '\n%d\t%f\t%f\t%f\t%f\t%f', o, vcpmgSim(o), R2effSim(o), eR2effSim(o), ...
                    ISim(o), eISim(o));
            end
        end
    end
end

close(handles.rd_simulator_export_gui)
fclose(FILE);
fprintf('Done!\n');    


% --- Executes on button press in button_export_GUARDD.
function button_export_GUARDD_Callback(hObject, eventdata, handles)
% hObject    handle to button_export_GUARDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Send simulations to GUARDD data structures
% Create NEW groups for each simulation set

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.rd_simulator_export_gui, 'handles_main');
simulationSession      = getappdata(handles_main.main_gui, 'simulationSession');
session         = getappdata(handles_main.main_gui, 'session');

% Simulation specifications from table
table_data  = get(handles.table_export, 'Data');
vcpmgMin    = table_data{1};
vcpmgMax    = table_data{2};
Nobs        = table_data{3};            
fNoise      = table_data{4} / 100;

% Iterate through each simulation
for cs = 1:simulationSession.Ncs
    curveset = simulationSession.curvesets{cs};
    
    if( curveset.Nc > 0 )
        % Simulation Set -> Group (which contains one curveset)
        new_group = Group(session);
        new_group.setName(curveset.name);   
        new_group.setIndex(cs);
        new_group.setNote(sprintf('Simulation Specs\nT0=%0.1fC, PA0=%0.1f%%, kex0=%0.1f/s\nEab=%0.1fkcal/mol, Pab=%0.1e/s\nEba=%0.1fkcal/mol, Pba=%0.1e/s\ndH=%0.1fkcal/mol, dS=%0.1fcal/mol/K\nAX=%s\ndwH=%0.1fppm, dwX=%0.1fppm', ...
                            curveset.T0-273, ...
                            curveset.PA0*100, ...
                            curveset.kex0, ...
                            curveset.Eab/1000, ...
                            curveset.Pab, ...
                            curveset.Eba/1000, ...
                            curveset.Pba, ...
                            curveset.dH/1000, ...
                            curveset.dS, ...
                            curveset.AX_String, ...
                            curveset.dwHppm, ...
                            curveset.dwXppm ));
                        
        session.addGroup(new_group);

        % Set -> Curveset
        new_curveset = Curveset;
        index   = cs;
        atom    = '';
        residue = 'Sim';
        new_curveset.setAssignment(index, atom, residue);
        new_curveset.setName(curveset.name);
        new_group.addCurveset(new_curveset);

        % Iterate through each simulation in the set
        for c = 1:curveset.Nc
            curve = curveset.curves{c};
            sim_params = curve.params;

            % Simulation set -> curve in the curveset
            new_curve = Curve;

            index   = c;
            atom    = '';
            residue = 'Sim';
            new_curve.setAssignment(index, atom, residue);
            
            % Set the dataset as the simulation that owns it
            new_curve.setDataset( curveset )

            B0      = curve.B0;
            Temp    = curve.Temp;
            TCPMG   = curve.TCPMG;
            SQX     = curve.SQX;
            
            AX_String = curveset.AX_String;
            new_curve.setConditions( B0, Temp, TCPMG, SQX, AX_String );
            
            % Simulate the data
            [vcpmgSim, R2effSim, eR2effSim] = ...
                curve.simulateData(vcpmgMin, vcpmgMax, Nobs, fNoise);
            
            new_curve.setData( vcpmgSim, R2effSim, eR2effSim );
            new_curve.setNotes( sprintf('T=%dC dwH=%0.0fr/s dwX=%0.0fr/s PA=%0.1f%%%% kex=%0.0f/s R20=%0.0fHz TCPMG=%0.1fms Noise=%0.1f%%%%', ...                                
                            curve.getTempC(), ...
                            curve.params(1), ...
                            curve.params(2), ...
                            curve.getParamValueForDisplay(3), ...
                            curve.getParamValueForDisplay(4), ...
                            curve.getParamValueForDisplay(5), ...
                            curve.TCPMG*1000, ...
                            fNoise*100 ...
                            ) );

            new_curveset.addCurve(new_curve);
        end
        
        % Update the fit parameters for this group
        CONSTRAIN_RATE_ANALYSIS = false;
        new_group.updateFitParams(CONSTRAIN_RATE_ANALYSIS); 
    end
end

%session.updateGroups();

close(handles.rd_simulator_export_gui)

% --- Executes when entered data in editable cell(s) in table_export.
function table_export_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_export (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


%% Check if values entered are reasonable

% TODO % Check if values entered in table are reasonable


% --- Executes when rd_simulator_export_gui is resized.
function rd_simulator_export_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to rd_simulator_export_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Alter settings in the program via session.param_info (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/06/07 Create code
%
% TO DO


function varargout = GUARDD_settings(varargin)
% GUARDD_SETTINGS M-file for GUARDD_settings.fig
%      GUARDD_SETTINGS, by itself, creates a new GUARDD_SETTINGS or raises the existing
%      singleton*.
%
%      H = GUARDD_SETTINGS returns the handle to a new GUARDD_SETTINGS or the handle to
%      the existing singleton*.
%
%      GUARDD_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_SETTINGS.M with the given input arguments.
%
%      GUARDD_SETTINGS('Property','Value',...) creates a new GUARDD_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_settings

% Last Modified by GUIDE v2.5 07-Jun-2011 15:11:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_settings_OutputFcn, ...
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


% --- Executes just before GUARDD_settings is made visible.
function GUARDD_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_settings (see VARARGIN)

% Choose default command line output for GUARDD_settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_settings wait for user response (see UIRESUME)
% uiwait(handles.settings_gui);



%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.settings_gui, 'handles_main', handles_main);

%% Initialize the table
session = getappdata(handles_main.main_gui, 'session');
set(handles.table_settings, 'Data', generateParamTable(session) );


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in button_ok.
function button_ok_Callback(hObject, eventdata, handles)
% hObject    handle to button_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Automatically generate groups via NMR probe, with constraints from table
handles_main        = getappdata(handles.settings_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% Get data from table
table_data = get(handles.table_settings, 'Data');

% Commit changes from table to session
session.saveDataFromTable( table_data );

% Done!
close(handles.settings_gui);


% --- Executes when settings_gui is resized.
function settings_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to settings_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in table_group_constraints.
function table_group_constraints_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_group_constraints (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_reset.
function button_reset_Callback(hObject, eventdata, handles)
% hObject    handle to button_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Reset the table
handles_main        = getappdata(handles.settings_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');
set(handles.table_group, 'Data', generateParamTable(session) )

% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Cancel button closes GUI
close(handles.settings_gui);

% --- Executes when entered data in editable cell(s) in table_settings.
function table_settings_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_settings (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% Table is changed
set(handles.button_ok, 'Enable', 'on');
set(handles.button_reset, 'Enable', 'on');


% --- Executes when selected cell(s) is changed in table_settings.
function table_settings_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_settings (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
%% Proceed only if there is something selected
if( isempty( eventdata.Indices ) )
    return
end

% Check what the name of the row that is selcted
table_data = get(handles.table_settings, 'Data');
row_selected = eventdata.Indices(1);
parameter_name = table_data{ row_selected, 1 };

% TODO % Intelligent want to determine if source is from a file
% If it is the output directory
if( strcmp(parameter_name, 'outputDir') )
    % Ask user to select file
    pathname = uigetdir('.', 'Select output directory');
    
    % Make directory relative to current
    relative_path = strrep(pathname, pwd, '.');

    % Commit the selection to the table
    if( ~isequal(pathname,0) )
        table_data{row_selected, 2} = relative_path;
        set(handles.table_settings, 'Data', table_data);
        
        set(handles.button_ok, 'Enable', 'on');
        set(handles.button_reset, 'Enable', 'on');
    end
end

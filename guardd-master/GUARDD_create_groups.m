% Export simulated data from RD simulator to ASCII file (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/06/06 Create code
%
% TO DO


function varargout = GUARDD_create_groups(varargin)
% GUARDD_CREATE_GROUPS M-file for GUARDD_create_groups.fig
%      GUARDD_CREATE_GROUPS, by itself, creates a new GUARDD_CREATE_GROUPS or raises the existing
%      singleton*.
%
%      H = GUARDD_CREATE_GROUPS returns the handle to a new GUARDD_CREATE_GROUPS or the handle to
%      the existing singleton*.
%
%      GUARDD_CREATE_GROUPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_CREATE_GROUPS.M with the given input arguments.
%
%      GUARDD_CREATE_GROUPS('Property','Value',...) creates a new GUARDD_CREATE_GROUPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_create_groups_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_create_groups_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_create_groups

% Last Modified by GUIDE v2.5 06-Jun-2011 14:57:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_create_groups_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_create_groups_OutputFcn, ...
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


% --- Executes just before GUARDD_create_groups is made visible.
function GUARDD_create_groups_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_create_groups (see VARARGIN)

% Choose default command line output for GUARDD_create_groups
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_create_groups wait for user response (see UIRESUME)
% uiwait(handles.create_groups_gui);



%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.create_groups_gui, 'handles_main', handles_main);


%% Initialize the table for export params
row_string      = {'Index', 'Atom', 'Residue', 'TempC', 'B0', 'SQX'};
units_string = {'Residue number (if assigned)', 'Atom on residue (if assigned)', ...
                'Residue type (if assigned)', 'Celcius', 'MHz', '1=Single Quantum, 0=Multiple Quantum'};

% Add HTML tag for italics
for r = 1:length(units_string)
    units_string{r} = sprintf('<html><i>%s</i>', units_string{r});
end
            
set(handles.table_group_constraints, 'RowName', row_string);
set(handles.table_group_constraints, 'ColumnWidth', {80,80,280});

column_string = {'Constrain?', 'Value', 'Description'};
set(handles.table_group_constraints, 'ColumnName', column_string );

table_data = { false, 1, units_string{1}; ...
               false, 'NH', units_string{2}; ...
               false, 'Ile', units_string{3}; ...
               false, 25, units_string{4}; ...
               false, 800, units_string{5}; ...
               false, 1, units_string{6} ...
                };

set(handles.table_group_constraints, 'Data', table_data);



% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_create_groups_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in button_generate_groups.
function button_generate_groups_Callback(hObject, eventdata, handles)
% hObject    handle to button_generate_groups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Automatically generate groups via NMR probe, with constraints from table
handles_main        = getappdata(handles.create_groups_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% Check the table data to make the proper specifications
table_data = get(handles.table_group_constraints, 'Data');
Nrows = size(table_data,1);

varargin = {};

param_names = get(handles.table_group_constraints, 'RowName');

for r = 1:Nrows
    param           = param_names{r};
    constrain_param = table_data{r,1};
    constrain_value = table_data{r,2};
    
    if( constrain_param )
        varargin{end+1} = upper(param);
        varargin{end+1} = constrain_value;
    end
end

session.generateGroups(varargin{:});
close(handles.create_groups_gui);


% --- Executes when create_groups_gui is resized.
function create_groups_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to create_groups_gui (see GCBO)
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

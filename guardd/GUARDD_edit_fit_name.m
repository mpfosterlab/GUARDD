% Rename the currently selected fitResult
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/06/10 Create code
%
% TO DO


function varargout = GUARDD_edit_fit_name(varargin)
% GUARDD_EDIT_FIT_NAME M-file for GUARDD_edit_fit_name.fig
%      GUARDD_EDIT_FIT_NAME, by itself, creates a new GUARDD_EDIT_FIT_NAME or raises the existing
%      singleton*.
%
%      H = GUARDD_EDIT_FIT_NAME returns the handle to a new GUARDD_EDIT_FIT_NAME or the handle to
%      the existing singleton*.
%
%      GUARDD_EDIT_FIT_NAME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_EDIT_FIT_NAME.M with the given input arguments.
%
%      GUARDD_EDIT_FIT_NAME('Property','Value',...) creates a new GUARDD_EDIT_FIT_NAME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_edit_fit_name_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_edit_fit_name_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_edit_fit_name

% Last Modified by GUIDE v2.5 10-Jun-2011 11:57:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_edit_fit_name_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_edit_fit_name_OutputFcn, ...
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


% --- Executes just before GUARDD_edit_fit_name is made visible.
function GUARDD_edit_fit_name_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_edit_fit_name (see VARARGIN)

% Choose default command line output for GUARDD_edit_fit_name
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_edit_fit_name wait for user response (see UIRESUME)
% uiwait(handles.edit_fit_name_gui);



%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

index_fitResult    = find(strcmp(varargin, 'fitResult'));
% The actual handle is the index after (+1) the name
fitResult            = varargin{index_fitResult+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.edit_fit_name_gui, 'handles_main', handles_main);
setappdata(handles.edit_fit_name_gui, 'fitResult', fitResult);


%% Initialize edit box
set(handles.edit_fitResult_name, 'String', fitResult.name);



% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_edit_fit_name_OutputFcn(hObject, eventdata, handles) 
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
fitResult        = getappdata(handles.edit_fit_name_gui, 'fitResult');

new_name  = get(handles.edit_fitResult_name, 'String');

fitResult.setName(new_name);
close(handles.edit_fit_name_gui);


% --- Executes when edit_fit_name_gui is resized.
function edit_fit_name_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to edit_fit_name_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function edit_fitResult_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fitResult_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.edit_fit_name_gui);

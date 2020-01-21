% Dialog to load NMRPipe file with CPMG RD data
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/06/07? Start coding

function varargout = GUARDD_load_curves(varargin)
% LOAD_CURVES_GUI M-file for load_curves_gui.fig
%      LOAD_CURVES_GUI, by itself, creates a new LOAD_CURVES_GUI or raises the existing
%      singleton*.
%
%      H = LOAD_CURVES_GUI returns the handle to a new LOAD_CURVES_GUI or the handle to
%      the existing singleton*.
%
%      LOAD_CURVES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_CURVES_GUI.M with the given input arguments.
%
%      LOAD_CURVES_GUI('Property','Value',...) creates a new LOAD_CURVES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_load_curves_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_load_curves_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help load_curves_gui

% Last Modified by GUIDE v2.5 01-Feb-2011 14:33:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_load_curves_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_load_curves_OutputFcn, ...
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


% --- Executes just before load_curves_gui is made visible.
function GUARDD_load_curves_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to load_curves_gui (see VARARGIN)

% Choose default command line output for load_curves_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes load_curves_gui wait for user response (see UIRESUME)
% uiwait(handles.load_curves_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));

% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.load_curves_gui, 'handles_main', handles_main);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_load_curves_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_go_nmrpipe.
function button_go_nmrpipe_Callback(hObject, eventdata, handles)
% hObject    handle to button_go_nmrpipe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Load the nmrPipe files into the selected dataset

nlin_filename   = get(handles.edit_nlin_filename, 'String');
vcpmg_filename  = get(handles.edit_vcpmg_filename, 'String');

% Check the filename
if( ~exist(nlin_filename,'file') )
    msgbox('Invalid nlin.tab file. Please reload', 'Invalid nlin.tab', 'Error');
    return

elseif( ~exist(vcpmg_filename,'file') )
    msgbox('Invalid vcpmg text file. Please reload', 'Invalid vcpmg.txt', 'Error');
    return
    
% Both required files exist
else
    % Save the loading specs to the main GUI
    LoadCurves.MODE             = 'NMRPIPE';
    LoadCurves.nlin_filename    = nlin_filename;
    LoadCurves.vcpmg_filename   = vcpmg_filename;
    
    % Call the load curves GUI using variables from main handle
    handles_main    = getappdata(handles.load_curves_gui, 'handles_main');
    setappdata(handles_main.main_gui, 'LoadCurves', LoadCurves);
    
    % Signal completion of the GUI by closing the window
    %  Note, a waitfor() function remains from the GUI caller, GUARDD_data_manager
    close(handles.load_curves_gui);    
end

% --- Executes during object creation, after setting all properties.
function edit_nlin_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nlin_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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



function edit_vcpmg_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vcpmg_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vcpmg_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_vcpmg_filename as a double


% --- Executes during object creation, after setting all properties.
function edit_vcpmg_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vcpmg_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_browse_nlin_filename.
function button_browse_nlin_filename_Callback(hObject, eventdata, handles)
% hObject    handle to button_browse_nlin_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Get the nlin.tab file
title = 'Select nlin.tab file';
file_extension = '*.tab';
[filename, pathname] = uigetfile(file_extension, title);

% Set the edit box string as filename, if the user selected one
if( ~isequal(filename,0) )
    set(handles.edit_nlin_filename, 'String', sprintf('%s/%s', pathname, filename));
end

% --- Executes on button press in button_browse_vcpmg_filename.
function button_browse_vcpmg_filename_Callback(hObject, eventdata, handles)
% hObject    handle to button_browse_vcpmg_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Get the nlin.tab file
title = 'Select vcpmg file';
file_extension = '*.txt';
[filename, pathname] = uigetfile(file_extension, title);

% Set the edit box string as filename, if the user selected one
if( ~isequal(filename,0) )
    set(handles.edit_vcpmg_filename, 'String', sprintf('%s/%s', pathname, filename));
end

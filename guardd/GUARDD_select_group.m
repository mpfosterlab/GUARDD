% List all the groups, and select one (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/06/01 Start coding
%
% TO DO


function varargout = GUARDD_select_group(varargin)
% GUARDD_SELECT_GROUP M-file for GUARDD_select_group.fig
%      GUARDD_SELECT_GROUP, by itself, creates a new GUARDD_SELECT_GROUP or raises the existing
%      singleton*.
%
%      H = GUARDD_SELECT_GROUP returns the handle to a new GUARDD_SELECT_GROUP or the handle to
%      the existing singleton*.
%
%      GUARDD_SELECT_GROUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_SELECT_GROUP.M with the given input arguments.
%
%      GUARDD_SELECT_GROUP('Property','Value',...) creates a new GUARDD_SELECT_GROUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_select_group_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_select_group_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_select_group

% Last Modified by GUIDE v2.5 01-Jun-2011 10:22:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_select_group_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_select_group_OutputFcn, ...
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


% --- Executes just before GUARDD_select_group is made visible.
function GUARDD_select_group_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_select_group (see VARARGIN)

% Choose default command line output for GUARDD_select_group
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_select_group wait for user response (see UIRESUME)
% uiwait(handles.select_group_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));

% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.select_group_gui, 'handles_main', handles_main);

refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_select_group_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% This is called to update all display elements
function refresh_display(handles)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.select_group_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

%% GUI stuff
set(handles.button_select_group, 'Enable', 'on');
set(handles.button_select_curveset, 'Enable', 'off');
set(handles.button_select_curve, 'Enable', 'off');


%% Groups
if( session.Ng > 0 )       
    % Populate list box
    listbox_string = cell(1,session.Ng);
    for g = 1:session.Ng
        listbox_string{g} = sprintf('<html>(%02d) <i>%s</i></html>', ...
            g, session.groups{g}.name);
    end
    
    % Make sure user does not select an item that does not exist
    if( get(handles.listbox_groups, 'Value') > length(listbox_string) )
        set(handles.listbox_groups, 'Value', 1);
    end
    set(handles.listbox_groups, 'String', listbox_string);
    
    % Get the selected item from the listbox
    g_selected      = get(handles.listbox_groups, 'Value');
    if( g_selected > session.Ng )
        g_selected = 1;
        set(handles.listbox_groups, 'Value', g_selected);
    end
    group_selected  = session.groups{g_selected};
    
    % Table
    set(handles.table_group, 'Data', generateParamTable(group_selected) )

% No groups to display
else
    set(handles.listbox_groups, 'String', '(No groups)');
    set(handles.listbox_groups, 'Value', 1);
    group_selected = NaN;
end

%% Curvesets in group
if( isobject(group_selected) && group_selected.Ncs > 0 )
    % Curvesets in group    
    listbox_string = cell(1,group_selected.Ncs);
    for cs = 1:group_selected.Ncs
        listbox_string{cs} = sprintf('<html>(%02d) <i>%s</i></html>', ...
            cs, group_selected.curvesets{cs}.name);
    end
    
    % Make sure user does not select an item that does not exist
    if( get(handles.listbox_curvesets, 'Value') > length(listbox_string) )
        set(handles.listbox_curvesets, 'Value', 1);
    end
    set(handles.listbox_curvesets, 'String', listbox_string);

    % Get the selected item from the listbox
    cs_selected         = get(handles.listbox_curvesets, 'Value');
    if( cs_selected > group_selected.Ncs )
        cs_selected = 1;
        set(handles.listbox_curvesets, 'Value', cs_selected);
    end
    curveset_selected   = group_selected.curvesets{cs_selected};

    % Table
    set(handles.table_curveset, 'Data', generateParamTable(curveset_selected) );
    
else
    set(handles.listbox_curvesets, 'String', '(No curvesets in this group)');   
    set(handles.listbox_curvesets, 'Value', 1);
    curveset_selected = NaN;
end

%% Curves in curveset
if( isobject(curveset_selected) && curveset_selected.Nc > 0 )
    % List of items
    listbox_string = cell(1,curveset_selected.Nc);
    for c = 1:curveset_selected.Nc
        listbox_string{c} = sprintf('<html>(%02d) <i>%s</i></html>', ...
            c, curveset_selected.curves{c}.name);
    end
    
    % Make sure user does not select an item that does not exist
    if( get(handles.listbox_curves_in_curveset, 'Value') > length(listbox_string) )
        set(handles.listbox_curves_in_curveset, 'Value', 1);
    end
    set(handles.listbox_curves_in_curveset, 'String', listbox_string);

    % Get the selected item from the listbox
    c_cs_selected               = get(handles.listbox_curves_in_curveset, 'Value');
    if( c_cs_selected > curveset_selected.Nc )
        c_cs_selected = 1;
        set(handles.listbox_curves_in_curveset, 'Value', c_cs_selected);
    end
    curve_curveset_selected     = curveset_selected.curves{c_cs_selected};

    % Table of properties for selected item
    set(handles.table_curve_in_curveset, 'Data', generateParamTable(curve_curveset_selected) );
    
else
    set(handles.listbox_curves_in_curveset, 'String', '(No curves in this curveset)');
    set(handles.listbox_curves_in_curveset, 'Value', 1);
    curve_curveset_selected = NaN;
end



% --- Executes when select_group_gui is resized.
function select_group_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to select_group_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox_groups.
function listbox_groups_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_groups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_groups contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_groups
%% User picked item on listbox
refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function listbox_groups_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_groups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_datasets.
function listbox_datasets_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_datasets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_datasets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_datasets

%% User picked item on listbox
refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function listbox_datasets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_datasets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_curves_in_dataset.
function listbox_curves_in_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_curves_in_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_curves_in_dataset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_curves_in_dataset
%% User picked item on listbox
refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function listbox_curves_in_dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_curves_in_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_curves_in_curveset.
function listbox_curves_in_curveset_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_curves_in_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_curves_in_curveset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_curves_in_curveset
%% User picked item on listbox
refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function listbox_curves_in_curveset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_curves_in_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_curvesets.
function listbox_curvesets_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_curvesets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_curvesets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_curvesets
%% User picked item on listbox
refresh_display(handles);

% --- Executes during object creation, after setting all properties.
function listbox_curvesets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_curvesets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_select_group.
function button_select_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_select_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Return selection
% Call the load curves GUI using variables from main handle
handles_main    = getappdata(handles.select_group_gui, 'handles_main');

g_SELECTED      = get(handles.listbox_groups, 'Value');
setappdata(handles_main.main_gui, 'g_SELECTED', g_SELECTED);

'GROUP'
g_SELECTED

% Signal completion of the GUI by closing the window
%  Note, a waitfor() function remains from the GUI caller, GUARDD_data_manager
close(handles.select_group_gui);


% --- Executes on button press in button_select_curveset.
function button_select_curveset_Callback(hObject, eventdata, handles)
% hObject    handle to button_select_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

'curveset'

% --- Executes on button press in button_select_curve.
function button_select_curve_Callback(hObject, eventdata, handles)
% hObject    handle to button_select_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

'curve'

% Add / remove data (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/08?/? Start coding
% 2011/01/11 Convert to GUARDD program
% 2011/01/28 Start coding revised manager GUI (for classes)
% 2011/05/00 Updated for GUARDD

function varargout = GUARDD_data_manager(varargin)
% GUARDD_DATA_MANAGER M-file for GUARDD_data_manager.fig
%      GUARDD_DATA_MANAGER, by itself, creates a new GUARDD_DATA_MANAGER or raises the existing
%      singleton*.
%
%      H = GUARDD_DATA_MANAGER returns the handle to a new GUARDD_DATA_MANAGER or the handle to
%      the existing singleton*.
%
%      GUARDD_DATA_MANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_DATA_MANAGER.M with the given input arguments.
%
%      GUARDD_DATA_MANAGER('Property','Value',...) creates a new GUARDD_DATA_MANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_data_manager_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_data_manager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_data_manager

% Last Modified by GUIDE v2.5 09-Jun-2011 15:53:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_data_manager_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_data_manager_OutputFcn, ...
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


% --- Executes just before GUARDD_data_manager is made visible.
function GUARDD_data_manager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_data_manager (see VARARGIN)

% Choose default command line output for GUARDD_data_manager
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_data_manager wait for user response (see UIRESUME)
% uiwait(handles.data_manager_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));

% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.data_manager_gui, 'handles_main', handles_main);

refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_data_manager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% This is called to update all display elements
function refresh_display(handles)

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

%% GUI stuff
set(handles.button_save_dataset, 'Enable', 'off');
set(handles.button_revert_dataset, 'Enable', 'off');
set(handles.button_save_group, 'Enable', 'off');
set(handles.button_revert_group, 'Enable', 'off');
set(handles.button_save_curve_in_dataset, 'Enable', 'off');
set(handles.button_revert_curve_in_dataset, 'Enable', 'off');
set(handles.button_save_curvesets, 'Enable', 'off');
set(handles.button_revert_curvesets, 'Enable', 'off');
set(handles.button_save_curve_in_curveset, 'Enable', 'off');
set(handles.button_revert_curve_in_curveset, 'Enable', 'off');

%% Datasets
if( session.Nds > 0 )
    % Items can be removed
    set(handles.button_remove_dataset, 'Enable', 'on');
    
    % Populate list box    
    listbox_string = cell(1,session.Nds);
    for ds = 1:session.Nds
        listbox_string{ds} = sprintf('<html>(%02d) <i>%s</i></html>', ...
            ds, session.datasets{ds}.name);
    end    
    
    % Make sure user does not select an item that does not exist
    if( get(handles.listbox_datasets, 'Value') > length(listbox_string) )
        set(handles.listbox_datasets, 'Value', 1);
    end
    set(handles.listbox_datasets, 'String', listbox_string);
    
    % Get the selected item from the listbox
    ds_selected         = get(handles.listbox_datasets, 'Value');
    if( ds_selected > session.Nds )
        ds_selected = 1;
        set(handles.listbox_datasets, 'Value', ds_selected);
    end
    dataset_selected    = session.datasets{ds_selected};

    % Table    
    set(handles.table_dataset, 'Data', generateParamTable(dataset_selected) )    
    
    
    %% Curves in dataset
    if( dataset_selected.Nc > 0 )
        % Items can be removed
        set(handles.button_remove_curve_from_dataset, 'Enable', 'on');
        
        listbox_string = cell(1,dataset_selected.Nc);
        for c = 1:dataset_selected.Nc
            listbox_string{c} = sprintf('<html>(%02d) <i>%s</i></html>', ...
                c, dataset_selected.curves{c}.name);
        end
        
        % Make sure user does not select an item that does not exist
        if( get(handles.listbox_curves_in_dataset, 'Value') > length(listbox_string) )
            set(handles.listbox_curves_in_dataset, 'Value', 1);
        end
        set(handles.listbox_curves_in_dataset, 'String', listbox_string);
        
        % Get the selected item from the listbox
        c_ds_selected           = get(handles.listbox_curves_in_dataset, 'Value');
        if( c_ds_selected > dataset_selected.Nc )
            c_ds_selected = 1;
            set(handles.listbox_curves_in_dataset, 'Value', c_ds_selected);
        end
        curve_dataset_selected  = dataset_selected.curves{c_ds_selected};
        
        % Table of properties for selected item
        set(handles.table_curve_in_dataset, 'Data', generateParamTable(curve_dataset_selected) )

        % Display the selected curve
        %h = subplot(1,1,1,'Parent', handles.panel_display_curve_in_dataset);
        delete( get(handles.panel_display_curve_in_dataset, 'Children') );
        h = axes('Parent', handles.panel_display_curve_in_dataset);
        hold(h, 'all');
        
        errorbar(h, curve_dataset_selected.vcpmg, curve_dataset_selected.R2eff, curve_dataset_selected.eR2eff, '.k');
        YLIM = get(h, 'YLim');
        set(h, 'YLim', [0, YLIM(2)]);
        set(h, 'XLim', [0, max(curve_dataset_selected.vcpmg)]);
        box(h, 'on');
        
    else
        % No curves in dataset
        set(handles.listbox_curves_in_dataset, 'Value', 1);
        set(handles.listbox_curves_in_dataset, 'String', '(No curves in this dataset)');
        set(handles.button_remove_curve_from_dataset, 'Enable', 'off');
        curve_dataset_selected = NaN;        
    end
else
    set(handles.listbox_datasets, 'String', '(No datasets in this session)');
    set(handles.listbox_datasets, 'Value', 1);
    set(handles.button_remove_dataset, 'Enable', 'off');        
    dataset_selected        = NaN;
    
    set(handles.listbox_curves_in_dataset, 'String', '(No curves in this dataset)');
    set(handles.button_remove_curve_from_dataset, 'Enable', 'off');
    curve_dataset_selected = NaN;    
end

%% Groups
if( session.Ng > 0 )
    % Items can be removed
    set(handles.button_remove_group, 'Enable', 'on');
    set(handles.button_new_curveset, 'Enable', 'on');
    set(handles.button_copy_group, 'Enable', 'on');
    
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
    
    set(handles.button_remove_group, 'Enable', 'off');
    set(handles.button_new_curveset, 'Enable', 'off');
    set(handles.button_copy_group, 'Enable', 'off');    
end

%% Curvesets in group
if( isobject(group_selected) && group_selected.Ncs > 0 )
    % Items can be removed
    set(handles.button_remove_curveset, 'Enable', 'on');
    set(handles.button_copy_curveset_to_group, 'Enable', 'on');
    
    % An active curve is required to add to the curveset
    if( isobject(curve_dataset_selected) )
        set(handles.button_add_curve_to_curveset, 'Enable', 'on');
    end
    
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
    
    set(handles.button_remove_curveset, 'Enable', 'off');
    set(handles.button_add_curve_to_curveset, 'Enable', 'off');
    set(handles.button_copy_curveset_to_group, 'Enable', 'off');
end

%% Curves in curveset
if( isobject(curveset_selected) && curveset_selected.Nc > 0 )
    % Items can be removed
    set(handles.button_remove_curve_from_curveset, 'Enable', 'on');
    
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
    
    % Display the selected curve
    %h = subplot(1,1,1,'Parent', handles.panel_display_curve_in_dataset);
    delete( get(handles.panel_display_curve_in_curveset, 'Children') );
    h = axes('Parent', handles.panel_display_curve_in_curveset);
    hold(h, 'all');

    errorbar(h, curve_curveset_selected.vcpmg, curve_curveset_selected.R2eff, curve_curveset_selected.eR2eff, '.k');
    YLIM = get(h, 'YLim');
    set(h, 'YLim', [0, YLIM(2)]);
    set(h, 'XLim', [0, max(curve_curveset_selected.vcpmg)]);
    box(h, 'on');
    
else
    set(handles.listbox_curves_in_curveset, 'String', '(No curves in this curveset)');
    set(handles.listbox_curves_in_curveset, 'Value', 1);
    delete( get(handles.panel_display_curve_in_curveset, 'Children') );
    curve_curveset_selected = NaN;
    
    set(handles.button_remove_curve_from_curveset, 'Enable', 'off');
end



% --- Executes when data_manager_gui is resized.
function data_manager_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to data_manager_gui (see GCBO)
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



% --- Executes when entered data in editable cell(s) in table_dataset.
function table_dataset_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_dataset (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% If the user alters a value in the table
set(handles.button_save_dataset, 'Enable', 'on');
set(handles.button_revert_dataset, 'Enable', 'on');


% --- Executes when selected cell(s) is changed in table_dataset.
function table_dataset_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_dataset (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

%% Proceed only if there is something selected
if( isempty( eventdata.Indices ) )
    return
end

%% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
data            = getappdata(handles_main.main_gui, 'data');
ds_selected     = getappdata(handles.data_manager_gui, 'ds_selected');


% Check what the name of the row that is selcted
table_data = get(handles.table_dataset, 'Data');
row_selected = eventdata.Indices(1);
parameter_name = table_data{ row_selected, 1 };

% TODO % Intelligent want to determine if source is from a file
%% If it is the output directory
if( strcmp(parameter_name, 'outputDir') )
    % Ask user to select file
    pathname = uigetdir('.', 'Select output directory');

    % Commit the selection to the table
    if( ~isequal(pathname,0) )
        table_data{row_selected, 2} = pathname;
        set(handles.table_dataset, 'Data', table_data);
        
        set(handles.button_save, 'Enable', 'on');
        set(handles.button_revert, 'Enable', 'on');
        
        % Add the updated parameter to the list
        parameters_updated = getappdata(handles.data_manager_gui, 'parameters_updated');
        parameters_updated = [parameters_updated, parameter_name];
        parameters_updated = unique( parameters_updated );
        setappdata(handles.data_manager_gui, 'parameters_updated', parameters_updated);        
    end
end


title = 'NONE';

%% Load the desired file
if( strcmp(parameter_name, 'tauFileName') )
    title = 'Load taufile';
    file_extension = '*.txt';  
    
elseif( strcmp(parameter_name, 'nlinFileName') )
    title = 'Load nlin.tab file';
    file_extension = '*.tab';
    
elseif( strcmp(parameter_name, 'seqFileName') )
    title = 'Load sequence file';
    file_extension = '*.txt';
    
elseif( strcmp(parameter_name, 'R2TableFile') )
    title = 'Load R2 table file';
    file_extension = '*.csv';
    
elseif( strcmp(parameter_name, 'rangefile') )
    title = 'Set range file';
    file_extension = '*.txt';
end

if( ~strcmp(title, 'NONE') )
    % Ask user to select file
    [filename, pathname] = uigetfile(file_extension, title);

    % Commit the selection to the table
    if( ~isequal(filename,0) )
        filename_absolute = sprintf('%s%s', pathname, filename);        
        
        table_data{row_selected, 2} = filename_absolute;
        set(handles.table_dataset, 'Data', table_data);
        
        set(handles.button_save, 'Enable', 'on');
        set(handles.button_revert, 'Enable', 'on');
        
        % Add the updated parameter to the list
        parameters_updated = getappdata(handles.data_manager_gui, 'parameters_updated');
        parameters_updated = [parameters_updated, parameter_name];
        parameters_updated = unique( parameters_updated );
        setappdata(handles.data_manager_gui, 'parameters_updated', parameters_updated);        
    end
end


% --- Executes on button press in button_new_group.
function button_new_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Add a group

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

group = Group(session);
session.addGroup(group);

set(handles.listbox_groups, 'Value', session.Ng);

setappdata(handles_main.main_gui, 'session', session);
refresh_display(handles);


% --- Executes on button press in button_remove_group.
function button_remove_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_remove_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Remove group from session

if( strcmp(questdlg('Do you want to remove this group?', ...
        'Remove group?', 'Remove', 'Cancel', 'Cancel'),'Cancel') )
    return
end

handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

g_selected          = get(handles.listbox_groups, 'Value');
group_selected      = session.groups{g_selected};

% Remove it
session.removeGroup(group_selected);
refresh_display(handles);

% --- Executes on button press in button_remove_dataset.
function button_remove_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to button_remove_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Remove dataset from session
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% Assume selected dataset is valid
ds_selected         = get(handles.listbox_datasets, 'Value');
dataset_selected    = session.datasets{ds_selected};

session.removeDataset(dataset_selected);
refresh_display(handles);


% --- Executes on button press in button_add_dataset.
function button_add_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to button_add_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Add a dataset

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

dataset = Dataset(session);
session.addDataset(dataset);

setappdata(handles_main.main_gui, 'session', session);
refresh_display(handles);



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


% --- Executes on button press in button_add_curve_to_curveset.
function button_add_curve_to_curveset_Callback(hObject, eventdata, handles)
% hObject    handle to button_add_curve_to_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Add a curve from the dataset to the curveset in the group

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% (1) Get curveset in the group
% Assume selected group is valid
g_selected          = get(handles.listbox_groups, 'Value');
group_selected      = session.groups{g_selected};

% Assume selected curveset is valid
cs_selected         = get(handles.listbox_curvesets, 'Value');
curveset_selected   = group_selected.curvesets{cs_selected};

% (2) Get curve in the dataset
% Assume selected dataset is valid
ds_selected         = get(handles.listbox_datasets, 'Value');
dataset_selected    = session.datasets{ds_selected};

% Assume selected curve is valid
c_selected          = get(handles.listbox_curves_in_dataset, 'Value');
curve_selected      = dataset_selected.curves{c_selected};

% Add it to the curveset ONLY if it is not there already
if( ~curveset_selected.containsCurve(curve_selected) )    
    curveset_selected.addCurve(curve_selected);
    %setappdata(handles_main.main_gui, 'session', session);
    
    CONSTRAIN_RATE_ANALYSIS = false;
    group_selected.updateFitParams(CONSTRAIN_RATE_ANALYSIS);
else
	msgbox('This curve is already in the set', 'Curve in set', 'Error');
end

refresh_display(handles);

% --- Executes on button press in button_remove_curve_from_curveset.
function button_remove_curve_from_curveset_Callback(hObject, eventdata, handles)
% hObject    handle to button_remove_curve_from_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Remove curve from curveset
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

g_selected          = get(handles.listbox_groups, 'Value');
group_selected      = session.groups{g_selected};

cs_selected         = get(handles.listbox_curvesets, 'Value');
curveset_selected   = group_selected.curvesets{cs_selected};

c_selected          = get(handles.listbox_curves_in_curveset, 'Value');
curve_selected      = curveset_selected.curves{c_selected};

% Remove it
curveset_selected.removeCurve(curve_selected);
CONSTRAIN_RATE_ANALYSIS = false;
group_selected.updateFitParams(CONSTRAIN_RATE_ANALYSIS);
refresh_display(handles);

% --- Executes on button press in button_remove_curveset.
function button_remove_curveset_Callback(hObject, eventdata, handles)
% hObject    handle to button_remove_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Remove curveset from group
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

g_selected          = get(handles.listbox_groups, 'Value');
group_selected      = session.groups{g_selected};

cs_selected         = get(handles.listbox_curvesets, 'Value');
curveset_selected   = group_selected.curvesets{cs_selected};

% Remove it
group_selected.removeCurveset(curveset_selected);
% AUTOMATIC -> group_selected.updateFitParams(CONSTRAIN_RATE_ANALYSIS);
refresh_display(handles);

% --- Executes on button press in button_new_curveset.
function button_new_curveset_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Add a curveset

% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Assume selected group is valid
g_selected      = get(handles.listbox_groups, 'Value');
group_selected  = session.groups{g_selected};

curveset = Curveset;
group_selected.addCurveset(curveset);
% AUTOMATIC -> group_selected.updateFitParams(CONSTRAIN_RATE_ANALYSIS);

setappdata(handles_main.main_gui, 'session', session);
refresh_display(handles);

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


% --------------------------------------------------------------------
function menu_load_sequence_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_sequence_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Load sequence file

% User selects a file
[filename, pathname]    = uigetfile('*.txt', 'Select sequence file');
sequence_filename       = sprintf('%s/%s', pathname, filename);

if( exist(sequence_filename, 'file') )

    % Extract the sequence_array, if possible
    FILE = fopen(sequence_filename);

    % Read the line of text from the file
    inputText       = textscan( FILE, '%s' );
    sequence_array  = inputText{1};

    fclose(FILE);

    % Access main window handle which is stored in this window's main handle
    handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
    session         = getappdata(handles_main.main_gui, 'session');

    % Commit these changes to the session
    session.setSequence( sequence_array );
    refresh_display(handles);

else
    % File does not exist (report message if desired)

end


% --------------------------------------------------------------------
function menu_input_Callback(hObject, eventdata, handles)
% hObject    handle to menu_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_remove_curve_from_dataset.
function button_remove_curve_from_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to button_remove_curve_from_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Remove curve from dataset (and hence all instances of the curve)
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Get the selected item from the listbox
ds_selected         = get(handles.listbox_datasets, 'Value');
dataset_selected    = session.datasets{ds_selected};

% Get the selected item from the listbox
c_ds_selected           = get(handles.listbox_curves_in_dataset, 'Value');
curve_dataset_selected  = dataset_selected.curves{c_ds_selected};

% Remove the curve
dataset_selected.removeCurve(curve_dataset_selected);
refresh_display(handles);

% --- Executes on button press in button_load_curves.
function button_load_curves_Callback(hObject, eventdata, handles)
% hObject    handle to button_load_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Load curves via GUI

% Call the load curves GUI using variables from main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');

% Pause GUI execution until the data are loaded
h = GUARDD_load_curves('GUARDD', handles_main);
waitfor(h);

% Get the specs to load the curves
LoadCurves      = getappdata(handles_main.main_gui, 'LoadCurves');
session         = getappdata(handles_main.main_gui, 'session');

% NMRPipe format
if( isstruct(LoadCurves) && strcmp(LoadCurves.MODE, 'NMRPIPE') )
    nlin_filename    = LoadCurves.nlin_filename;
    vcpmg_filename   = LoadCurves.vcpmg_filename;
    
    % Assume selected dataset is valid
    ds_selected         = get(handles.listbox_datasets, 'Value');
    dataset_selected    = session.datasets{ds_selected};

    dataset_selected.readNlin(nlin_filename, vcpmg_filename);

    %setappdata(handles_main.main_gui, 'session', session);
    refresh_display(handles);

else
    
end


% --- Executes on button press in button_save_group.
function button_save_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Save data from the table
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');
g_selected          = get(handles.listbox_groups, 'Value');
group_selected      = session.groups{g_selected};

% Read and save data from table
table_data          = get(handles.table_group, 'Data');
group_selected.saveDataFromTable(table_data);

refresh_display(handles);

% --- Executes on button press in button_revert_group.
function button_revert_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_revert_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Abort changes to table
refresh_display(handles);

% --- Executes on button press in button_save_dataset.
function button_save_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Save data from the table
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');
ds_selected         = get(handles.listbox_datasets, 'Value');
dataset_selected    = session.datasets{ds_selected};


table_data          = get(handles.table_dataset, 'Data');
dataset_selected.saveDataFromTable(table_data);

%setappdata(handles_main.main_gui, 'session', session);
refresh_display(handles);

% --- Executes on button press in button_revert_dataset.
function button_revert_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to button_revert_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Abort changes to table
refresh_display(handles);

% --- Executes on button press in button_save_curve_in_dataset.
function button_save_curve_in_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_curve_in_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Save data from the table
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% Get the selected item from the listbox
ds_selected         = get(handles.listbox_datasets, 'Value');
dataset_selected    = session.datasets{ds_selected};

% Get the selected item from the listbox
c_ds_selected           = get(handles.listbox_curves_in_dataset, 'Value');
curve_dataset_selected  = dataset_selected.curves{c_ds_selected};


table_data          = get(handles.table_curve_in_dataset, 'Data');
curve_dataset_selected.saveDataFromTable(table_data);

CONSTRAIN_RATE_ANALYSIS = false;
group_selected.updateFitParams(CONSTRAIN_RATE_ANALYSIS);
%setappdata(handles_main.main_gui, 'session', session);
refresh_display(handles);

% --- Executes on button press in button_revert_curve_in_dataset.
function button_revert_curve_in_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to button_revert_curve_in_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Abort changes to table
refresh_display(handles);

% --- Executes on button press in button_save_curve_in_curveset.
function button_save_curve_in_curveset_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_curve_in_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Save data from the table
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');
g_selected          = get(handles.listbox_groups, 'Value');
group_selected      = session.groups{g_selected};

cs_selected         = get(handles.listbox_curvesets, 'Value');
curveset_selected   = group_selected.curvesets{cs_selected};

c_selected          = get(handles.listbox_curves_in_curveset, 'Value');
curve_selected      = curveset_selected.curves{c_selected};

% Read and save data from table
table_data          = get(handles.table_curve_in_curveset, 'Data');
curve_selected.saveDataFromTable(table_data);

CONSTRAIN_RATE_ANALYSIS = false;
group_selected.updateFitParams(CONSTRAIN_RATE_ANALYSIS);

refresh_display(handles);

% --- Executes on button press in button_revert_curve_in_curveset.
function button_revert_curve_in_curveset_Callback(hObject, eventdata, handles)
% hObject    handle to button_revert_curve_in_curveset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Abort changes to table
refresh_display(handles);

% --- Executes on button press in button_save_curvesets.
function button_save_curvesets_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_curvesets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Save data from the table
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');
g_selected          = get(handles.listbox_groups, 'Value');
group_selected      = session.groups{g_selected};

cs_selected         = get(handles.listbox_curvesets, 'Value');
curveset_selected   = group_selected.curvesets{cs_selected};

% Read and save data from table
table_data          = get(handles.table_curveset, 'Data');
curveset_selected.saveDataFromTable(table_data);

refresh_display(handles);

% --- Executes on button press in button_revert_curvesets.
function button_revert_curvesets_Callback(hObject, eventdata, handles)
% hObject    handle to button_revert_curvesets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Abort changes to table
refresh_display(handles);

% --- Executes when entered data in editable cell(s) in table_curve_in_dataset.
function table_curve_in_dataset_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_curve_in_dataset (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% If the user alters a value in the table
set(handles.button_save_curve_in_dataset, 'Enable', 'on');
set(handles.button_revert_curve_in_dataset, 'Enable', 'on');


% --- Executes when entered data in editable cell(s) in table_group.
function table_group_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_group (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% If the user alters a value in the table
set(handles.button_save_group, 'Enable', 'on');
set(handles.button_revert_group, 'Enable', 'on');


% --- Executes when entered data in editable cell(s) in table_curveset.
function table_curveset_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_curveset (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% If the user alters a value in the table
set(handles.button_save_curvesets, 'Enable', 'on');
set(handles.button_revert_curvesets, 'Enable', 'on');


% --- Executes when entered data in editable cell(s) in table_curve_in_curveset.
function table_curve_in_curveset_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_curve_in_curveset (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%% If the user alters a value in the table
set(handles.button_save_curve_in_curveset, 'Enable', 'on');
set(handles.button_revert_curve_in_curveset, 'Enable', 'on');


% --------------------------------------------------------------------
function menu_generate_groups_Callback(hObject, eventdata, handles)
% hObject    handle to menu_generate_groups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Automatically generate groups via NMR probe
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

session.generateGroups();
refresh_display(handles);


% --------------------------------------------------------------------
function menu_load_data_script_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_data_script (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Load data via script
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% User selects a file
[filename, pathname]    = uigetfile('*.txt', 'Select sequence file');
script_filename         = sprintf('%s/%s', pathname, filename);

%script_filename='/home/ian/foster_lab/00-working/GUARDD/GUARDD/data//TRAP-GUARDD-Load_Data.txt'

if( exist(script_filename, 'file') )
    session.loadDatasets(script_filename);    
else    
    % Report error message if desired    
end

refresh_display(handles);


% --------------------------------------------------------------------
function menu_group_Callback(hObject, eventdata, handles)
% hObject    handle to menu_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_group_sort_Callback(hObject, eventdata, handles)
% hObject    handle to menu_group_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Load data via script
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

session.sortGroups();
refresh_display(handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_dataset_sort_curves_Callback(hObject, eventdata, handles)
% hObject    handle to menu_dataset_sort_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Load data via script
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% Get selected dataset
ds = get(handles.listbox_datasets, 'Value');
dataset = session.datasets{ds};

dataset.sortCurves();
refresh_display(handles);

% --------------------------------------------------------------------
function menu_group_autogenerate_Callback(hObject, eventdata, handles)
% hObject    handle to menu_group_autogenerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Automatically generate groups via NMR probe
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

session.generateGroups();
refresh_display(handles);


% --------------------------------------------------------------------
function menu_dataset_sort_curves_all_Callback(hObject, eventdata, handles)
% hObject    handle to menu_dataset_sort_curves_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Load data via script
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% Get selected dataset
for ds = 1:session.Nds
    dataset = session.datasets{ds};
    dataset.sortCurves();
end

refresh_display(handles);


% --- Executes on button press in button_copy_group.
function button_copy_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_copy_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Copy current group
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% Get the selected item from the listbox
g_selected      = get(handles.listbox_groups, 'Value');
group_selected  = session.groups{g_selected};

% Copy it and add it to the session
group_copy = group_selected.copy();
session.addGroup(group_copy);

% Pick the new group
set(handles.listbox_groups, 'Value', session.Ng);
refresh_display(handles);


% --------------------------------------------------------------------
function menu_curvesets_sort_Callback(hObject, eventdata, handles)
% hObject    handle to menu_curvesets_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Automatically generate curvesets in the group via NMR probe
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

if( session.Ng == 0 )
    errordlg('Must create/select group first');
end

% Get the selected item from the listbox
g_selected      = get(handles.listbox_groups, 'Value');
group_selected  = session.groups{g_selected};

if( group_selected.Nf > 0 )
    errordlg('GUARDD will not sort curvesets and/or curves because this will invalidate prior fits (puts them out of order). Solution: copy this group, then sort, then fit.');

    
else
    % Sort them
    group_selected.sortCurvesets();
    refresh_display(handles);
end


% --------------------------------------------------------------------
function menu_curveset_autogenerate_Callback(hObject, eventdata, handles)
% hObject    handle to menu_curveset_autogenerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Automatically generate curvesets in the group via NMR probe
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

if( session.Ng == 0 )
    errordlg('Must create/select group first');
end

% Get the selected item from the listbox
g_selected      = get(handles.listbox_groups, 'Value');
group_selected  = session.groups{g_selected};

session.generateCurvesetsForGroup(group_selected);
refresh_display(handles);


% --- Executes on button press in button_generate_groups_specs.
function button_generate_groups_specs_Callback(hObject, eventdata, handles)
% hObject    handle to button_generate_groups_specs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Automatically generate groups via NMR probe
handles_main        = getappdata(handles.data_manager_gui, 'handles_main');
session             = getappdata(handles_main.main_gui, 'session');

% Check the table data to make the proper specifications
%table_data = get(handles.table_generate_group_specs, 'Data');

session.generateGroups('TEMP', 310);
refresh_display(handles);


% --------------------------------------------------------------------
function menu_curveset_copy_to_group_Callback(hObject, eventdata, handles)
% hObject    handle to menu_curveset_copy_to_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Copy the selected curveset to a group

% Call the load curves GUI using variables from main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Pause GUI execution until the data are loaded
h = GUARDD_select_group('GUARDD', handles_main);
waitfor(h);

% Get the specs to load the curves
g_SELECTED          = getappdata(handles_main.main_gui, 'g_SELECTED');
destination_group   = session.groups{g_SELECTED};    

% Get the selected CURVESET item from the listbox
g           = get(handles.listbox_groups, 'Value');
group       = session.groups{g};    
cs          = get(handles.listbox_curvesets, 'Value');
curveset    = group.curvesets{cs};

% Now copy the curveset to that group
curveset_copy       = curveset.copy();
destination_group.addCurveset(curveset_copy);

% Choose the destimation group from the listbox
set(handles.listbox_groups, 'Value', g_SELECTED);
refresh_display(handles);


% --- Executes on button press in button_copy_curveset_to_group.
function button_copy_curveset_to_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_copy_curveset_to_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Copy the selected curveset to a group

% Call the load curves GUI using variables from main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Pause GUI execution until the data are loaded
h = GUARDD_select_group('GUARDD', handles_main);
waitfor(h);

% Get the specs to load the curves
g_SELECTED          = getappdata(handles_main.main_gui, 'g_SELECTED');
destination_group   = session.groups{g_SELECTED};    

% Get the selected CURVESET item from the listbox
g           = get(handles.listbox_groups, 'Value');
group       = session.groups{g};    
cs          = get(handles.listbox_curvesets, 'Value');
curveset    = group.curvesets{cs};

% Now copy the curveset to that group
curveset_copy       = curveset.copy();
destination_group.addCurveset(curveset_copy);

% Choose the destimation group from the listbox
set(handles.listbox_groups, 'Value', g_SELECTED);
refresh_display(handles);


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_output_datasets_Callback(hObject, eventdata, handles)
% hObject    handle to menu_output_datasets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Output datasets to a file that can be loaded as a script

% Call the load curves GUI using variables from main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');
VERSION         = getappdata(handles_main.main_gui, 'VERSION');

% Pick filename
default_filename    = sprintf('output-GUARDD-Datasets--%s.txt', datestr(now, 'yyyy.mm.dd-HH-MM'));

if( ~exist(session.outputDir, 'dir') )
    mkdir(session.outputDir);
end

[filename, filepath] = uiputfile('*.txt', 'Save results to file', ...
    sprintf('./%s/%s', session.outputDir, default_filename) );

% Save it
if( isequal(filename,0) )
    return
end
FILE = fopen(sprintf('%s/%s', filepath, filename), 'w');

fprintf('\nWriting %s...', filename);
t0 = tic();

fprintf(FILE, '# GUARDD %s', VERSION);
fprintf(FILE, '\n# (C) Ian Kleckner 2010-2011');
fprintf(FILE, '\n# Exported RD Data from datasets');
fprintf(FILE, '\n# This file can be loaded as a script via Data Manager');
fprintf(FILE, '\n# Created on %s', datestr(now, 'yyyy.mm.dd-HH-MM'));

session.exportDatasets(FILE);
fprintf('Done! (%0.1f sec)', toc(t0));

% --------------------------------------------------------------------
function menu_output_groups_Callback(hObject, eventdata, handles)
% hObject    handle to menu_output_groups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Output groups to a text file that can be loaded as a script

% Call the load curves GUI using variables from main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');
VERSION         = getappdata(handles_main.main_gui, 'VERSION');

% Pick filename
default_filename    = sprintf('output-GUARDD-Groups--%s.txt', datestr(now, 'yyyy.mm.dd-HH-MM'));

if( ~exist(session.outputDir, 'dir') )
    mkdir(session.outputDir);
end

[filename, filepath] = uiputfile('*.txt', 'Save results to file', ...
    sprintf('./%s/%s', session.outputDir, default_filename) );

% Save it
if( isequal(filename,0) )
    return
end
FILE = fopen(sprintf('%s/%s', filepath, filename), 'w');

fprintf('\nWriting %s...', filename);
t0 = tic();

fprintf(FILE, 'GUARDD %s', VERSION);
fprintf(FILE, '\n(C) Ian Kleckner 2010-2011');
fprintf(FILE, '\nExported groups made from datasets');
fprintf(FILE, '\nCreated on %s', datestr(now, 'yyyy.mm.dd-HH-MM'));

session.outputGroups(FILE);
fprintf('Done! (%0.1f sec)\n', toc(t0));


% --------------------------------------------------------------------
function menu_generate_groups_constraints_Callback(hObject, eventdata, handles)
% hObject    handle to menu_generate_groups_constraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Generate groups using constraints

% Call the load curves GUI using variables from main handle
handles_main    = getappdata(handles.data_manager_gui, 'handles_main');

% Pause GUI execution until the data are loaded
h = GUARDD_create_groups('GUARDD', handles_main);
waitfor(h);
refresh_display(handles);

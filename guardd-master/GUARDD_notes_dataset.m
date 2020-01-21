% View / eidt notes on entire dataset via dedicated window (requires .fig file)
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/02/01 Start coding
% 2011/01/11 Convert to GUARDD program
% 
% TO DO
%  Nothing

function varargout = GUARDD_notes_dataset(varargin)
% GUARDD_NOTES_DATASET M-file for GUARDD_notes_dataset.fig
%      GUARDD_NOTES_DATASET, by itself, creates a new GUARDD_NOTES_DATASET or raises the existing
%      singleton*.
%
%      H = GUARDD_NOTES_DATASET returns the handle to a new GUARDD_NOTES_DATASET or the handle to
%      the existing singleton*.
%
%      GUARDD_NOTES_DATASET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_NOTES_DATASET.M with the given input arguments.
%
%      GUARDD_NOTES_DATASET('Property','Value',...) creates a new GUARDD_NOTES_DATASET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_notes_dataset_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_notes_dataset_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_notes_dataset

% Last Modified by GUIDE v2.5 11-Jan-2011 16:50:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_notes_dataset_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_notes_dataset_OutputFcn, ...
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


% --- Executes just before GUARDD_notes_dataset is made visible.
function GUARDD_notes_dataset_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_notes_dataset (see VARARGIN)

% Choose default command line output for GUARDD_notes_dataset
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_notes_dataset wait for user response (see UIRESUME)
% uiwait(handles.notes_dataset_gui);

%% Save handle from GUI that called this GUI
% Find main GUI string in list of input arguments
index_main_gui_input    = find(strcmp(varargin, 'GUARDD'));
% The actual handle is the index after (+1) the name
handles_main            = varargin{index_main_gui_input+1};

% Store the main window's handle in this window's data
% Now this window can access all variables, etc. from main window
setappdata(handles.notes_dataset_gui, 'handles_main', handles_main);

%% Initialize notes if needed
data = getappdata(handles_main.main_gui, 'data');
if( ~isfield(data, 'notes_dataset') )
    data.Nnotes         = 0;
    data.notes_dataset  = [];
    setappdata(handles_main.main_gui, 'data', data);
end

% Diable save and restore buttons unless text is altered
set(handles.button_save_note, 'Enable', 'off');
set(handles.button_restore_note, 'Enable', 'off');

refresh_display(handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_notes_dataset_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% This is called to update all display elements
function refresh_display(handles)
    
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.notes_dataset_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

%{
initialize notes in opening function
displya nots in refresh_display

save notes, tec.

%}

Nnotes = session.Nnotes;

if( Nnotes > 0 || session.Ng > 0 )
    %% Populate list with note titles
    notes_title{1}          = 'View all';
    notes_title{Nnotes+1}   = ' ';
    
    for n = 1:Nnotes
        note = session.notes{n};
        notes_title{n+1} = note.name;
    end    
    set(handles.popup_notes,'String',notes_title);    
     
    % Get selected note from list
    n = get(handles.popup_notes, 'Value')-1;
    
    % Display information on current note    
    % This will show all notes
    if( n == 0 )
        
        % Concatenate all notes for simultaneous display
        all_notes   = {};
        line        = 0;
        
        for n = 1:Nnotes
            note = session.notes{n};
            
            % Vertical string array yields clean multiple lines of text
            % Other formatting messes up multiline text for some reason
%             all_notes = strvcat( all_notes, ...
%                 sprintf('--== %s [%s]', ...
%                     session.note_names{n}, ...
%                     datestr(session.note_times{n}, 'mmm dd, yyyy at HH:MM:SS')), ...                    
%                     session.notes{n}, ...
%                 ' ' );

            if( line > 0 )
                line = line+1;
                all_notes{line}     = '';
            end
            
            line = line+1;
            all_notes{line}     = sprintf('------- %s [%s] -------', note.name, ...
                                    datestr(note.time, 'mmm dd, yyyy at HH:MM:SS'));
            %line = line+1;
            %all_notes{line}     = session.notes{n};
            
            % Must add each line from the session to a new line in all_notes
            Nlines_n        = max(size(note.text));
            for l_n = 1:Nlines_n
                line = line+1;
                all_notes{line}     = note.text{l_n};
            end
        end
        
        % Concatenate notes for each group
        for g = 1:session.Ng
            group = session.groups{g};
            
            if( line > 0 )
                line = line+1;
                all_notes{line}     = '';
            end
            
            line = line+1;
            all_notes{line}     = sprintf('------ %s ------', group.name );
                
            % Must add each line from the group to a new line in all_notes
            Nlines_g = max(size(group.note.text));
            for l_g = 1:Nlines_g
                line = line+1;
                all_notes{line}     = group.note.text{l_g};
            end
        end
        
        % Display the final result
        set(handles.edit_notes, 'String', all_notes);
        
        % Disable time of note
        set( handles.text_time, 'String', ''); 

        % Title
        set(handles.panel_notes, 'Title', sprintf('All notes [read only]'));
        set(handles.edit_name, 'String', '');
        set(handles.edit_name, 'Enable', 'off');
        
        % Buttons
        set(handles.button_save_note, 'Enable', 'off');
        set(handles.button_restore_note, 'Enable', 'off');
        set(handles.button_delete_note, 'Enable', 'off');
        set(handles.text_time, 'Enable', 'off');
        
    % Show the nth note
    else
        % Show content of note
        note = session.notes{n};
        
        %set(handles.edit_notes, 'String', session.displayNote(n));
        set(handles.edit_notes, 'String', note.text );
        
        % Show time of note
        set( handles.text_time, 'String', sprintf('%s', ...
             datestr(note.time, 'mmm dd, yyyy\nHH:MM:SS')) );

        % Title
        set(handles.panel_notes, 'Title', sprintf('Note: %s', note.name));
        set(handles.edit_name, 'String', note.name);
        set(handles.edit_name, 'Enable', 'on');
        
        % Buttons
        set(handles.button_delete_note, 'Enable', 'on');
        set(handles.text_time, 'Enable', 'on');
    end
    
    %% Enable GUI elements    
    set(handles.popup_notes, 'Enable', 'on');    
    set(handles.button_save_note, 'Enable', 'off');
    set(handles.button_restore_note, 'Enable', 'off');    
    %set(handles.button_timestamp, 'Enable', 'on');
    

% No notes
else
    % Disable list and buttons
    set(handles.popup_notes, 'Enable', 'off');    
    set(handles.button_save_note, 'Enable', 'off');
    set(handles.button_restore_note, 'Enable', 'off');
    set(handles.button_delete_note, 'Enable', 'off');
    set(handles.button_timestamp, 'Enable', 'off');
    set(handles.edit_name, 'Enable', 'off');
    set(handles.text_time, 'Enable', 'off');
    set(handles.panel_notes, 'Title', sprintf('No notes to show'));    
end
    

% --- Executes when notes_dataset_gui is resized.
function notes_dataset_gui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to notes_dataset_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Resize each element of the GUI

% Get size of the GUI window which is just resized
fP = get(handles.notes_dataset_gui, 'Position');

% Main panel (scale size, fix relative position)
left    = 0.5;
bottom  = 0.1;
width   = fP(3)-left*2;
height  = fP(4)-8;
set(handles.panel_notes, 'Position', [left, bottom, width, height]);

% Notes within main panel (scale size, fix relative position)
left    = 0.1;
bottom  = 0.1;
width   = fP(3)-(0.5+left)*4;
height  = fP(4)-8-2.2;
set(handles.edit_notes, 'Position', [left, bottom, width, height]);

% Adjust position of button (scale position, fix size)
width   = 59;
height  = 1.5;
left    = 0.5;
bottom  = fP(4)-2;
set(handles.popup_notes, 'Position', [left, bottom, width, height]);

%% Buttons

% Adjust position of button (scale position, fix size)
width   = 10;
height  = 2;
left    = 0.1;
bottom  = fP(4)-4.2;
set(handles.button_new_note, 'Position', [left, bottom, width, height]);

% Adjust position of button (scale position, fix size)
width   = 10;
height  = 2;
left    = 0.1+10+0.5;
bottom  = fP(4)-4.2;
set(handles.button_save_note, 'Position', [left, bottom, width, height]);

% Adjust position of button (scale position, fix size)
width   = 10;
height  = 2;
left    = 0.1+10+0.5+10+0.5;
bottom  = fP(4)-4.2;
set(handles.button_restore_note, 'Position', [left, bottom, width, height]);

% Adjust position of button (scale position, fix size)
width   = 10;
height  = 2;
left    = 0.1+10+0.5+10+0.5+10+0.5;
bottom  = fP(4)-4.2;
set(handles.button_delete_note, 'Position', [left, bottom, width, height]);

% Adjust position of button (scale position, fix size)
width   = 17;
height  = 2;
left    = 0.1+10+0.5+10+0.5+10+0.5+10+0.5;
bottom  = fP(4)-4.2;
set(handles.button_timestamp, 'Position', [left, bottom, width, height]);

%% Text boxes

% Adjust position of button (scale position, fix size)
width   = 32;
height  = 2;
left    = 0.1;
bottom  = fP(4)-6.5;
set(handles.edit_name, 'Position', [left, bottom, width, height]);

% Adjust position of button (scale position, fix size)
width   = 25;
height  = 2.2;
left    = 0.1+32+1;
bottom  = fP(4)-7;
set(handles.text_time, 'Position', [left, bottom, width, height]);

%{
% Adjust position of checkbox (scale position, fix size)
width   = 15;
height  = 1.5;
left    = 29 - width/2;
bottom  = fP(4)-2.2;
set(handles.checkbox_best_fit_is_ok, 'Position', [left, bottom, width, height]);

% Adjust position of checkbox (scale position, fix size)
width   = 25;
height  = 1.5;
left    = 48 - width/2;
bottom  = fP(4)-2.2;
set(handles.checkbox_exhibits_exchange, 'Position', [left, bottom, width, height]);
%}


function edit_notes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_notes as text
%        str2double(get(hObject,'String')) returns contents of edit_notes as a double

%% Enable save and restore buttons
n = get(handles.popup_notes, 'Value')-1;

% Cannot save or edit the "all notes"
if( n == 0 )
    refresh_display(handles);
    
else
    set(handles.button_save_note, 'Enable', 'on');
    set(handles.button_restore_note, 'Enable', 'on');
end


% --- Executes during object creation, after setting all properties.
function edit_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_save_note.
function button_save_note_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_note (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Save the new text notes and title the user entered
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.notes_dataset_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

n = get(handles.popup_notes, 'Value')-1;

% Cannot save the "View all" note
if( n > 0 )
    text = get(handles.edit_notes, 'String');
    name = get(handles.edit_name, 'String');

    session.setNote(n, text, name);

    set(handles.button_save_note, 'Enable', 'off');
    set(handles.button_restore_note, 'Enable', 'off');

    refresh_display(handles);
end


% --- Executes on button press in button_restore_note.
function button_restore_note_Callback(hObject, eventdata, handles)
% hObject    handle to button_restore_note (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% This will erase what the user typed before saving it

set(handles.button_save_note, 'Enable', 'off');
set(handles.button_restore_note, 'Enable', 'off');

refresh_display(handles);

% --- Executes on key press with focus on edit_notes and none of its controls.
function edit_notes_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%% If user starts typing in notes panel

% Enable save and notes buttons
set(handles.button_save_note, 'Enable', 'on');
set(handles.button_restore_note, 'Enable', 'on');


% --- Executes on button press in button_timestamp.
function button_timestamp_Callback(hObject, eventdata, handles)
% hObject    handle to button_timestamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Add timestamp to notes
current_note = get(handles.edit_notes, 'String');

set( handles.edit_notes, 'String', sprintf('%s%s', ...
     current_note, datestr(now, '[mmm dd,yyyy at HH:MM:SS]')) );
 
 % Don't necessarily save the notes though


% --- Executes on selection change in popup_notes.
function popup_notes_Callback(hObject, eventdata, handles)
% hObject    handle to popup_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_notes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_notes

%% Change user selection

refresh_display(handles);


% --- Executes during object creation, after setting all properties.
function popup_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_new_note.
function button_new_note_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_note (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Add a new note to the end of the list
% Access main window handle which is stored in this window's main handle
handles_main    = getappdata(handles.notes_dataset_gui, 'handles_main');
session         = getappdata(handles_main.main_gui, 'session');

% Set up default new notes
session.addNote( '[Type notes here]', '[Default name]');

% Update GUI
set(handles.popup_notes, 'Value', session.Nnotes+1);
refresh_display(handles);


% --- Executes on button press in button_delete_note.
function button_delete_note_Callback(hObject, eventdata, handles)
% hObject    handle to button_delete_note (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Remove the selected note
n = get(handles.popup_notes, 'Value')-1;

% Ask user if they want to delete fit
if( n ~= 0 && strcmp(questdlg('Are you sure you want to delete this note?', ...
    'Delete note?', 'Delete note', 'Cancel','Cancel'),'Delete note') )

    % Access main window handle which is stored in this window's main handle
    handles_main    = getappdata(handles.notes_dataset_gui, 'handles_main');
    session         = getappdata(handles_main.main_gui, 'session');
    
    session.removeNote( n );
    set(handles.popup_notes, 'Value', n);

    refresh_display(handles);
end

% --- Executes on key press with focus on edit_name and none of its controls.
function edit_name_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_name (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
%% If user starts typing in notes panel

% Enable save and notes buttons
set(handles.button_save_note, 'Enable', 'on');
set(handles.button_restore_note, 'Enable', 'on');

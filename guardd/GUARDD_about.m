% Describes GUARDD software and display license
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/07/01 Start coding

function varargout = GUARDD_about(varargin)
% GUARDD_ABOUT MATLAB code for GUARDD_about.fig
%      GUARDD_ABOUT, by itself, creates a new GUARDD_ABOUT or raises the existing
%      singleton*.
%
%      H = GUARDD_ABOUT returns the handle to a new GUARDD_ABOUT or the handle to
%      the existing singleton*.
%
%      GUARDD_ABOUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUARDD_ABOUT.M with the given input arguments.
%
%      GUARDD_ABOUT('Property','Value',...) creates a new GUARDD_ABOUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUARDD_about_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUARDD_about_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUARDD_about

% Last Modified by GUIDE v2.5 01-Jul-2011 12:05:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUARDD_about_OpeningFcn, ...
                   'gui_OutputFcn',  @GUARDD_about_OutputFcn, ...
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


% --- Executes just before GUARDD_about is made visible.
function GUARDD_about_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUARDD_about (see VARARGIN)

% Choose default command line output for GUARDD_about
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUARDD_about wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% Draw the GUI elements

% License
LICENSE_FILE = '00-LICENSE-GPL3.txt';
FILE = fopen(LICENSE_FILE, 'r');
text = textscan(FILE,'%s', 'delimiter', '\n', 'whitespace', '');
fclose(FILE);
set(handles.edit_license, 'String', text{1});


% Version
% Find main GUI string in list of input arguments
index_version   = find(strcmp(varargin, 'VERSION'));

% The actual handle is the index after (+1) the name
VERSION         = varargin{index_version+1};
set(handles.text_version, 'String', sprintf('%s',VERSION));

% References
CITATION_HTML = '<html>Kleckner, I. R., & Foster, M. P. (submitted). <b>GUARDD: User-friendly MATLAB software for rigorous analysis of CPMG RD NMR data.</b< <i>J. Biomol. NMR.</i>';
CITATION = 'Kleckner, I. R., & Foster, M. P. (submitted). GUARDD: User-friendly MATLAB software for rigorous analysis of CPMG RD NMR data. J. Biomol. NMR.';
set(handles.text_citation, 'String', CITATION);



% --- Outputs from this function are returned to the command line.
function varargout = GUARDD_about_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_license_Callback(hObject, eventdata, handles)
% hObject    handle to edit_license (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_license as text
%        str2double(get(hObject,'String')) returns contents of edit_license as a double


% --- Executes during object creation, after setting all properties.
function edit_license_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_license (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_citation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_citation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_citation as text
%        str2double(get(hObject,'String')) returns contents of edit_citation as a double


% --- Executes during object creation, after setting all properties.
function edit_citation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_citation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

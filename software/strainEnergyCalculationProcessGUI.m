function varargout = strainEnergyCalculationProcessGUI(varargin)
% strainEnergyCalculationProcessGUI M-file for strainEnergyCalculationProcessGUI.fig
%      strainEnergyCalculationProcessGUI, by itself, creates a new strainEnergyCalculationProcessGUI or raises the existing
%      singleton*.
%
%      H = strainEnergyCalculationProcessGUI returns the handle to a new strainEnergyCalculationProcessGUI or the handle to
%      the existing singleton*.
%
%      strainEnergyCalculationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in strainEnergyCalculationProcessGUI.M with the given input arguments.
%
%      strainEnergyCalculationProcessGUI('Property','Value',...) creates a new strainEnergyCalculationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before strainEnergyCalculationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to strainEnergyCalculationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
%
% This file is part of TFM_Package.
% 
% TFM_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TFM_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TFM_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Edit the above text to modify the response to help strainEnergyCalculationProcessGUI

% Last Modified by GUIDE v2.5 15-May-2017 10:00:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @strainEnergyCalculationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @strainEnergyCalculationProcessGUI_OutputFcn, ...
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


% --- Executes just before strainEnergyCalculationProcessGUI is made visible.
function strainEnergyCalculationProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

set(handles.checkbox_useFOV,'Value',funParams.useFOV);
set(handles.checkbox_cellMask,'Value',funParams.useCellMask);
set(handles.checkbox_performForceBlobAnal,'Value',funParams.performForceBlobAnalysis);
set(handles.checkbox_exportCSV,'Value',funParams.exportCSV);

% Choose default command line output for strainEnergyCalculationProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = strainEnergyCalculationProcessGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Check user input
userData = get(handles.figure1, 'UserData');

funParams.useFOV=get(handles.checkbox_useFOV,'Value');
funParams.useCellMask=get(handles.checkbox_cellMask,'Value');
funParams.performForceBlobAnalysis=get(handles.checkbox_performForceBlobAnal,'Value');
funParams.exportCSV=get(handles.checkbox_exportCSV,'Value');

% Process Sanity check ( only check underlying data )
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);



% Hint: get(hObject,'Value') returns toggle state of checkbox_cellMask


% --- Executes on button press in checkbox_useFOV.
function checkbox_useFOV_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_useFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_useFOV


% --- Executes on button press in checkbox_cellMask.
function checkbox_cellMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cellMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cellMask


% --- Executes on button press in checkbox_performForceBlobAnal.
function checkbox_performForceBlobAnal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_performForceBlobAnal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_performForceBlobAnal

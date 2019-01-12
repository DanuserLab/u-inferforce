function varargout = displacementFieldCorrectionProcessGUI(varargin)
% displacementFieldCorrectionProcessGUI M-file for displacementFieldCorrectionProcessGUI.fig
%      displacementFieldCorrectionProcessGUI, by itself, creates a new displacementFieldCorrectionProcessGUI or raises the existing
%      singleton*.
%
%      H = displacementFieldCorrectionProcessGUI returns the handle to a new displacementFieldCorrectionProcessGUI or the handle to
%      the existing singleton*.
%
%      displacementFieldCorrectionProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in displacementFieldCorrectionProcessGUI.M with the given input arguments.
%
%      displacementFieldCorrectionProcessGUI('Property','Value',...) creates a new displacementFieldCorrectionProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before displacementFieldCorrectionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to displacementFieldCorrectionProcessGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help displacementFieldCorrectionProcessGUI

% Last Modified by GUIDE v2.5 10-Mar-2017 13:02:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @displacementFieldCorrectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @displacementFieldCorrectionProcessGUI_OutputFcn, ...
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


% --- Executes just before displacementFieldCorrectionProcessGUI is made visible.
function displacementFieldCorrectionProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

if ~isempty(funParams.outlierThreshold),
    set(handles.checkbox_outlierThreshold,'Value',1);
    set(handles.edit_outlierThreshold,'String',funParams.outlierThreshold);
else
    set(handles.checkbox_filterOutliers,'Value',0);
end
set(handles.checkbox_doRotReg,'Value',funParams.doRotReg);
set(handles.checkbox_fill,'Value',funParams.fillVectors);

% Choose default command line output for displacementFieldCorrectionProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = displacementFieldCorrectionProcessGUI_OutputFcn(~, ~, handles) 
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

if get(handles.checkbox_outlierThreshold,'Value')
    value = get(handles.edit_outlierThreshold,'String');
    if isnan(str2double(value))
         errordlg('Please enter a valid value for the outlier threshold.','Setting Error','modal')
        return;
    end
    funParams.outlierThreshold=str2double(value);
else
    funParams.outlierThreshold=[];
end

funParams.doRotReg=get(handles.checkbox_doRotReg,'Value');
funParams.fillVectors=get(handles.checkbox_fill,'Value');

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


% --- Executes on button press in checkbox_outlierThreshold.
function checkbox_outlierThreshold_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.edit_outlierThreshold,'Enable','on');
else
    enableState=set(handles.edit_outlierThreshold,'Enable','off');    
end


% --- Executes on button press in checkbox_fill.
function checkbox_fill_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fill

function varargout = EfficientSubpixelRegistrationProcessGUI(varargin)
% EfficientSubpixelRegistrationProcessGUI M-file for EfficientSubpixelRegistrationProcessGUI.fig
%      EfficientSubpixelRegistrationProcessGUI, by itself, creates a new EfficientSubpixelRegistrationProcessGUI or raises the existing
%      singleton*.
%
%      H = EfficientSubpixelRegistrationProcessGUI returns the handle to a new EfficientSubpixelRegistrationProcessGUI or the handle to
%      the existing singleton*.
%
%      EfficientSubpixelRegistrationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EfficientSubpixelRegistrationProcessGUI.M with the given input arguments.
%
%      EfficientSubpixelRegistrationProcessGUI('Property','Value',...) creates a new EfficientSubpixelRegistrationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EfficientSubpixelRegistrationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EfficientSubpixelRegistrationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
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

% Edit the above text to modify the response to help EfficientSubpixelRegistrationProcessGUI

% Last Modified by GUIDE v2.5 20-Apr-2017 09:30:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EfficientSubpixelRegistrationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @EfficientSubpixelRegistrationProcessGUI_OutputFcn, ...
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


% --- Executes just before EfficientSubpixelRegistrationProcessGUI is made visible.
function EfficientSubpixelRegistrationProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

set(handles.edit_referenceFramePath,'String',funParams.referenceFramePath);

nframes = userData.MD.nFrames_;
frameSelStr = arrayfun(@(x) ['Frame #: ' num2str(x)], 1:nframes, 'UniformOutput',false);
frameSelStr = ['Image Path', frameSelStr];
set(handles.referenceFrame_popupmenu,'String',frameSelStr);

% set(handles.referenceFrame_popupmenu,'Callback',@frameSelStr);

userData.numParams = {'usfac'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),userData.numParams);


% Read the first image and update the sliders max value and steps  
set(handles.listbox_selectedChannels,'Callback',@(h,event) update_data(h,event,guidata(h)));
    
% Override default channels callback function
set(handles.checkbox_all,'Callback',@(hObject,eventdata)...
    checkallChannels(hObject,eventdata,guidata(hObject)));
set(handles.pushbutton_select,'Callback',@(hObject,eventdata)...
    selectChannel(hObject,eventdata,guidata(hObject)));
set(handles.pushbutton_delete,'Callback',@(hObject,eventdata)...
    deleteChannel(hObject,eventdata,guidata(hObject)));

% Choose default command line output for EfficientSubpixelRegistrationProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = EfficientSubpixelRegistrationProcessGUI_OutputFcn(~, ~, handles) 
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

if isfield(userData, 'previewFig') && ishandle(userData.previewFig) 
    delete(userData.previewFig) 
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
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
else
    channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
    funParams.ChannelIndex = channelIndex;
end

if handles.referenceFrame_popupmenu.Value == 1 %,== 'Select Image Path')
    % Retrieve reference frame path
    funParams.referenceFramePath=get(handles.edit_referenceFramePath,'String');    
    funParams.referenceFrameNum = 0;
else
    funParams.referenceFramePath='';
    funParams.referenceFrameNum = handles.referenceFrame_popupmenu.Value + 1; % Assumes first one is path selection
    
end

if handles.referenceFrame_popupmenu.Value == 1 && isempty(funParams.referenceFramePath)
    errordlg('No reference frame selected, please select path or choose frame #.')
    return;
end

% Read numeric information
for i = 1:numel(userData.numParams)
    value = get(handles.(['edit_' userData.numParams{i}]),'String');
    if isempty(value)
        errordlg('Please enter a valid value.','Setting Error','modal')
        return;
    end
    funParams.(userData.numParams{i})=str2double(value); 
end


if logical(mod(funParams.usfac,1)) || (funParams.usfac <= 0)
    errordlg('Please enter a valid value.','Setting Error','modal')
    return;
end

% Process Sanity check ( only check underlying data )
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles, funParams);

% --- Executes on button press in pushbutton_selectReferenceFrame.
function pushbutton_selectReferenceFrame_Callback(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
[file path]=uigetfile({'*.tif;*.tiff;*.ome.tiff;*.TIF;*.stk;*.STK;*.bmp;*.BMP;*.jpg;*.JPG',...
    'Image files (*.tif,*.tiff,*.stk,*.bmp,*.jpg,*.ome.tiff)'},...
    'Select the reference frame',userData.MD.outputDirectory_);
if ~isequal(file,0) && ~isequal(path,0)
    set(handles.edit_referenceFramePath,'String',[path file]);
end

 % --- Executes on button press in checkbox_crop.
function update_data(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

% Retrieve the channel index
props=get(handles.listbox_selectedChannels,{'UserData','Value'});
if isempty(props{1}), return; end
chanIndx = props{1}(props{2});

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


function selectChannel(hObject, eventdata, handles)

selectChannel_Callback(hObject, eventdata, handles);
update_data(hObject, eventdata, handles);

function deleteChannel(hObject, eventdata, handles)

deleteChannel_Callback(hObject, eventdata, handles);
update_data(hObject, eventdata, handles);

function checkallChannels(hObject, eventdata, handles)

checkallChannels_Callback(hObject, eventdata, handles);
update_data(hObject, eventdata, handles);


% --------------------------------------------------------------------
function uipanel_1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in referenceFrame_popupmenu.
function referenceFrame_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to referenceFrame_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns referenceFrame_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from referenceFrame_popupmenu

if hObject.Value == 1
    set(handles.edit_referenceFramePath, 'Enable', 'on');
    set(handles.pushbutton_selectReferenceFrame, 'Enable', 'on');
elseif hObject.Value >= 2
    set(handles.edit_referenceFramePath, 'Enable', 'off');    
    set(handles.pushbutton_selectReferenceFrame, 'Enable', 'off'); 
%     set(handles.text39, 'ForegroundColor', [0 0 0]); 
end


% --- Executes during object creation, after setting all properties.
function referenceFrame_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to referenceFrame_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

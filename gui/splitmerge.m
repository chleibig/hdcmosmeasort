function varargout = splitmerge(varargin)
% SPLITMERGE M-file for splitmerge.fig
%      SPLITMERGE, by itself, creates a new SPLITMERGE or raises the existing
%      singleton*.
%
%      H = SPLITMERGE returns the handle to a new SPLITMERGE or the handle to
%      the existing singleton*.
%
%      SPLITMERGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPLITMERGE.M with the given input arguments.
%
%      SPLITMERGE('Property','Value',...) creates a new SPLITMERGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before splitmerge_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to splitmerge_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help splitmerge

% Last Modified by GUIDE v2.5 07-Oct-2013 15:14:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @splitmerge_OpeningFcn, ...
                   'gui_OutputFcn',  @splitmerge_OutputFcn, ...
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

% --- Executes just before splitmerge is made visible.
function splitmerge_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to splitmerge (see VARARGIN)

% Choose default command line output for splitmerge
handles.output = hObject;

%Load preliminary sorting results into guidata.
dataimport(hObject, handles);
%update handles structure that now contains the output of the automatic
%sorting stage.
handles = guidata(hObject);

%Set initially all unit states to unchecked.
%1 - unchecked, 2 - threshold ok, 3 - to be saved, % 4 - to be deleted
if ~isfield([handles.units],'state')
    [handles.units.state] = deal(1);
end
% {'unchecked', 'threshold ok', 'to be saved', 'to be deleted'} <-> 
handles.stateColor = {[0 0 1],[0 1 1],[0 1 0],[1 0 0]}; 


%Set labels in listbox1 accordingly.
set(handles.listbox1, 'String', ...
    cellfun(@num2str,num2cell(1:size(handles.S,1)),'UniformOutput',false));

%Initial plots.
drawrois(hObject,handles);

% UIWAIT makes splitmerge wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%--------------------------------------------------------------------------


% --- Outputs from this function are returned to the command line.---------
function varargout = splitmerge_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%--------------------------------------------------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACK FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% -------------------------------------------------------------------------


function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end
% -------------------------------------------------------------------------



function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)
% -------------------------------------------------------------------------



function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.figure1)
% -------------------------------------------------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choice of IC component 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in listbox1-----------------------------
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

popup_sel_index = get(handles.listbox1, 'Value');

selectcomponent(hObject, eventdata, handles, popup_sel_index);

%--------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude adjustment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on slider movement and on IC selection
function thresholdslider_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

popup_sel_index = get(handles.listbox1,'Value');

%set threshold.
handles.threshold = get(hObject, 'Value');

% update unit
updateunitdata(hObject, handles, popup_sel_index, handles.threshold);
% update handles structure
handles = guidata(hObject);

% update plots
updateunitplots(hObject, handles, popup_sel_index, handles.threshold);
% -------------------------------------------------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit/IC state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in listbox2.----------------------------
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2

popup_sel_index = get(handles.listbox1,'Value');
handles.units(popup_sel_index).state = get(hObject,'Value');
set(handles.sPlot,...
    'Color',handles.stateColor{handles.units(popup_sel_index).state});
N_units = length(handles.units);
fractionUnchecked = 100*nnz([handles.units.state] == 1)/N_units;
fractionThresholdOk = 100*(nnz([handles.units.state] == 2) + ... 
                           nnz([handles.units.state] == 3) + ...
                           nnz([handles.units.state] == 4))/N_units;
fractionToBeSaved = 100*nnz([handles.units.state] == 3)/N_units;
fractionToBeDeleted = 100*nnz([handles.units.state] == 4)/N_units;
set(handles.textFractionUnchecked,'String',num2str(fractionUnchecked,'%3.1f'));
set(handles.textFractionThresholdOk,'String',num2str(fractionThresholdOk,'%3.1f'));
set(handles.textFractionToBeSaved,'String',num2str(fractionToBeSaved,'%3.1f'));
set(handles.textFractionToBeDeleted,'String',num2str(fractionToBeDeleted,'%3.1f'));

guidata(hObject,handles);



% -------------------------------------------------------------------------

% --- Executes on button press in pushbuttonGetduplicate.------------------
function pushbuttonGetduplicate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonGetduplicate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%collect information of currently selected unit.

% currentUnitIndex = get(handles.listbox1, 'Value');
% handles.marked(currentUnitIndex) = true;
% guidata(hObject,handles);

currentUnit = get(handles.listbox1,'Value');
distMax = str2double(get(handles.editDistMax,'String'));
if isnan(distMax)
    errordlg('Maximum distance has to be a numeric value','Bad Input','modal');
    return
end

duplicatecandidates(currentUnit,...
    handles.units, handles.ROIs,handles.params,...
    distMax, handles.data);
% -------------------------------------------------------------------------



% --- Executes on button press in pushbuttonSTA.---------------------------
function pushbuttonSTA_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSTA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

popup_sel_index = get(handles.listbox1, 'Value');
GetSTA(handles.data,handles.units(popup_sel_index).time,handles.params.sr,1);
% -------------------------------------------------------------------------



function editDistMax_Callback(hObject, eventdata, handles)
% hObject    handle to editDistMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDistMax as text
%        str2double(get(hObject,'String')) returns contents of editDistMax as a double
% -------------------------------------------------------------------------


% --- Executes on button press in pushbuttonRedrawSpatial.
function pushbuttonRedrawSpatial_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRedrawSpatial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
drawrois(hObject,handles);
% -------------------------------------------------------------------------



% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if nnz([handles.units.state] == 3) + nnz([handles.units.state] == 4) ~= ...
        length(handles.units)
    msgbox(['You are not ready with manual intervention, saving '...
        'intermediate results only!']);
end
dataexport(hObject, eventdata, handles);
delete(handles.figure1);
% -------------------------------------------------------------------------



% --- Executes on button press in pushbuttonShowRawData.
function pushbuttonShowRawData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShowRawData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popup_sel_index = get(handles.listbox1, 'Value');

sensorRow = handles.units(popup_sel_index).boss_row;
sensorCol = handles.units(popup_sel_index).boss_col;
dataRow = find(handles.params.sensor_rows == sensorRow);
dataCol = find(handles.params.sensor_cols == sensorCol);

N_ROW = length(handles.params.sensor_rows);
N_COL = length(handles.params.sensor_cols);

figure;suptitle(['Raw data region around unit position centered on ( ',...
    num2str(sensorRow),', ', num2str(sensorCol), ')']);
for r = 1:5
    for c = 1:5
        subplot(5,5,sub2ind([5 5],c,r));
        if ((dataRow-3+r)>0) && ((dataCol-3+c)>0) && ...
                ((dataRow-3+r)<=N_ROW) && ((dataCol-3+c)<=N_COL)
            plot(squeeze(handles.data(dataRow-3+r,dataCol-3+c,:)));
        end
    end
end

%rawDataTrace = resample(handles.data(dataRow,dataCol,:),...
%    handles.params.upsample,1);
% figure;plot(rawDataTrace);
% title(['Row: ' num2str(sensorRow) ' Column: ' num2str(sensorCol)]);
% ylabel('mV');
% xlabel('samples');

% -------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%Multiple selection could be enabled via Max-Min > 1
set(hObject,'Max',1);
set(hObject,'Min',0);
% -------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Set labels for unit state accordingly.
set(hObject, 'String', {'unchecked', 'threshold ok', 'to be saved', 'to be deleted'});
% -------------------------------------------------------------------------



% --- Executes during object creation, after setting all properties.
function thresholdslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min',-100, 'Max', 0);
% -------------------------------------------------------------------------



% --- Executes during object creation, after setting all properties.
function editDistMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDistMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------



% --- Executes during object creation, after setting all properties.
function textFractionUnchecked_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textFractionUnchecked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textFractionToBeSaved_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textFractionToBeSaved (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function textFractionToBeDeleted_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textFractionToBeDeleted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function textFractionThresholdOk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textFractionThresholdOk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRAPPERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
function dataimport(hObject, handles)
%Currently get data from workspace and push it to handles structure;
%Other solutions might be thought of...
ROIs = evalin('base','ROIs');
handles.ROIs = ROIs;
handles.ROIsAsCC = evalin('base','ROIsAsCC');
handles.params = evalin('base','params');
handles.data = evalin('base','data');
%all sources.
S = cat(1,ROIs.S);
S(skewness(S') > 0,:) = -1*S(skewness(S') > 0,:);
handles.S = S;
clear S
%all units.
units = [];
for i = 1:length(ROIs);
    if ~isempty(ROIs(i).units);
        [ROIs(i).units.k] = deal(i);
        units = [units ROIs(i).units];
    end
end
handles.units = units;
clear units
%Update guidata
guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function dataexport(hObject, eventdata, handles)
%Currently push data changes back to workspace;
%Other solutions might be thought of...

%The manual intervention either changes units or removes them, these
%changes we have to store accordingly.
%go through handels.ROIs, overwrite previous data with changes and push
%the results to the workspace.

dbstop if error

for i = 1:length(handles.ROIs);
    %to be deleted <-> state '4'
    toDelete = [handles.units.k] == i & [handles.units.state] == 4;
    handles.ROIs(i).units_del = handles.units(toDelete);
    handles.ROIs(i).S_del = handles.S(toDelete,:);
    %map toDelete to index range in ROI.
    toDelete = toDelete([handles.units.k] == i);
    if ~isfield(handles.ROIs(i),'A_del'); handles.ROIs(i).A_del = [];end
    %handles.ROIs(i).A_del = cat(2,handles.ROIs(i).A_del,...
     %                             handles.ROIs(i).A_tau(:,toDelete,:));    
    clear toDelete
    %to be saved <-> state '3'
    toSave = [handles.units.k] == i & [handles.units.state] <= 3;
    handles.ROIs(i).units = handles.units(toSave);
    handles.ROIs(i).S = handles.S(toSave,:);
    %map toSave to index range in ROI.
    toSave = toSave([handles.units.k] == i);
    %handles.ROIs(i).A_tau = handles.ROIs(i).A_tau(:,toSave,:);
    clear toSave                                                          
end

assignin('base','ROIs',handles.ROIs);
assignin('base','params',handles.params);

%Update guidata
guidata(hObject, handles);

%Save results to basic event file.
SaveResults(handles.ROIs, handles.params);

%--------------------------------------------------------------------------
function selectcomponent(hObject, eventdata, handles, whichone)

unit = handles.units(whichone);

%Get data to plot.
sr = handles.params.sr;
upsample = handles.params.upsample;
S = resample(handles.S(whichone,:),upsample,1);
N_SAMPLES = length(S);
% indices = round(handles.units(whichone).time * ...
%     sr * upsample);
noise_std = median(abs(S)/0.6745);
%thrFactor = handles.params.thrFactor;
thrFactor = handles.params.thrFactor;
handles.threshold = -thrFactor*noise_std;
[indices,pos_amplitudes] = find_peaks(-S,thrFactor*noise_std,ceil(sr*upsample));
amplitudes = -1*pos_amplitudes;

if ~isempty(amplitudes)
    %Get waveforms.
    pre = round(0.5 * sr*upsample);
    post = round(0.5 * sr*upsample);
    f_tot = pre + post + 1;
    %take only those peaks for which the desired window is entirely contained
    %in the data:
    amplitudes = amplitudes( ((indices - pre) >= 1) & ((indices + post) <= N_SAMPLES) );
    indices = indices( ((indices - pre) >= 1) & ((indices + post) <= N_SAMPLES) );
    N_pks = length(indices);
    pks = zeros(f_tot,N_pks);
    for j = 1:N_pks
        pks(:,j) = S(indices(j)-pre:indices(j)+post);
    end
    [cts, bin_ctrs] = hist(amplitudes,floor(sqrt(length(amplitudes))));
end

%Plotting

%entire IC.
axes(handles.axes1);
cla reset;
handles.sPlot = plot(S,'Color', handles.stateColor{unit.state});
hold on;
title(strcat('amplSD = ',num2str(unit.amplitudeSD),'; RSTD = ',...
             num2str(unit.RSTD),'; sep. = ',num2str(unit.separability)));

         
hold on;
%all threshold crossings.
plot(indices, S(indices),'ro','LineStyle','none');
handles.axes1thr = plot([0 N_SAMPLES], [handles.threshold handles.threshold],'r');
%those from units.
% handles.axes1amp = plot(round(unit.time*upsample*sr),...
%                       unit.amplitude,'go','LineStyle','none');

if ~isempty(amplitudes)
    %IC waveforms.
    axes(handles.axes2);
    cla reset;
    plot(pks,'Color',[0.75 0.75 0.75]);
    hold on;
    handles.axes2thr = plot([0 f_tot], [handles.threshold handles.threshold],'r');

    %amplitude histogram.
    axes(handles.axes3);
    cla reset;
    bar(bin_ctrs,cts); 
    hold on;
    handles.axes3thr = plot([handles.threshold handles.threshold], ylim,'r');
end


%Spatial visualization.
axes(handles.axes4);hold on;
if isfield(handles,'unitMarker') && ishandle(handles.unitMarker)
    delete(handles.unitMarker);
end
handles.unitMarker = plot(unit.boss_col,unit.boss_row,...
    'o','MarkerSize',25,'Color','black','LineWidth',3);



%To be used by slider:
if ~isempty(amplitudes)
    handles.amplitudes = amplitudes;
    handles.indices = indices;
    handles.pks = pks;
    handles.time = handles.indices/(upsample*sr);
else
    handles.amplitudes = [];
    handles.indices = [];
    handles.pks = [];
    handles.time = [];
end


handles.threshold = max(unit.amplitude);

%Adjust slider.
set(handles.thresholdslider, 'Min',min(amplitudes),...
                             'Value',handles.threshold);
                         
%Update unit state.
set(handles.listbox2,'Value',unit.state);

%update guidata.
guidata(hObject,handles);

thresholdslider_Callback(handles.thresholdslider, eventdata, handles);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function updateunitdata(hObject, handles, whichone, newThreshold)
% updateunit(hObject, handles, whichone, newThreshold)

%valid threshold crossings.
valid = handles.amplitudes <= newThreshold;

handles.units(whichone).time = handles.time(valid);
handles.units(whichone).amplitude = handles.amplitudes(valid);
handles.units(whichone).amplitudeSD = std(handles.units(whichone).amplitude);
handles.units(whichone).RSTD = ...
                            handles.units(whichone).amplitudeSD/...
                       mean(abs(handles.units(whichone).amplitude));
if nnz(~valid) > 0
handles.units(whichone).separability = ...
    mean(abs(handles.units(whichone).amplitude)) - ...
    mean(abs(handles.amplitudes(~valid)));
else
    handles.units(whichone).separability = ...
    mean(abs(handles.units(whichone).amplitude));
end
    
handles.units(whichone).SDscore = [];

%which ROI <-> k.
k = handles.units(whichone).k;
% STA.
dataTmp = reshape(handles.ROIs(k).X,...
    [length(handles.ROIs(k).sensor_rows)...
    length(handles.ROIs(k).sensor_cols) size(handles.ROIs(k).X,2)]);
handles.units(whichone).STA = GetSTA(dataTmp,handles.time(valid),...
                            handles.params.sr,0);
clear dataTmp
% unit position.
[row_max,col_max] = find(max(abs(handles.units(whichone).STA),[],3)...
    == max(max(max(abs(handles.units(whichone).STA)))));
handles.units(whichone).boss_row = handles.ROIs(k).sensor_rows(row_max);
handles.units(whichone).boss_col = handles.ROIs(k).sensor_cols(col_max);
clear k
% update guidata
guidata(hObject,handles);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function updateunitplots(hObject, handles, whichone, newThreshold)
%updateunitplots(hObject, handles, whichone, newThreshold)

%valid threshold crossings.
valid = handles.amplitudes <= newThreshold;

axes(handles.axes1);
if isfield(handles,'axes1amp');
    delete(handles.axes1amp(ishandle(handles.axes1amp)));
end
handles.axes1amp = plot(handles.indices(valid),...
                      handles.amplitudes(valid),'go','LineStyle','none');
if isfield(handles,'axes1thr');
    delete(handles.axes1thr(ishandle(handles.axes1thr)));
end
handles.axes1thr = plot(xlim, [newThreshold newThreshold],'r');

axes(handles.axes2);
if isfield(handles,'axes2pks');
    delete(handles.axes2pks(ishandle(handles.axes2pks)));
end
handles.axes2pks = plot(handles.pks(:,valid),'Color','green');
if isfield(handles,'axes2thr');
    delete(handles.axes2thr(ishandle(handles.axes2thr)));
end
handles.axes2thr = plot(xlim, [newThreshold newThreshold],'r');

axes(handles.axes3);
if isfield(handles,'axes3thr');
    delete(handles.axes3thr(ishandle(handles.axes3thr)));
end
handles.axes3thr = plot([newThreshold newThreshold], ylim,'r');

axes(handles.axes4);
if isfield(handles,'unitLabels')
    delete(handles.unitLabels(ishandle(handles.unitLabels)));
end
cols = [handles.units.boss_col];
rows = [handles.units.boss_row];
notDeleted = [handles.units.state] <= 3;
cols = cols(notDeleted);
rows = rows(notDeleted);
idx = 1:length(handles.units);
idx = idx(notDeleted);
handles.unitLabels = text(cols, rows,...
         arrayfun(@num2str,idx,'UniformOutput',false));

axes(handles.axes4);hold on;
if isfield(handles,'unitMarker') && ishandle(handles.unitMarker)
    delete(handles.unitMarker);
end
handles.unitMarker = plot(handles.units(whichone).boss_col,...
                          handles.units(whichone).boss_row,...
                        'o','MarkerSize',25,'Color','black','LineWidth',3);

% update guidata
guidata(hObject,handles);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function drawrois(hObject,handles)

axes(handles.axes4);
showrois(handles.ROIsAsCC,...
    'fillColumnFlag',logical(handles.params.d_col/handles.params.pitch-1));
cols = [handles.units.boss_col];
rows = [handles.units.boss_row];
notDeleted = [handles.units.state] <= 3;
cols = cols(notDeleted);
rows = rows(notDeleted);
idx = 1:length(handles.units);
idx = idx(notDeleted);
handles.unitLabels = text(cols, rows,...
         arrayfun(@num2str,idx,'UniformOutput',false));


% Update handles structure
guidata(hObject, handles);
%--------------------------------------------------------------------------






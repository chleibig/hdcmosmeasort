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

% Last Modified by GUIDE v2.5 15-Apr-2014 16:28:08

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

% Background colors
bgcolor = [1 1 1];
set(hObject, 'Color', bgcolor);
set(findobj(hObject,'-property', 'BackgroundColor'), 'BackgroundColor', bgcolor);


% Choose default command line output for splitmerge
handles.output = hObject;

%Load preliminary sorting results into guidata.
dataimport(hObject, handles);
%update handles structure that now contains the output of the automatic
%sorting stage.
handles = guidata(hObject);

%Set initially all unit states to unchecked.
%1 - unchecked, 2 - mixture, 3 - single, % 4 - delete
if ~isfield([handles.units],'state')
    [handles.units.state] = deal(1);
end


%Set labels in listbox1 accordingly.
set(handles.listbox1, 'String', ...
    cellfun(@num2str,num2cell(1:size(handles.S,1)),'UniformOutput',false));


%Set values of edit text fields
set(handles.editDistMax,'String',num2str(handles.params.d_max));
set(handles.editCoinThr,'String',num2str(handles.params.coin_thr));
set(handles.editSimThr,'String',num2str(handles.params.sim_thr));

%Initial plots.
drawrois(hObject,handles);
handles.unitIDsAsStrings = ...
    arrayfun(@num2str,1:length(handles.units),'UniformOutput',false);
handles.unitIDsAsStringsSorted = handles.unitIDsAsStrings;
guidata(hObject,handles);
drawunitlabels(hObject,handles);


%set background colors by state in unit list
bgColors = handles.stateColor([handles.units.state]);
unitIDsAsColoredStrings = arrayfun(@(x) ...
                    setbgcolor(bgColors{x},handles.unitIDsAsStrings{x}),...
                                 1:length(bgColors),'UniformOutput',false);
clear bgColors
set(handles.listbox1,'String',unitIDsAsColoredStrings);

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
% Choice of Source component 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in listbox1-----------------------------
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

currentUnit = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});

selectcomponent(hObject, eventdata, handles, currentUnit);

%--------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude adjustment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on slider movement and on Source selection
function thresholdslider_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

currentUnit = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});

%set threshold.
handles.threshold = get(hObject, 'Value');

% update unit
updateunitdata(hObject, handles, currentUnit,...
                                        handles.threshold);
% update handles structure
handles = guidata(hObject);

% update plots
updateunitplots(hObject, handles, currentUnit, handles.threshold);
% -------------------------------------------------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit/Source state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in listboxUnitState.----------------------------
function listboxUnitState_Callback(hObject, eventdata, handles)
% hObject    handle to listboxUnitState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listboxUnitState contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxUnitState

newState = get(hObject,'Value');

unitListIdx = get(handles.listbox1,'Value');
currentUnit = str2double(handles.unitIDsAsStrings{unitListIdx});


currentUnitNewBgColorAsString = setbgcolor(handles.stateColor{newState},...
                                           num2str(currentUnit));
tmp = get(handles.listbox1,'String');
tmp{unitListIdx} = currentUnitNewBgColorAsString;
set(handles.listbox1,'String',tmp);
clear tmp

handles.units(currentUnit).state = newState;
guidata(hObject,handles);

updateunitstatefractions(hObject, handles, currentUnit);

%draw spatial positions of units
axes(handles.axes4);
if isfield(handles,'unitLabels')
    delete(handles.unitLabels(ishandle(handles.unitLabels)));
end
drawunitlabels(hObject, handles);
% -------------------------------------------------------------------------



% --- Executes on selection change in listboxSortCriteria.-----------------
function listboxSortCriteria_Callback(hObject, eventdata, handles)
% hObject    handle to listboxSortCriteria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listboxSortCriteria contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        listboxSortCriteria

% The sorting is encoded in handles.unitIDsAsStrings -> whenever a unit
% is selected via the unit listbox make sure to use the string value as
% an identifier and not the position in the list.

%reset index to default order first
handles.unitIDsAsStrings = handles.unitIDsAsStringsSorted;
%color backgrounds according to unit states
bgColors = handles.stateColor([handles.units.state]);
unitIDsAsColoredStrings = arrayfun(@(x) ...
                    setbgcolor(bgColors{x},handles.unitIDsAsStrings{x}),...
                                 1:length(bgColors),'UniformOutput',false);
clear bgColors

popup_sel_index = get(handles.listboxSortCriteria,'Value');
switch handles.sortCriteria{popup_sel_index}
    case 'index'
         
    case 'separability'
        [unused,idx] = sort([handles.units.separability],'descend');
        handles.unitIDsAsStrings = handles.unitIDsAsStrings(idx);
        unitIDsAsColoredStrings = unitIDsAsColoredStrings(idx);
    
    case 'IsoIBg'
        [unused,idx] = sort([handles.units.IsoIBg],'descend');
        handles.unitIDsAsStrings = handles.unitIDsAsStrings(idx);
        unitIDsAsColoredStrings = unitIDsAsColoredStrings(idx);
        
    case 'IsoINN'
        [unused,idx] = sort([handles.units.IsoINN],'descend');
        handles.unitIDsAsStrings = handles.unitIDsAsStrings(idx);
        unitIDsAsColoredStrings = unitIDsAsColoredStrings(idx);

    case 'RSTD'
        [unused,idx] = sort([handles.units.RSTD],'ascend');
        handles.unitIDsAsStrings = handles.unitIDsAsStrings(idx);
        unitIDsAsColoredStrings = unitIDsAsColoredStrings(idx);
         
    case 'skewness'
        [unused,idx] = sort(abs([handles.skewn]),'descend');
        handles.unitIDsAsStrings = handles.unitIDsAsStrings(idx);
        unitIDsAsColoredStrings = unitIDsAsColoredStrings(idx);
         
    case 'kurtosis'
        [unused,idx] = sort([handles.kurtosis],'descend');
        handles.unitIDsAsStrings = handles.unitIDsAsStrings(idx);
        unitIDsAsColoredStrings = unitIDsAsColoredStrings(idx);
        
    case 'SNR'
        try
            [unused,idx] = sort([handles.units.snr],'descend');
            handles.unitIDsAsStrings = handles.unitIDsAsStrings(idx);
            unitIDsAsColoredStrings = unitIDsAsColoredStrings(idx);
        catch
            msgbox('SNR is not available. Sort by index instead.')
        end
        
    otherwise
end

set(handles.listbox1,'String',unitIDsAsColoredStrings);
clear unitIDsAsColoredStrings

guidata(hObject, handles);

% -------------------------------------------------------------------------



% --- Executes on button press in pushbuttonAcceptThresholds. -------------
function pushbuttonAcceptThresholds_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAcceptThresholds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.units([handles.units.state] == 1).state] = deal(2);
guidata(hObject,handles);

%Update visualization
currentUnit = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});
updateunitstatefractions(hObject, handles, currentUnit);
%update background colors in unit list
bgColors = handles.stateColor([handles.units.state]);
unitIDsAsColoredStrings = arrayfun(@(x) ...
                    setbgcolor(bgColors{x},handles.unitIDsAsStrings{x}),...
                                 1:length(bgColors),'UniformOutput',false);
clear bgColors
set(handles.listbox1,'String',unitIDsAsColoredStrings);


% -------------------------------------------------------------------------



% --- Executes on button press in pushbuttonShowLocalRedundancy.-----------
function pushbuttonShowLocalRedundancy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShowLocalRedundancy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

maxDist = str2double(get(handles.editDistMax,'String'));
minCoin = str2double(get(handles.editCoinThr,'String'));
minSim = str2double(get(handles.editSimThr,'String'));
currentUnit = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});

[numDistinct, sizes, members] = ...
                       redundantcandidates(handles.units, handles.ROIs, ...
                                handles.params, handles.data,...
                                maxDist, minCoin, minSim, currentUnit);
fprintf(['Found %g putative different underlying unit(s)\n'],numDistinct);
figure;
for i = 1:numDistinct
   subplot(numDistinct,1,i);plot(handles.S(members{i},:)');
   legend(...
       arrayfun(@(x) ...
       [num2str(x) ' sep. = ' num2str(handles.units(x).separability,2)], ...
       members{i},'UniformOutput',false)...
       )
end
suplabel('samples','x');
suplabel('source activation','y');

% -------------------------------------------------------------------------



% --- Executes on button press in pushbuttonRemoveAllRedundancy.
function pushbuttonRemoveAllRedundancy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemoveAllRedundancy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

maxDist = str2double(get(handles.editDistMax,'String'));
minCoin = str2double(get(handles.editCoinThr,'String'));
minSim = str2double(get(handles.editSimThr,'String'));

[numDistinct, sizes, members] = redundantcandidates(handles.units,...
                             handles.ROIs, handles.params, handles.data,...
                             maxDist, minCoin, minSim);

N_deleted = 0;                           
for i = 1:numDistinct
    if sizes(i) > 1
        %keep only the most separable one.
        markAsDelete = ([handles.units(members{i}).separability] ~= ...
                            max([handles.units(members{i}).separability]));
        fprintf('Changing state of unit(s) %g to delete.\n',...
            members{i}(markAsDelete));
        [handles.units(members{i}(markAsDelete)).state] = deal(4);
        N_deleted = N_deleted + nnz(markAsDelete);
    end
end

%update params struct if threshold values were more aggressive in the sense
%of enabling a larger fraction of units to be deleted 
if maxDist > handles.params.d_max
    handles.params.d_max = maxDist;
end
if minCoin < handles.params.coin_thr;
    handles.params.coin_thr = minCoin;
end
if minSim < handles.params.sim_thr;
    handles.params.sim_thr = minSim;
end

guidata(hObject,handles);

currentUnit = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});
updateunitstatefractions(hObject, handles, currentUnit);
%update background colors in unit list
listboxSortCriteria_Callback(hObject, eventdata, handles);

%redraw spatial positions of units to visualize removed redundancy
drawunitlabels(hObject, handles);
fprintf('Changed %g unit states to delete.\n Remove all done.\n',N_deleted);

% -------------------------------------------------------------------------




% --- Executes on button press in pushbuttonSTA.---------------------------
function pushbuttonSTA_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSTA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentUnit = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});
computetemplate(handles.data,...
                handles.units(currentUnit).time,handles.params.sr,1);
% -------------------------------------------------------------------------



function editDistMax_Callback(hObject, eventdata, handles)
% hObject    handle to editDistMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDistMax as text
%        str2double(get(hObject,'String')) returns contents of editDistMax as a double

maxDist = str2double(get(hObject,'String'));
if isnan(maxDist)
    errordlg('Entered value must be numeric!','Bad Input','modal');
    return
end
guidata(hObject, handles);

% -------------------------------------------------------------------------


function editCoinThr_Callback(hObject, eventdata, handles)
% hObject    handle to editCoinThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCoinThr as text
%        str2double(get(hObject,'String')) returns contents of editCoinThr as a double

minCoin = str2double(get(hObject,'String'));
if isnan(minCoin)
    errordlg('Entered value must be numeric!','Bad Input','modal');
    return
end

guidata(hObject, handles);
% -------------------------------------------------------------------------



function editSimThr_Callback(hObject, eventdata, handles)
% hObject    handle to editSimThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSimThr as text
%        str2double(get(hObject,'String')) returns contents of editSimThr
%        as a double

minSim = str2double(get(hObject,'String'));
if isnan(minSim)
    errordlg('Entered value must be numeric!','Bad Input','modal');
    return
end

guidata(hObject, handles);
% -------------------------------------------------------------------------



function editMin_Callback(hObject, eventdata, handles)
% hObject    handle to editMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMin as text
%        str2double(get(hObject,'String')) returns contents of editMin as a double

featMin = str2double(get(hObject,'String'));
if isnan(featMin)
    errordlg('Entered value must be numeric!','Bad Input','modal');
    return
end
handles.featMin = featMin;
guidata(hObject, handles);
%--------------------------------------------------------------------------



function editMax_Callback(hObject, eventdata, handles)
% hObject    handle to editMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMax as text
%        str2double(get(hObject,'String')) returns contents of editMax as a
%        double

featMax = str2double(get(hObject,'String'));
if isnan(featMax)
    errordlg('Entered value must be numeric!','Bad Input','modal');
    return
end
handles.featMax = featMax;
guidata(hObject, handles);
%--------------------------------------------------------------------------



% --- Executes on button press in pushbuttonRedrawSpatial.
function pushbuttonRedrawSpatial_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRedrawSpatial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
drawrois(hObject,handles);
drawunitlabels(hObject, handles);
hold(handles.axes4,'on');
whichone = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});
%Only generate handle for the first time
handles.unitMarker = plot(handles.axes4,...
    handles.units(whichone).boss_col,...
    handles.units(whichone).boss_row,...
    'o','MarkerSize',25,'Color','black','LineWidth',3);
guidata(hObject,handles);
% -------------------------------------------------------------------------



% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if nnz([handles.units.state] == 3) + nnz([handles.units.state] == 4) ~= ...
%         length(handles.units)
%     msgbox(['You are not ready with manual intervention, saving '...
%         'intermediate results only!']);
% end
choice = questdlg(['Do you want to save the current state?'...
                  '(overwrite ROIs and params in workspace'...
                  ' and ' handles.params.filenameResults ' on disk)'],...
                  '','Continue','Cancel','Cancel');
if strcmp(choice,'Cancel')
   return
end
dataexport(hObject, eventdata, handles);
delete(handles.figure1);
% -------------------------------------------------------------------------



% --- Executes on button press in pushbuttonShowRawData.
function pushbuttonShowRawData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShowRawData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentUnit = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});

sensorRow = handles.units(currentUnit).boss_row;
sensorCol = handles.units(currentUnit).boss_col;
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



% --- Executes on button press in pushbuttonFeatureHist.
function pushbuttonFeatureHist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFeatureHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selected = get(handles.listboxSortCriteria, 'Value');
feature = handles.sortCriteria{selected};

N_UNITS = length(handles.units);

switch feature
    case 'index'
         msgbox('This would result in a flat distribution;)')
    case 'separability'
        figure;hist([handles.units.separability],sqrt(N_UNITS));
        xlabel(feature);ylabel('counts');
        set(findall(gcf,'type','text'),'fontSize',12);
        
    case 'IsoIBg'
        figure;hist([handles.units.IsoIBg],sqrt(N_UNITS));
        xlabel(feature);ylabel('counts');
        set(findall(gcf,'type','text'),'fontSize',12);
        
    case 'IsoINN'
        figure;hist([handles.units.IsoINN],sqrt(N_UNITS));
        xlabel(feature);ylabel('counts');
        set(findall(gcf,'type','text'),'fontSize',12);
        
    case 'RSTD'
        figure;hist([handles.units.RSTD],sqrt(N_UNITS));
        xlabel(feature);ylabel('counts');
        set(findall(gcf,'type','text'),'fontSize',12);
         
    case 'skewness'
        figure;hist([handles.skewn],sqrt(length(handles.skewn)));
        xlabel(feature);ylabel('counts');
        set(findall(gcf,'type','text'),'fontSize',12);
         
    case 'kurtosis'
        figure;hist([handles.kurtosis],sqrt(length(handles.kurtosis)));
        xlabel(feature);ylabel('counts');
        set(findall(gcf,'type','text'),'fontSize',12);
        
    case 'SNR'
        msgbox('SNR is not a selective feature!')

    otherwise
end
%--------------------------------------------------------------------------



% --- Executes on button press in pushbuttonCropUnits.
function pushbuttonCropUnits_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCropUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


selected = get(handles.listboxSortCriteria, 'Value');
feature = handles.sortCriteria{selected};
featMin = handles.featMin;
featMax = handles.featMax;

choice = questdlg(['Do you want to delete units with ' feature ...
                  ' < ' num2str(featMin) ' or > ' num2str(featMax) '?'],...
                  '','Continue','Cancel','Cancel');

if strcmp(choice,'Cancel')
   return
end

toDelete = false(size([handles.units.state]));
%update params struct if threshold values were more aggressive in the sense
%of enabling a larger fraction of units to be deleted
switch feature
    case 'index'
        msgbox('This is not a good idea;)')
    case 'separability'
        if ~isinf(featMax)
            featMax = Inf;
            msgbox(['Units with separability exceeding max. value will '...
                    'not be deleted!']);
        end
        toDelete = ( [handles.units.separability] < featMin ) | ...
                   ( [handles.units.separability] > featMax );
        
        if featMin > handles.params.minSeparability
            handles.params.minSeparability = featMin;
        end
        
    case 'IsoIBg'
        if ~isinf(featMax)
            featMax = Inf;
            msgbox(['Units with IsoIBg exceeding max. value will '...
                'not be deleted!']);
        end
        toDelete = ( [handles.units.IsoIBg] < featMin ) | ...
            ( [handles.units.IsoIBg] > featMax );

        if featMin > handles.params.minIsoIBg
            handles.params.minIsoIBg = featMin;
        end
        
    case 'IsoINN'
        if ~isinf(featMax)
            featMax = Inf;
            msgbox(['Units with IsoINN exceeding max. value will '...
                'not be deleted!']);
        end
        toDelete = ( [handles.units.IsoINN] < featMin ) | ...
            ( [handles.units.IsoINN] > featMax );

        if featMin > handles.params.minIsoINN
            handles.params.minIsoINN = featMin;
        end

    case 'RSTD'
        if ~isinf(abs(featMin))
            featMin = -Inf;
            msgbox(['Units with RSTD less than the min. value will '...
                    'not be deleted!']);
        end
        toDelete = ( [handles.units.RSTD] < featMin ) | ...
                   ( [handles.units.RSTD] > featMax );
        if featMax < handles.params.maxRSTD
            handles.params.maxRSTD = featMax;
        end
        
    case 'skewness'
        if ~isinf(abs(featMin))
            featMin = -Inf;
            msgbox(['Units with skewness less than the min. value will '...
                    'not be deleted!']);
        end
        toDelete = ( [handles.skewn] < featMin ) | ...
                   ( [handles.skewn] > featMax );
        if featMin > handles.params.min_skewness
            handles.params.min_skewness = featMin;
        end
               
    case 'kurtosis'
        if ~isinf(featMax)
            featMax = Inf;
            msgbox(['Units with kurtosis exceeding max. value will '...
                    'not be deleted!']);
        end
        toDelete = ( [handles.kurtosis] < featMin ) | ...
            ( [handles.kurtosis] > featMax );
        if featMin > handles.params.minKurtosis
            handles.params.minKurtosis = featMin;
        end

    case 'SNR'
        msgbox(['SNR is not a selective feature! ' ...
            'No unit will be set to delete'])

    otherwise
end

[handles.units(toDelete).state] = deal(4);

guidata(hObject,handles);

currentUnit = str2double(handles.unitIDsAsStrings{get(handles.listbox1,'Value')});
updateunitstatefractions(hObject, handles, currentUnit);
%update background colors in unit list
listboxSortCriteria_Callback(hObject, eventdata, handles);

%redraw spatial positions of units to visualize removed redundancy
drawunitlabels(hObject, handles);
fprintf('Changed %g unit states to delete.\n',nnz(toDelete));



%--------------------------------------------------------------------------



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
function listboxUnitState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxUnitState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% {'unchecked', 'mixture', 'single', 'delete'} <-> 
%handles.stateColor = {[0 0 1],[0 1 1],[0 1 0],[1 0 0]}; 
handles.stateColor = {'white','aqua','lime','red'}; 
guidata(hObject, handles);
stateLabels = {'unchecked', 'mixture', 'single', 'delete'};
stateLabelsColored = arrayfun(@(x) ...
                    setbgcolor(handles.stateColor{x},stateLabels{x}),...
                                 1:length(stateLabels),'UniformOutput',false);
set(hObject, 'String', stateLabelsColored);
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
function textFractionSingle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textFractionSingle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function textFractionDelete_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textFractionDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function textFractionMixture_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textFractionMixture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function editCoinThr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCoinThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function editSimThr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSimThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function listboxSortCriteria_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxSortCriteria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%sort criteria.
handles.sortCriteria = {'index', 'separability', 'IsoIBg', 'IsoINN',...
                        'RSTD','skewness','kurtosis','SNR'};
guidata(hObject, handles);
%Set labels for unit state accordingly.
set(hObject, 'String', handles.sortCriteria);



% --- Executes during object creation, after setting all properties.
function editMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
featMin = -Inf;
set(hObject,'String',num2str(featMin));
handles.featMin = featMin;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function editMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
featMax = +Inf;
set(hObject,'String',num2str(featMax));
handles.featMax = featMax;
guidata(hObject, handles);




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
params = evalin('base','params');
%handles.data = evalin('base','data');
handles.data = readdatablock(params.filename,...
                    1,length(params.sensor_rows),...
                    1,length(params.sensor_cols),...
                    1,length(params.frameStartTimes));
                
handles.params = params;
%all sources.
S = cat(1,ROIs.S);
handles.S = S;

%all units.
units = unitsfromrois(ROIs);

if ~isfield(units,'noise_std')
    for i = 1:length(units)
        noise_std = median(abs(handles.S(i,:))/0.6745);
        units(i).noise_std = noise_std;
    end
end

if ~isfield(units,'skewn')
    fprintf('Computing skewness...\n');
    skewn = num2cell(skewness(S'));
    [units.skewn] = deal(skewn{:});
    clear skewn
end
handles.skewn = [units.skewn];

if ~isfield(units,'kurt')
    fprintf('Computing kurtosis...\n');
    kurt = num2cell(kurtosis(S'));
    [units.kurt] = deal(kurt{:});
end
handles.kurtosis = [units.kurt];

clear S


%Calculate IsoI for legacy results

if ~isfield(units,'IsoIBg') || ~isfield(units,'IsoINN')
    t1 = clock;
    fprintf('Compute isolation information measures for all units...\n');
    sr = handles.params.sr;
    upsample = handles.params.upsample;
    thrFactor = handles.params.thrFactor;
    for i = 1:length(units)
        noise_std = units(i).noise_std;
        tmpS = resample(handles.S(i,:),upsample,1);
        [unused,pos_amplitudes] = searchpeaks(-tmpS,...
            thrFactor*noise_std,ceil(sr*upsample));
        valid = pos_amplitudes > min(abs(units(i).amplitude));
        IsoI = isoi(pos_amplitudes(valid),pos_amplitudes(~valid));
        % no distinction between Bg and NN here:
        units(i).IsoIBg = IsoI;
        units(i).IsoINN = IsoI;
    end
    t2 = clock;
    fprintf('...done in %g seconds.\n',etime(t2,t1));
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

[handles.ROIs] = units2rois(handles.units, handles.ROIs);

% for i = 1:length(handles.ROIs);
%     %delete <-> state '4'
%     toDelete = [handles.units.k] == i & [handles.units.state] == 4;
    
%     if ~isfield(handles.ROIs(i),'units_del');
%         handles.ROIs(i).units_del = [];
%     end
%     handles.ROIs(i).units_del = ...
%                        [handles.ROIs(i).units_del handles.units(toDelete)];
    
%     if ~isfield(handles.ROIs(i),'S_del');
%         handles.ROIs(i).S_del = [];
%     end
%     handles.ROIs(i).S_del = [handles.ROIs(i).S_del; handles.S(toDelete,:)];
    
%     %map toDelete to index range in ROI.
%     toDelete = toDelete([handles.units.k] == i);
    
%     if ~isfield(handles.ROIs(i),'A_del'); 
%         handles.ROIs(i).A_del = [];
%     end
%     handles.ROIs(i).A_del = cat(2,handles.ROIs(i).A_del,...
%                                   handles.ROIs(i).A_tau(:,toDelete,:));    
%     clear toDelete

%     %single <-> state '3'
%     toSave = [handles.units.k] == i & [handles.units.state] <= 3;
%     handles.ROIs(i).units = handles.units(toSave);
%     handles.ROIs(i).S = handles.S(toSave,:);
%     %map toSave to index range in ROI.
%     toSave = toSave([handles.units.k] == i);
%     handles.ROIs(i).A_tau = handles.ROIs(i).A_tau(:,toSave,:);
%     clear toSave                                                          
% end

assignin('base','ROIs',handles.ROIs);
assignin('base','params',handles.params);

%Update guidata
guidata(hObject, handles);

%Save results to basic event file.
saveresults(handles.ROIs, handles.params);

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
% noise_std = median(abs(S)/0.6745);
%noise_std = std(S);
noise_std = unit.noise_std;
%thrF;actor = handles.params.thrFactor;
thrFactor = handles.params.thrFactor;
handles.threshold = -thrFactor*noise_std;

[indices,pos_amplitudes] = searchpeaks(-S,thrFactor*noise_std,ceil(sr*upsample));

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

%entire Source.
%axes(handles.axes1);
cla(handles.axes1,'reset');
hold(handles.axes1,'on');
%handles.sPlot = plot(handles.axes1,S,'Color',...
%                     handles.stateColor{unit.state});
handles.sPlot = plot(handles.axes1,S,'Color',[0 0 0]);

title(handles.axes1, strcat('amplSD = ',num2str(unit.amplitudeSD,2),...
      '; RSTD = ', num2str(unit.RSTD,2),'; sep. = ',...
      num2str(unit.separability,2),...
      '; IsoIBg = ', num2str(unit.IsoIBg,2),...
      '; IsoINN = ', num2str(unit.IsoINN,2),... 
      '; skewn. = ', num2str(handles.skewn(whichone),2),...
      '; kurt. = ', num2str(handles.kurtosis(whichone),2)));%,...
      %'; SNR = ', num2str(unit.snr)));

%all threshold crossings.
plot(handles.axes1, indices, S(indices),'ro','LineStyle','none');
handles.axes1thr = plot(handles.axes1, ...
    [0 N_SAMPLES], [handles.threshold handles.threshold],'r');
%those from units.
% handles.axes1amp = plot(round(unit.time*upsample*sr),...
%                       unit.amplitude,'go','LineStyle','none');

if ~isempty(amplitudes)
    %Source waveforms.
    %axes(handles.axes2);
    cla(handles.axes2,'reset');
    %Setup axes.
    hold(handles.axes2,'on');
    plot(handles.axes2, pks,'Color',[0.75 0.75 0.75]);
    handles.axes2thr = plot(handles.axes2,...
    [0 f_tot], [handles.threshold handles.threshold],'r');

    %amplitude histogram.
    %axes(handles.axes3);
    cla(handles.axes3,'reset');
    hold(handles.axes3,'on');
    bar(handles.axes3, bin_ctrs,cts,'k'); 
    handles.axes3thr = plot(handles.axes3,...
                            [handles.threshold handles.threshold],...
                            get(handles.axes3,'ylim'),'r');
end


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
set(handles.listboxUnitState,'Value',unit.state);

%update guidata.
guidata(hObject,handles);

thresholdslider_Callback(handles.thresholdslider, eventdata, handles);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function updateunitdata(hObject, handles, whichone, newThreshold)
%updateunitdata(hObject, handles,whichone, newThreshold)

%valid threshold crossings.
valid = handles.amplitudes <= newThreshold;

handles.units(whichone).time = handles.time(valid);
handles.units(whichone).amplitude = handles.amplitudes(valid);
handles.units(whichone).amplitudeSD = std(handles.units(whichone).amplitude);
handles.units(whichone).RSTD = ...
                            handles.units(whichone).amplitudeSD/...
                       mean(abs(handles.units(whichone).amplitude));
if nnz(~valid) > 0
    handles.units(whichone).separability = (...
        mean(abs(handles.units(whichone).amplitude)) - ...
        mean(abs(handles.amplitudes(~valid))) ) /...
        handles.units(whichone).noise_std;
    IsoI = isoi(handles.amplitudes(valid),handles.amplitudes(~valid));
    %IsoIBg and IsoINN are equivalent, because Bg is not clustered.
    handles.units(whichone).IsoIBg = IsoI;
    handles.units(whichone).IsoINN = IsoI;
else
    handles.units(whichone).separability = ...
        mean(abs(handles.units(whichone).amplitude)) / ...
        handles.units(whichone).noise_std;
    %IsoI measures are not defined because there are no other amplitudes
    handles.units(whichone).IsoIBg = NaN;
    handles.units(whichone).IsoINN = NaN;
end
    
handles.units(whichone).SDscore = [];

%which ROI <-> k.
k = handles.units(whichone).k;
% STA.
% dataTmp = reshape(handles.ROIs(k).X,...
%     [length(handles.ROIs(k).sensor_rows)...
%     length(handles.ROIs(k).sensor_cols) size(handles.ROIs(k).X,2)]);
dataTmp = handles.data(...
       (handles.ROIs(k).sensor_rows(1) <= handles.params.sensor_rows) & ...
       (handles.params.sensor_rows <= handles.ROIs(k).sensor_rows(end)),...
       (handles.ROIs(k).sensor_cols(1) <= handles.params.sensor_cols) & ...
       (handles.params.sensor_cols <= handles.ROIs(k).sensor_cols(end)),:);
handles.units(whichone).STA = computetemplate(dataTmp,...
                                              handles.time(valid),...
                                              handles.params.sr,0);
clear dataTmp
% unit position.
extrSTA = max(max(max(abs(handles.units(whichone).STA))));
[row_max,col_max] = find(max(abs(handles.units(whichone).STA),[],3) == extrSTA);
boss_row = handles.ROIs(k).sensor_rows(row_max);
boss_col = handles.ROIs(k).sensor_cols(col_max);
handles.units(whichone).boss_row = boss_row;
handles.units(whichone).boss_col = boss_col;
%handles.units(whichone).snr = extrSTA/...
%                            handles.params.sigma(...
%                            handles.params.sensor_rows == boss_row,...
%                            handles.params.sensor_cols == boss_col);
clear k
% update guidata
guidata(hObject,handles);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function updateunitplots(hObject, handles, whichone, newThreshold)
%updateunitplots(hObject, handles, whichone, newThreshold)

%valid threshold crossings.
valid = handles.amplitudes <= newThreshold;


%axes(handles.axes1);%DIRECT AXES CALL IS SLOW!!
if isfield(handles,'axes1amp');
    delete(handles.axes1amp(ishandle(handles.axes1amp)));
end
handles.axes1amp = plot(handles.axes1,handles.indices(valid),...
                      handles.amplitudes(valid),'go','LineStyle','none');
if isfield(handles,'axes1thr');
    delete(handles.axes1thr(ishandle(handles.axes1thr)));
end
handles.axes1thr = plot(handles.axes1, get(handles.axes1,'xlim'), ...
                        [newThreshold newThreshold],'r');

%axes(handles.axes2);
if isfield(handles,'axes2pks');
    delete(handles.axes2pks(ishandle(handles.axes2pks)));
end
handles.axes2pks = plot(handles.axes2,...
                        handles.pks(:,valid),'Color','black');
if isfield(handles,'axes2thr');
    delete(handles.axes2thr(ishandle(handles.axes2thr)));
end
handles.axes2thr = plot(handles.axes2,get(handles.axes2,'xlim'),...
                        [newThreshold newThreshold],'r');

%axes(handles.axes3);
if isfield(handles,'axes3thr');
    delete(handles.axes3thr(ishandle(handles.axes3thr)));
end
handles.axes3thr = plot(handles.axes3, [newThreshold newThreshold],...
                        get(handles.axes3,'ylim'),'r');

%draw spatial positions of units
axes(handles.axes4);
if isfield(handles,'unitLabels')
    delete(handles.unitLabels(ishandle(handles.unitLabels)));
end
drawunitlabels(hObject, handles);

% hold(handles.axes4,'on');
% if isfield(handles,'unitMarker') && ishandle(handles.unitMarker)
%     delete(handles.unitMarker);
% end
% handles.unitMarker = plot(handles.axes4,...
%                          handles.units(whichone).boss_col,handles.units(whichone).boss_row,...
%                         'o','MarkerSize',25,'Color','black','LineWidth',3);
if isfield(handles,'unitMarker') && ishandle(handles.unitMarker)
    hold(handles.axes4,'on');
    set(handles.unitMarker,...
        'XData',handles.units(whichone).boss_col,...
        'YData',handles.units(whichone).boss_row);
% set(handles.unitMarker,'XData',50,...
%         'YData',50);
else
hold(handles.axes4,'on');
%Only generate handle for the first time
handles.unitMarker = plot(handles.axes4,...
                          handles.units(whichone).boss_col,...
                          handles.units(whichone).boss_row,...
                        'o','MarkerSize',25,'Color','black','LineWidth',3);
end

% update guidata
guidata(hObject,handles);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function drawrois(hObject,handles)

axes(handles.axes4);
showrois(handles.ROIsAsCC,...
    'fillColumnFlag',logical(handles.params.d_col/handles.params.pitch-1),...
    'shuffle',1);

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function drawunitlabels(hObject, handles)

set(handles.axes4,'DrawMode','fast');
cols = [handles.units.boss_col];
rows = [handles.units.boss_row];
notDeleted = [handles.units.state] <= 3;
% cols = cols(notDeleted);
% rows = rows(notDeleted);
%idx = 1:length(handles.units);
%idx = idx(notDeleted);
% handles.unitLabels = text(cols, rows,...
%     arrayfun(@num2str,idx,'UniformOutput',false),...
%     'Parent',handles.axes4);
handles.unitLabels = text(cols(notDeleted), rows(notDeleted),...
    handles.unitIDsAsStringsSorted(notDeleted),...
    'Parent',handles.axes4);

guidata(hObject, handles);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function updateunitstatefractions(hObject, handles, currentUnit)

%state color Source activation plot
% set(handles.sPlot,...
%    'Color',handles.stateColor{handles.units(currentUnit).state});

%fraction of Source / unit states
N_units = length(handles.units);
fractionUnchecked = 100*nnz([handles.units.state] == 1)/N_units;
fractionMixture = 100*(nnz([handles.units.state] == 2) + ... 
                           nnz([handles.units.state] == 3) + ...
                           nnz([handles.units.state] == 4))/N_units;
fractionSingle = 100*nnz([handles.units.state] == 3)/N_units;
fractionDelete = 100*nnz([handles.units.state] == 4)/N_units;
set(handles.textFractionUnchecked,'String',num2str(fractionUnchecked,'%3.1f'));
set(handles.textFractionMixture,'String',num2str(fractionMixture,'%3.1f'));
set(handles.textFractionSingle,'String',num2str(fractionSingle,'%3.1f'));
set(handles.textFractionDelete,'String',num2str(fractionDelete,'%3.1f'));

%state of currently selected unit
set(handles.listboxUnitState,'Value',handles.units(currentUnit).state);

guidata(hObject, handles);
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
function bgColoredString = setbgcolor(color,string)
%bgColoredString = setbgcolor(color,string)
bgColoredString = ['<html><DIV bgcolor="' color '"><pre>' string ...
                   '                                 </pre></DIV></html>'];

%--------------------------------------------------------------------------









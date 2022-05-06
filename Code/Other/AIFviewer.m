function varargout = AIFviewer(varargin)
% AIFVIEWER MATLAB code for AIFviewer.fig
%      AIFVIEWER, by itself, creates a new AIFVIEWER or raises the existing
%      singleton*.
%
%      H = AIFVIEWER returns the handle to a new AIFVIEWER or the handle to
%      the existing singleton*.
%
%      AIFVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AIFVIEWER.M with the given input arguments.
%
%      AIFVIEWER('Property','Value',...) creates a new AIFVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AIFviewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AIFviewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AIFviewer

% Last Modified by GUIDE v2.5 25-Jul-2019 16:40:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AIFviewer_OpeningFcn, ...
                   'gui_OutputFcn',  @AIFviewer_OutputFcn, ...
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


% --- Executes just before AIFviewer is made visible.
function AIFviewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AIFviewer (see VARARGIN)

% if strcmpi(varargin{1},'dsc')
%     dsc = varargin{2};
% else
%     error('invalid input');
% end
%questdlg( );

% DSC is expected to be [row, col, time, slice] for perfusion
if isempty(varargin) || ischar(varargin{1})
    if isempty(varargin)
        selpath = uigetdir;
    end
    if ischar(varargin{1})
        selpath = varargin{1};
    end
    %datadir = dir(selpath);
    %dscpath = fullfile(selpath,'P001','ep2d_perf');
    dscpath = selpath;
    if ~exist(fullfile(dscpath,'1.dcm'))
        dscpath = selpath;
        if isempty(dir(fullfile(dscpath,'*.dcm')))
            f = errordlg(['Could not find dicoms in dsc path: ' dscpath],'File Error');
            error(['Could not find dicoms in dsc path: ' dscpath]);
        else
            dsc = loadDSCv2(dscpath);
        end
    else
        % Use this for 1.dcm, 2.dcm, ... format
        dsc = loadDSC(dscpath);
    end
    dsc = double(dsc);
%     if ~exist(dscpath,'dir')
%         f = errordlg(['Could not find dsc path: ' dscpath],'File Error');
%         error(['Could not find dsc path: ' dscpath]);
%     end
    %dsc = loadDSC(dscpath);
    handles.selpath = selpath;
    handles.dscpath = dscpath;
    handles.dsc = dsc;
    % Check for WM .mat and Results .mat
    if exist(fullfile(selpath,'WM_Mask_P001GE_M.mat'),'file') && ...
            exist(fullfile(selpath,'Result_MSwcf2','P001GE_M.mat'),'file')
        handles.calcCBF = 1;        
        handles = prepCBF(handles);       
    else
        handles.calcCBF = 0;
        set(handles.txt_curCBF,'Visible','off');
        set(handles.txt_curCBF_return,'Visible','off');
        set(handles.panel_perf,'Visible','off');
    end
    handles.debug_perf = 0;
    set(handles.button_CBFmap,'Enable','off');
    set(handles.button_CBFmap,'Visible','off');
else
    dsc = varargin{1};
    
    %dsc = smoothDSC(dsc);
    
    set(handles.panel_perf,'Visible','off');
    handles.calcCBF = 0;
    
    if length(varargin) > 1
        handles.images = varargin{2};
        %handles.image_names = varargin{3};
        handles.debug_perf = 1;
        set(handles.button_CBFmap,'Enable','off');
    else
        handles.debug_perf = 0;
        set(handles.button_CBFmap,'Enable','off');
        set(handles.toggle_conc,'Enable','off');
    end
end
set(handles.button_saveROI,'Enable','off');

% Initialize some vars
cslc = 1;
ctp = 1;
nslcs = size(dsc,4);
ntp = size(dsc,3);
handles.ntp = ntp;
handles.nslcs = nslcs;
handles.cslc = cslc;
handles.ctp = ctp;
handles.dsc = dsc;
handles.coord = [1 1];
handles.overlay = zeros(size(dsc,1),size(dsc,2),nslcs);
handles.dr = [];
handles.aifcoord = [1 1];
handles.savecoord = 0;
handles.aifloclist = [];

% Show initial dsc image
handles = update_imshow(handles);

% Plot initial aif
axes(handles.axes_aif);
aif = reshape(dsc(handles.aifcoord(2),handles.aifcoord(1),:,handles.cslc),1,handles.ntp);
plot(aif);

% Update text display
handles = update_textDisp(handles);

% Set slice slider range and increment step
set(handles.slider_slc,'Value',1);
set(handles.slider_slc,'Min',1);
set(handles.slider_slc,'Max',nslcs);
% Increment by one
if nslcs == 1
    set(handles.slider_slc,'Enable','off');
    %set(handles.slider_slc,'SliderStep',[0 0]);
else
    set(handles.slider_slc,'SliderStep',[1/(nslcs-1) 1/(nslcs-1)]);
end

% Set time slider range and increment step
set(handles.slider_time,'Value',1);
set(handles.slider_time,'Min',1);
set(handles.slider_time,'Max',ntp);
% Increment by one
if ntp == 1
    set(handles.slider_time,'Enable','off');
    %set(handles.slider_slc,'SliderStep',[0 0]);
else
    set(handles.slider_time,'SliderStep',[1/(ntp-1) 1/(ntp-1)]);
end

% Choose default command line output for AIFviewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AIFviewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AIFviewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_slc_Callback(hObject, eventdata, handles)
% hObject    handle to slider_slc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.cslc = round(get(hObject,'Value'));
handles = update_imshow(handles);
handles = update_textDisp(handles);
if handles.calcCBF
handles = update_curCBF(handles);
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_slc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_slc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_time_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.ctp = round(get(hObject,'Value'));
handles = update_imshow(handles);
handles = update_textDisp(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on mouse press over axes background.
function axes_dsc_ButtonDownFcn(hObject, eventdata)
% hObject    handle to axes_dsc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get coordinate of click (x,y)
handles = guidata(hObject);
coord = get(hObject,'CurrentPoint');
coord = round([coord(1,1) coord(1,2)]);
handles.aifcoord = coord;

% Plot aif in that coordinate
dsc = handles.dsc;
axes(handles.axes_aif);
aif = reshape(dsc(coord(2),coord(1),:,handles.cslc),1,handles.ntp);
plot(aif);
% hold on;
% plot(smooth(aif,.05,'sgolay',1));
% plot(smooth(aif,5,'sgolay',3));
% hold off;

handles.coord = coord;
handles = update_imshow(handles);
handles = update_textDisp(handles);
if handles.calcCBF
handles = update_curCBF(handles);
end
guidata(hObject,handles);


% Updates dsc img with current slice and timepoint
function handles = update_imshow(handles)

dsc = handles.dsc;
overlay = handles.overlay;
overlay = logical(overlay);
%overlay = bwperim(overlay);
dr = handles.dr;
if isempty(dr)
    dr(1) = min(dsc(:));
    dr(2) = max(dsc(:));
end

% Convert image to int8 and color
image = dsc(:,:,handles.ctp,handles.cslc);
image(image > dr(2)) = dr(2);
image(image < dr(1)) = dr(1);
image = image-dr(1);
scale = 255/(dr(2)-dr(1));%double(max(image(:)));
image = uint8(double(image).*scale);
imager = image;
imageg = image;
imageb = image;
% Show saved overlay in blue
imager(overlay(:,:,handles.cslc)) = 0;
imageb(overlay(:,:,handles.cslc)) = 255;
imageg(overlay(:,:,handles.cslc)) = 0;

% Show current coord in red
imager(handles.aifcoord(2),handles.aifcoord(1)) = 255;
imageg(handles.aifcoord(2),handles.aifcoord(1)) = 0;
imageb(handles.aifcoord(2),handles.aifcoord(1)) = 0;

% Show aif roi in green
if isfield(handles,'roi')
    imager(handles.roi(:,:,handles.cslc)) = 0;
    imageb(handles.roi(:,:,handles.cslc)) = 0;
    imageg(handles.roi(:,:,handles.cslc)) = image(handles.roi(:,:,handles.cslc));
end

% Save coordinates as overlay
if handles.savecoord
    overlay(handles.aifcoord(2),handles.aifcoord(1),handles.cslc) = 1;
    handles.overlay = overlay;
    handles.savecoord = 0;
end
olimg = cat(3,imager,imageg,imageb);

% Show image and make sure parent axes of image is "pickable"
handles.img_dsc = imshow(olimg,[],'Parent',handles.axes_dsc);
handles.axes_dsc = handles.img_dsc.Parent;
set(handles.axes_dsc,'Tag','axes_dsc');
set(handles.axes_dsc,'PickableParts','all');
set(handles.axes_dsc,'ButtonDownFcn',@axes_dsc_ButtonDownFcn);
set(handles.img_dsc,'PickableParts','none');


% Show and update image info in text display 
function handles = update_textDisp(handles)

dsc = handles.dsc;
imgVals = sprintf('Slc: %d \t Time: %d \t X(col): %d (%d) \t Y(row): %d (%d) \t Value: %f',...
    handles.cslc,handles.ctp,handles.coord(1),handles.aifcoord(1),handles.coord(2),handles.aifcoord(2),dsc(handles.coord(2),handles.coord(1),handles.ctp,handles.cslc));
set(handles.txt_imgVals,'String',imgVals);

if handles.debug_perf

    oaif = double(reshape(dsc(handles.aifcoord(2),handles.aifcoord(1),:,handles.cslc),1,handles.ntp));
    %aif = smooth(aif,.05,'sgolay',1)';
    
%     tmpfit = fit([1:length(oaif)]',oaif','smoothingspline','smoothingparam',.8);
%             aif = reshape(tmpfit(1:length(oaif)),[1 length(oaif)]);
%             aif(aif < 0) = 0;
    aif = oaif;
    [tmpibat,tmpbat,tmpbrt] = findBAT_test(aif);
    
%     ibat = handles.images{strcmpi(handles.image_names,'IBATmap')};
%     bat = handles.images{strcmpi(handles.image_names,'BATmap')};
%     brt = handles.images{strcmpi(handles.image_names,'BRTmap')};
%     atd = handles.images{strcmpi(handles.image_names,'ATDmap')};
%     cbf = handles.images{strcmpi(handles.image_names,'qCBF_nSVD')};
%     conc = handles.images{strcmpi(handles.image_names,'dBATmap')};
%     tmax = handles.images{strcmpi(handles.image_names,'TmaxMap_nSVD')};
%     aifVals = sprintf('iBAT: (%d) %d \t BAT: (%d) %d \t BRT: (%d) %d \t ATD: %d \t CBF: %.2f \t Conc: %.2f \t Tmax: %.2f',...
%         tmpibat,...
%         ibat(handles.aifcoord(2),handles.aifcoord(1),handles.cslc),...
%         tmpbat,...
%         bat(handles.aifcoord(2),handles.aifcoord(1),handles.cslc),...
%         tmpbrt,...
%         brt(handles.aifcoord(2),handles.aifcoord(1),handles.cslc),...
%         atd(handles.aifcoord(2),handles.aifcoord(1),handles.cslc),...
%         cbf(handles.aifcoord(2),handles.aifcoord(1),handles.cslc),...
%         conc(handles.aifcoord(2),handles.aifcoord(1),handles.cslc),...
%         tmax(handles.aifcoord(2),handles.aifcoord(1),handles.cslc));

    ibat = handles.images.etc.IBATmap;
    bat = handles.images.etc.BATmap;
    brt = handles.images.etc.BRTmap;
    aifVals = sprintf('iBAT: %d \t BAT: %d \t BRT: %d \t (%d, %d, %d)',tmpibat,tmpbat,tmpbrt,...
        ibat(handles.aifcoord(2),handles.aifcoord(1),handles.cslc),...
        bat(handles.aifcoord(2),handles.aifcoord(1),handles.cslc),...
        brt(handles.aifcoord(2),handles.aifcoord(1),handles.cslc));

    set(handles.txt_aifVals,'String',aifVals,'FontSize',13);
else
    oaif = double(reshape(dsc(handles.aifcoord(2),handles.aifcoord(1),:,handles.cslc),1,handles.ntp));
    %aif = smooth(oaif,max(.05*handles.ntp,5),'sgolay',1)';
    
%     tmpfit = fit([1:length(oaif)]',oaif','smoothingspline','smoothingparam',.8);
%             aif = reshape(tmpfit(1:length(oaif)),[1 length(oaif)]);
%             aif(aif < 0) = 0;
    aif = oaif;
            
    [tmpibat,tmpbat,tmpbrt] = findBAT_test(aif);
    [a,b,c] = findBAT(aif);
    aifVals = sprintf('iBAT: %d \t BAT: %d \t BRT: %d \t (%d, %d, %d)',tmpibat,tmpbat,tmpbrt,a,b,c);
    set(handles.txt_aifVals,'String',aifVals,'FontSize',13);
end


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get current coordinates of mouse cursor and display it
dsc = handles.dsc;
coord = get(handles.axes_dsc,'CurrentPoint');
coord = round([coord(1,1) coord(1,2)]);

if coord(1) >= 1 && coord(2) >= 1 &&...
        coord(1) <= size(dsc,2) && coord(2) <= size(dsc,1)
    handles.coord = coord;
    handles = update_textDisp(handles);
end

guidata(hObject,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

cobj = gco;
if ~sum(strcmpi(cobj.Tag,{'slider_slc','slider_time','edit_drlow','edit_drhigh'}))
%set(handles.slider_slc,'Enable','off');
% Move aif coordinates by arrow keys and save by hitting return (enter key)
key = eventdata.Key;
switch key
    case 'escape'
        close;
        return;
    case 'rightarrow'
        if handles.aifcoord(1) < size(handles.dsc,2)
            handles.aifcoord = handles.aifcoord + [1 0];
        end
    case 'leftarrow'
        if handles.aifcoord(1) > 1
            handles.aifcoord = handles.aifcoord + [-1 0];
        end
    case 'uparrow'
        if handles.aifcoord(2) > 1
            handles.aifcoord = handles.aifcoord + [0 -1];
        end
    case 'downarrow'
        if handles.aifcoord(2) < size(handles.dsc,1)
            handles.aifcoord = handles.aifcoord + [0 1];
        end
    case 'return'
        handles.savecoord = 1;
        handles.aifloclist = [handles.aifloclist; handles.aifcoord(2) handles.aifcoord(1) handles.cslc];
        if handles.calcCBF
            handles = update_avgCBF(handles);
        end
end

dsc = double(handles.dsc);
axes(handles.axes_aif);
oaif = reshape(dsc(handles.aifcoord(2),handles.aifcoord(1),:,handles.cslc),1,handles.ntp);
% plot(aif);
% hold on;
% plot(smooth(aif,.05,'sgolay',1));
% %plot(smooth(aif,.05,'sgolay',3));
% plot(smooth(aif,5,'sgolay',3));


if handles.debug_perf
    %ibat = handles.images{strcmpi(handles.image_names,'IBATmap')};
    %bat = handles.images{strcmpi(handles.image_names,'BATmap')};
    %brt = handles.images{strcmpi(handles.image_names,'BRTmap')};
    
    ibat = handles.images.etc.IBATmap;
    bat = handles.images.etc.BATmap;
    brt = handles.images.etc.BRTmap;
    
    ibat = ibat(handles.aifcoord(2),handles.aifcoord(1),handles.cslc);
    bat = bat(handles.aifcoord(2),handles.aifcoord(1),handles.cslc);
    brt = brt(handles.aifcoord(2),handles.aifcoord(1),handles.cslc);
    
    if handles.toggle_conc.Value
        aif = smooth(oaif,5,'sgolay',3);
        conc = Sig2Conct_AIF(aif, [ibat bat brt], 30);
        plot(conc.Ct);
        hold on;
        %plot(smooth(conc.Ct,.05,'sgolay',1));
        aif = smooth(oaif,.05,'sgolay',3);
        conc = Sig2Conct_AIF(aif, [ibat bat brt], 30);
        plot(conc.Ct);
        xline(tmpbat,'r--');
        xline(tmpbrt,'r--');
        hold off;
    else
        
%         tmpfit = fit([1:length(oaif)]',oaif','smoothingspline','smoothingparam',.8);
%             aif = reshape(tmpfit(1:length(oaif)),[1 length(oaif)]);
%             aif(aif < 0) = 0;
            
        aif = oaif;
        %aif = smooth(oaif,5,'sgolay',3);    
        [tmpibat,tmpbat,tmpbrt] = findBAT_test(aif);
        %aif = smooth(oaif,.1,'sgolay',1);
        plot(oaif);
        hold on;
        plot(aif,'linewidth',2);
        %aif = smooth(oaif,.05,'sgolay',3);
        %plot(aif,'linewidth',2);
        xline(tmpbat,'r--');
        xline(tmpbrt,'r--');
        xline(bat,'k--');
        xline(brt,'k--');
        hold off;
    end
else
    %plot(oaif);
        %hold on;
        
%         tmpfit = fit([1:length(oaif)]',oaif','smoothingspline','smoothingparam',.8);
%             aif = reshape(tmpfit(1:length(oaif)),[1 length(oaif)]);
%             aif(aif < 0) = 0;
            
        aif = oaif;
        %aif = smooth(oaif,max(.05*handles.ntp,5),'sgolay',3)';        
        %plot(aif,'linewidth',2);
%         aif = smooth(oaif,5,'sgolay',3);
%         plot(aif);
        [tmpibat,tmpbat,tmpbrt] = findBAT_test(aif);
        
        %conc = Sig2Conct_AIF(aif, [tmpibat tmpbat tmpbrt], 30);
        %aif = conc.Ct;
        plot(oaif);
        hold on;
        plot(aif,'linewidth',2);
        xline(tmpbat,'r--');
        xline(tmpbrt,'r--');
        hold off;
        
end

handles.coord = handles.aifcoord;

handles = update_textDisp(handles);
if handles.calcCBF
handles = update_curCBF(handles);
end
handles = update_imshow(handles);
%set(handles.slider_slc,'Enable','on');
guidata(hObject,handles);
end


function edit_drlow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_drlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_drlow as text
%        str2double(get(hObject,'String')) returns contents of edit_drlow as a double

if ~isnan(str2double(get(handles.edit_drhigh,'String'))) &&...
        ~isnan(str2double(get(handles.edit_drlow,'String')))
    handles.dr = [str2double(get(handles.edit_drlow,'String')) str2double(get(handles.edit_drhigh,'String'))];
    handles = update_textDisp(handles);
    handles = update_imshow(handles);
    guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function edit_drlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_drlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_drhigh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_drhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_drhigh as text
%        str2double(get(hObject,'String')) returns contents of edit_drhigh as a double

if ~isnan(str2double(get(handles.edit_drhigh,'String'))) &&...
        ~isnan(str2double(get(handles.edit_drlow,'String')))
    handles.dr = [str2double(get(handles.edit_drlow,'String')) str2double(get(handles.edit_drhigh,'String'))];
    handles = update_textDisp(handles);
    handles = update_imshow(handles);
    guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function edit_drhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_drhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_saveAIFloc.
function button_saveAIFloc_Callback(hObject, eventdata, handles)
% hObject    handle to button_saveAIFloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','aifloclist',handles.aifloclist);


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

if eventdata.VerticalScrollCount > 0
    if handles.cslc < handles.nslcs
        handles.cslc = handles.cslc + 1;
    end
elseif eventdata.VerticalScrollCount < 0
    if handles.cslc > 1
        handles.cslc = handles.cslc - 1;
    end
end

set(handles.slider_slc,'Value',handles.cslc);
handles = update_textDisp(handles);
if handles.calcCBF
handles = update_curCBF(handles);
end
handles = update_imshow(handles);
guidata(hObject,handles);
    

% --- Executes on button press in button_autoAIF.
function button_autoAIF_Callback(hObject, eventdata, handles)
% hObject    handle to button_autoAIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

overlay = handles.overlay;
dscpath = handles.dscpath;
% if sum(overlay(:))
%     f = warndlg(

% Create dialog to wait for Auto AIF to finish
cf = handles.figure1;
fpos = get(cf,'Position');
d = dialog('Name','Auto AIF',...
            'Position',[fpos(1)+fpos(3)/2-fpos(3)/16 fpos(2)+fpos(4)/2-fpos(4)/16 fpos(3)/8 fpos(4)/8],...
            'WindowStyle','modal');
dpos = get(d,'Position');
txt = uicontrol('Parent',d,...
            'Position',[dpos(3)/2-dpos(3)/4 2*dpos(4)/3-dpos(4)/8 dpos(3)/2 dpos(4)/4],...
            'Style','text',...
            'String','Running Auto AIF...');
btn = uicontrol('Parent',d,...   
            'Position',[dpos(3)/2-dpos(3)/8 dpos(4)/3-dpos(4)/8 dpos(3)/4 dpos(4)/4],...
            'String','Done',...
            'Enable','off',...
            'Callback','delete(gcf)');
drawnow;
global sampdelayfix;
sampdelayfix = 1;
% Run autoAIF code
[Signal_AIF,cutoffs_AIF,position,Conct_AIF,Images,n_AIFslice,AIF,positionAIF] = autoaif_wy_Philips_v2(dscpath,handles.nslcs,handles.ntp,[]);
%positionAIF = {[1 1 1]};
sampdelayfix = 0;
%positionAIF;
overlay = zeros(size(overlay));
aifloclist = [];
for ii = 1:length(positionAIF)
    tmppos = positionAIF{ii};
    overlay(tmppos(1),tmppos(2),tmppos(3)) = 1;
    aifloclist = [aifloclist; tmppos];
end
handles.overlay = overlay;
handles.aifloclist = aifloclist;

% Update after dialog closes
set(btn,'Enable','on');
uiwait(d);

handles = update_imshow(handles);
handles = update_avgCBF(handles);
guidata(hObject,handles);


% Load DSC dicoms
function dsc = loadDSC(dscpath)

dscpath = checkFileSep(dscpath);
dcmdir = dir(fullfile(dscpath,'*dcm'));

tmpdsc = [];
tmpdsc_stack = [];
f = waitbar(0,'Loading DSC');
tmphdr = dicominfo(fullfile(dscpath,'1.dcm'));

if isfield(tmphdr,'NumberOfTemporalPositions')
    ntp = tmphdr.NumberOfTemporalPositions;
    N_slices = length(dcmdir)/ntp;
else
    for ii = 1:length(dcmdir)
        tmphdr = dicominfo(fullfile(dscpath, [num2str(ii) '.dcm']));
        slclist(ii) = tmphdr.SliceLocation;
    end
    N_slices = length(unique(slclist));
    ntp = length(dcmdir)/N_slices;
end

%ntp = tmphdr.NumberOfTemporalPositions;
for ii = 1:length(dcmdir)
    tmpdsc = cat(3,tmpdsc,dicomread(fullfile(dscpath,[num2str(ii) '.dcm'])));
    if mod(ii,ntp) == 0
        tmpdsc_stack = cat(4,tmpdsc_stack,tmpdsc);
        tmpdsc = [];
    end
    waitbar(ii/length(dcmdir),f);
end
dsc = tmpdsc_stack;
close(f);

% Load DSC dicoms
function dsc = loadDSCv2(dscpath)

dscpath = checkFileSep(dscpath);
dcmdir = dir(fullfile(dscpath,'*dcm'));

tmpdsc = [];
tmpdsc_stack = [];
f = waitbar(0,'Loading DSC');
tmphdr = dicominfo(fullfile(dscpath,dcmdir(1).name));

if isfield(tmphdr,'NumberOfTemporalPositions')
    ntp = tmphdr.NumberOfTemporalPositions;
    N_slices = length(dcmdir)/ntp;
else
    for ii = 1:length(dcmdir)
        tmphdr = dicominfo(fullfile(dscpath, dcmdir(ii).name));
        slclist(ii) = tmphdr.SliceLocation;
    end
    N_slices = length(unique(slclist));
    ntp = length(dcmdir)/N_slices;
end

%ntp = tmphdr.NumberOfTemporalPositions;
for ii = 1:length(dcmdir)
    tmpdsc = cat(3,tmpdsc,dicomread(fullfile(dscpath,dcmdir(ii).name)));
    if mod(ii,ntp) == 0
        tmpdsc_stack = cat(4,tmpdsc_stack,tmpdsc);
        tmpdsc = [];
    end
    waitbar(ii/length(dcmdir),f);
end
dsc = tmpdsc_stack;
close(f);

% --- Executes on button press in button_clearAIF.
function button_clearAIF_Callback(hObject, eventdata, handles)
% hObject    handle to button_clearAIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dsc = handles.dsc;
handles.overlay = zeros(size(dsc,1),size(dsc,2),size(dsc,4));
handles.aifloclist = [];
handles = update_imshow(handles);
if handles.calcCBF
    handles = update_avgCBF(handles);
end

guidata(hObject,handles);


% --- Executes on button press in button_viewAIF.
function button_viewAIF_Callback(hObject, eventdata, handles)
% hObject    handle to button_viewAIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dsc = handles.dsc;
aifloclist = handles.aifloclist;

if isempty(aifloclist)
    f = warndlg('No AIF data stored.');
else
    favg = figure;
    fall = figure;
    avgsignal = zeros(1,size(dsc,3));
    for ii = 1:size(aifloclist,1)
        tmppos = aifloclist(ii,:);
        tmpsignal = double(reshape(dsc(tmppos(1),tmppos(2),:,tmppos(3)),[1 size(dsc,3)]));
        
        figure(fall);subplot(round(sqrt(length(aifloclist))),ceil(sqrt(length(aifloclist))),ii);
        plot(tmpsignal);title(['[' num2str(tmppos(1)) ',' num2str(tmppos(2)) ',' num2str(tmppos(3)) ']']);
        
        avgsignal = avgsignal + tmpsignal;
    end   
    avgsignal = avgsignal/length(aifloclist);
    figure(favg);
    plot(avgsignal);title('Average');
end


% --- Executes on button press in button_roiAIF.
function button_roiAIF_Callback(hObject, eventdata, handles)
% hObject    handle to button_roiAIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dsc = handles.dsc;
axes(handles.axes_dsc);
roi = roipoly;
if ~isempty(roi)
handles.roi = false([size(roi) handles.nslcs]);
handles.roi(:,:,handles.cslc) = roi;
tmpdsc = reshape(dsc(:,:,:,handles.cslc),[size(dsc,1)*size(dsc,2) handles.ntp]);
roi = reshape(roi,[size(dsc,1)*size(dsc,2) 1]);
meansig = mean(tmpdsc(roi,:),1);
figure;
plot(meansig);

handles.meansig = meansig;
handles = update_imshow(handles);
set(handles.button_saveROI,'Enable','on');
end
guidata(hObject,handles);


% --- Executes on button press in button_saveROI.
function button_saveROI_Callback(hObject, eventdata, handles)
% hObject    handle to button_saveROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','aif_curve',handles.meansig);
% aif_roi = zeros([size(handles.roi) handles.nslcs]);
% aif_roi(:,:,handles.cslc) = handles.roi;
assignin('base','aif_roi',handles.roi);


% Update text display for CBF
function handles = update_curCBF(handles)

dsc = handles.dsc;
CmtWM = handles.CmtWM;
tmpsignal = double(reshape(dsc(handles.aifcoord(2),handles.aifcoord(1),:,handles.cslc),[1 handles.ntp]));
[ibat, bat, brt] = findBAT(tmpsignal);
if ibat > 0 && bat > 0 && brt < 200 && bat < CmtWM(2)+5 && tmpsignal(1) > 400
    AIF = Sig2Conct_AIF(tmpsignal, [ibat bat brt], handles.echoTime);
    fResidue = fResidue_nSVD(AIF,handles.Cmt,handles.Dt,handles.Wij);
    CBF = 100*60*handles.Kh/handles.psi*max(fResidue)*handles.WCF*handles.ratio_CBV;
    set(handles.txt_curCBF_return,'String',num2str(CBF));
    set(handles.txt_curBATs_return,'String',['[' num2str(ibat) ',' num2str(bat) ',' num2str(brt) ']']);
else
    set(handles.txt_curCBF_return,'String','bad aif');
    set(handles.txt_curBATs_return,'String',['[' num2str(ibat) ',' num2str(bat) ',' num2str(brt) ']']);
end


% Update text display for all average CBF
function handles = update_avgCBF(handles)

dsc = handles.dsc;
CmtWM = handles.CmtWM;
aifloclist = handles.aifloclist;

avgsignal = zeros(1,handles.ntp);
for ii = 1:size(aifloclist,1)
    tmploc = aifloclist(ii,:);
    tmpsignal = double(reshape(dsc(tmploc(1),tmploc(2),:,tmploc(3)),[1 handles.ntp]));
    avgsignal = avgsignal + tmpsignal;
end
avgsignal = avgsignal/size(aifloclist,1);

[ibat, bat, brt] = findBAT(avgsignal);
if ibat > 0 && bat > 0 && brt < 200 && bat < CmtWM(2)+5 && avgsignal(1) > 400
    AIF = Sig2Conct_AIF(avgsignal, [ibat bat brt], handles.echoTime);
    fResidue = fResidue_nSVD(AIF,handles.Cmt,handles.Dt,handles.Wij);
    CBF = 100*60*handles.Kh/handles.psi*max(fResidue)*handles.WCF*handles.ratio_CBV;
    set(handles.txt_allavgCBF_return,'String',num2str(CBF));
    set(handles.txt_allavgBATs_return,'String',['[' num2str(ibat) ',' num2str(bat) ',' num2str(brt) ']']);
else
    set(handles.txt_allavgCBF_return,'String','bad aif');
    set(handles.txt_allavgBATs_return,'String',['[' num2str(ibat) ',' num2str(bat) ',' num2str(brt) ']']);
end


% Prepare data for CBF calculation
function handles = prepCBF(handles)

% Get WM signal
global glblTargetPath;
global injectionNum;
global seqType;
global sampdelayfix;
sampdelayfix = 0;
glblTargetPath = handles.selpath;
injectionNum = 1;
seqType = 'GE_M';
path_DSC = handles.dscpath;
[wm_signal,CmtWM] = getWMCtBAT(path_DSC,size(handles.dsc,4),size(handles.dsc,3));

% Header info
header = dicominfo([path_DSC '\1.dcm']);
Dt = (header.RepetitionTime)/1000;
echoTime = header.EchoTime;
Kh = (1-0.45)/(1-0.25);
psi = 1.04;
Wij = 0.2;

resultdata = load(fullfile(handles.selpath,'Result_MSwcf2','P001GE_M.mat'));
ratio_CBV = resultdata.ROIs.values.CF.ratio_CBV;
WCF = resultdata.ROIs.values.CF.WCF;

% Get WM conc signal
Cmt = Sig2Conct_AIF(wm_signal, [CmtWM(1:3)], echoTime);
set(handles.txt_wmBATs_return,'String',['[' num2str(CmtWM(1)) ',' num2str(CmtWM(2)) ',' num2str(CmtWM(3)) ']']);

handles.Dt = Dt;
handles.echoTime = echoTime;
handles.Kh = Kh;
handles.psi = psi;
handles.Wij = Wij;
handles.ratio_CBV = ratio_CBV;
handles.WCF = WCF;
handles.Cmt = Cmt;
handles.CmtWM = CmtWM;
handles.wm_signal = wm_signal;


% --- Executes on button press in button_showWMAIF.
function button_showWMAIF_Callback(hObject, eventdata, handles)
% hObject    handle to button_showWMAIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
plot(handles.wm_signal);title('WM Signal');


% --- Executes on button press in button_CBFmap.
function button_CBFmap_Callback(hObject, eventdata, handles)
% hObject    handle to button_CBFmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
%cbf = handles.images{strcmpi(handles.image_names,'qCBF_nSVD')};
cbf = handles.images.DD.qCBF_nSVD;
imshow(cbf(:,:,handles.cslc),[0 240]);colormap(gca,'jet');colorbar;


% --- Executes on button press in toggle_conc.
function toggle_conc_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_conc

isconc = eventdata.Source.Value;
dsc = handles.dsc;
axes(handles.axes_aif);
oaif = double(reshape(dsc(handles.aifcoord(2),handles.aifcoord(1),:,handles.cslc),1,handles.ntp));
if handles.debug_perf
    %ibat = handles.images{strcmpi(handles.image_names,'IBATmap')};
    %bat = handles.images{strcmpi(handles.image_names,'BATmap')};
    %brt = handles.images{strcmpi(handles.image_names,'BRTmap')};
    
    ibat = handles.images.etc.IBATmap;
    bat = handles.images.etc.BATmap;
    brt = handles.images.etc.BRTmap;
    
    tmpibat = ibat(handles.aifcoord(2),handles.aifcoord(1),handles.cslc);
    tmpbat = bat(handles.aifcoord(2),handles.aifcoord(1),handles.cslc);
    tmpbrt = brt(handles.aifcoord(2),handles.aifcoord(1),handles.cslc);
    
    if isconc
        aif = smooth(oaif,5,'sgolay',3);
        conc = Sig2Conct_AIF(aif, [tmpibat tmpbat tmpbrt], 30);
        plot(conc.Ct);
        hold on;
        %plot(smooth(conc.Ct,.05,'sgolay',1));
        aif = smooth(oaif,.05,'sgolay',1);
        conc = Sig2Conct_AIF(aif, [tmpibat tmpbat tmpbrt], 30);
        plot(conc.Ct);
        xline(tmpbat,'r--');
        xline(tmpbrt,'r--');
        hold off;
    else
        plot(oaif);
        hold on;
        aif = smooth(oaif,.05,'sgolay',1);
        plot(aif);
        aif = smooth(oaif,5,'sgolay',3);
        plot(aif);
        xline(tmpbat,'r--');
        xline(tmpbrt,'r--');
        hold off;
    end
end
guidata(hObject,handles);


function ndsc = smoothDSC(odsc)

odsc = double(odsc);
for kk = 1:size(odsc,4)
    tmpmask = automaskns(imfilter(odsc(:,:,:,kk),fspecial('gaussian',[3 3],2)));
    ndsc(:,:,:,kk) = fSpatialFilter(odsc(:,:,:,kk),[2 3],tmpmask);
end
    

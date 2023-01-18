function varargout = digitzer(varargin)
% DIGITZER Graph Digitizing Application Beta 0.1 (Use at your own risk.)
%           Accuracy is not garenteed or implied.
%
%   Released by Andrew March
%       Massachusetts Institute of Technology
%
%   Helps to digitize scanned data, or plots contained in image files.
%     Current LIMITATIONs in Beta 0.1:
%         - Rotation of plots has not been accounted for.
%         - Initial window size is not based on screen resolution (resize
%         does work)
%         - Points list is not sorted, and that is not an option
%         - When point selection is begun, previous list of points is
%         purged and new list is begun fresh.
%         - Only works with linear scales, log and semi-log to be added.
%
% For help with the code, or to report bugs, email: amarch@mit.edu
%       Version 0.1 released 10/13/2006

%      DIGITZER, by itself, creates a new DIGITZER or raises the existing
%      singleton*.
%
%      H = DIGITZER returns the handle to a new DIGITZER or the handle to
%      the existing singleton*.
%
%      DIGITZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIGITZER.M with the given input arguments.
%
%      DIGITZER('Property','Value',...) creates a new DIGITZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before digitzer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to digitzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help digitzer

% Last Modified by GUIDE v2.5 11-Oct-2006 17:43:01

% Begin initialization code - DO NOT EDIT
global X;

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @digitzer_OpeningFcn, ...
                   'gui_OutputFcn',  @digitzer_OutputFcn, ...
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

% --- Executes just before digitzer is made visible.
function digitzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to digitzer (see VARARGIN)

% Choose default command line output for digitzer
global X;
X=zeros(4,2);
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% hObject
% handles
% scrsz = get(0,'ScreenSize');
% % %set(fig,'Position',[scrsz(1)+.1*scrsz(3),scrsz(2)+.1*scrsz(4),.8*scrsz(3
% % ),.8*scrsz(4)]);
% set(handles.figure1,'Position',[scrsz(1)+.1*scrsz(3),scrsz(2)+.1*scrsz(4)
% ,.8*scrsz(3),.8*scrsz(4)]);


% --- Outputs from this function are returned to the command line.
function varargout = digitzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function file_name_Callback(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_name as text
%        str2double(get(hObject,'String')) returns contents of file_name as a double


% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_file.
function open_file_Callback(hObject, eventdata, handles)
% hObject    handle to open_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
       {'*.bmp;*.jpg;*.jpeg;*.gif;*.png', 'Image Files (*.bmp, *.jpg, *.jpeg, *.gif,*.png)'; ...
        '*.*',                   'All Files (*.*)'}, ...
        'Pick an Image file');
while(isequal(filename,0) || isequal(pathname,0))
   [filename, pathname] = uigetfile( ...
       {'*.bmp;*.jpg;*.jpeg;*.gif;*.png', 'Image Files (*.bmp, *.jpg, *.jpeg, *.gif,*.png)'; ...
        '*.*',                   'All Files (*.*)'}, ...
        'Pick an Image file');
end
set(handles.file_name,'String',[pathname,filename])
fid=imread([pathname,filename]);
axes(handles.imag1);
image(fid);




% --- Executes on button press in x0.
function x0_Callback(hObject, eventdata, handles)
% hObject    handle to x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X;

[x,y]=ginput(1);
X(1,1)=x;
X(1,2)=y;
set(handles.x0x,'String',x);
set(handles.x0y,'String',y);
axes(handles.imag1);
hold on;
plot(x,y,'*r');
hold off;


function x0x_Callback(hObject, eventdata, handles)
% hObject    handle to x0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x0x as text
%        str2double(get(hObject,'String')) returns contents of x0x as a double


% --- Executes during object creation, after setting all properties.
function x0x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x0y_Callback(hObject, eventdata, handles)
% hObject    handle to x0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x0y as text
%        str2double(get(hObject,'String')) returns contents of x0y as a double


% --- Executes during object creation, after setting all properties.
function x0y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in y0.
function y0_Callback(hObject, eventdata, handles)
% hObject    handle to y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X;
[x,y]=ginput(1);
X(2,1)=x;
X(2,2)=y;
set(handles.y0x,'String',x);
set(handles.y0y,'String',y);
axes(handles.imag1);
hold on;
plot(x,y,'*r');
hold off;


function y0x_Callback(hObject, eventdata, handles)
% hObject    handle to y0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y0x as text
%        str2double(get(hObject,'String')) returns contents of y0x as a double


% --- Executes during object creation, after setting all properties.
function y0x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y0y_Callback(hObject, eventdata, handles)
% hObject    handle to y0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y0y as text
%        str2double(get(hObject,'String')) returns contents of y0y as a double


% --- Executes during object creation, after setting all properties.
function y0y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y0y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in xe.
function xe_Callback(hObject, eventdata, handles)
% hObject    handle to xe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X;
[x,y]=ginput(1);
X(3,1)=x;
X(3,2)=y;
set(handles.xex,'String',x);
set(handles.xey,'String',y);
axes(handles.imag1);
hold on;
plot(x,y,'*r');
hold off;


function xex_Callback(hObject, eventdata, handles)
% hObject    handle to xex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xex as text
%        str2double(get(hObject,'String')) returns contents of xex as a double


% --- Executes during object creation, after setting all properties.
function xex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xey_Callback(hObject, eventdata, handles)
% hObject    handle to xey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xey as text
%        str2double(get(hObject,'String')) returns contents of xey as a double


% --- Executes during object creation, after setting all properties.
function xey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ye.
function ye_Callback(hObject, eventdata, handles)
% hObject    handle to ye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X;
[x,y]=ginput(1);
X(4,1)=x;
X(4,2)=y;
set(handles.yex,'String',x);
set(handles.yey,'String',y);
axes(handles.imag1);
hold on;
plot(x,y,'*r');
hold off;

function yex_Callback(hObject, eventdata, handles)
% hObject    handle to yex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yex as text
%        str2double(get(hObject,'String')) returns contents of yex as a double


% --- Executes during object creation, after setting all properties.
function yex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yey_Callback(hObject, eventdata, handles)
% hObject    handle to yey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yey as text
%        str2double(get(hObject,'String')) returns contents of yey as a double


% --- Executes during object creation, after setting all properties.
function yey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in points.
function points_Callback(hObject, eventdata, handles)
% hObject    handle to points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns points contents as cell array
%        contents{get(hObject,'Value')} returns selected item from points


% --- Executes during object creation, after setting all properties.
function points_CreateFcn(hObject, eventdata, handles)
% hObject    handle to points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in output.
function output_Callback(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vals=[];
[filename, pathname] = uiputfile( ...
       {        '*.*',                   'All Files (*.*)'}, ...
        'Create and ASCII file');
if(~isequal(filename,0) && ~isequal(pathname,0))
    fid=fopen([pathname,filename],'w');
    points=get(handles.points,'String');
    for i=1:length(points(:,1))
        fprintf(fid,[points(i,:),'\n']);
        if(i>=2)
            vals=[vals;str2double(points(i,1:10)),str2double(points(i,11:end))];
        end
    end
    fclose(fid);
    assignin('base','Digitizer_Points',vals);
else
    warndlg('Cancelled by user.')
end

% --- Executes on button press in digitize.
function digitize_Callback(hObject, eventdata, handles)
% hObject    handle to digitize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global X;
i=1;
bad=false;
while i<=4 && ~bad
    if(X(i,1)<=0 || X(i,2)<=0)
        bad=true;
    else
        i=i+1;
    end
end
x0=str2num(get(handles.x0v,'String'));
y0=str2num(get(handles.y0v,'String'));
xe=str2num(get(handles.xev,'String'));
ye=str2num(get(handles.yev,'String'));
if(isnan(x0) || isnan(y0) || isnan(xe) || isnan(ye))
    bad=true;
    errordlg('Values assocaited with axis endpoints must be set.','Set Endpoint Values');
end
if(~bad)
    points=[];
    rot=(atan2((X(3,2)-X(1,2)),(X(3,1)-X(1,1)))+atan2((X(4,1)-X(2,1)),abs((X(4,2)-X(2,2)))))/2;
    rmat=[cos(rot),-sin(rot);sin(rot),cos(rot)];
    %atan2((X(3,2)-X(1,2)),(X(3,1)-X(1,1)))
    %atan2((X(4,1)-X(2,1)),abs((X(4,2)-X(2,2))))
    % Program does not account for rotated plots
    set(handles.rotation,'String',[num2str(rot*180/pi),' degrees'])
    %prog=gcf;
    uiwait(msgbox('To stop selecting points, press any key.','Any key to stop.'));

    axes(handles.imag1);
    hold on;
    [x,y,button]=ginput(1);
    while length(button)==1 && (button==1 | button==2 | button==3)
        x_cur=(x-X(1,1))*(xe-x0)/(X(3,1)-X(1,1))+x0;
        y_cur=(y-X(2,2))*(ye-y0)/(X(4,2)-X(2,2))+y0;
        point=eye(2)*[x_cur;y_cur];
        %point=rmat*[x_cur;y_cur];
        points=[points;point'];
        plot(x,y,'og')
        str=['        X            Y    '];
        for k=1:length(points(:,1))
            s=sprintf('%10.5f %10.5f',points(k,1),points(k,2));
            str=strvcat(str,s);
        end
        set(handles.points,'String',str);
        [x,y,button]=ginput(1);
    end
    hold off;
end

function rotation_Callback(hObject, eventdata, handles)
% hObject    handle to rotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotation as text
%        str2double(get(hObject,'String')) returns contents of rotation as a double


% --- Executes during object creation, after setting all properties.
function rotation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function x0v_Callback(hObject, eventdata, handles)
% hObject    handle to x0v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x0v as text
%        str2double(get(hObject,'String')) returns contents of x0v as a double


% --- Executes during object creation, after setting all properties.
function x0v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x0v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y0v_Callback(hObject, eventdata, handles)
% hObject    handle to y0v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y0v as text
%        str2double(get(hObject,'String')) returns contents of y0v as a double


% --- Executes during object creation, after setting all properties.
function y0v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y0v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xev_Callback(hObject, eventdata, handles)
% hObject    handle to xev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xev as text
%        str2double(get(hObject,'String')) returns contents of xev as a double


% --- Executes during object creation, after setting all properties.
function xev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yev_Callback(hObject, eventdata, handles)
% hObject    handle to yev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yev as text
%        str2double(get(hObject,'String')) returns contents of yev as a double


% --- Executes during object creation, after setting all properties.
function yev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


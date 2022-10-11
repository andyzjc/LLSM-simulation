function varargout = LSSimulator(varargin)
% LSSIMULATOR MATLAB code for LSSimulator.fig
%      LSSIMULATOR, by itself, creates a new LSSIMULATOR or raises the existing
%      singleton*.
%
%      H = LSSIMULATOR returns the handle to a new LSSIMULATOR or the handle to
%      the existing singleton*.
%
%      LSSIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LSSIMULATOR.M with the given input arguments.
%
%      LSSIMULATOR('Property','Value',...) creates a new LSSIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LSSimulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LSSimulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSSimulator

% Last Modified by GUIDE v2.5 10-Oct-2022 12:03:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSSimulator_OpeningFcn, ...
                   'gui_OutputFcn',  @LSSimulator_OutputFcn, ...
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


% --- Executes just before LSSimulator is made visible.
function LSSimulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LSSimulator (see VARARGIN)

% Choose default command line output for LSSimulator
    handles.output = hObject;

% Update handles structure
    guidata(hObject, handles);

    addpath('functions')

% set default to hex
    setDefault(handles);

% Calculate Physical Parameters
    CalculatePhysics(handles);

% UIWAIT makes LSSimulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSSimulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_wavelength_exc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_wavelength_exc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_wavelength_exc as text
%        str2num(get(hObject,'String')) returns contents of edit_wavelength_exc as a double
    wavelength = str2double(handles.edit_wavelength_exc.String);
    if wavelength <= 0 || isnan(wavelength)
        wavelength = 488;
        handles.edit_wavelength_exc.String = num2str(wavelength);
    end
    CalculatePhysics(handles)

% --- Executes during object creation, after setting all properties.
function edit_wavelength_exc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_wavelength_exc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_wavelength_det_Callback(hObject, eventdata, handles)
% hObject    handle to edit_wavelength_det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_wavelength_det as text
%        str2num(get(hObject,'String')) returns contents of edit_wavelength_det as a double
    wavelength = str2double(handles.edit_wavelength_det.String);
    if wavelength <= 0 || isnan(wavelength)
        wavelength = 488;
        handles.edit_wavelength_det.String = num2str(wavelength);
    end
    CalculatePhysics(handles)

% --- Executes during object creation, after setting all properties.
function edit_wavelength_det_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_wavelength_det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_refractive_index_Callback(hObject, eventdata, handles)
% hObject    handle to edit_refractive_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_refractive_index as text
%        str2num(get(hObject,'String')) returns contents of edit_refractive_index as a double
    index = str2double(handles.edit_refractive_index.String);
    if index <= 0 || isnan(index)
        index = 1.33;
        handles.edit_refractive_index.String = num2str(index);
    end
    CalculatePhysics(handles)


% --- Executes during object creation, after setting all properties.
function edit_refractive_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_refractive_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Image_Size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Image_Size_Nxz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Image_Size_Nxz as text
%        str2num(get(hObject,'String')) returns contents of edit_Image_Size_Nxz as a double
    size = str2double(handles.edit_Image_Size.String);
    if mod(size,2) == 0 
        size = size+1;
        handles.edit_Image_Size.String = num2str(size);
    end
    if isnan(size)
        size = 257;
        handles.edit_Image_Size.String = num2str(size);
    end
    CalculatePhysics(handles)


% --- Executes during object creation, after setting all properties.
function edit_Image_Size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Image_Size_Nxz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NAdet_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NAdet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NAdet as text
%        str2num(get(hObject,'String')) returns contents of edit_NAdet as a double
    NAdet = str2double(handles.edit_NAdet.String);
    if NAdet <= 0 || isnan(NAdet)
        NAdet = 1.0;
        handles.edit_NAdet.String = num2str(NAdet);
    end
    CalculatePhysics(handles)


% --- Executes during object creation, after setting all properties.
function edit_NAdet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NAdet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_XZ_sampling_Callback(hObject, eventdata, handles)
% hObject    handle to edit_XZ_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_XZ_sampling as text
%        str2num(get(hObject,'String')) returns contents of edit_XZ_sampling as a double
    factor = str2double(handles.edit_XZ_sampling.String);
    if factor < 2 || isnan(factor)
        factor = 2;
        handles.edit_XZ_sampling.String = num2str(factor);
    end
    CalculatePhysics(handles)

% --- Executes during object creation, after setting all properties.
function edit_XZ_sampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_XZ_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_propagation_step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_propagation_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_propagation_step as text
%        str2num(get(hObject,'String')) returns contents of edit_propagation_step as a double
    factor = str2double(handles.edit_propagation_step.String);
    if factor < 0 || isnan(factor)
        factor = 1;
        handles.edit_propagation_step.String = num2str(factor);
    end
    CalculatePhysics(handles)

% --- Executes during object creation, after setting all properties.
function edit_propagation_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_propagation_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_standingwave.
function radiobutton_standingwave_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_standingwave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_standingwave
global Data
handles.radiobutton_standingwave.Value = 1;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 0;
Data.sti_case = 'sw';
Data.angle = [90,270];
Data.weighting = [0,0];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);

% --- Executes on button press in `Square_lattice.
function radiobutton_Square_lattice_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Square_lattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Square_lattice
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 1;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 0;
Data.sti_case = 'square';
Data.angle = [0,90,180,270];
Data.weighting = [0,0,0,0];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);



% --- Executes on button press in pushbutton_simulate.
function pushbutton_simulate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    ClosePlots;
    CalculatePhysics(handles);
    getParameters(handles);
    Simulate(handles);
    disp("Ploting")
    PrettyPlots;
    if handles.checkbox_save2dir.Value == 1
        SaveResults;
    end


% --- Executes on button press in radiobutton_hex_lattice.
function radiobutton_hex_lattice_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_hex_lattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of radiobutton_hex_lattice
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 1;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 0;
Data.sti_case = 'hex';
Data.angle = [30, 90, 150, 210, 270, 330];
Data.weighting = [0,0,0,0,0,0];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


% --- Executes on button press in radiobutton_2dgaussian.
function radiobutton_2dgaussian_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_2dgaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of radiobutton_2dgaussian
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 1;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 0;
Data.sti_case = '2dgaussian';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


% --- Executes on button press in radiobutton_bessel.
function radiobutton_bessel_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_bessel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of radiobutton_bessel
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 1;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 0;
Data.sti_case = 'bessel';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


% --- Executes on button press in radiobutton_1dgaussian.
function radiobutton_1dgaussian_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_1dgaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_1dgaussian
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 1;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 0;
Data.sti_case = '1dgaussian';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


function edit_NAmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NAmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NAmax as text
%        str2num(get(hObject,'String')) returns contents of edit_NAmax as a double
    NAmax = str2double(handles.edit_NAmax.String);
    if NAmax <= 0 || isnan(NAmax)
        NAmax = 0.6;
        handles.edit_NAmax.String = num2str(NAmax);
    end


% --- Executes during object creation, after setting all properties.
function edit_NAmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NAmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NAmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NAmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NAmin as text
%        str2num(get(hObject,'String')) returns contents of edit_NAmin as a double
    NAmin = str2double(handles.edit_NAmin.String);
    if NAmin < 0 || isnan(NAmin)
        NAmin = 0.5;
        handles.edit_NAmin.String = num2str(NAmin);
    end


% --- Executes during object creation, after setting all properties.
function edit_NAmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NAmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dither_period_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dither_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dither_period as text
%        str2num(get(hObject,'String')) returns contents of edit_dither_period as a double
    Dither_period = str2double(handles.edit_dither_period.String);
    if Dither_period < 0 || isnan(Dither_period)
        Dither_period = 3;
        handles.edit_dither_period.String = num2str(Dither_period);
    end


% --- Executes during object creation, after setting all properties.
function edit_dither_period_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dither_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dither_step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dither_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dither_step as text
%        str2num(get(hObject,'String')) returns contents of edit_dither_step as a double
    Dither_step = str2double(handles.edit_dither_step.String);
    if Dither_step <= 0 || isnan(Dither_step)
        Dither_step = 201;
        handles.edit_dither_step.String = num2str(Dither_step);
    end


% --- Executes during object creation, after setting all properties.
function edit_dither_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dither_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_gaussian_boundwidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gauss_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gauss_bound as text
%        str2num(get(hObject,'String')) returns contents of edit_gauss_bound as a double
    boundwidth = str2double(handles.edit_gaussian_boundwidth.String);
    if boundwidth <= 0 || isnan(boundwidth)
        boundwidth = 3;
        handles.edit_gaussian_boundwidth.String = num2str(boundwidth);
    end

% --- Executes during object creation, after setting all properties.
function edit_gauss_bound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gauss_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_other_lattice.
function radiobutton_other_lattice_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_other_lattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_other_lattice
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 1;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 0;
Data.sti_case = 'other';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


function edit_beam_angle_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beam_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beam_angle as text
%        str2num(get(hObject,'String')) returns contents of edit_beam_angle as a double
    global Data
    types = ['sw','hex','square','other'];
    if contains(types,Data.sti_case)
        angle = str2num(handles.edit_beam_angle.String);
        if isempty(angle) || sum(angle) == 0
            handles.edit_beam_angle.String = [];
            error("Invalid beam angle for lattice.")
        end
    else
        handles.edit_beam_angle.String = [];
    end

% --- Executes during object creation, after setting all properties.
function edit_beam_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beam_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_beam_weighting_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beam_weighting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beam_weighting as text
%        str2num(get(hObject,'String')) returns contents of edit_beam_weighting as a double
    global Data
    types = ['sw','hex','square','other'];
    if contains(types,Data.sti_case)
        angle = str2num(handles.edit_beam_weighting.String);
        if isempty(angle) 
            handles.edit_beam_weighting.String = [];
            error("Invalid beam weighting for lattice.")
        end
    else
        handles.edit_beam_weighting.String = [];
    end


% --- Executes during object creation, after setting all properties.
function edit_beam_weighting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beam_weighting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gauss_circular_NA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gauss_circular_NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gauss_circular_NA as text
%        str2num(get(hObject,'String')) returns contents of edit_gauss_circular_NA as a double
    circularNA = str2double(handles.edit_gauss_circular_NA.String);
    if circularNA <= 0 || isnan(circularNA)
        circularNA = 0.3;
        handles.edit_gauss_circular_NA.String = num2str(circularNA);
    end


% --- Executes during object creation, after setting all properties.
function edit_gauss_circular_NA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gauss_circular_NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_1dairy.
function radiobutton_1dairy_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_1dairy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_1dairy
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 1;
Data.sti_case = '1dairy';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


% --- Executes on button press in radiobutton_2dairy.
function radiobutton_2dairy_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_2dairy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_2dairy
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 1;
handles.radiobutton_1dairy.Value = 0;
Data.sti_case = '2dairy';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


function edit_gaussianwidth_NA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gaussianwidth_NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gaussianwidth_NA as text
%        str2num(get(hObject,'String')) returns contents of edit_gaussianwidth_NA as a double
    gausswidthNA = str2double(handles.edit_gaussianwidth_NA.String);
    if gausswidthNA <= 0 || isnan(gausswidthNA)
        gausswidthNA = 0.3;
        handles.edit_gaussianwidth_NA.String = num2str(gausswidthNA);
    end



% --- Executes during object creation, after setting all properties.
function edit_gaussianwidth_NA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gaussianwidth_NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_airy_phase_Callback(hObject, eventdata, handles)
% hObject    handle to edit_airy_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_airy_phase as text
%        str2num(get(hObject,'String')) returns contents of edit_airy_phase as a double
    airyphase = str2double(handles.edit_airy_phase.String);
    if airyphase < 0 || isnan(airyphase)
        airyphase = pi/2;
        handles.edit_airy_phase.String = "pi/2";
    end


% --- Executes during object creation, after setting all properties.
function edit_airy_phase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_airy_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_select_dir.
function pushbutton_select_dir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global Data

    Data.dir = uigetdir("Select an output directory");
    handles.text_directory.String = Data.dir;
    

% --- Executes during object creation, after setting all properties.
function text_directory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    global Data

    Data.dir = pwd;
    hObject.String = pwd;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global Data
    global Lattice
    ClosePlots;
    handles.text_directory = 'Directory';
    Data = [];
    Lattice = [];
    clear Data;

% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global Data
    ClosePlots;
    disp("Ploting")
    PrettyPlots;


% --- Executes on button press in checkbox_fastdither.
function checkbox_fastdither_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fastdither (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fastdither


% --- Executes on button press in checkbox_realspaceGaussianBound.
function checkbox_realspaceGaussianBound_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_realspaceGaussianBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_realspaceGaussianBound


% --- Executes on button press in checkbox_tophatbeam.
function checkbox_tophatbeam_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tophatbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tophatbeam


% --- Executes on button press in checkbox_gaussianbeam.
function checkbox_gaussianbeam_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_gaussianbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_gaussianbeam


% --- Executes on button press in checkbox_UniformTophatProfile.
function checkbox_UniformTophatProfile_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_UniformTophatProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_UniformTophatProfile


% --- Executes on button press in checkbox_gaussianBeamProfile.
function checkbox_gaussianBeamProfile_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_gaussianBeamProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_gaussianBeamProfile


% --- Executes on button press in checkbox_TophatProfile.
function checkbox_TophatProfile_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_TophatProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_TophatProfile


% --- Executes on button press in checkbox_Gaussian_profile.
function checkbox_Gaussian_profile_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Gaussian_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Gaussian_profile


% --- Executes on button press in radiobutton_SW_Lattice.
function radiobutton_SW_Lattice_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_SW_Lattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_SW_Lattice
global Data
handles.radiobutton_standingwave.Value = 0;
handles.radiobutton_Square_lattice.Value = 0;
handles.radiobutton_hex_lattice.Value = 0;
handles.radiobutton_bessel.Value = 0;
handles.radiobutton_other_lattice.Value = 0;
handles.radiobutton_2dgaussian.Value = 0;
handles.radiobutton_1dgaussian.Value = 0;
handles.radiobutton_2dairy.Value = 0;
handles.radiobutton_1dairy.Value = 0;
handles.radiobutton_SW_Lattice.Value = 1;
Data.sti_case = 'SWLattice';
Data.angle = [90,90,270,270];
Data.weighting = [0,0,0,0];

handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


% --- Executes on button press in radiobutton_SW_Lattice_Hex.
function radiobutton_SW_Lattice_Hex_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_SW_Lattice_Hex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_SW_Lattice_Hex


% --- Executes on button press in radiobutton_SW_Lattice_Square.
function radiobutton_SW_Lattice_Square_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_SW_Lattice_Square (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_SW_Lattice_Square



% --- Executes on button press in radiobutton_SW_Incoherent.
function radiobutton_SW_Incoherent_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_SW_Incoherent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_SW_Incoherent


% --- Executes on button press in radiobutton_SW_Coherent.
function radiobutton_SW_Coherent_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_SW_Coherent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_SW_Coherent

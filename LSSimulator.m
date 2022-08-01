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

% Last Modified by GUIDE v2.5 01-Aug-2022 12:36:19

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

% Calculate Physical Parameters
Calculate_physics(handles);

% set default to hex
handles.radiobutton_hex_lattice.Value = 1;
setDefault(handles)
getParameters(handles)

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
%        str2double(get(hObject,'String')) returns contents of edit_wavelength_exc as a double
    wavelength = str2double(handles.edit_wavelength_exc.String);
    if wavelength <= 0 || isnan(wavelength)
        wavelength = 488;
        handles.edit_wavelength_exc.String = num2str(wavelength);
    end
    Calculate_physics(handles)

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
%        str2double(get(hObject,'String')) returns contents of edit_wavelength_det as a double
    wavelength = str2double(handles.edit_wavelength_det.String);
    if wavelength <= 0 || isnan(wavelength)
        wavelength = 488;
        handles.edit_wavelength_det.String = num2str(wavelength);
    end
    Calculate_physics(handles)

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
%        str2double(get(hObject,'String')) returns contents of edit_refractive_index as a double
    index = str2double(handles.edit_refractive_index.String);
    if index <= 0 || isnan(index)
        index = 1.33;
        handles.edit_refractive_index.String = num2str(index);
    end
    Calculate_physics(handles)


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
% hObject    handle to edit_Image_Size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Image_Size as text
%        str2double(get(hObject,'String')) returns contents of edit_Image_Size as a double
    size = str2double(handles.edit_Image_Size.String);
    if mod(size,2) == 0 
        size = size+1;
        handles.edit_Image_Size.String = num2str(size);
    end
    if isnan(size)
        size = 513;
        handles.edit_Image_Size.String = num2str(size);
    end
    Calculate_physics(handles)


% --- Executes during object creation, after setting all properties.
function edit_Image_Size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Image_Size (see GCBO)
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
%        str2double(get(hObject,'String')) returns contents of edit_NAdet as a double
    NAdet = str2double(handles.edit_NAdet.String);
    if NAdet <= 0 || isnan(NAdet)
        NAdet = 1.0;
        handles.edit_NAdet.String = num2str(NAdet);
    end
    Calculate_physics(handles)


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
%        str2double(get(hObject,'String')) returns contents of edit_XZ_sampling as a double
    factor = str2double(handles.edit_XZ_sampling.String);
    if factor < 2 || isnan(factor)
        factor = 2;
        handles.edit_XZ_sampling.String = num2str(factor);
    end
    Calculate_physics(handles)

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
%        str2double(get(hObject,'String')) returns contents of edit_propagation_step as a double
    factor = str2double(handles.edit_propagation_step.String);
    if factor < 0 || isnan(factor)
        factor = 1;
        handles.edit_propagation_step.String = num2str(factor);
    end
    Calculate_physics(handles)

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
Data.case = 'sw';
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
Data.case = 'square';
Data.angle = [0,90,180,270];
Data.weighting = [0,0,0,0];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);



% --- Executes on button press in pushbutton_simulate.
function pushbutton_simulate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_terminate.
function pushbutton_terminate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_terminate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
Data.case = 'hex';
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
Data.case = '2dgaussian';
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
Data.case = 'bessel';
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
Data.case = '1dgaussian';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


function edit_NAmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NAmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NAmax as text
%        str2double(get(hObject,'String')) returns contents of edit_NAmax as a double
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
%        str2double(get(hObject,'String')) returns contents of edit_NAmin as a double
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
%        str2double(get(hObject,'String')) returns contents of edit_dither_period as a double
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
%        str2double(get(hObject,'String')) returns contents of edit_dither_step as a double
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


function edit_gauss_bound_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gauss_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gauss_bound as text
%        str2double(get(hObject,'String')) returns contents of edit_gauss_bound as a double
    gauss_boundwidth = str2double(handles.edit_gaussian_boundwidth.String);
    if gauss_boundwidth <= 0 || isnan(gauss_boundwidth)
        gauss_boundwidth = 3;
        handles.edit_gaussian_boundwidth.String = num2str(gauss_boundwidth);
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
Data.case = 'other';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


function edit_beam_angle_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beam_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beam_angle as text
%        str2double(get(hObject,'String')) returns contents of edit_beam_angle as a double


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
%        str2double(get(hObject,'String')) returns contents of edit_beam_weighting as a double


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
%        str2double(get(hObject,'String')) returns contents of edit_gauss_circular_NA as a double
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



function edit_slit_na_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slit_na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slit_na as text
%        str2double(get(hObject,'String')) returns contents of edit_slit_na as a double
    slitNA = str2double(handles.edit_slit_na.String);
    if slitNA <= 0 || isnan(slitNA)
        slitNA = 0.1;
        handles.edit_slit_na.String = num2str(slitNA);
    end


% --- Executes during object creation, after setting all properties.
function edit_slit_na_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slit_na (see GCBO)
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
Data.case = '1dairy';
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
Data.case = '2dairy';
Data.angle = [ ];
Data.weighting = [ ];
handles.edit_beam_angle.String = num2str(Data.angle);
handles.edit_beam_weighting.String = num2str(Data.weighting);


function edit_gaussianwidth_NA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gaussianwidth_NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gaussianwidth_NA as text
%        str2double(get(hObject,'String')) returns contents of edit_gaussianwidth_NA as a double
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
%        str2double(get(hObject,'String')) returns contents of edit_airy_phase as a double
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


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


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
    handles.text_directory = 'Directory';

function Calculate_physics(handles)
    global Data

    % Physical Parameter 
    Data.N = str2double(handles.edit_Image_Size.String); % pixels
    Data.xz_scale = str2double(handles.edit_XZ_sampling.String);
    Data.y_scale = str2double(handles.edit_propagation_step.String);
    Data.n = str2double(handles.edit_refractive_index.String);
    Data.lambda_exc = str2double(handles.edit_wavelength_exc.String) / 1000; % um 
    Data.lambda_det = str2double(handles.edit_wavelength_det.String) / 1000;
    Data.wavelength_exc =  Data.lambda_exc / Data.n;
    Data.wavelength_det = Data.lambda_det / Data.n;
    Data.k_wave_exc = 1/Data.wavelength_exc;
    Data.k_wave_det = 1/Data.wavelength_det;
    Data.NAdet = str2double(handles.edit_NAdet.String);
    Data.k_det = Data.k_wave_det * Data.NAdet / Data.n;
    Data.k_bound = Data.xz_scale * Data.k_wave_exc;    
    Data.deltak = 2 * Data.k_bound / Data.N;
    Data.deltax = 1/(2 * Data.k_bound);

    % excitation
    [ax, ~] = meshgrid(  -(Data.N-1)/2 : (Data.N-1)/2 ) ; 
    Data.kx_exc = Data.deltak * ax;  %in unit wavelength
    Data.kz_exc = Data.kx_exc';
    Data.ky_exc = sqrt(Data.k_wave_exc^2 - Data.kx_exc.^2 - Data.kz_exc.^2);
    Data.ky_exc(Data.kx_exc.^2 + Data.kz_exc.^2 > Data.k_wave_exc.^2 ) = 0;
    Data.x_exc = Data.deltax * ax; 
    Data.z_exc = Data.x_exc'; 
    Data.y_exc = (-(Data.N+1)/2+1 : (Data.N+1)/2-1) * Data.deltax * Data.y_scale; 
    
    % detection
    Data.kx_det = Data.deltak * ax;
    Data.ky_det = Data.kx_det';
    Data.kz_det = sqrt(Data.k_wave_det^2 - Data.kx_det.^2 - Data.ky_det.^2);
    Data.x_det = Data.deltax * ax;
    Data.y_det = Data.x_det';
    Data.z_det = (-(Data.N+1)/2+1 : (Data.N+1)/2-1) * Data.deltax * Data.y_scale;
    
    % for displaying
    Data.KX_exc = Data.kx_exc(1,:) / Data.k_wave_exc;
    Data.KZ_exc = Data.KX_exc';
    Data.X_exc = Data.x_exc(1,:)  / Data.wavelength_exc; % value * wavelength = physical value (um)
    Data.Z_exc = Data.X_exc'; 
    Data.Y_exc = Data.y_exc / Data.wavelength_exc; 
    
    Data.KX_det = Data.kx_det(1,:) / Data.k_wave_det;
    Data.KY_det = Data.KX_det';
    Data.X_det = Data.x_det(1,:)  / Data.wavelength_det;
    Data.Y_det = Data.X_det';
    Data.Z_det = Data.y_det / Data.wavelength_det;
    
    % Calculate FOV
    handles.text_step_size.String = num2str(Data.deltax,'%.2f') + " x " + ...
                                    num2str(Data.deltax,'%.2f') + " x " + ...
                                    num2str(Data.deltax * Data.y_scale,'%.2f');

    handles.text_FOV.String = num2str(Data.deltax*Data.N,'%.1f') + " x " + ...
                              num2str(Data.deltax* Data.N,'%.1f') + " x " + ...
                              num2str(Data.deltax * Data.y_scale * (Data.N-1),'%.1f');



function setDefault(handles)
    global Data
    Data.case = 'hex';

    handles.edit_NAmax.String = num2str(0.6);
    handles.edit_NAmin.String = num2str(0.5);
    handles.edit_dither_period.String = num2str(3);
    handles.edit_dither_step.String = num2str(201);
    handles.edit_gaussian_boundwidth.String = num2str(3);
    handles.edit_beam_angle.String = num2str([30, 90, 150, 210, 270, 330]);
    handles.edit_beam_weighting.String = num2str([0,0,0,0,0,0]);
    handles.edit_gauss_circular_NA.String = num2str(0.5);
    handles.edit_airy_phase.String = "pi/2";
    handles.edit_slit_na.String = num2str(0.1);
    handles.edit_gaussianwidth_NA.String = num2str(0.3);


function getParameters(handles)
    global Data
    Data.NAmax = str2double(handles.edit_NAmax);
    Data.k_NAmax = Data.NAmax /Data.n * Data.k_wave_exc; % k
    Data.NAmin = str2double(handles.edit_NAmin);  
    Data.k_NAmin = Data.NAmin /Data.n * Data.k_wave_exc;
    Data.NAideal = (Data.NAmin + Data.NAmax)/2;
    Data.k_ideal = Data.k_wave_exc * Data.NAideal / Data.n;
    Data.dither_period = str2double(handles.edit_dither_period.String); % um
    Data.dither_step = str2double(handles.edit_dither_step.String); % number of s teps per dither period 
    Data.gauss_bound_width = str2double(handles.edit_gaussian_boundwidth.String); % Gaussian Bounding, um
    
    Data.angle = str2double(handles.edit_beam_angle.String);
    Data.weighting =  str2double(handles.edit_beam_weighting.String);
    Data.apertureNA = str2double(handles.edit_gauss_circular_NA.String);
    Data.airyphase = str2double(handles.edit_airy_phase.String);
    Data.slitNA = str2double(handles.edit_slit_na.String);
    Data.gaussian_width_NA = str2double(handles.edit_gaussianwidth_NA.String);
    

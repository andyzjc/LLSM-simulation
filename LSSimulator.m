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

% Last Modified by GUIDE v2.5 30-Jul-2022 13:18:42

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


% --- Executes on button press in radiobutton_Square_lattice.
function radiobutton_Square_lattice_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Square_lattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Square_lattice


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton_1dAiry.
function radiobutton_1dAiry_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_1dAiry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_1dAiry


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes on button press in pushbutton_start_simulation.
function pushbutton_start_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_terminate_simulation.
function pushbutton_terminate_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_terminate_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_save2dir.
function checkbox_save2dir_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_save2dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_save2dir


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


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


% --- Executes on button press in radiobutton_2dGaussian.
function radiobutton_2dGaussian_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_2dGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_2dGaussian


% --- Executes on button press in radiobutton_bessel.
function radiobutton_bessel_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_bessel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_bessel


% --- Executes on button press in radiobutton_1dgaussian.
function radiobutton_1dgaussian_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_1dgaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_1dgaussian


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton9



function edit_NAmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NAmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NAmax as text
%        str2double(get(hObject,'String')) returns contents of edit_NAmax as a double


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


% --- Executes on button press in checkbox_gauss_bound.
function checkbox_gauss_bound_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_gauss_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_gauss_bound



function edit_gauss_bound_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gauss_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gauss_bound as text
%        str2double(get(hObject,'String')) returns contents of edit_gauss_bound as a double


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



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
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


% --- Executes on button press in radiobutton_2d_airy.
function radiobutton_2d_airy_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_2d_airy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_2d_airy



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

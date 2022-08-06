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

% Last Modified by GUIDE v2.5 05-Aug-2022 14:58:31

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
% getParameters(handles)

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
%        str2num(get(hObject,'String')) returns contents of edit_wavelength_det as a double
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
%        str2num(get(hObject,'String')) returns contents of edit_refractive_index as a double
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
%        str2num(get(hObject,'String')) returns contents of edit_Image_Size as a double
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
%        str2num(get(hObject,'String')) returns contents of edit_NAdet as a double
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
%        str2num(get(hObject,'String')) returns contents of edit_XZ_sampling as a double
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
%        str2num(get(hObject,'String')) returns contents of edit_propagation_step as a double
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
    global Data
    ClosePlots
    getParameters(handles);
    Simulate(handles);
    disp("Ploting")
    if Data.Dither == 1
        PrettyPlotsLattice;
    else
        PrettyPlots;
    end

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



function edit_slit_na_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slit_na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slit_na as text
%        str2num(get(hObject,'String')) returns contents of edit_slit_na as a double
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
    if Data.Dither == 1
        PrettyPlotsLattice;
    else
        PrettyPlots;
    end


function Calculate_physics(handles)
    global Data

    % Physical Parameter 
    Data.N = str2num(handles.edit_Image_Size.String); % pixels
    Data.xz_scale = str2num(handles.edit_XZ_sampling.String);
    Data.y_scale = str2num(handles.edit_propagation_step.String);
    Data.n = str2num(handles.edit_refractive_index.String);
    Data.lambda_exc = str2num(handles.edit_wavelength_exc.String) / 1000; % um 
    Data.lambda_det = str2num(handles.edit_wavelength_det.String) / 1000;
    Data.wavelength_exc =  Data.lambda_exc / Data.n;
    Data.wavelength_det = Data.lambda_det / Data.n;
    Data.k_wave_exc = 1/Data.wavelength_exc;
    Data.k_wave_det = 1/Data.wavelength_det;
    Data.NAdet = str2num(handles.edit_NAdet.String);
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
    Data.kz_det = Data.ky_exc;
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
    Data.Z_det = Data.z_det / Data.wavelength_det;
    
    % Calculate FOV
    handles.text_step_size.String = num2str(Data.deltax,'%.2f') + " x " + ...
                                    num2str(Data.deltax,'%.2f') + " x " + ...
                                    num2str(Data.deltax * Data.y_scale,'%.2f');

    handles.text_FOV.String = num2str(Data.deltax*Data.N,'%.1f') + " x " + ...
                              num2str(Data.deltax* Data.N,'%.1f') + " x " + ...
                              num2str(Data.deltax * Data.y_scale * (Data.N-1),'%.1f');


function setDefault(handles)
    global Data
    Data.sti_case = 'hex';

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
    Data.NAmax = str2num(handles.edit_NAmax.String);
    Data.k_NAmax = Data.NAmax /Data.n * Data.k_wave_exc; % k
    Data.NAmin = str2num(handles.edit_NAmin.String);  
    Data.k_NAmin = Data.NAmin /Data.n * Data.k_wave_exc;
    Data.NAideal = (Data.NAmin + Data.NAmax)/2;
    Data.k_ideal = Data.k_wave_exc * Data.NAideal / Data.n;
    Data.dither_period = str2num(handles.edit_dither_period.String); % um
    Data.dither_step = str2num(handles.edit_dither_step.String); % number of s teps per dither period 
    Data.gauss_bound_width = str2num(handles.edit_gaussian_boundwidth.String); % Gaussian Bounding, um
    
    Data.angle = str2num(handles.edit_beam_angle.String);
    Data.weighting =  exp(str2num(handles.edit_beam_weighting.String));
    Data.apertureNA = str2num(handles.edit_gauss_circular_NA.String);
    Data.k_apertureNA = Data.apertureNA /Data.n * Data.k_wave_exc;
    Data.airyphase = exp(str2num(handles.edit_airy_phase.String));
    Data.slitNA = str2num(handles.edit_slit_na.String);
    Data.k_slitNA = Data.slitNA /Data.n * Data.k_wave_exc;
    

function Simulate(handles)
    global Data
    tic
    % Define pupil function
    switch Data.sti_case
        case 'sw'
            CreateLattice;
            Data.Dither = 1;
        case 'square'
            CreateLattice;
            Data.Dither = 1;
        case 'hex'
            CreateLattice;
            Data.Dither = 1;
        case 'bessel'
            Data.Pupil_fun_exc = ((Data.k_NAmax > sqrt(Data.kx_exc.^2 + Data.kz_exc.^2)) ...
                             .* (Data.k_NAmin < sqrt(Data.kx_exc.^2 + Data.kz_exc.^2)));
            Data.Dither = 0;
        case 'other'
            CreateLattice;
            Data.Dither = 1;
        case '2dgaussian'
            Data.Pupil_fun_exc = exp( -(Data.kx_exc.^2 + Data.kz_exc.^2)/ (2 * Data.k_apertureNA) );
            Data.Dither = 0;
        case '1dgaussian'
            Data.Pupil_fun_exc = (Data.kx_exc >= -Data.k_slitNA/2 &...
                             Data.kx_exc <= Data.k_slitNA/2) .* ...
                             exp( -(Data.kx_exc.^2 + Data.kz_exc.^2)/ (2 * Data.k_apertureNA) );
            Data.Dither = 0;
        case '2dairy'
            Data.Pupil_fun_exc = exp(1i * Data.airyphase) .* (Data.k_apertureNA >= (Data.kx_exc.^2 + Data.kz_exc.^2));
            Data.Dither = 0;
        case '1dairy'
            Data.Pupil_fun_exc = exp(1i * Data.airyphase) .* (Data.kx_exc >= -Data.k_slitNA/2 &...
                             Data.kx_exc <= Data.k_slitNA/2) .* ...
                             Data.k_apertureNA >= (Data.kx_exc.^2 + Data.kz_exc.^2) ;
            Data.Dither = 0;
    end

    disp("Case: " + Data.sti_case)

    % Detection
    Data.Pupil_fun_det = Data.k_det > sqrt(Data.kx_det.^2 + Data.ky_det.^2);

    % Simulation start here
    PSF_exc_3d = zeros(Data.N,Data.N,Data.N);
    OTF_exc_3d = zeros(Data.N,Data.N,Data.N);
    OTF_exc_3d_phase = zeros(Data.N,Data.N,Data.N);
    PSF_det_3d = zeros(Data.N,Data.N,Data.N);

    disp("Propagating")
    % propagation
    for i = 1:length(Data.y_exc)
        propagator_exc = exp(2*pi * 1i * Data.ky_exc * Data.y_exc(i));
        PSF_exc_3d(:,:,i) = abs( fftshift( ifft2(Data.Pupil_fun_exc .* propagator_exc) ) ).^2;
        OTF_exc_3d(:,:,i) = abs(fftshift(fft2(PSF_exc_3d(:,:,i))));
        OTF_exc_3d_phase(:,:,i) = angle(OTF_exc_3d(:,:,i));
    end  

    % detection propagation
    for ii = 1:length(Data.z_det)
        propagator_det = exp(2*pi * 1i * Data.kz_det * Data.z_det(ii));
        PSF_det_3d(:,:,ii) = abs( fftshift( ifft2(Data.Pupil_fun_det .* propagator_det) ) ).^2;
    end
    
    % dithering for lattice
    if Data.Dither == 1
        disp("Dithering")
        PSF_exc_3d_dither = zeros(Data.N,Data.N,Data.N);
        OTF_exc_3d_dither = PSF_exc_3d_dither;
        OTF_exc_3d_dither_phase = PSF_exc_3d_dither;

        % dithering along x exc
        for j = 1:Data.dither_step
            PSF_exc_3d_dither = PSF_exc_3d_dither + ...
                circshift(PSF_exc_3d,round(j * Data.dither_period / Data.deltax / Data.dither_step),2);
        end
        
        for k = 1:length(Data.y_exc)
            OTF_exc_3d_dither(:,:,k) = abs(fftshift(fft2(PSF_exc_3d_dither(:,:,k))));
            OTF_exc_3d_dither_phase(:,:,k) = angle( OTF_exc_3d_dither(:,:,k) );
        end

        Data.PSF_exc_3d_dither = PSF_exc_3d_dither/max(max(max(PSF_exc_3d_dither)));
        Data.OTF_exc_3d_dither = OTF_exc_3d_dither/max(max(max(OTF_exc_3d_dither)));
        Data.OTF_exc_3d_dither_phase = OTF_exc_3d_dither_phase;
    end
    
    % Overall 
    Overall_PSF_axial = squeeze(PSF_exc_3d(:,:,(Data.N+1)/2)) .* squeeze(PSF_det_3d(:,(Data.N+1)/2,:))' ; 
    Overall_PSF_lateral = squeeze(PSF_exc_3d(:,(Data.N+1)/2,:)) .* squeeze(PSF_det_3d(:,:,(Data.N+1)/2));
    Overall_OTF_axial = abs(fftshift(fft2(Overall_PSF_axial)));
    Overall_OTF_lateral =  abs(fftshift(fft2(Overall_PSF_lateral)));
    
    % Normalize
    Data.PSF_exc_3d = PSF_exc_3d/max(max(max(PSF_exc_3d)));
    Data.OTF_exc_3d = OTF_exc_3d/max(max(max(OTF_exc_3d)));
    Data.OTF_exc_3d_phase = OTF_exc_3d_phase;
    Data.PSF_det_3d = PSF_det_3d/max(max(max(PSF_det_3d)));
    Data.Overall_PSF_axial = Overall_PSF_axial/max(max(Overall_PSF_axial));
    Data.Overall_PSF_lateral = Overall_PSF_lateral/max(max(Overall_PSF_lateral));
    Data.Overall_OTF_axial = Overall_OTF_axial/max(max(Overall_OTF_axial));
    Data.Overall_OTF_lateral = Overall_OTF_lateral/max(max(Overall_OTF_lateral));  
    toc
    

function CreateLattice
    global Lattice
    global Data
    
    theta = Data.angle;
    weighting = Data.weighting;
    
    Lattice.Illumi_ideal = zeros(size(Data.x_exc));
    kxposition = Data.k_ideal * cosd(theta) /Data.deltak; % pixel
    kzposition = Data.k_ideal * sind(theta) /Data.deltak; % pixel
    
    for j = 1:length(kxposition)
    
        Lattice.Illumi_ideal( ...
            (Data.N+1)/2 + round(kzposition(j)) ,...
            (Data.N+1)/2 + round(kxposition(j)) ) = 1 * weighting(j);
    end
    Lattice.E_ideal = ifft2(Lattice.Illumi_ideal); 
    
    % bounded lattice 
    gauss_bound = exp(-2 * Data.z_exc.^2 / (Data.gauss_bound_width)^2);
    Lattice.E_bound = gauss_bound .* Lattice.E_ideal;
    
    % bounded back pupil
    Lattice.Illum_bound = abs(fft2(fftshift(Lattice.E_bound))).^2;
    Lattice.Illum_bound = Lattice.Illum_bound/max(max(Lattice.Illum_bound));
    
    % Generate mask 
    Lattice.A_mask = ((Data.k_NAmax > sqrt(Data.kx_exc.^2 + Data.kz_exc.^2)) .* (Data.k_NAmin < sqrt(Data.kx_exc.^2 + Data.kz_exc.^2)));
    
    % Pupil functions
    Data.Pupil_fun_exc = Lattice.Illum_bound .* Lattice.A_mask;

function ClosePlots
    global Figures
    if isfield(Figures,'fig1') 
        delete(Figures.fig1)
    end
    if isfield(Figures,'fig2') 
        delete(Figures.fig2)
    end
    if isfield(Figures,'fig3') 
        delete(Figures.fig3)
    end

function PrettyPlotsLattice
    global Data
    global Lattice
    global Figures

    KX_exc = Data.KX_exc;
    KZ_exc = Data.KZ_exc;
    KX_det = Data.KX_det;
    X_exc = Data.X_exc;
    Z_exc = Data.Z_exc;
    Y_exc = Data.Y_exc;
    X_det = Data.X_det;
    Z_det = Data.Z_det;

    N = Data.N;
    n = Data.n;
    
    % Figure 1 - Rear Pupil 
    Figures.fig1 = figure(1);
    Figures.fig1.Name = "XZ-excitation, Y = 0";
    Figures.fig1.WindowState = 'maximized';
    colormap(hot(256))

     subplot(3,4,1);
image11 = imagesc(KX_exc,KZ_exc, real(Lattice.Illumi_ideal) );
    title("Ideal Lattice, " +...
          "N_{xz} = " + num2str(N) + ...
          ", K_{bound} = " + num2str(Data.k_bound) + " (1/um)")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image11.Parent.XLim = [-1,1];
    image11.Parent.YLim = [-1,1];
    colorbar;    

    subplot(3,4,2);
image12 = imagesc(X_exc,Z_exc, abs(Lattice.E_ideal));
    title("Ideal Lattice," + ...
          "\lambda_{exc}/n = " + num2str(Data.lambda_exc, '%.3f') + "um / " + ...
           num2str(n))
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    axis image
    colorbar;

    subplot(3,4,3);
image13 = imagesc(X_exc, Z_exc, abs(Lattice.E_bound));
    title("Bounded Ideal Lattice, " + ...
          "Bound width = " + num2str(Data.gauss_bound_width) + "um")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    axis image
    colorbar;
    
    subplot(3,4,4)
image14 = imagesc( KX_exc, KZ_exc,...
                  Lattice.Illum_bound );
    title("Rear pupil illumination after bounding")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image14.Parent.XLim = [-1,1];
    image14.Parent.YLim = [-1,1];
    colorbar;

Illum_mask = imfuse(Lattice.Illum_bound,Lattice.A_mask,"falsecolor","ColorChannels","green-magenta");
    subplot(3,4,5)
image15 = imagesc( KX_exc, KZ_exc,...
                  Illum_mask);
    title("Masking, " +...
          "NA_{max} = " + num2str(Data.NAmax) +...
          ", NA_{min} = " + num2str(Data.NAmin) )
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image15.Parent.XLim = [-1,1];
    image15.Parent.YLim = [-1,1];

    subplot(3,4,6)
image16 = imagesc( KX_exc, KZ_exc,...
                  Data.Pupil_fun_exc );
    title("Bounded Rear Pupil")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image16.Parent.XLim = [-1,1];
    image16.Parent.YLim = [-1,1];
    colorbar;

    subplot(3,4,7)
image17 = imagesc(X_exc, Z_exc, Data.PSF_exc_3d(:,:,(N+1)/2));
    title("XZ-Excitation PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(3,4,8)
image18 = imagesc( KX_exc,...
                  KZ_exc,...
                  Data.OTF_exc_3d(:,:,(N+1)/2) ) ;
    title("XZ-Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image
    image18.Parent.XLim = [-1,1];
    image18.Parent.YLim = [-1,1];

    subplot(3,4,9)
    hold on;
image19 = imagesc(X_exc, Z_exc ,Data.PSF_exc_3d_dither(:,:,(N+1)/2));
    title("Dithered XZ-Excitation PSF, " + ...
          "T_d = " + num2str(Data.dither_period) + "um, " +...
          "Dither Step = " + num2str(Data.dither_step))
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(3,4,10)
    zPSF = squeeze(Data.PSF_exc_3d_dither(:,(N+1)/2,(N+1)/2))/max(squeeze(Data.PSF_exc_3d_dither(:,(N+1)/2,(N+1)/2)));
image110 = plot( Z_exc, zPSF);
    title("Dithered Z-Excitation PSF")
    ylabel("Normalized a.u. ")
    xlabel("z/\lambda")
    image110.LineWidth = 2;
    image110.Color = 'r';
    image110.Parent.YAxis.TickValues = linspace(0,1,11);
    image110.Parent.XLim = [-15,15];
    grid on

    subplot(3,4,11)
image111 = imagesc( KX_exc,...
                   KZ_exc,...
                   Data.OTF_exc_3d_dither(:,:,(N+1)/2)  );
    title("Dithered XZ-Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image
    image111.Parent.XLim = [-1,1];
    image111.Parent.YLim = [-1,1];

    subplot(3,4,12)
image112 = plot( KZ_exc, squeeze(Data.OTF_exc_3d_dither(:,(N+1)/2,(N+1)/2)));
    title("Dithered XZ-Excitation OTF, " + ...
        "K_x = " + num2str(0) )
    ylabel("Normalized a.u. ")
    xlabel("kz * \lambda")
    image112.LineWidth = 2;
    image112.Color = 'r';
    image112.Parent.XLim = [-2,2];
    image112.Parent.YAxis.TickValues = linspace(0,1,11);
    grid on
    drawnow

% Figure 2 - Excitation
    Figures.fig2 = figure(2);
    Figures.fig2.Name = "Focal PSF/OTF";
    Figures.fig2.WindowState = 'maximized';
    colormap(hot(256))

     subplot(2,3,1);
image21 = imagesc(X_exc,Z_exc,Data.PSF_exc_3d(:,:,(N+1)/2) );
    title("XZ-Excitation PSF, "  + "Y = 0")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(2,3,2);
image22 = imagesc(X_exc,Z_exc,Data.PSF_exc_3d_dither(:,:,(N+1)/2));
     title("Dithered XZ-Excitation PSF, " + ...
          "T_d = " + num2str(Data.dither_period) + "um, " +...
          "Step = " + num2str(Data.dither_step) + ...
          ", Y = 0")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    axis image
    colorbar;
    
    subplot(2,3,3);
    hold on
image23 = plot( KZ_exc, Data.OTF_exc_3d(:,(N+1)/2,(N+1)/2));
phase1 = plot( KZ_exc, Data.OTF_exc_3d_phase(:,(N+1)/2,(N+1)/2));
    title("Z-Excitation-OTF, " + "K_X = 0, " + "K_Y = 0")
    ylabel("Normalized a.u. ")
    xlabel("kz * \lambda")
    image23.Color = 'r';
    image23.LineWidth = 2;
    image23.Parent.XLim = [-2,2];
    phase1.Color = 'g';
    phase1.LineWidth = 2;
    lgd = legend("Amplitude","Phase");
    colorbar;
    grid on
    hold off

    subplot(2,3,4);
image24 = imagesc(Y_exc, Z_exc, squeeze(Data.PSF_exc_3d(:,(N+1)/2,:)));
    title("YZ-Excitation-PSF, " + "X = 0" )
    axis image
    xlabel("y/\lambda")
    ylabel("z/\lambda")
    colorbar

    subplot(2,3,5);
image25 = imagesc(Y_exc,Z_exc, squeeze(Data.PSF_exc_3d_dither(:,(N+1)/2,:)) );
    title("Dithered YZ-Excitation-PSF, " + "X = 0" )
    xlabel("y/\lambda")
    ylabel("z/\lambda")
    axis image
    colorbar;

    subplot(2,3,6);
    yPSF_exc = squeeze(Data.PSF_exc_3d((N+1)/2,(N+1)/2,:));
    % Calculate yFWHM
    index = find(yPSF_exc >= 0.5);
    Data.yFWHM = Y_exc(index(end)) - Y_exc(index(1)); % pixels
image26 = plot(Y_exc, yPSF_exc );
    title("Y-Excitation-PSF, " + "X = 0, Z = 0, " + ...
          "yFWHM = " + num2str(Data.yFWHM) + "\lambda")
    xlabel("y/\lambda")
    ylabel("Normalized a.u. ")
    image26.Color = 'r';
    image26.LineWidth = 2;
    colorbar;
    grid on
    drawnow
    
%% Figure 3 - Overall PSF/OTF
    Figures.fig3 = figure(3);  
    Figures.fig3.Name = "Overall XZ-Axial PSF/OTF, focal plane";
    Figures.fig3.WindowState = 'maximized';
    colormap(hot(256))
        
    subplot(2,4,1)
image31 = imagesc(X_exc,Z_exc,Data.PSF_exc_3d_dither(:,:,(N+1)/2) );
    title("Dithered XZ-Excitation PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(2,4,2)
 image32 = imagesc(X_det,Z_det,squeeze(Data.PSF_det_3d(:,(N+1)/2,:))');
    title("XZ-Detection PSF ")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;  
    image32.Parent.XLim = [-5,5];
    image32.Parent.YLim = [-5,5];

    subplot(2,4,3)
 image33 = imagesc(X_exc,Z_exc,Data.Overall_PSF_axial);
    title("Overall PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;  
    image33.Parent.XLim = [-5,5];
    image33.Parent.YLim = [-5,5];

    subplot(2,4,5)
image35 = imagesc(KX_exc,...
                  KZ_exc,...
                  Data.OTF_exc_3d(:,:,(N+1)/2) );
    title("XZ-Excitation OTF ")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;
    image35.Parent.XLim = [-2,2];
    image35.Parent.YLim = [-2,2];

    subplot(2,4,6)
 image36 = imagesc(KX_det,...
                  KZ_exc,...
                  abs(fftshift(fft2(squeeze(Data.PSF_det_3d(:,(N+1)/2,:))'))));
    title("XZ-Detection OTF ")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;  
    image36.Parent.XLim = [-2,2];
    image36.Parent.YLim = [-2,2];

    subplot(2,4,7)
 image37 = imagesc(X_exc,Z_exc,Data.Overall_OTF_axial);
    title("Overall OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;  

    h1 = subplot(2,4,[4,8]);
    hold on
line_exc = plot(zPSF,Z_exc);
    line_exc.Color = 'g';
    line_exc.LineWidth = 2;
line_det = plot(squeeze(Data.PSF_det_3d((N+1)/2,(N+1)/2,:)),Z_exc);
    line_det.Color = 'b';
    line_det.LineWidth = 2;
line_overall = plot(squeeze(Data.Overall_PSF_axial(:,(N+1)/2)),Z_exc);
    line_overall.Color = 'r';
    line_overall.LineWidth = 2;
% line_lateral = plot(squeeze(PSF_det_3d(:,(N+1)/2,(N+1)/2)),Z_exc);
%     line_lateral.Color = 'k';
%     line_lateral.LineWidth = 2;
    title("Overall Axial-PSF")
    ylabel("z/\lambda")
    xlabel("Normalized a.u. ")
    lgd = legend("Excitation", "Detection","Overall");  
        lgd.FontWeight = 'bold';
        lgd.FontSize = 7;
        lgd.LineWidth = 1;
    h1.YLim = [-6,6];
    h1.YTick = linspace(-6,6,13);
    grid on
    hold off
    drawnow
    
function PrettyPlots
    global Data
    global Figures

    KX_exc = Data.KX_exc;
    KZ_exc = Data.KZ_exc;
    KX_det = Data.KX_det;
    X_exc = Data.X_exc;
    Z_exc = Data.Z_exc;
    Y_exc = Data.Y_exc;
    X_det = Data.X_det;
    Z_det = Data.Z_det;

    N = Data.N;
    n = Data.n;

    Figures.fig1 = figure(1);
    Figures.fig1.Name = "XZ-excitation, Y = 0";
    Figures.fig1.WindowState = 'maximized';
    colormap(hot(256))

    subplot(1,3,1)
image16 = imagesc( KX_exc, KZ_exc,...
                  real(Data.Pupil_fun_exc) );
    title("Rear Pupil")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image16.Parent.XLim = [-1,1];
    image16.Parent.YLim = [-1,1];
    colorbar;

    subplot(1,3,2)
image17 = imagesc(X_exc, Z_exc, Data.PSF_exc_3d(:,:,(N+1)/2));
    title("XZ-Excitation PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(1,3,3)
image18 = imagesc( KX_exc,...
                  KZ_exc,...
                  Data.OTF_exc_3d(:,:,(N+1)/2) ) ;
    title("XZ-Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image
    image18.Parent.XLim = [-1,1];
    image18.Parent.YLim = [-1,1];
    drawnow

% Figure 2 - Excitation
    Figures.fig2 = figure(2);
    Figures.fig2.Name = "Focal PSF/OTF";
    Figures.fig2.WindowState = 'maximized';
    colormap(hot(256))

     subplot(2,2,1);
image21 = imagesc(X_exc,Z_exc,Data.PSF_exc_3d(:,:,(N+1)/2) );
    title("XZ-Excitation PSF, "  + "Y = 0")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;
    
    subplot(2,2,2);
image23 = plot( KZ_exc, Data.OTF_exc_3d(:,(N+1)/2,(N+1)/2));
    title("Z-Excitation-OTF, " + "K_X = 0, " + "K_Y = 0")
    ylabel("Normalized a.u. ")
    xlabel("kz * \lambda")
    image23.Color = 'r';
    image23.LineWidth = 2;
    image23.Parent.XLim = [-2,2];
    colorbar;
    grid on

    subplot(2,2,3);
image24 = imagesc(Y_exc, Z_exc, squeeze(Data.PSF_exc_3d(:,(N+1)/2,:)));
    title("YZ-Excitation-PSF, " + "X = 0" )
    axis image
    xlabel("y/\lambda")
    ylabel("z/\lambda")
    colorbar

    subplot(2,2,4);
    yPSF_exc = squeeze(Data.PSF_exc_3d((N+1)/2,(N+1)/2,:));
    % Calculate yFWHM
    index = find(yPSF_exc >= 0.5);
    Data.yFWHM = Y_exc(index(end)) - Y_exc(index(1)); % pixels
image26 = plot(Y_exc, yPSF_exc );
    title("Y-Excitation-PSF, " + "X = 0, Z = 0, " + ...
          "yFWHM = " + num2str(Data.yFWHM) + "\lambda")
    xlabel("y/\lambda")
    ylabel("Normalized a.u. ")
    image26.Color = 'r';
    image26.LineWidth = 2;
    colorbar;
    grid on
    drawnow
    
%% Figure 3 - Overall PSF/OTF
    Figures.fig3 = figure(3);  
    Figures.fig3.Name = "Overall XZ-Axial PSF/OTF, focal plane";
    Figures.fig3.WindowState = 'maximized';
    colormap(hot(256))
        
    subplot(2,4,1)
image31 = imagesc(X_exc,Z_exc,Data.PSF_exc_3d(:,:,(N+1)/2) );
    title("XZ-Excitation PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(2,4,2)
 image32 = imagesc(X_det,Z_det,squeeze(Data.PSF_det_3d(:,(N+1)/2,:))');
    title("XZ-Detection PSF ")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;  
    image32.Parent.XLim = [-5,5];
    image32.Parent.YLim = [-5,5];

    subplot(2,4,3)
 image33 = imagesc(X_exc,Z_exc,Data.Overall_PSF_axial);
    title("Overall PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;  
    image33.Parent.XLim = [-5,5];
    image33.Parent.YLim = [-5,5];

    subplot(2,4,5)
image35 = imagesc(KX_exc,...
                  KZ_exc,...
                  Data.OTF_exc_3d(:,:,(N+1)/2) );
    title("XZ-Excitation OTF ")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;
    image35.Parent.XLim = [-2,2];
    image35.Parent.YLim = [-2,2];

    subplot(2,4,6)
 image36 = imagesc(KX_det,...
                  KZ_exc,...
                  abs(fftshift(fft2(squeeze(Data.PSF_det_3d(:,(N+1)/2,:))'))));
    title("XZ-Detection OTF ")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;  
    image36.Parent.XLim = [-2,2];
    image36.Parent.YLim = [-2,2];

    subplot(2,4,7)
 image37 = imagesc(X_exc,Z_exc,Data.Overall_OTF_axial);
    title("Overall OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;  

    zPSF = squeeze(Data.PSF_exc_3d(:,(N+1)/2,(N+1)/2))/max(squeeze(Data.PSF_exc_3d(:,(N+1)/2,(N+1)/2)));
    h1 = subplot(2,4,[4,8]);
    hold on
line_exc = plot(zPSF,Z_exc);
    line_exc.Color = 'g';
    line_exc.LineWidth = 2;
line_det = plot(squeeze(Data.PSF_det_3d((N+1)/2,(N+1)/2,:)),Z_exc);
    line_det.Color = 'b';
    line_det.LineWidth = 2;
line_overall = plot(squeeze(Data.Overall_PSF_axial(:,(N+1)/2)),Z_exc);
    line_overall.Color = 'r';
    line_overall.LineWidth = 2;
% line_lateral = plot(squeeze(PSF_det_3d(:,(N+1)/2,(N+1)/2)),Z_exc);
%     line_lateral.Color = 'k';
%     line_lateral.LineWidth = 2;
    title("Overall Axial-PSF")
    ylabel("z/\lambda")
    xlabel("Normalized a.u. ")
    lgd = legend("Excitation", "Detection","Overall");  
        lgd.FontWeight = 'bold';
        lgd.FontSize = 7;
        lgd.LineWidth = 1;
    h1.YLim = [-6,6];
    h1.YTick = linspace(-6,6,13);
    grid on
    hold off
    drawnow

function SaveResults
    global Data
    global Figures

    disp("Saving")
    
    if Data.Dither == 1
         foldername = datestr(now,'mm-dd-yyyy') + "-" +Data.sti_case + "-NAmax=" + num2str(Data.NAmax) + "-NAmin=" +...
                      num2str(Data.NAmin) + "-boundWidth=" + num2str(Data.gauss_bound_width) + "-n=" + num2str(Data.n) ...
                      + "-Exclambda=" + num2str(Data.lambda_exc) +  "-Prop=" + num2str(2*max(Data.Y_exc));
    else
        foldername = datestr(now,'mm-dd-yyyy') + "-" +Data.sti_case + "-ApertureNA=" + num2str(Data.apertureNA) +...
                     "-slitNA=" + num2str(Data.slitNA) + "-n=" + num2str(Data.n) ...
                      + "-Exclambda=" + num2str(Data.lambda_exc) + "-Prop=" + num2str(2*max(Data.y_exc)) + "um";
    end

    savingdir = Data.dir + "/" + foldername;
    if ~exist(savingdir, 'dir')
       mkdir(savingdir)
    end
    
    save(savingdir + '/results.mat','-struct','Data','-v7.3');
    savefig(Figures.fig1, savingdir + '/XZ-Excitation.fig')
    savefig(Figures.fig2, savingdir + '/PSF-OTF.fig')
    savefig(Figures.fig3,savingdir + '/Overall.fig')

    saveas(Figures.fig1, savingdir + '/XZ-Excitation.png')
    saveas(Figures.fig2, savingdir + '/PSF-OTF.png')
    saveas(Figures.fig3, savingdir + '/Overall.png')
    
    tic

    for i = 1:size(Data.Y_exc,2)
        Figures.fig4 = figure("Visible","off",'WindowState','maximized');
        Figures.fig4.Name = "PSF/OTF";
        h1 = subplot(2,2,1);
        imagesc(Data.X_exc,Data.Z_exc,squeeze(Data.PSF_exc_3d(:,:,i)));
        title("XZ-Excitation PSF - Y=" + num2str(Data.Y_det(i)) + "\lambda")
        xlabel("x/\lambda")
        ylabel("z/\lambda")
        h1.XAxis.FontSize = 15;
        h1.YAxis.FontSize = 15;
        h1.XAxis.FontWeight = 'bold';
        h1.YAxis.FontWeight = 'bold';
        colormap(hot(256))
        colorbar;
        axis image;

        h2 = subplot(2,2,2);
        imagesc(Data.KX_exc, Data.KZ_exc, squeeze(Data.OTF_exc_3d(:,:,i)));
        title("XZ-Excitation OTF - Y=" + num2str(Data.Y_det(i)) + "\lambda")
        xlabel("kx * \lambda")
        ylabel("kz * \lambda")
        h2.XAxis.FontSize = 15; 
        h2.YAxis.FontSize = 15; 
        h2.XAxis.FontWeight = 'bold'; 
        h2.YAxis.FontWeight = 'bold'; 
        colormap(hot(256))
        colorbar; 
        axis image; 

        h3 = subplot(2,2,3);
        if Data.Dither == 1
            line = plot(Data.Z_exc, squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,i)));
            title("Dithered Z-Excitation PSF - Y=" + num2str(Data.Y_det(i)) + "\lambda")
        else
            line = plot(Data.Z_exc, squeeze(Data.PSF_exc_3d(:,(Data.N+1)/2,i)));
            title("Z-Excitation PSF - Y=" + num2str(Data.Y_det(i)) + "\lambda")
        end
        xlabel("z/\lambda")
        ylabel("Normalized a.u. ")
        h3.XAxis.FontSize = 15; 
        h3.YAxis.FontSize = 15; 
        h3.XAxis.FontWeight = 'bold'; 
        h3.YAxis.FontWeight = 'bold'; 
        line.Color = 'r';
        line.LineWidth = 2;
        colorbar; 
        grid on;

        h4 = subplot(2,2,4);
        hold on
        if Data.Dither == 1
            amp = plot(Data.KZ_exc, squeeze(Data.OTF_exc_3d_dither(:,(Data.N+1)/2,i)));
            phase = plot(Data.KZ_exc, squeeze(Data.OTF_exc_3d_dither_phase(:,(Data.N+1)/2,i)));
            title("Dithered Z-Excitation OTF - Y=" + num2str(Data.Y_det(i)) + "\lambda")
        else
            amp = plot(Data.KZ_exc, squeeze(Data.OTF_exc_3d(:,(Data.N+1)/2,i)));
            phase = plot(Data.KZ_exc, squeeze(Data.OTF_exc_3d_phase(:,(Data.N+1)/2,i)));
            title("Z-Excitation OTF - Y=" + num2str(Data.Y_det(i)) + "\lambda")
        end
        xlabel("kz * \lambda")
        ylabel("Normalized a.u. ")
        h4.XAxis.FontSize = 15; 
        h4.YAxis.FontSize = 15; 
        h4.XAxis.FontWeight = 'bold'; 
        h4.YAxis.FontWeight = 'bold'; 
        amp.Color = 'r';
        amp.LineWidth = 2;
        phase.Color = 'g';
        phase.LineWidth = 2;
        colorbar; 
        hold off
        grid on;

        frame = getframe(Figures.fig4);
        imwrite(frame2im(frame),hot, savingdir + '/XZ-PSF-OTF.tif','WriteMode','append')
    end
    toc

    disp("Saved to " + savingdir)
    


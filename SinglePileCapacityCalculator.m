function varargout = SinglePileCapacityCalculator(varargin)
% SinglePileCapacityCalculator MATLAB code for SinglePileCapacityCalculator.fig
%      SinglePileCapacityCalculator, by itself, creates a new SinglePileCapacityCalculator or raises the existing
%      singleton*.
%
%      H = SinglePileCapacityCalculator returns the handle to a new SinglePileCapacityCalculator or the handle to
%      the existing singleton*.
%
%      SinglePileCapacityCalculator('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SinglePileCapacityCalculator.M with the given input arguments.
%
%      SinglePileCapacityCalculator('Property','Value',...) creates a new SinglePileCapacityCalculator or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SinglePileCapacityCalculator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SinglePileCapacityCalculator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SinglePileCapacityCalculator

% Last Modified by GUIDE v2.5 11-Jun-2021 14:17:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SinglePileCapacityCalculator_OpeningFcn, ...
                   'gui_OutputFcn',  @SinglePileCapacityCalculator_OutputFcn, ...
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


% --- Executes just before SinglePileCapacityCalculator is made visible.
function SinglePileCapacityCalculator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SinglePileCapacityCalculator (see VARARGIN)
set(handles.EDT_Qult,'enable','off');
set(handles.EDT_Qs,'enable','off');
set(handles.EDT_Qp,'enable','off');
set(handles.EDT_now,'enable','off');
set(handles.EDT_pls,'enable','off');
% Choose default command line output for SinglePileCapacityCalculator
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SinglePileCapacityCalculator wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EDT_layer_number_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_layer_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_layer_number as text
%        str2double(get(hObject,'String')) returns contents of EDT_layer_number as a double
%layer number inserted as a double to handles.layer_number
handles.layer_number=str2double(get(hObject,'String'))  ;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function EDT_layer_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_layer_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PUSH_layer_number.
function PUSH_layer_number_Callback(hObject, eventdata, handles)
% hObject    handle to PUSH_layer_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%a cell of row names are created and a cell of layer number to 8 are
%created
for i=1:handles.layer_number;
    layer_number{i}=['Layer ' num2str(i)];
    empty{i,8}=[];
end
%layer names are assigned to the row names of the table
set(handles.TABLE_input,'RowName',layer_number);
%to show a path to the user, a message is printed
set(handles.EDT_now,'string','Now, using the table on the right');
set(handles.EDT_pls,'string','please enter your data for soil profile');
%empty cells are introduced to the table
set(handles.TABLE_input,'Data',empty);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PUSH_draw.
function PUSH_draw_Callback(hObject, eventdata, handles)
% hObject    handle to PUSH_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%data is inserted into the table
handles.tableData=get(handles.TABLE_input,'Data');
%according to the data entered table, handles are created for each column
handles.soil_type=handles.tableData(:,1);
handles.start_depth=handles.tableData(:,2);
handles.end_depth=handles.tableData(:,3);
handles.unit_weight=handles.tableData(:,4);
handles.Cu=handles.tableData(:,5);
handles.Adhesion=handles.tableData(:,6);
handles.friction=handles.tableData(:,7);
handles.Nq=handles.tableData(:,8);
%critical diameter is computed
Dcr=15*handles.diameter_pile;
%last end depth is assigned to length of pile for easy computation
%this value will be corrected before drawing the graph
handles.end_depth{end}=handles.length_pile;

%a for loop is created to compute skin friction and bearing for each layer
for j=1:handles.layer_number;
    
    %eff stress of all layers are calculated by the function below:
    eff_stress=eff_stress_calculation([handles.end_depth{:}],[handles.unit_weight{:}],handles.diameter_pile);
        
    %if critical depth occurs in the first layer, this condition is
    %computed
        if (j==1)&(handles.end_depth{j}<Dcr);
            eff_stress_avg=(eff_stress(j))/2;
     
      %if a depth is smaller or equal to critical depth, direct average is
      %taken
        elseif handles.end_depth{j}<=Dcr
               eff_stress_avg=(eff_stress(j)+eff_stress(j-1))/2;
               
      %if a depth is larger than critical depth
        elseif handles.end_depth{j}>Dcr
            
            %if the previous layer was smaller than dcr, this means Dcr is
            %first occuring here
            %this is inserted to account for cases where Dcr did not
            %occured in the last layer
            if handles.end_depth{j-1}<=Dcr
                %area based average is taken here to have more accurate
                %results
                eff_stress_avg=((eff_stress(j-1)+eff_stress(j))*(Dcr-handles.end_depth{j-1})/2+(eff_stress(j)*(handles.end_depth{j}-Dcr)))/(handles.end_depth{j}-handles.end_depth{j-1});
            %else case eff stress remains same
            else
                eff_stress_avg=eff_stress(j+1);
            end         
    
        end
        %for cohesive soils skin friction is computed as a cell for each
        %layer
    if handles.soil_type{j}=="Cohesive"
        handles.Qs{j}=pi*(handles.end_depth{j}-handles.start_depth{j})*handles.Adhesion{j}*handles.Cu{j}*handles.diameter_pile;
        %for cohesionless soils skin friction is computed as a cell for
        %each layer
    elseif handles.soil_type{j}=="Cohesionless"    
     handles.Qs{j}=pi*handles.diameter_pile*0.5*eff_stress_avg*tand(0.75*handles.friction{j})*(handles.end_depth{j}-handles.start_depth{j})
    end
end
%if the last layer is cohesive, tip is computed accordingly
    if handles.soil_type{handles.layer_number}=="Cohesive"
        handles.Qp=9*handles.Cu{handles.layer_number}*pi*((handles.diameter_pile)^2)/4;
      %if the last layer is cohesionless, tip is computed accordingly.
    elseif handles.soil_type{j}=="Cohesionless"
        eff_stress_tip=eff_stress(end);
        handles.Qp=handles.Nq{handles.layer_number}*eff_stress_tip*pi*((handles.diameter_pile)^2)/4;
    end
    %elements in the handles.Qs are summed to find total skin friction
handles.Qs_last=sum([handles.Qs{:}]);100    
%skin and tip are summed to have total capacity
handles.Qtotal=handles.Qs_last+handles.Qp;
%end depth is assigned as last value of the table again
handles.end_depth{end}=handles.tableData{handles.layer_number,3};
%elements of handles.Qs cell are deleted so that they wont remain in after
%computations
handles.Qs(:)=[];

%D values are created so that Dcr will be introduced as an important depth
for i=1:handles.layer_number
    %layers with end depth smaller than and equal to Dcr goes into this loop
    if handles.end_depth{i}<=Dcr
        D(i)=handles.end_depth{i};
    %layers with end depths bigger than Dcr goes into this loop
    elseif handles.end_depth{i}>Dcr
        %in cases where Dcr exists in a depth before than the current layer, this if
        %condition is computed so that code will not repeat itself and add
        %an extra Dcr
        if handles.end_depth{i-1}>Dcr
            %if Dcr is inserted beforehand to the array, just end depth of
            %this array is inserted to the i+1'th row.
            D(i+1)=handles.end_depth{i};
        else    
            if handles.end_depth{i-1}<Dcr
                %if Dcr is faced first time, this loop creates a row for Dcr in
                %D array.
                D(i+1)=handles.end_depth{i};
                D(i)=Dcr;
            else
                %if the previous layer had an end_depth exactly same with
                %Dcr, this loop is computed
            D(i)=handles.end_depth{i};
            end
         
        end
    end
end

%This will give non-zero unique values of D
[r,c,v]=find(D); 
%example situtation, if Dcr exists in layer 2, above loop will succesfully
%work until layer 3. However, if there is a layer 4, above loop will
%introduce a 0 value in the D equation. To overcome this issue, non-zero
%values of D is equated to Depth array itself.
D=v;


%effective stres are computed again to account for the end of last layer
%which might not be the end of pile
%first computation was done with the last point end of the pile length to
%compute average stresses properly.
eff_stress=eff_stress_calculation([handles.end_depth{:}],[handles.unit_weight{:}],handles.diameter_pile);



%by the rows below, required figure were drawn
%Here I inserted the initial values of depth and effective stress.
%Since I computed the effective stress function according to end depths,
%they were not inserted at the beginning
for i=(length(D)+1):-1:2;
    D(i)=D(i-1);
    eff_stress(i)=eff_stress(i-1);
end
D(1)=handles.start_depth{1};
eff_stress(1)=handles.start_depth{1}*(handles.unit_weight{1}-9.81);
%D values multiplied by minus to draw in a descending manner
D=-D;
%eff_stress is scaled by 1/20
eff_stress_scaled=eff_stress./20;
%eff stress vs depth is drawn
plot(eff_stress_scaled,D,'LineWidth',1,'Color',[0 0 0]);
%limitations on the x axis were put considering the largest
%eff_stress_scaled value
xlim([-5 eff_stress_scaled(end)+5]);
%pile was drawn at the position of -3
rectangle( 'Position' , [-3,-handles.length_pile,handles.diameter_pile,handles.length_pile] , 'FaceColor' , [0.859 0.855 0.843] );

%figure is divided into layers
for i=1:handles.layer_number
    line([-5 eff_stress_scaled(end)+5], [-handles.end_depth{i} -handles.end_depth{i}],'LineWidth',1,'Color',[0 0 0]);
    LayerName = ['Layer ', num2str(i)] ;
    text(eff_stress_scaled(end)+2,-(handles.end_depth{i}+handles.start_depth{i})/2,LayerName);          
end

%inserting a reference line to the effective stress
line([0 0],[-handles.start_depth{1} -handles.end_depth{handles.layer_number}],'LineWidth',2,'Color',[0 0 0]);
%inserting a dashed line for the end of pile
line([-5 eff_stress_scaled(end)],[-handles.length_pile -handles.length_pile],'LineWidth',0.5,'Color',[0 0 0],'LineStyle','--');
%inserting values of effective stresses
text(0.2,handles.start_depth{1}-0.1,num2str(eff_stress(1)));  
text(-0.5,handles.start_depth{1}+0.5,'Effective Stress(kPa)');  
for i=2:length(D);
    text(eff_stress_scaled(i)+0.2,D(i)+0.6,num2str(eff_stress(i))); 
    if D(i)==-Dcr;
     text_dcr=['(Dcr=' num2str(Dcr) 'm)'];
     text(eff_stress_scaled(i)+eff_stress_scaled(i)/3,D(i)+0.6,text_dcr);  
     %eff_stress_scacled(i)/4 is added to the location of eff_stress of Dcr
     %so that a nice space will be left between value and this will be
     %changed by the number largeness.
end
end
% space left at the bottom of the graph
D(end+1)=D(end)-2;
%axes set by the required values
rev_D=fliplr(D);
set(gca , 'ytick', rev_D);
set(gca,'xtick',[]);

%setting values as RESULTS
set(handles.EDT_Qs,'string',num2str(handles.Qs_last));
set(handles.EDT_Qp,'string',num2str(handles.Qp));
set(handles.EDT_Qult,'string',num2str(handles.Qtotal));

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in PUSH_export_pdf.
function PUSH_export_pdf_Callback(hObject, eventdata, handles)
% hObject    handle to PUSH_export_pdf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%below functions are used if whole figure is desired to be exported to the
%pdf
%orient(gcf,'landscape');
%print(gcf, '-dpdf', 'test.pdf');

%creating global variables to print in another script
global TableResult
global diameter
global Qs
global length
global Qp
global Qtotal

%for the first header of the table property is assigned to an string array
tableheader(1)="Property";
%Then row names( layer names) are assigned to the remaning headers
tableheader(2:handles.layer_number+1)=handles.TABLE_input.RowName;

%a for loop is conducted to take table data to cell structure result and
%then it converted into a string in the inner for loop
for i=2:handles.layer_number+1
    result{1,1}(i,:)=handles.tableData(i-1,:);
    for j=2:8
        result{1,1}{i,j}=num2str(result{1,1}{i,j});
    end
end

%transpose of the column names is assigned to result cell
result{1,1}(1,:)=transpose(handles.TABLE_input.ColumnName);
%some of the column names are shorted so that 3 layer and property
%information would fit in the pdf nicely.
result{1,1}{1,1}='Type';
result{1,1}{1,2}='Start';
result{1,1}{1,3}='End';
result{1,1}{1,6}='AdhFact';
result{1,1}{1,7}='phiAngle';

%variables are taken from handles.
diameter=handles.diameter_pile;
length=handles.length_pile;
Qs=handles.Qs_last;
Qp=handles.Qp;
Qtotal=handles.Qtotal;

%transpose of the result cell is converted to a table
TableResult=cell2table(transpose(result{1,1}));
%table was assigned with headers stated
TableResult.Properties.VariableNames=tableheader;

%to exclude code information this row is stated
export.showCode = false;
%export format was decided as pdf
export.format = 'pdf';
%specific script was published
publish('PileCapacityReport.m',export)



function EDT_Qs_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_Qs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_Qs as text
%        str2double(get(hObject,'String')) returns contents of EDT_Qs as a double


% --- Executes during object creation, after setting all properties.
function EDT_Qs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_Qs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDT_Qp_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_Qp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_Qp as text
%        str2double(get(hObject,'String')) returns contents of EDT_Qp as a double


% --- Executes during object creation, after setting all properties.
function EDT_Qp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_Qp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDT_Qult_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_Qult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: get(hObject,'String') returns contents of EDT_Qult as text
%        str2double(get(hObject,'String')) returns contents of EDT_Qult as a double


% --- Executes during object creation, after setting all properties.
function EDT_Qult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_Qult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDT_length_pile_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_length_pile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_length_pile as text
%        str2double(get(hObject,'String')) returns contents of EDT_length_pile as a double
handles.length_pile=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function EDT_length_pile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_length_pile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDT_diameter_pile_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_diameter_pile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_diameter_pile as text
%        str2double(get(hObject,'String')) returns contents of EDT_diameter_pile as a double
handles.diameter_pile=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function EDT_diameter_pile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_diameter_pile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Outputs from this function are returned to the command line.
function varargout = SinglePileCapacityCalculator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EDT_now_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_now (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_now as text
%        str2double(get(hObject,'String')) returns contents of EDT_now as a double


% --- Executes during object creation, after setting all properties.
function EDT_now_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_now (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDT_pls_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_pls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_pls as text
%        str2double(get(hObject,'String')) returns contents of EDT_pls as a double


% --- Executes during object creation, after setting all properties.
function EDT_pls_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_pls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% ������ %%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = test1(varargin)
% TEST1 MATLAB code for test1.fig
%      TEST1, by itself, creates a new TEST1 or raises the existing
%      singleton*.
%
%      H = TEST1 returns the handle to a new TEST1 or the handle to
%      the existing singleton*.
%
%      TEST1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST1.M with the given input arguments.
%
%      TEST1('Property','Value',...) creates a new TEST1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test1

% Last Modified by GUIDE v2.5 28-Jun-2023 22:20:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test1_OpeningFcn, ...
                   'gui_OutputFcn',  @test1_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% openingfucntion %%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before test1 is made visible.
function test1_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test1 (see VARARGIN)

% Choose default command line output for test1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes test1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%% output function %%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = test1_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%% �༭�ı�1 %%%%%%%%%%%%%%%%%%%%%%%%%%
function edit1_Callback(~, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
% 
% data=get(handles.edit1,'string');
% data=strsplit(data,' ');
% data=str2double(data);
% handles.Yraw_new= [(1:308)' data'];
% axes(handles.axes1);
% plot(handles.Yraw_new(:,2),'k','linewidth',1,'Color',[0.7 0.7 0.7]);
% guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% �༭�ı�2 %%%%%%%%%%%%%%%%%%%%%%%%%%
function edit2_Callback(~, ~, ~)

function edit2_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(~, ~, ~)

function edit3_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(~, ~, ~)

function edit4_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu1_Callback(hObject, ~, handles)

val = get(handles.popupmenu1,'value');
list = get(handles.popupmenu1,'string');
[currentFolder, ~, ~] = fileparts(mfilename('fullpath'));
[parent_path, ~, ~] = fileparts(currentFolder);
handles.atlas_img=[parent_path,'/rawdata/brain_atlas/',list{val}];%% recepter matrix of DK308
handles.aal = list{val};
guidata(hObject,handles);

function popupmenu2_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu2_Callback(hObject, ~, handles)

val = get(handles.popupmenu2,'value');
list = get(handles.popupmenu2,'string');
handles.brain_sphere=list{val};%% recepter matrix of DK308
guidata(hObject,handles);

function popupmenu3_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu3_Callback(hObject, ~, handles)

val = get(handles.popupmenu3,'value');
list = get(handles.popupmenu3,'string');
handles.interpolation=list{val};%% recepter matrix of DK308
guidata(hObject,handles);

% --- Executes during object deletion, before destroying properties.
function edit1_DeleteFcn(~, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function edit2_DeleteFcn(~, ~, ~)

function edit3_DeleteFcn(~, ~, ~)

function edit4_DeleteFcn(~, ~, ~)

%% %%%%%%%%%%%%%%%%%%%%%%%%% ��ť1 %%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(~, ~, ~)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Filename, Pathname]=uigetfile('*.txt;*.xlsx','Select Input Ground Motion'); % Raw download PEER data
file=[Pathname,Filename];data=importdata(file);handles.Yraw_new= data;handles.trait= data(:,2:size(data,2));
axes(handles.axes1);plot(handles.trait(:,1),'k','linewidth',1,'Color',[0.7 0.7 0.7]);
handles.dim=str2double(get(handles.edit10,'string')); handles.rep=str2double(get(handles.edit11,'string')); handles.boot=str2double(get(handles.edit12,'string'));

guidata(hObject,handles)

%% %%%%%%%%%%%%%%%%%%%%%%%%% ��ť2 %%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton2_CreateFcn(~, ~, ~)

function pushbutton2_Callback(hObject, eventdata, handles)

folder_name = uigetdir;set(handles.edit4,'String',folder_name);
handles.out_path = [folder_name,'/',datestr(datetime('now'),'yyyy_mm_dd_HH_MM'),'/'];mkdir(handles.out_path);
parpool;
guidata(hObject,handles)
 
function pushbutton3_CreateFcn(~, ~, ~)

function pushbutton3_Callback(hObject, eventdata, handles)

disp('Loading genes ...');
% gene--trait
handles.gene_matrix = Method_Gene_matrix(handles.atlas_img,handles.brain_sphere,handles.interpolation,handles.out_path);

[~,ia,ib] = intersect(handles.gene_matrix.ROIs,handles.Yraw_new(:,1));% "deblank" to remove the blank of the string-end.
X = handles.gene_matrix.expr(ia,:);% X = zscore(gene_matrix.expr(ia,:));
Y = zscore(handles.trait(ib,:));% Y = zscore([Yraw_new(ib,2),protein_matrix.expr(ib,:)]);
fig_name = 'gene';PLSperm_gene = Method_PLSpermutation(X,Y,handles.dim,handles.rep,handles.out_path,fig_name);
temp=[2 1];max_dim=temp(((PLSperm_gene.PCVAR(1) > PLSperm_gene.PCVAR(2)) & (PLSperm_gene.Pvalue_perm(1)<0.05))+1);
handles.PLSboot_gene = Method_PLSbootstrap(X,Y,handles.dim,max_dim,handles.boot,handles.gene_matrix.symbols,handles.out_path,fig_name);
% gene--Protein
handles.protein_matrix = Method_Protein_matrix(handles.atlas_img,handles.brain_sphere,handles.out_path);

[~,ia,ib] = intersect(handles.gene_matrix.ROIs,handles.protein_matrix.ROIs);% "deblank" to remove the blank of the string-end.
PLSboot_gp =[];X = handles.gene_matrix.expr(ia,:);% X = zscore(gene_matrix.expr(ia,:));
for idx=1:size(handles.protein_matrix.symbols,1)
    Y = zscore(handles.protein_matrix.expr(ib,idx));% Y = zscore([Yraw_new(ib,2),protein_matrix.expr(ib,:)]);
    fig_name = [handles.protein_matrix.symbols{idx}];PLSperm_gene = Method_PLSpermutation(X,Y,handles.dim,handles.rep,handles.out_path,fig_name);
    temp=[2 1];max_dim=temp(((PLSperm_gene.PCVAR(1) > PLSperm_gene.PCVAR(2)) & (PLSperm_gene.Pvalue_perm(1)<0.05))+1);
    tempboot = Method_PLSbootstrap(X,Y,handles.dim,max_dim,handles.boot,handles.gene_matrix.symbols,handles.out_path,fig_name);
    PLSboot_gp =horzcat(PLSboot_gp,tempboot.Z);
    disp(['Loading PETs: ',num2str(idx),'/',num2str(size(handles.protein_matrix.symbols,1))]);
end
handles.PLSboot_gp_Z = PLSboot_gp;
% Protein--trait
[~,ia,ib] = intersect(handles.protein_matrix.ROIs,handles.Yraw_new(:,1));% "deblank" to remove the blank of the string-end.
X = zscore(handles.protein_matrix.expr(ia,:));Y = zscore(handles.trait(ib,:));
fig_name = 'protein';% SNP, gene, protein
PLSperm_protein = Method_PLSpermutation(X,Y,handles.dim,handles.rep,handles.out_path,fig_name);
temp=[2 1];max_dim=temp(((PLSperm_protein.PCVAR(1) > PLSperm_protein.PCVAR(2)) & (PLSperm_protein.Pvalue_perm(1)<0.05))+1);
handles.PLSboot_protein = Method_PLSbootstrap(X,Y,handles.dim,max_dim,handles.boot,handles.protein_matrix.symbols,handles.out_path,fig_name);
% save expression Matrix and PLS-Zscore
gene_matrix=handles.gene_matrix; protein_matrix=handles.protein_matrix;
PLSboot_gene=handles.PLSboot_gene;PLSboot_gp_Z=handles.PLSboot_gp_Z;PLSboot_protein=handles.PLSboot_protein;
temp_file = [handles.out_path,'PLS_of_gpt_',handles.aal,'_',handles.brain_sphere,'_',handles.interpolation,'_ROIs.mat'];
save(temp_file,'gene_matrix','protein_matrix','PLSboot_gene','PLSboot_gp_Z','PLSboot_protein');
disp(['Loaded',num2str(size(handles.gene_matrix.symbols,1)),' genes and ',num2str(size(handles.protein_matrix.symbols,1)),'PET maps!']);
guidata(hObject,handles)

function pushbutton4_Callback(hObject, eventdata, handles)

disp('Loading PPI ...');current_file = mfilename('fullpath');   % ��ȡ��ǰ���д���ĵ�ַ
[current_path, ~, ~] = fileparts(current_file);[parent_path, ~, ~] = fileparts(current_path);  %��ȡ��һ���ַ
PPI = readtable([parent_path,'/rawdata/9606.protein.links.v11.5_new.txt'],'Delimiter', '\t');
proteins_40 = {'HTR1A','HTR1A','HTR1B','HTR1B','HTR2A','HTR2A','HTR4','HTR6','SLC6A4','CNR1','CNR1','DRD1',...
    'DRD2','DRD2','DRD2','DRD2','DRD2','SLC6A3','SLC6A3','GABRA1','GABRA1','HTR1B','CHRM1','OPRM1','OPRM1',...
    'SLC6A2','SLC6A2','SLC6A4','SLC6A4','SLC6A4','SV2A','VAT1L','VAT1L','VAT1L','VAT1L','GRM5','GRM5','GRM5','GRM5','CHRNB4'}';
if ~exist('handles.PLSboot_gp', 'var')
    load('E:\Project\XomicsEnrich\enrichment_out\PLS_of_gpt_left34_ROIs.mat');% PLSboot_gene40_map2whole308,geneSymbol
    handles.gene_matrix=gene_matrix; handles.protein_matrix=protein_matrix;
    handles.PLSboot_gene=PLSboot_gene;handles.PLSboot_gp_Z=PLSboot_gp_Z;handles.PLSboot_protein=PLSboot_protein;
end
thres=str2double(get(handles.edit3,'string'));PLSboot_gene40_map2whole308=[handles.PLSboot_gene.Z,handles.PLSboot_gp_Z];
PLSboot_gene_PPI = [];influence_mat = zeros(size(PLSboot_gene40_map2whole308));% ����Ӱ�쿪��Ϊȫ�ر�
for idx = 1: length(handles.protein_matrix.symbols)
    % �ҵ���ǰ���׷������Ӷȵı�
    idx_ppi = intersect(union(find(strcmp(PPI.protein1,proteins_40(idx))),find(strcmp(PPI.protein2,proteins_40(idx)))),find(PPI.combined_score>=quantile(PPI.combined_score,thres,1)));
    PPI_string = unique(union(PPI.protein1(idx_ppi),PPI.protein2(idx_ppi)));% ����STRING���ҵ��뵱ǰ���׹����Ļ���
    [~,act_idx,~,] = intersect(string(handles.gene_matrix.symbols),string(PPI_string));
    influence_mat(act_idx,idx) = 1;% �򿪵��׶���Щ�����Ӱ�쿪��
    PLSboot_gene_PPI(:,idx) = influence_mat(:,idx).*PLSboot_gene40_map2whole308(:,idx+1);% д��Ӱ���ЧӦֵ
end
handles.PLSboot_gene_PPI=PLSboot_gene_PPI;handles.influence_mat=influence_mat;
[~,ia,ib] = intersect(handles.gene_matrix.ROIs,handles.Yraw_new(:,1));candidate_proteins = proteins_40(abs(handles.PLSboot_protein.Z) >= quantile(abs(handles.PLSboot_protein.Z),thres,1));
% candidate_proteins = proteins_40;
thresh_combined_score = quantile(PPI.combined_score,thres,1);% 0.5
PLSboot_gene_PPI2 = zeros(length(handles.gene_matrix.symbols),length(candidate_proteins));
for idx = 1: length(candidate_proteins)
    % �ҵ���ǰ���׷������Ӷȵı�
    idx_ppi = intersect(union(find(strcmp(PPI.protein1,candidate_proteins(idx))),find(strcmp(PPI.protein2,candidate_proteins(idx)))),find(PPI.combined_score>=thresh_combined_score));
    PPI_string = unique(union(PPI.protein1(idx_ppi),PPI.protein2(idx_ppi)));% ����STRING���ҵ��뵱ǰ����first-stage�����Ļ���
    [~,act_idx,~,] = intersect(string(handles.gene_matrix.symbols),string(PPI_string));
%     geneSymbol_PPI = handles.gene_matrix.symbols(act_idx);
    X_PPI = handles.gene_matrix.expr(ia,act_idx);Y = zscore(handles.trait(ib,:));
    max_dim = 1;geneSymbol_PPI = handles.gene_matrix.symbols(act_idx);
    fig_name = ['PLSslim_protein_',num2str(idx)];
    temp = Method_PLSbootstrap(X_PPI,Y,handles.dim,max_dim,handles.boot,geneSymbol_PPI,handles.out_path,fig_name);% �򿪵��׶���Щ�����Ӱ�쿪��
    PLSboot_gene_PPI2(act_idx,idx) = temp.Z;% д��Ӱ���ЧӦֵ
end
handles.PLSboot_gene_PPI2=PLSboot_gene_PPI2;
disp('Loading PPI Done!');
guidata(hObject,handles)

function pushbutton4_CreateFcn(~, ~, ~)

function pushbutton5_Callback(hObject, eventdata, handles)
disp('Running ...');
[C,ia,ib] = intersect(handles.Yraw_new(:,1),handles.protein_matrix.ROIs);% "deblank" to remove the blank of the string-end.
Y = zscore([handles.trait(ia,:),handles.protein_matrix.expr(ib,:)]);
[~,ia,ib] = intersect(handles.gene_matrix.ROIs,C);% "deblank" to remove the blank of the string-end.
X = handles.gene_matrix.expr(ia,:);Y = Y(ib,:);% X = zscore(gene_matrix.expr(ia,:));
fig_name = ['41'];PLSperm_41 = Method_PLSpermutation(X,Y,handles.dim,handles.rep,handles.out_path,fig_name);
temp=[2 1];max_dim=temp(((PLSperm_41.PCVAR(1) > PLSperm_41.PCVAR(2)) & (PLSperm_41.Pvalue_perm(1)<0.05))+1);
PLSboot_41 = Method_PLSbootstrap(X,Y,handles.dim,max_dim,handles.boot,handles.gene_matrix.symbols,handles.out_path,fig_name);

% [PLSboot_41.Z,add_Z,add_Z_abs,max_Z,combine_Z_fisher,combine_Z_brown]
PLSboot_gene40_map2whole308 = horzcat(handles.PLSboot_gene.Z,handles.PLSboot_gp_Z);
add_Z = PLSboot_gene40_map2whole308(:,1) + sum(handles.PLSboot_protein.Z'.*PLSboot_gene40_map2whole308(:,2:41),2)/40;
add_Z_abs = abs(PLSboot_gene40_map2whole308(:,1)) + sum(abs(handles.PLSboot_protein.Z'.*PLSboot_gene40_map2whole308(:,2:41)),2)/40;
temp = [PLSboot_gene40_map2whole308(:,1)';(handles.PLSboot_protein.Z'.*PLSboot_gene40_map2whole308(:,2:41))'];max_Z = temp((abs(temp)==max(abs(temp))));
combine_Z_fisher = sum(PLSboot_gene40_map2whole308,2);
% v = 1./(size(PLSboot_gene40_map2whole308,2)-3+PLSboot_gene40_map2whole308.^2);combine_Z_brown = sum(PLSboot_gene40_map2whole308 .* v,2) ./ sum(v,2);

% [R(:,1),score1,score2,max_R]
[R,R_P] = corr(X,Y);
[weig,weig_P] = corr(Y(:,2:41),Y(:,1));
% R = abs(R);weig = abs(weig);
R_sort =[];R_sort_idx =[];
for trait=1:size(Y,2)
    [temp_R,temp_idx] = sort(abs(R(:,trait)));% ascend,����
    [~,temp_idx] = sort(temp_idx); % rank of origin r-value
    R_sort =[R_sort,temp_R];% �洢����Ծ���
    R_sort_idx =[R_sort_idx,temp_idx];% �洢���������
end
R_sort_idx(R_P>0.01) = -Inf;% ����Բ�������ָ������Ϊ����ֵInf
score1 = 1./(length(handles.gene_matrix.symbols)+1 - R_sort_idx(:,1)) + sum(weig'./(length(handles.gene_matrix.symbols)+1 - R_sort_idx(:,2:41)),2);
R_sort_idx(R_P>0.01) = 0;% ����Բ�������ָ������Ϊ0
score2 = R_sort_idx(:,1) + sum(weig'.*(R_sort_idx(:,2:41)),2);% ʹ�ã�����-������Ϊ�÷�
temp = [R(:,1)';(weig'.*R(:,2:41))'];[~,index]=max(abs(temp),[],1);% max for each column
max_R = [];for idx=1:size(temp,2);max_R = [max_R;temp(index(idx),idx)];end %�ֺŵ����кϲ�


% [add_Z_PPI,max_Z_PPI,score1_PPI,score2_PPI,max_R_PPI,combine_Z_fisher_PPI]
PLSboot_protein=handles.PLSboot_protein;PLSboot_gene_PPI=handles.PLSboot_gene_PPI;geneSymbol=0;
add_Z_PPI = PLSboot_gene40_map2whole308(:,1) + sum(PLSboot_protein.Z'.*PLSboot_gene_PPI,2)/40;
temp = [PLSboot_gene40_map2whole308(:,1)';(PLSboot_protein.Z'.*PLSboot_gene_PPI)'];[~,index]=max(abs(temp),[],1);% max for each column
max_Z_PPI = [];for idx=1:size(temp,2); max_Z_PPI = [max_Z_PPI;temp(index(idx),idx)];end %�ֺŵ����кϲ�
new_Z_abs_PPI = abs(PLSboot_gene40_map2whole308(:,1)) + sum(abs(PLSboot_protein.Z'.*PLSboot_gene_PPI),2)/40;
R_sort_idx(R_P>0.01) = -Inf;% ����Բ�������ָ������Ϊ����ֵInf
score1_PPI = 1./(length(handles.gene_matrix.symbols)+1 - R_sort_idx(:,1)) + sum((weig'./(length(handles.gene_matrix.symbols)+1 - R_sort_idx(:,2:41))).*PLSboot_gene_PPI,2);
R_sort_idx(R_P>0.01) = 0;% ����Բ�������ָ������Ϊ0
score2_PPI = R_sort_idx(:,1) + sum((weig'.*(R_sort_idx(:,2:41))).*PLSboot_gene_PPI,2);% ʹ�ã�����-������Ϊ�÷�
temp = [R(:,1)';(weig'.*R(:,2:41).*PLSboot_gene_PPI)'];max_R_PPI = temp((abs(temp)==max(abs(temp))));
influence_mat=handles.influence_mat;influence_mat(:,1)=1;combine_Z_fisher_PPI = sum(influence_mat.*PLSboot_gene40_map2whole308,2);

%[PLSg,PLSg_p,PLSg_p_max,PLSg_p_IVW,PLSg_p_Z]
PLSboot_gene=handles.PLSboot_gene;PLSboot_gene_PPI2=handles.PLSboot_gene_PPI2;
PLSg = mean(PLSboot_gene_PPI2,2);
PLSg_p = PLSboot_gene.Z(:,1)+mean(PLSboot_gene_PPI2,2);
temp = [PLSboot_gene.Z(:,1)';PLSboot_gene_PPI2'];PLSg_p_max = temp((abs(temp)==max(abs(temp))));
PLSg_p_IVW = PLSboot_gene.Z(:,1)/var(PLSboot_gene.Z(:,1))+mean(PLSboot_gene_PPI2,2)/var(mean(PLSboot_gene_PPI2,2));
PLSg_p_Z = zscore(PLSboot_gene.Z(:,1))+zscore(mean(PLSboot_gene_PPI2,2));

res =table(handles.gene_matrix.symbols,PLSboot_gene.Z,PLSboot_41.Z,add_Z,max_Z,combine_Z_fisher,R(:,1),score2,max_R,max_Z_PPI,max_R_PPI,combine_Z_fisher_PPI,PLSg,PLSg_p,PLSg_p_max);
res.Properties.VariableNames(1:15) = [{'Symbol'},{'PLSboot_gene'},{'PLSboot_41'},{'add_Z'},{'max_Z'},{'combine_Z_fisher'},{'R'},{'score2'},{'max_R'},{'max_Z_PPI'},{'max_R_PPI'},{'combine_Z_fisher_PPI'},{'PLSg'},{'PLSg_p'},{'PLSg_p_max'}];
writetable(res,[handles.out_path,'map2whole308_PLS_boot_total.txt'],'Delimiter','\t','Encoding','UTF-8');

disp('Running Done!');delete(gcp('nocreate'))
% set(hObject,'ColumnName',{'Symbol','PLSboot_gene','PLSboot_41','add_Z','max_Z','combine_Z_fisher','R','score2','max_R','max_Z_PPI','max_R_PPI','combine_Z_fisher_PPI','PLSg','PLSg_p','PLSg_p_max'});
% reset(handles.uitable1);
% set(handles.uitable1,'ColumnEditable',logical(ones(1,15)));set(handles.uitable1,'ColumnName',res.Properties.VariableNames);
% set(handles.uitable1,'Data',res{1:50,:});

guidata(hObject,handles)

function pushbutton5_CreateFcn(~, ~, ~)



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

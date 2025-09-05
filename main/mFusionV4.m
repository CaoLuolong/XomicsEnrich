%%%%%%%%%%%%%%%%%%%%%%%%%%% main func %%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = mFusionV4(varargin)
% MFUSIONV4 MATLAB code for mFusionV4.fig
%      MFUSIONV4, by itself, creates a new MFUSIONV4 or raises the existing
%      singleton*.
%
%      H = MFUSIONV4 returns the handle to a new MFUSIONV4 or the handle to
%      the existing singleton*.
%
%      MFUSIONV4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MFUSIONV4.M with the given input arguments.
%
%      MFUSIONV4('Property','Value',...) creates a new MFUSIONV4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mFusionV4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mFusionV4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mFusionV4

% Last Modified by GUIDE v2.5 04-Sep-2025 11:06:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mFusionV4_OpeningFcn, ...
                   'gui_OutputFcn',  @mFusionV4_OutputFcn, ...
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
% --- Executes just before mFusionV4 is made visible.
function mFusionV4_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%% output function %%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = mFusionV4_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%% edit1 %%%%%%%%%%%%%%%%%%%%%%%%%%
function edit1_Callback(~, ~, ~)
function edit1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% edit2 %%%%%%%%%%%%%%%%%%%%%%%%%%
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
[parent_path, ~, ~] = fileparts(currentFolder);%work_dir of this code file
% parent_path=handles.parent_path;
handles.atlas_img=[parent_path,'/rawdata/brain_atlas/',list{val}];%% which atlas.
handles.aal = list{val};handles.parent_path = parent_path;addpath([handles.parent_path,'/utils']);

folder_name = [parent_path,'/Output22']; if exist(folder_name)==0 mkdir(folder_name); end
handles.out_path = [folder_name,'/',datestr(datetime('now'),'yyyy_mm_dd_HH_MM'),'/'];mkdir(handles.out_path);
set(handles.edit4,'String',handles.out_path);
guidata(hObject,handles);

function popupmenu2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu2_Callback(hObject, ~, handles)
val = get(handles.popupmenu2,'value');
list = get(handles.popupmenu2,'string');
handles.brain_sphere=list{val};%% left or right sphere.
guidata(hObject,handles);

function popupmenu3_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu3_Callback(hObject, ~, handles)
val = get(handles.popupmenu3,'value');
list = get(handles.popupmenu3,'string');
handles.interpolation=list{val};%% interpolation or not.
guidata(hObject,handles);

function popupmenu5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu5_Callback(hObject, eventdata, handles)
val = get(handles.popupmenu5,'value');list = get(handles.popupmenu5,'string');
% handles.PPIdb = '/rawdata/9606.protein.links.v11.5_new.txt';
switch val
    case 1
        handles.PPIdb = '/rawdata/9606.protein.links.v11.5_new.txt';
    case 2
        handles.PPIdb = '/rawdata/9606.protein.links.v11.5_new.txt';
    case 3
        handles.PPIdb = '/rawdata/9606.protein.physical.links.v11.5_new.txt';
    case 4
        handles.PPIdb = '/rawdata/9606.protein.links.v12.0_new.txt';
    case 5
        handles.PPIdb = '/rawdata/9606.protein.physical.links.v12.0_new.txt';
end
guidata(hObject,handles);

% --- Executes during object deletion, before destroying properties.
function edit1_DeleteFcn(~, ~, ~)
function edit2_DeleteFcn(~, ~, ~)
function edit3_DeleteFcn(~, ~, ~)
function edit4_DeleteFcn(~, ~, ~)
%% %%%%%%%%%%%%%%%%%%%%%%%%% pushbutton1 %%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(~, ~, ~)
function pushbutton1_Callback(hObject, eventdata, handles)
[Filename, Pathname]=uigetfile('*.txt;*.csv;*.xlsx','Select Input Ground Motion'); % Raw download PEER data
file=[Pathname,Filename];data=importdata(file);
% 检查 data 是否为结构体
if isstruct(data)
    % 如果是结构体，尝试提取 data 字段
    if isfield(data, 'data')
        data = data.data;
    else
        error('The imported data is a structure, but it does not contain a "data" field.');
    end
end
% data = data.data;
nanRows = any(isnan(data), 2);
data = data(~nanRows, :);% 删除含 NaN 的行
handles.Yraw_new= data;handles.trait= data(:,2:end);%there must not have missing value or NaN in 'file'.
handles.input=file;
set(handles.edit1,'String',file);

guidata(hObject,handles)
%% %%%%%%%%%%%%%%%%%%%%%%%%% pushbutton2 %%%%%%%%%%%%%%%%%%%%%%%%%%
% function pushbutton2_CreateFcn(~, ~, ~)
% 
% function pushbutton2_Callback(hObject, eventdata, handles)
% 
% guidata(hObject,handles)
 
function pushbutton3_CreateFcn(~, ~, ~)
function pushbutton3_Callback(hObject, eventdata, handles)
% [currentFolder, ~, ~] = fileparts(mfilename('fullpath'));
% [currentFolder, ~, ~] = fileparts(currentFolder)
% currentFolder = 'E:/Project/mFusionV4/';
handles.dim=str2double(get(handles.edit10,'string')); 
handles.rep=str2double(get(handles.edit11,'string')); 
handles.boot=str2double(get(handles.edit12,'string'));
diary([handles.out_path,'Run.log']);
clc;disp(['mFusion Version: 1.2.1']);disp(['Start Time: ',datestr(datetime('now'),'yyyy_mm_dd_HH_MM')]);
disp(['Trait: ',handles.input]);disp(['Out Path: ',handles.out_path]);
handles.depth=str2double(get(handles.edit2,'string'));handles.conf_score=str2double(get(handles.edit3,'string'));% set PPI dpeth and threshold.
handles.temp_file = [handles.out_path,'PLS_',num2str(handles.depth),'_',num2str(handles.conf_score),'_',handles.aal,'_',handles.brain_sphere,'_',handles.interpolation];
disp(['Atlas atlas::::sphere::::interpolation','=',handles.atlas_img,'::::',handles.brain_sphere,'::::',handles.interpolation]);
disp(['PLS dimension::::permutation::::bootstrap','=',num2str(handles.dim),'::::',num2str(handles.rep),'::::',num2str(handles.boot)]);
disp(['PPI data-base::::node-strth::::node-depth::::edge-scores','=',handles.PPIdb,'::::',get(handles.edit13,'string'),'::::',get(handles.edit2,'string'),'::::',get(handles.edit3,'string')]);
disp('Loading genes ...');
if isempty(gcp('nocreate')); parpool; end % open the parpool
rng(123); % set the random seed.
% gene--trait
handles.gene_matrix = Method_Gene_matrix(handles.atlas_img,handles.brain_sphere,handles.interpolation,handles.out_path);

[~,ia,ib] = intersect(handles.gene_matrix.ROIs,handles.Yraw_new(:,1));% "deblank" to remove the blank of the string-end.
X = handles.gene_matrix.expr(ia,:);% X = zscore(gene_matrix.expr(ia,:));
Y = zscore(handles.trait(ib,:));% Y = zscore([Yraw_new(ib,2),protein_matrix.expr(ib,:)]);
fig_name = 'gene';PLSperm_gene = Method_PLSpermutation(X,Y,handles.dim,handles.rep,handles.out_path,fig_name);
% temp=[2 1];max_dim=temp(((PLSperm_gene.PCTVAR(1) > PLSperm_gene.PCTVAR(2)) & (PLSperm_gene.Pvalue_perm(1)<0.05))+1);
temp=[find(PLSperm_gene.PCTVAR==max(PLSperm_gene.PCTVAR)) 1];max_dim=temp(((PLSperm_gene.PCTVAR(1) > PLSperm_gene.PCTVAR(2)) & (PLSperm_gene.Pvalue_perm(1)<0.05))+1);
handles.PLSboot_gene = Method_PLSbootstrap(X,Y,handles.dim,max_dim,handles.boot,handles.gene_matrix.symbols,handles.out_path,fig_name);
axes(handles.axes1);
plot(1:handles.dim,cumsum(100*PLSperm_gene.PCTVAR),'-o','LineWidth',1.5,'Color',[140/255,0,0]);ylim([0 100]);xlabel('gene dimension','FontSize',14);ylabel(strcat('Variance(%)'),'FontSize',14);
% gene--Protein
handles.protein_matrix = Method_Protein_matrix(handles.atlas_img,handles.brain_sphere,handles.out_path);
% [~,ia,ib] = intersect(handles.gene_matrix.ROIs,handles.protein_matrix.ROIs);% "deblank" to remove the blank of the string-end.
% PLSboot_gp_Z =zeros(size(handles.gene_matrix.expr,2),length(handles.protein_matrix.symbols));X = handles.gene_matrix.expr(ia,:);% X = zscore(gene_matrix.expr(ia,:));
% % PLSboot_gp_Z =rand(size(handles.gene_matrix.expr,2),length(handles.protein_matrix.symbols));X = handles.gene_matrix.expr(ia,:);% X = zscore(gene_matrix.expr(ia,:));
% for idx=1:length(handles.protein_matrix.symbols)
%     Y = zscore(handles.protein_matrix.expr(ib,idx));% Y = zscore([Yraw_new(ib,2),protein_matrix.expr(ib,:)]);
%     fig_name = [handles.protein_matrix.symbols{idx}];
% %     PLSperm_gene = Method_PLSpermutation(X,Y,handles.dim,handles.rep,handles.out_path,fig_name);temp=[2 1];max_dim=temp(((PLSperm_gene.PCTVAR(1) > PLSperm_gene.PCTVAR(2)) & (PLSperm_gene.Pvalue_perm(1)<0.05))+1);
%     [~,~,~,~,~,PCTVAR,~,~]=plsregress(X,Y,handles.dim);max_dim = find(PCTVAR(2,:)==max(PCTVAR(2,:)));% Logically, here max_dim should be calculated using permutation and then specified.
%     tempboot = Method_PLSbootstrap(X,Y,handles.dim,max_dim,handles.boot,handles.gene_matrix.symbols,handles.out_path,fig_name);
%     PLSboot_gp_Z(:,idx) =tempboot.Z;
%     disp(['Loading PETs: ',num2str(idx),'/',num2str(size(handles.protein_matrix.symbols,1))]);
% end
% load([handles.parent_path,'/rawdata/GenePetsPLS.mat']); % PLS from DK68 atlas
load([handles.parent_path,'/rawdata/GenePetsPLS308.mat']); % PLS from DK308 atlas
handles.PLSboot_gp_Z = PLSboot_gp_Z';
% Protein--trait
[~,ia,ib] = intersect(handles.protein_matrix.ROIs,handles.Yraw_new(:,1));% "deblank" to remove the blank of the string-end.
X = zscore(handles.protein_matrix.expr(ia,:));Y = zscore(handles.trait(ib,:));
fig_name = 'protein';% SNP, gene, protein
PLSperm_protein = Method_PLSpermutation(X,Y,handles.dim,handles.rep,handles.out_path,fig_name);
temp=[find(PLSperm_protein.PCTVAR==max(PLSperm_protein.PCTVAR)) 1];max_dim=temp(((PLSperm_protein.PCTVAR(1) > PLSperm_protein.PCTVAR(2)) & (PLSperm_protein.Pvalue_perm(1)<0.05))+1);
handles.PLSboot_protein = Method_PLSbootstrap(X,Y,handles.dim,max_dim,handles.boot,handles.protein_matrix.symbols,handles.out_path,fig_name);
axes(handles.axes2);
plot(1:handles.dim,cumsum(100*PLSperm_protein.PCTVAR),'-o','LineWidth',1.5,'Color',[140/255,0,0]);ylim([0 100]);xlabel('PET dimension','FontSize',14);ylabel(strcat('Variance(%)'),'FontSize',14);
% if ~isempty(gcp('nocreate')); delete(gcp('nocreate')); end % close the parpool
% save the two images
hsave = figure;% temp_fig = [handles.out_path,'PLS_of_gpt_',handles.aal,'_',handles.brain_sphere,'_',handles.interpolation,'_ROIs.fig'];
subplot(121);plot(1:handles.dim,cumsum(100*PLSperm_gene.PCTVAR),'-o','LineWidth',1.5,'Color',[140/255,0,0]);ylim([0 100]);xlabel('gene dimension','FontSize',14);ylabel(strcat('Variance(%)'),'FontSize',14);
subplot(122);plot(1:handles.dim,cumsum(100*PLSperm_protein.PCTVAR),'-o','LineWidth',1.5,'Color',[140/255,0,0]);ylim([0 100]);xlabel('PET dimension','FontSize',14);ylabel(strcat('Variance(%)'),'FontSize',14);
savefig(hsave,[handles.temp_file,'_ROIs.fig']);close(hsave);  
% save expression Matrix and PLS-Zscore
gene_matrix=handles.gene_matrix; protein_matrix=handles.protein_matrix;
PLSboot_gene=handles.PLSboot_gene;PLSboot_protein=handles.PLSboot_protein;
% depth=str2double(get(handles.edit2,'string'));conf_score=str2double(get(handles.edit3,'string'));% set PPI dpeth and threshold.
% temp_file = [handles.out_path,'PLS_of_gpt_',handles.depth,'_',handles.conf_score,'_',handles.aal,'_',handles.brain_sphere,'_',handles.interpolation,'_ROIs.mat'];
save([handles.temp_file,'_ROIs.mat'],'gene_matrix','protein_matrix','PLSperm_gene','PLSboot_gene','PLSperm_protein','PLSboot_protein','PLSboot_gp_Z');
disp(['Loaded ',num2str(size(handles.gene_matrix.symbols,1)),' genes and ',num2str(size(handles.protein_matrix.symbols,1)),' PET maps!']);
if (max(PLSperm_gene.PCTVAR)<=0.1)
    disp('These genes cannot express trait, please check the trait again!');
%     diary off;
%     error('quit!');
end
if (max(PLSperm_protein.PCTVAR)<=0.1)
    disp('These molecules cannot express trait, , please check the trait again!');
%     diary off;
%     error('quit!');
end
clear PLSperm_gene PLSperm_protein hsave gene_matrix protein_matrix PLSboot_gene tempboot PLSboot_gp_Z;

disp('Loading PPI ...');
% PPI = readtable([handles.parent_path,'/rawdata/9606.protein.links.v11.5_new.txt'],'Delimiter', '\t');
PPI = readtable([handles.parent_path,handles.PPIdb],'Delimiter', '\t');
PPI_protein1=PPI.protein1;PPI_protein2=PPI.protein2;PPI_combined_score=PPI.combined_score;clear PPI;
% proteins_40 = {'HTR1A','HTR1A','HTR1B','HTR1B','HTR2A','HTR2A','HTR4','HTR6','SLC6A4','CNR1','CNR1','DRD1',...
%     'DRD2','DRD2','DRD2','DRD2','DRD2','SLC6A3','SLC6A3','GABRA1','GABRA1','HTR1B','CHRM1','OPRM1','OPRM1',...
%     'SLC6A2','SLC6A2','SLC6A4','SLC6A4','SLC6A4','SV2A','VAT1L','VAT1L','VAT1L','VAT1L','GRM5','GRM5','GRM5','GRM5','CHRNB4'}';
proteins_40 = {'HTR1A','HTR1A','HTR1B','HTR1B','HTR1B','HTR2A','HTR2A','HTR2A','HTR4','HTR6',...
    'SLC6A4','SLC6A4','SLC6A4','SLC6A4','CNR1','CNR1','DRD1','DRD2','DRD2','DRD2','DRD2','DRD2',...
    'SLC6A3','SLC6A3','SLC6A3','GABRA1','GABRA1','HRH3','OPRM1','OPRM1','SLC6A2','SLC6A2','KIF17','SV2A',...
    'VAT1L','VAT1L','VAT1L','VAT1L','VAT1L','CHRM1','GRM5','GRM5','GRM5','GRM5','CHRNB4'}';

% depth=str2double(get(handles.edit2,'string'));conf_score=str2double(get(handles.edit3,'string'));% set PPI dpeth and threshold.
PLSboot_gene40_map2whole308=[handles.PLSboot_gene.Z,handles.PLSboot_gp_Z];
PLSboot_gene_PPI = zeros(size(handles.PLSboot_gp_Z));influence_mat = PLSboot_gene_PPI;% Set the influence switch to all off
strong_edges=find(PPI_combined_score>=quantile(PPI_combined_score,handles.conf_score,1));steps=handles.depth;

for idx = 1: length(handles.protein_matrix.symbols)
    % Find the edge where the current protein matches the connectivity
    idx_ppi = intersect(union(find(strcmp(PPI_protein1,proteins_40(idx))),find(strcmp(PPI_protein2,proteins_40(idx)))),strong_edges);
    PPI_string_temp1 = unique(union(PPI_protein1(idx_ppi),PPI_protein2(idx_ppi)));% Find the gene associated with the current protein from the STRING library
    disp([num2str(idx),'/',num2str(size(handles.protein_matrix.symbols,1)),':',num2str(length(PPI_string_temp1))]);
    if steps==1
        PPI_string = PPI_string_temp1;
    elseif steps==2
        idx_ppi2 =[];
        for idx2 = 1: length(PPI_string_temp1)
            % Find the edge where the current protein matches the connectivity
            idx_ppi = intersect(union(find(strcmp(PPI_protein1,PPI_string_temp1(idx2))),find(strcmp(PPI_protein2,PPI_string_temp1(idx2)))),strong_edges);   
            idx_ppi2 = [idx_ppi2;idx_ppi];
        end; clear idx_ppi;
        PPI_string = unique(union(PPI_protein1(idx_ppi2),PPI_protein2(idx_ppi2)));% Find the gene associated with the current protein from the STRING library
    elseif steps==3
        idx_ppi2 =[];idx_ppi3 =[];
        for idx2 = 1: length(PPI_string_temp1)
            % Find the edge where the current protein matches the connectivity
            idx_ppi = intersect(union(find(strcmp(PPI_protein1,PPI_string_temp1(idx2))),find(strcmp(PPI_protein2,PPI_string_temp1(idx2)))),strong_edges);   
            idx_ppi2 = [idx_ppi2;idx_ppi];
        end; clear idx_ppi;
        PPI_string_temp2 = unique(union(PPI_protein1(idx_ppi2),PPI_protein2(idx_ppi2)));% Find the gene associated with the current protein from the STRING library
        for idx2 = 1: length(PPI_string_temp2)
            % Find the edge where the current protein matches the connectivity
            idx_ppi = intersect(union(find(strcmp(PPI_protein1,PPI_string_temp2(idx2))),find(strcmp(PPI_protein2,PPI_string_temp2(idx2)))),strong_edges);   
            idx_ppi3 = [idx_ppi3;idx_ppi];
        end; clear idx_ppi;
        PPI_string = unique(union(PPI_protein1(idx_ppi3),PPI_protein2(idx_ppi3)));% Find the gene associated with the current protein from the STRING library
    else
        disp('=================There is no meaning for depth > 3, change into 1 automately.');
        PPI_string = PPI_string_temp1;
    end
    [~,act_idx,~,] = intersect(string(handles.gene_matrix.symbols),string(PPI_string));
%     if(size(act_idx,1)<10)% Prevent finding fewer than 10 genes.
%         idx_ppi = union(find(strcmp(PPI_protein1,proteins_40(idx))),find(strcmp(PPI_protein2,proteins_40(idx))));
%         temp=PPI_combined_score(idx_ppi);t=sort(temp(:));
%         [temp,~]=find(temp<=t(20),20);temp=idx_ppi(temp);
%         PPI_string = unique(union(PPI_protein1(temp),PPI_protein2(temp)));
%         [~,act_idx,~,] = intersect(string(handles.gene_matrix.symbols),string(PPI_string));
%     end
    influence_mat(act_idx,idx) = 1;% Turn on the protein's effect on these genes

    PLSboot_gene_PPI(:,idx) = influence_mat(:,idx).*PLSboot_gene40_map2whole308(:,idx+1);% write the protein's effect into matrix
end
handles.PLSboot_gene_PPI=PLSboot_gene_PPI;handles.influence_mat=influence_mat;handles.proteins_40=proteins_40;
% % candidate_proteins = proteins_40;
% [~,ia,ib] = intersect(handles.gene_matrix.ROIs,handles.Yraw_new(:,1));candidate_proteins = find(abs(handles.PLSboot_protein.Z) >= quantile(abs(handles.PLSboot_protein.Z),handles.conf_score,1));
% PLSboot_gene_PPI2 = zeros(length(handles.gene_matrix.symbols),length(candidate_proteins));
% for idx = 1: length(candidate_proteins)
%     % Find the edge where the current protein matches the connectivity
%     act_idx = find(influence_mat(:,candidate_proteins(idx)) == 1);
%     X_PPI = handles.gene_matrix.expr(ia,act_idx);Y = zscore(handles.trait(ib,:));
%     [~,~,~,~,~,PCTVAR,~,~]=plsregress(X_PPI,Y,handles.dim);max_dim = find(PCTVAR(2,:)==max(PCTVAR(2,:)));% Logically, here max_dim should be calculated using permutation and then specified.
%     geneSymbol_PPI = handles.gene_matrix.symbols(act_idx);
%     fig_name = ['PLSslim_protein_',num2str(idx)];
%     temp = Method_PLSbootstrap(X_PPI,Y,handles.dim,max_dim,handles.boot,geneSymbol_PPI,handles.out_path,fig_name);% Turn on the protein's effect on these genes
%     PLSboot_gene_PPI2(act_idx,idx) = temp.Z;% write the protein's effect into matrix
% end
% handles.PLSboot_gene_PPI2=PLSboot_gene_PPI2;
disp('Loading PPI Done!');
disp('Running pathway analysis...');
path_closeness = influence_mat.* handles.PLSboot_gp_Z;%pathway的紧密度
[~, ~, val_b] = find(path_closeness);

score_gene = repmat(handles.PLSboot_gene.Z,1,45);%复制一列到多列
score_PET = repmat(handles.PLSboot_protein.Z',15408,1);%复制一行到多行
score_path = (abs(score_gene)+abs(score_PET) ) /sqrt(2);% Z值相加，结果作为边的打分
%     score_path = (abs(score_gene)+abs(score_PET) ) /sqrt(2).*sign(score_gene).*sign(score_PET);% Z值相加并带符号，结果作为边的打分
%     score_path = (abs(PLSboot_gp_Z)+abs(score_gene)+abs(score_PET) ) /sqrt(3);% Z值相加,并考虑GP的紧密程度，结果作为边的打分
%     score_path = (abs(PLSboot_gp_Z)+abs(score_gene)+abs(score_PET) ) /sqrt(3).*sign(score_gene).*sign(score_PET);% Z值相加并带符号并考虑GP的紧密程度，结果作为边的打分
score_path = influence_mat.* score_path;
[row, col, val_a] = find(score_path);
results = table(val_a, handles.PLSboot_gene.gene(row), proteins_40(col),val_b,handles.PLSboot_protein.Z(col),handles.PLSboot_gene.Z(row),'VariableNames', {'Value','Gene', 'PETs', 'GPz','PTz','GTz'});
results.Edge= string(results.Gene) + "_" + string(results.PETs);
writetable(results, [handles.temp_file,'_PathwayScores.csv',], 'Delimiter', ',', 'WriteVariableNames', true);

% 有重复的edge，进行合并
uniqueKeys = unique(results.Edge);
mergedT = cell(length(uniqueKeys),4);
for idx = 1:length(uniqueKeys)        
    idx_rep = find(results.Edge == uniqueKeys{idx});% 找出当前键的所有索引
    % 如果当前键不止出现一次，则计算平均值；否则直接取值
    if length(idx_rep) > 1
        mergedT{idx,1} = mean(results.Value(idx_rep));
        mergedT{idx,2} = results.Gene{idx_rep(1)};
        mergedT{idx,3} = results.PETs{idx_rep(1)};
        mergedT{idx,4} = mean(results.GPz(idx_rep));
        mergedT{idx,5} = mean(results.PTz(idx_rep));
        mergedT{idx,6} = mean(results.GTz(idx_rep));
        mergedT{idx,7} = results.Edge{idx_rep(1)};
    else
        mergedT{idx,1} = results.Value(idx_rep);
        mergedT{idx,2} = results.Gene{idx_rep};
        mergedT{idx,3} = results.PETs{idx_rep};
        mergedT{idx,4} = results.GPz(idx_rep);
        mergedT{idx,5} = results.PTz(idx_rep);
        mergedT{idx,6} = results.GTz(idx_rep);
        mergedT{idx,7} = results.Edge{idx_rep};
    end
end
mergedT = table(mergedT(:,1),mergedT(:,2),mergedT(:,3),mergedT(:,4),mergedT(:,5),mergedT(:,6),mergedT(:,7),'VariableNames', {'Value','Gene', 'PETs', 'GPz','PTz','GTz','Edge'});
writetable(mergedT, [handles.temp_file,'_PathwayScores_unique.csv'], 'Delimiter', ',', 'WriteVariableNames', true);

disp('Running gene analysis...');pet_num=size(handles.PLSboot_gp_Z,2);
PLSboot_protein=handles.PLSboot_protein;PLSboot_gene=handles.PLSboot_gene;
PLSboot_gp_Z=handles.PLSboot_gp_Z;
PLSboot_gene_PPI=handles.PLSboot_gene_PPI;% PLSboot_gene_PPI=influence_mat.*PLSboot_gp_Z;
% PLSboot_gene_PPI2=handles.PLSboot_gene_PPI2;
% % [PLSboot_41.Z]
% [C,ia,ib] = intersect(handles.Yraw_new(:,1),handles.protein_matrix.ROIs);% "deblank" to remove the blank of the string-end.
% Y = zscore([handles.trait(ia,:),handles.protein_matrix.expr(ib,:)]);
% [~,ia,ib] = intersect(handles.gene_matrix.ROIs,C);% "deblank" to remove the blank of the string-end.
% X = handles.gene_matrix.expr(ia,:);Y = Y(ib,:);% X = zscore(gene_matrix.expr(ia,:));
% fig_name = [num2str(pet_num+1)];PLSperm_41 = Method_PLSpermutation(X,Y,handles.dim,handles.rep,handles.out_path,fig_name);
% temp=[2 1];max_dim=temp(((PLSperm_41.PCTVAR(1) > PLSperm_41.PCTVAR(2)) & (PLSperm_41.Pvalue_perm(1)<0.05))+1);
% PLSboot_41 = Method_PLSbootstrap(X,Y,handles.dim,max_dim,handles.boot,handles.gene_matrix.symbols,handles.out_path,fig_name);

% [add_Z,add_Z_abs,max_Z,combine_Z_fisher,combine_Z_brown]
PLSboot_gene40_map2whole308 = horzcat(PLSboot_gene.Z,handles.PLSboot_gp_Z);
PLS_multip = (repmat(handles.PLSboot_protein.Z',size(handles.PLSboot_gp_Z,1),1)+PLSboot_gene40_map2whole308(:,2:(pet_num+1)))/sqrt(2);% size: 15408*40
add_Z = (PLSboot_gene40_map2whole308(:,1) + sum(PLS_multip,2))/sqrt((pet_num+1));
% add_Z_abs = (abs(PLSboot_gene40_map2whole308(:,1)) + sum(abs(PLS_multip),2))/sqrt((pet_num+1));
temp = [PLSboot_gene40_map2whole308(:,1)';(PLS_multip)'];[~,index]=max(abs(temp),[],1);% max for each column?????????????
% temp = [PLSboot_gene40_map2whole308(:,1)';(sign(PLS_multip).*sqrt(abs(PLS_multip)))'];[~,index]=max(abs(temp),[],1);% max for each column?????????????
max_Z = [];for idx=1:size(temp,2);max_Z = [max_Z;temp(index(idx),idx)];end % ";" means row combine in MATLAB
% max_Z = temp((abs(temp)==max(abs(temp))));

combine_Z_fisher = sum(PLSboot_gene40_map2whole308,2)/sqrt((pet_num+1));
% % [R,max_R,max_R_PPI,score1,score1_PPI,score2,score2_PPI,]
% if size(handles.trait,2)==1
%     % [R(:,1),score1,score2,max_R]
%     [R,R_P] = corr(X,Y);
%     [weig,weig_P] = corr(Y(:,2:(pet_num+1)),Y(:,1));
%     % R = abs(R);weig = abs(weig);
%     R_sort =[];R_sort_idx =[];
%     for trait=1:size(Y,2)
%         [temp_R,temp_idx] = sort(abs(R(:,trait)));% ascend sort
%         [~,temp_idx] = sort(temp_idx); % rank of origin r-value
%         R_sort =[R_sort,temp_R];% corr matrix
%         R_sort_idx =[R_sort_idx,temp_idx];% corr order matrix
%     end
%     R_sort_idx(R_P>0.01) = -Inf;% The specified order for which the correlation is not significant is the maximum Inf
%     score1 = 1./(length(handles.gene_matrix.symbols)+1 - R_sort_idx(:,1)) + sum(weig'./(length(handles.gene_matrix.symbols)+1 - R_sort_idx(:,2:(pet_num+1))),2);
%     R_sort_idx(R_P>0.01) = 0;% The specified order for which the correlation is not significant is the 0
%     score2 = R_sort_idx(:,1) + sum(weig'.*(R_sort_idx(:,2:(pet_num+1))),2);% Use (total - sort) as the score
%     temp = [R(:,1)';(weig'.*R(:,2:(pet_num+1)))'];[~,index]=max(abs(temp),[],1);% max for each column?????????????
%     max_R = [];for idx=1:size(temp,2);max_R = [max_R;temp(index(idx),idx)];end % ";" means row combine in MATLAB
% % [add_Z_PPI,max_Z_PPI,score1_PPI,score2_PPI,max_R_PPI,combine_Z_fisher_PPI]
%     R_sort_idx(R_P>0.01) = -Inf;% The specified order for which the correlation is not significant is the maximum Inf
%     score1_PPI = 1./(length(handles.gene_matrix.symbols)+1 - R_sort_idx(:,1)) + sum((weig'./(length(handles.gene_matrix.symbols)+1 - R_sort_idx(:,2:(pet_num+1)))).*PLSboot_gene_PPI,2);
%     R_sort_idx(R_P>0.01) = 0;% The specified order for which the correlation is not significant is the 0
%     score2_PPI = R_sort_idx(:,1) + sum((weig'.*(R_sort_idx(:,2:(pet_num+1)))).*PLSboot_gene_PPI,2);% Use (total - sort) as the score
%     temp = [R(:,1)';(weig'.*R(:,2:(pet_num+1)).*PLSboot_gene_PPI)'];[~,index]=max(abs(temp),[],1);% max for each column?????????????
%     max_R_PPI = [];for idx=1:size(temp,2);max_R_PPI = [max_R_PPI;temp(index(idx),idx)];end 
%     % max_R_PPI = temp((abs(temp)==max(abs(temp))));
% else
%     R=zeros(size(max_Z));score1=zeros(size(max_Z));score2=score1;max_R=score1;score1_PPI=score1;score2_PPI=score1;max_R_PPI=score1;
% end

% conf_score=str2double(get(handles.edit3,'string'));
% cand_idx = find(abs(PLSboot_protein.Z) <= quantile(abs(PLSboot_protein.Z),handles.conf_score,1)); 
cand_idx=find(abs(PLSboot_protein.Z) <= str2double(get(handles.edit13,'string')));
PLSboot_gene_PPI(:,cand_idx) =0;PLSboot_protein.Z(cand_idx) =0;%PLSboot_gene_PPI=influence_mat.*PLSboot_gp_Z;%remove proteins which have small effect on trait.
PLS_multip_PPI = PLSboot_protein.Z'.*PLSboot_gene_PPI;% size: 15408*40
% add_Z_PPI = PLSboot_gene40_map2whole308(:,1) + sum(PLS_multip_PPI,2)/pet_num;
temp = [PLSboot_gene.Z,PLS_multip_PPI];temp(temp==0) = NaN;add_Z_PPI = sum(temp,2,'omitnan')./sqrt(sum(~isnan(temp), 2));
temp = [PLSboot_gene.Z';(PLS_multip_PPI)'];[~,index]=max(abs(temp),[],1);max_Z_PPI = [];for idx=1:size(temp,2); max_Z_PPI = [max_Z_PPI;temp(index(idx),idx)];end % max for each column

PPIweight=1;PLSboot_proteinREP=repmat(PLSboot_protein.Z',size(PLSboot_gp_Z,1),1);PLSboot_proteinREP(PLSboot_gene_PPI==0)=0;
PLS_multip_PPI2 = PPIweight*(abs(PLSboot_proteinREP)+abs(PLSboot_gene_PPI)).*sign(PLSboot_proteinREP).*sign(PLSboot_gene_PPI)/sqrt(2);% size: 15408*40
% add_Z_PPI = PLSboot_gene40_map2whole308(:,1) + sum(PLS_multip_PPI,2)/pet_num;
temp = [PLSboot_gene.Z,PLS_multip_PPI2];add_Z_PPI2_matrix=temp;
temp(temp==0) = NaN;add_Z_PPI2 = sum(temp,2,'omitnan')./sqrt(sum(~isnan(temp), 2));
temp = [PLSboot_gene.Z';(PLS_multip_PPI2)'];max_Z_PPI2_matrix=temp;
[~,index]=max(abs(temp),[],1);max_Z_PPI2 = [];for idx=1:size(temp,2); max_Z_PPI2 = [max_Z_PPI2;temp(index(idx),idx)];end % max for each column

% % [PLSg,PLSg_mean,PLSg_p,PLSg_p_max,PLSg_p_IVW,PLSg_p_Z]
% temp=PLSboot_gene_PPI2;temp(temp==0) = NaN;PLSg = sum(temp,2,'omitnan')./sqrt(sum(temp~=0,2));%Find the mean by the non-0 elements of the row
% PLSg_p = (PLSboot_gene.Z(:,1)+sum(temp,2,'omitnan'))./sqrt(sum(temp~=0,2)+1);
% temp = [PLSboot_gene.Z(:,1)';PLSboot_gene_PPI2'];[~,index]=max(abs(temp),[],1);% max for each column?????????????
% PLSg_p_max = [];for idx=1:size(temp,2);PLSg_p_max = [PLSg_p_max;temp(index(idx),idx)];end 
% % PLSg_p_IVW = PLSboot_gene.Z(:,1)/var(PLSboot_gene.Z(:,1))+mean(PLSboot_gene_PPI2,2)/var(mean(PLSboot_gene_PPI2,2));
% % PLSg_p_Z = zscore(PLSboot_gene.Z(:,1))+zscore(mean(PLSboot_gene_PPI2,2));

% if ~isempty(gcp('nocreate')); delete(gcp('nocreate')); end % close the parpool
influence_mat=handles.influence_mat;influence_mat(:,cand_idx) =0;
temp=[ones(size(influence_mat,1),1),influence_mat].*PLSboot_gene40_map2whole308;temp(temp==0) = NaN;combine_Z_fisher_PPI = sum(temp,2,'omitnan')./sqrt(sum(temp~=0,2));
geneSymbols = handles.gene_matrix.symbols;proteins_40=handles.proteins_40;
path_pt=[0;PLSboot_protein.Z];path_gp=[zeros(size(PLSboot_gene_PPI,1),1),PLSboot_gene_PPI];path_width=([ones(size(influence_mat,1),1),influence_mat].*PLSboot_gene40_map2whole308)';
save([handles.temp_file,'PPI_influence.mat'],'path_width','path_pt','path_gp','influence_mat','geneSymbols','proteins_40','add_Z_PPI2_matrix','max_Z_PPI2_matrix');

% res =table(handles.gene_matrix.symbols,PLSboot_gene.Z,PLSboot_41.Z,add_Z,max_Z,combine_Z_fisher,R(:,1),score1,score2,max_R,add_Z_PPI,max_Z_PPI,add_Z_PPI2,max_Z_PPI2,combine_Z_fisher_PPI,score1_PPI,score2_PPI,max_R_PPI,PLSg,PLSg_p,PLSg_p_max);
% res.Properties.VariableNames(1:21) = [{'Symbol'},{'PLS'},{'PLS_40'},{'add_Z'},{'max_Z'},{'combineZ'},{'R'},{'score1'},{'score2'},{'max_R'},{'add_Z_PPI'},{'max_Z_PPI'},{'addPPI'},{'maxPPI'},{'combine_Z_fisher_PPI'},{'score1_PPI'},{'score2_PPI'},{'max_R_PPI'},{'PLSg'},{'PLSg_p'},{'PLSg_p_max'}];
res =table(handles.gene_matrix.symbols,PLSboot_gene.Z,add_Z,max_Z,combine_Z_fisher,add_Z_PPI,max_Z_PPI,add_Z_PPI2,max_Z_PPI2,combine_Z_fisher_PPI);
% res.Properties.VariableNames(1:10) = [{'Symbol'},{'PLS'},{'add_Z'},{'max_Z'},{'combineZ'},{'add_Z_PPI'},{'max_Z_PPI'},{'addPPI'},{'maxPPI'},{'combine_Z_fisher_PPI'}];
res.Properties.VariableNames(1:10) = [{'Symbol'},{'PLS'},{'meanGPT'},{'maxGPT'},{'meanGP'},{'add_Z_PPI'},{'max_Z_PPI'},{'meanPPI'},{'maxPPI'},{'combine_Z_fisher_PPI'}];
writetable(res,[handles.temp_file,'_ROIs_total.txt'],'Delimiter','\t','Encoding','UTF-8');

res =table(handles.gene_matrix.symbols,PLSboot_gene.Z,combine_Z_fisher,add_Z,add_Z_PPI2,max_Z,max_Z_PPI2);
res.Properties.VariableNames(1:7) = [{'Symbol'},{'PLS'},{'meanGP'},{'meanGPT'},{'meanPPI'},{'maxGPT'},{'maxPPI'}];
writetable(res,[handles.temp_file,'_ROIs_slim_total.txt'],'Delimiter','\t','Encoding','UTF-8');

res =table(handles.gene_matrix.symbols,add_Z_PPI2);
res.Properties.VariableNames(1:2) = [{'Symbol'},{'meanPPI'}];
% temp_file = [handles.out_path,'PLS_',depth,'_',conf_score,'_',handles.aal,'_',handles.brain_sphere,'_',handles.interpolation,'_',num2str(size(Y,1)),'_ROIs_enrich_total.txt'];
writetable(res,[handles.temp_file,'_ROIs_enrich_total.txt'],'Delimiter','\t','Encoding','UTF-8');

% % pathway enrichment analysis.
% Rpath = 'E:\software\R\R-4.1.0\bin';
% [currentFolder, ~, ~] = fileparts(mfilename('fullpath'));
% [parent_path, ~, ~] = fileparts(currentFolder);%work_dir of this code file
% RunRcode([parent_path,'/utils/pathway_enrichment.R'],Rpath);
disp('Running Done!');
disp('Please do Gene Enrichment Analysis in R !');
disp(['End Time: ',datestr(datetime('now'),'yyyy_mm_dd_HH_MM')]);
diary off;
guidata(hObject,handles)

function edit10_Callback(hObject, eventdata, handles)
function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit11_Callback(hObject, eventdata, handles)
function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit12_Callback(hObject, eventdata, handles)
function edit12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_KeyPressFcn(hObject, eventdata, handles)
function edit2_KeyPressFcn(hObject, eventdata, handles)
function edit3_KeyPressFcn(hObject, eventdata, handles)

function popupmenu5_DeleteFcn(hObject, eventdata, handles)

function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function edit13_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

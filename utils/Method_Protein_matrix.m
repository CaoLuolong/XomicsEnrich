function Xmatrix = Method_Protein_matrix(atlas_img,brain_sphere,out_path)
% Function: Preprocession of receptor PET images.
% Updata date: 2022.12.16
% Email:luolongcao@163.com
% addpath('E:\Project\AHBAenrich\utils1026');
    %% recepter matrix of AAL116
%     atlas_img='D:\work_dir\AHBAenrich\rawdata\brain_atlas\AAL1_ROI_MNI_V4';%% recepter matrix of AAL116
% 	  atlas_img='D:\work_dir\AHBAenrich\rawdata\brain_atlas\AAL3v1_1mm';%% recepter matrix of AAL3
%     atlas_img='D:\work_dir\AHBAenrich\rawdata\brain_atlas\DK68_aparcaseg';% recepter matrix of DKT
%     atlas_img='D:\work_dir\AHBAenrich\rawdata\brain_atlas\500.aparc';%% recepter matrix of DK308
    v_DKT = spm_vol([atlas_img,'.Resliced.nii']);
    DKT_roi_name = readtable([atlas_img,'.nii.csv'],'Delimiter', ',');
    ROI_index = DKT_roi_name{:,2};
    [v_aal_img, XYZ] = spm_read_vols(v_DKT);    
    if strcmp(brain_sphere,'left')||strcmp(brain_sphere,'map2left')
%         delete right atlas
    v_aal_img(XYZ(1,:)>0) =0;
    end    
    ROI = unique(v_aal_img(:));
    [~,~,ib]=intersect(ROI,ROI_index);%there may be missing ROI labels in atlas or xmls.
    ROI = [table2cell(DKT_roi_name(ib,2)),table2cell(DKT_roi_name(ib,3))];

    clear location;
    for id = 1:size(ROI,1)
        location{id,1} = find(v_aal_img == ROI{id,1});
    end
%%     file_name = dir('F:\work_dir\PET\JuSpace_PETatlas\*.nii');
%     file_name = dir('D:\work_dir\AHBAenrich\rawdata\PET_Resliced\*.nii');
    file_name = dir([out_path,'../../rawdata/PET_Resliced_dup/*.nii']);
    result = [];
    colname = [];
    for id_PET = 1:size(file_name,1)
        file = strcat(file_name(id_PET).folder, '\', file_name(id_PET).name);
        temp = strsplit(file_name(id_PET).name, '.');
%         colname = [colname; strcat('PET_',temp(1))];
        colname = [colname; temp(1)];
        v_temp = spm_vol(file);
        [v_temp_img, ~] = spm_read_vols(v_temp);
        for id_ROI = 1:size(ROI,1)
            result(id_ROI,id_PET) = mean(v_temp_img(location{id_ROI,1}));
        end
    end

%%  output
%     out_result = table(['ROI_index', 'ROI_name',colname';ROI(:,1), ROIROI(:,2), num2cell(result)]);
    clear Xmatrix
    Xmatrix.parameter=[atlas_img];
    Xmatrix.expr=result; % expression matrix: #ROIs * #proteins.
    Xmatrix.ROIs=[ROI_index(ib)]; % ROIs names
    Xmatrix.ROIinfo=DKT_roi_name(ib,:);
    Xmatrix.symbols=colname; % genes symbols
     
%     temp_name = [out_path,datestr(datetime('now'),'yyyy_mm_dd_HH_MM'),'_PET_receptor_of_',num2str(size(ROI,1)),'_ROIs.mat'];
%     save(temp_name,'Xmatrix')
%     writetable(out_result, [temp_name, '.csv'], 'WriteVariableNames',false);
end
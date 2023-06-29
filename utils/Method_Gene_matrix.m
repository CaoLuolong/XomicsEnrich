function Xmatrix = Method_Gene_matrix(atlas_img,brain_sphere,interpolation,out_path)
% Function: Pre-procession of AHBA gene expression
% Updata date: 2022.09.20; this code contains Combination/Normalization. Get geneExpression matrix from [AAL or coordinates corresponding Y].
% Email:luolongcao@163.com
% atlas_img: an .nii image to segment the brain into ROIs.
% brain_sphere(1/3): left, map2left, whole
% interpolation(1/3): nearest, linear, no
    %%
%     addpath('E:/Project\AHBAenrich\utils1026');
    cd(out_path);cd('../..');
    load ./rawdata/AHBA_probe_Mean_reanote_1027.mat;%contain:expressionMean,probeInformation,sampleInfo
    load ./rawdata/AHBA_probe_Mean_reanote_expMeanRaw_1027.mat;%contain:Index_cortL,Index_subL_use,cort_expMeanRaw,sub_expMeanRaw
    load ./rawdata/AHBA_probe_Mean_reanote_expMeanRaw_scaled_1027.mat;%contain:genes,geneSymbol,cort_expMeanScaled,cort_l_expMS,sub_expMeanScaled,sub_l_expMS
    load ./rawdata/AHBA_ROI_index_50_1027.mat;%contain:sample_roi,Ncort_l_80,Nsub_l_80,Coordinatesall
    %%
%     atlas_img='.\rawdata\brain_atlas\brodmann';%% recepter matrix of brodmann
%     atlas_img='.\rawdata\brain_atlas\AAL1_ROI_MNI_V4';%% recepter matrix of AAL116
% 	  atlas_img='.\rawdata\brain_atlas\AAL3v1_1mm';%% recepter matrix of AAL3
%     atlas_img='.\rawdata\brain_atlas\DK68_aparcaseg';% recepter matrix of DKT,https://neurovault.org/images/23262/
%     atlas_img='.\rawdata\brain_atlas\500.aparc';%% recepter matrix of DK308
    v1 = spm_vol([atlas_img,'.Resliced.nii']);
    ROI_name = readtable([atlas_img,'.nii.csv'],'Delimiter', ',');
    [AAL_image, XYZ] = spm_read_vols(v1); %256*256*256
%     AAL_image(AAL_image<=1000)=0;spm_write_vol(v1,AAL_image); % delete non-cortex lables.
%     XYZ_lables = AAL_image(:);
    ROI_radius =5;%set a radius of a sphere to coordinata probes and ROIs.
    %% classify samples into ROIs by Euclidean distance.
    sum_point = [Ncort_l_80;Nsub_l_80];
    sum_exp = [vertcat(cort_l_expMS,sub_l_expMS)];
%     sum_point = [Ncort_l_80];
%     sum_exp = cort_l_expMS;
    if strcmp(brain_sphere,'left')
%         delete right atlas
        AAL_image(XYZ(1,:)>0) =0;
        ROI_name=ROI_name(table2array(ROI_name(:,4))<0,:);
    end
    if strcmp(brain_sphere,'map2left')
%         map right to left
        sum_point(:,3)=-abs(sum_point(:,3));
%         delete right atlas
        AAL_image(XYZ(1,:)>0) =0;
        ROI_name=ROI_name(table2array(ROI_name(:,4))<0,:);
    end
    if strcmp(brain_sphere,'whole')
%         continue;
    end
    if strcmp(brain_sphere,'map2whole')
%         map to both sphere
        temp = [-abs(sum_point(:,3));abs(sum_point(:,3))];
        sum_point = [sum_point;sum_point];sum_point(:,3) = temp;
        sum_exp = [sum_exp;sum_exp];
    end
    len = size(sum_point,1);
    sample_roi = zeros(len,2);
    B = XYZ';
    N = size(B,1);
    for i = 1:len % a method to get sample label
        A = sum_point(i,3:5);
        distance = sqrt(sum((B-repmat(A,N,1)).^2,2));% calculate distance of sample A and all point in B
        sample_roi(i,1) = i;
        temp_tab=tabulate(AAL_image(distance <= ROI_radius));
        if size(temp_tab,1) >1
            test=temp_tab(temp_tab(:,3)==max(temp_tab(2:end,3)),1);
            sample_roi(i,2) = test(1);
        end        
    end
    %% out put the gene expression Matrix
%     % delete right samples
    sample_roi_left = sample_roi;
%     sample_roi_left(find(sample_roi_left(:,2)<1000),2) = 0;

    % expression matrix: rows means ROI, columns means genes.
    AAL3_ROIs=unique(sample_roi_left(:,2));
    AAL3_ROIs=AAL3_ROIs(2:end);
    AAL3_ROIs_expr = zeros(length(AAL3_ROIs),size(sum_exp,2));
    for idx = 1:length(AAL3_ROIs)
%         temp_sample = find(sample_roi_left(:,2)==AAL3_ROIs(idx));
        AAL3_ROIs_expr(idx,:)=mean(sum_exp(sample_roi_left(:,2)==AAL3_ROIs(idx),:));
    end
    %% interpolation for some ROIs
    miss_rois=setdiff(table2array(ROI_name(:,2)),AAL3_ROIs);
    clear AAL3_ROIs_expr_miss
    if strcmp(interpolation ,'nearest') && ~isempty(miss_rois)
        B=sum_point(:,3:5);
        N=size(B,1);
        for idx =1:length(miss_rois)
            A=table2array(ROI_name(table2array(ROI_name(:,2))==miss_rois(idx),4:6));
            distance = sqrt(sum((B-repmat(A,N,1)).^2,2));
            AAL3_ROIs_expr_miss(idx,:)=mean(sum_exp(distance==min(distance),:),1);
        end
        AAL3_ROIs=vertcat(AAL3_ROIs,miss_rois);
        AAL3_ROIs_expr=vertcat(AAL3_ROIs_expr,AAL3_ROIs_expr_miss);
        [C,ia,ib]=intersect(AAL3_ROIs,sort(AAL3_ROIs));
        AAL3_ROIs_expr=AAL3_ROIs_expr(ia,:);
    end
    if strcmp(interpolation ,'linear') && ~isempty(miss_rois)
%         ad;------undone!!
        AAL3_ROIs=vertcat(AAL3_ROIs,miss_rois);
        AAL3_ROIs_expr=vertcat(AAL3_ROIs_expr,AAL3_ROIs_expr_miss);
        [C,ia,ib]=intersect(AAL3_ROIs,sort(AAL3_ROIs));
        AAL3_ROIs_expr=AAL3_ROIs_expr(ia,:);
    end
    if strcmp(interpolation ,'no') && ~isempty(miss_rois)
%         continue;
    end
    [C,ia,ib]=intersect(AAL3_ROIs,table2array(ROI_name(:,2)));
    AAL3_ROIs_expr=AAL3_ROIs_expr(ia,:);
    AAL3_ROIs=AAL3_ROIs(ia,:);
    %% output
    clear Xmatrix
    Xmatrix.parameter=[atlas_img,'-',brain_sphere,'-',interpolation];
    Xmatrix.sample=sample_roi; % each sample belongs to which ROI.
    Xmatrix.expr=AAL3_ROIs_expr; % expression matrix: #ROIs * #genes.
    Xmatrix.ROIs=AAL3_ROIs; % ROIs names
    Xmatrix.ROIinfo=ROI_name;
    Xmatrix.genes=genes; % genes names
    Xmatrix.symbols=geneSymbol; % genes symbols
%     temp_file = [out_path,datestr(datetime('now'),'yyyy_mm_dd_HH_MM'),'_genes_expr_of_',brain_sphere,num2str(length(AAL3_ROIs)),'_ROIs.mat'];
%     save(temp_file,'Xmatrix');
end

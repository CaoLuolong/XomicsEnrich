atlas_img = 'E:\Project\AHBAenrich\rawdata\brain_atlas\brodmann.nii';%% recepter matrix of brodmann
v1 = spm_vol([atlas_img]);
[AAL_image, XYZ] = spm_read_vols(v1); %256*256*256
XYZ = XYZ';
roi=unique(AAL_image(:));roi(1)=[];
v_test = v1;v_test.fname = atlas_img;img_test = zeros(size(AAL_image));
mean_MNI = zeros(length(roi),6);
for idx=1:length(roi)
index_left = find(AAL_image(:)==roi(idx) & XYZ(:,1)<0);index_right = find(AAL_image(:)==roi(idx) & XYZ(:,1)>=0);
img_test(index_left)=roi(idx);img_test(index_right)=roi(idx)+52;% total 52 ROIs in half sphere.
mean_MNI(idx,:) = [mean(XYZ(index_left,:),1),mean(XYZ(index_right,:),1)];
disp([num2str(idx),':',num2str(length(index_left)),':',num2str(length(index_right))])
end
unique(img_test(:));
spm_write_vol(v_test,img_test);
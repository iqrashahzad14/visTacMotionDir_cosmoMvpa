function [opt] = chooseMaskNEW(opt)

%     for iSub = 1:numel(opt.subjects)
% 
        subID = opt.subjects{iSub};
%         % get subject folder name
        subFolder = ['sub-', subID];
% %         radius=opt.radius;
% %         radiusFolder = ['radius-', radius];

        opt.maskPath = fullfile(fileparts(mfilename('fullpath')), '..','..', '..','derivatives' ,'cpp_spm-roi',subFolder);%,radiusFolder, radius);

        % masks to decode/use
        opt.maskName = {strcat('sub-',num2str(subID),'_hemi-L_space-MNI_label-lhMT_vox-500', '.nii'), ...
            strcat('sub-',num2str(subID),'_hemi-R_space-MNI_label-rhMT_vox-500', '.nii'), ...
            strcat('sub-',num2str(subID),'_hemi-L_space-MNI_label-lS1_vox-500', '.nii'), ...
            strcat('sub-',num2str(subID),'_hemi-L_space-MNI_label-lPC_vox-500', '.nii'), ...
            strcat('sub-',num2str(subID),'_hemi-R_space-MNI_label-rPC_vox-500', '.nii'), ...
            strcat('sub-',num2str(subID),'_hemi-L_space-MNI_label-lMTt_vox-500', '.nii'), ...
            strcat('sub-',num2str(subID),'_hemi-R_space-MNI_label-rMTt_vox-500', '.nii'),...
            };

        % use in output roi name
        opt.maskLabel = {'lhMT', 'rhMT', 'lS1', 'lPC', 'rPC','lMTt', 'rMTt'};
%     end
    
end
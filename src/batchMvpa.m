clear ;
clc;

%% set paths
  % spm
  warning('off');
  addpath(genpath('/Users/shahzad/Documents/MATLAB/spm12'));
  % cosmo
  cosmo = '~/Documents/MATLAB/CoSMoMVPA';
  addpath(genpath(cosmo));
  cosmo_warning('once');

  % libsvm
  libsvm = '~/Documents/MATLAB/libsvm';
  addpath(genpath(libsvm));
  % verify it worked.
  cosmo_check_external('libsvm'); % should not give an error
  
  % add cpp repo
  run ../../visTacMotionDir_fMRI_analysis/src/../lib/CPP_SPM/initCppSpm.m;
     
  % load your options
  opt = getOptionMvpa();

  %% run mvpa 
  
  % use parcels or NS masks?
%   roiSource = 'hmat'; % 'freesurfer', 'neurosynth', ...
%   accuracy = calculateMvpa(opt, roiSource);
%     maskVoxel = calculateMaskSize(opt);
    
    % take the most responsive xx nb of voxels
%   opt.mvpa.ratioToKeep = 100;% % 120, 224, 368, 100 150 250 350 420
%   %sub005, rad-8mm, 133; rad-10, 249; rad-12 414;  
%   
% %   accuracy = calculateMvpa(opt);
%   
% %   accuracyWithinModality = calculateMvpaWithinModality(opt);
%     accuracyWithinModality = calculateMvpaWithinModalityNEW(opt);

% %   accuracyCrossModal = calculateMvpaCrossModal(opt);
%       accuracyCrossModal = calculateMvpaCrossModalNEW(opt);
      
% %   accuracyModalityDecoding = calculateMvpaModalityDecoding(opt);
%         accuracyModalityDecoding = calculateMvpaModalityDecodingNEW(opt);


  opt.mvpa.ratioToKeep = 100;
  accuracyCrossModal = calculateMvpaCrossModalNEW(opt);
  
  opt.mvpa.ratioToKeep = 150;
  accuracyCrossModal = calculateMvpaCrossModalNEW(opt);
  
  opt.mvpa.ratioToKeep = 200;
  accuracyCrossModal = calculateMvpaCrossModalNEW(opt);
  
  opt.mvpa.ratioToKeep = 250;
  accuracyCrossModal = calculateMvpaCrossModalNEW(opt);
  
  opt.mvpa.ratioToKeep = 300;
  accuracyCrossModal = calculateMvpaCrossModalNEW(opt);
  
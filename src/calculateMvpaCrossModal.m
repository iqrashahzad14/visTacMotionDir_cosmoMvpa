function accu = calculateMvpaCrossModal(opt)

% main function which loops through masks and subjects to calculate the
% decoding accuracies for given conditions.
% dependant on SPM + CPP_SPM and CosMoMvpa toolboxes

  % get the smoothing parameter for 4D map
  funcFWHM = opt.funcFWHM;

  % choose masks to be used
  opt = chooseMask(opt);
    

  % set output folder/name
  savefileMat = fullfile(opt.pathOutput, ...
                         [opt.taskName, ...
                         '_CrossModal',...
                         '_radius', num2str(opt.radius),...
                          '_smoothing', num2str(funcFWHM), ...
                          '_ratio', num2str(opt.mvpa.ratioToKeep), ...
                          '_', datestr(now, 'yyyymmddHHMM'), '.mat']);


  %% MVPA options

  % set cosmo mvpa structure
  condLabelNb = [1 2 ];
  condLabelName = {'vertical', 'horizontal'};
  modalityLabelName = {'visual','tactile'};
  decodingConditionList = {'trainVisual_testTactile','trainTactile_testVision','both'};
  modalityLabelNb = [1 2];
  

  %% let's get going!

  % set structure array for keeping the results
  accu = struct( ...
                'subID', [], ...
                'mask', [], ...
                'accuracy', [], ...
                'radius', [], ...
                'prediction', [], ...
                'maskVoxNb', [], ...
                'choosenVoxNb', [], ...
                'image', [], ...
                'ffxSmooth', [], ...
                'roiSource', [], ...
                'decodingCondition', [], ...
                'permutation', [], ...
                'imagePath', []);

  count = 1;

  for iSub = 1:numel(opt.subjects)

    % get FFX path
    subID = opt.subjects{iSub};
    ffxDir = getFFXdir(subID, funcFWHM, opt);

    % get subject folder name
    subFolder = ['sub-', subID];

    for iImage = 1:length(opt.mvpa.map4D)

      for iMask = 1:length(opt.maskName)

        % choose the mask
        mask = fullfile(opt.maskPath, opt.maskName{iMask});

        % display the used mask
        disp(opt.maskName{iMask});
        
        % 4D image
        imageName = ['4D_', opt.mvpa.map4D{iImage}, '_', num2str(funcFWHM), '.nii'];
        image = fullfile(ffxDir, imageName);
%         idx = [];
%         ds = [];
%         targets = [];
%         chunks = [];
        
        for iModality=1:3 %see the types in decoding conditionlist
            if iModality==1
                test=1;
                decodingCondition=decodingConditionList(1);
            elseif iModality==2
                test = 2;
                decodingCondition=decodingConditionList(2);
            elseif iModality==3
                test = [];
                decodingCondition=decodingConditionList(3);
            end
            
                        
            % load cosmo input
            ds = cosmo_fmri_dataset(image, 'mask', mask);
            
            ds = cosmo_remove_useless_data(ds);

%             % Getting rid off zeros
%             zeroMask = all(ds.samples == 0, 1);
%             ds = cosmo_slice(ds, ~zeroMask, 2);

            % set cosmo structure
            ds = setCosmoStructure(opt, ds, condLabelNb, condLabelName, modalityLabelNb, modalityLabelName);

            % Demean  every pattern to remove univariate effect differences
            meanPattern = mean(ds.samples,2);  % get the mean for every pattern
            meanPattern = repmat(meanPattern,1,size(ds.samples,2)); % make a matrix with repmat
            ds.samples  = ds.samples - meanPattern; % remove the mean from every every point in each pattern

            % Slice the dataset accroding to modality
            modIdx = (ds.sa.modality==1) | (ds.sa.modality==2) ;
            ds = cosmo_slice(ds,modIdx) ;
            
            % Slice the dataset accroding to motion direction  
            %slice according to directions
            ds = cosmo_slice(ds,strcmp(ds.sa.labels,'vertical') | strcmp(ds.sa.labels,'horizontal')) ;
            
            % partitioning, for test and training : cross validation
            %partitions = cosmo_nfold_partitioner(ds);
            partitions=cosmo_nchoosek_partitioner(ds, 1, 'modality',test);

            
            % remove constant features
            %ds = cosmo_remove_useless_data(ds);

            % calculate the mask size
            maskVoxel = size(ds.samples, 2);

            

            % use the ratios, instead of the voxel number:
            opt.mvpa.feature_selection_ratio_to_keep = opt.mvpa.ratioToKeep;

            % ROI mvpa analysis
            [pred, accuracy] = cosmo_crossvalidate(ds, ...
                                       @cosmo_classify_meta_feature_selection, ...
                                       partitions, opt.mvpa);

            %% store output
            accu(count).subID = subID;
            accu(count).mask = opt.maskLabel{iMask};
            accu(count).maskVoxNb = maskVoxel;
            accu(count).choosenVoxNb = opt.mvpa.feature_selection_ratio_to_keep;
           % accu(count).choosenVoxNb = round(maskVoxel * maxRatio);
           % accu(count).maxRatio = maxRatio;
            accu(count).image = opt.mvpa.map4D{iImage};
            accu(count).ffxSmooth = funcFWHM;
            accu(count).accuracy = accuracy;
            accu(count).radius = opt.radius;
            accu(count).prediction = pred;
            accu(count).imagePath = image;
    %         accu(count).roiSource = roiSource;
            accu(count).decodingCondition = decodingCondition;

            %% PERMUTATION PART
            if opt.mvpa.permutate  == 1
              % number of iterations
              nbIter = 100;

              % allocate space for permuted accuracies
              acc0 = zeros(nbIter, 1);

              % make a copy of the dataset
              ds0 = ds;

              % for _niter_ iterations, reshuffle the labels and compute accuracy
              % Use the helper function cosmo_randomize_targets
              for k = 1:nbIter
                      % manaully randomize the targets (because of cross-modal error)
                        % In every modality seperatly and in every chunk ,
                        % randomize the directions
                    for iChunk=1:max(ds.sa.chunks)
                        for iTestModality = 1:max(ds.sa.modality)
                            ds0.sa.targets(ds.sa.chunks==iChunk & ds.sa.modality==iTestModality) = Shuffle(ds.sa.targets(ds.sa.chunks==iChunk & ds.sa.modality==iTestModality));
                        end
                    end
%                 ds0.sa.targets = cosmo_randomize_targets(ds);
                [~, acc0(k)] = cosmo_crossvalidate(ds0, ...
                                                   @cosmo_meta_feature_selection_classifier, ...
                                                   partitions, opt.mvpa);
              end

              p = sum(accuracy < acc0) / nbIter;
              fprintf('%d permutations: accuracy=%.3f, p=%.4f\n', nbIter, accuracy, p);

              % save permuted accuracies
              accu(count).permutation = acc0';
            end

            % increase the counter and allons y!
            count = count + 1;

            fprintf(['Sub'  subID ' - area: ' opt.maskLabel{iMask} ...
                     ', accuracy: ' num2str(accuracy) '\n\n\n']);

        end
      end
    end
  end
  %% save output

  % mat file
  save(savefileMat, 'accu');

  % csv but with important info for plotting
%   csvAccu = rmfield(accu, 'permutation');
%   csvAccu = rmfield(csvAccu, 'prediction');
%   csvAccu = rmfield(csvAccu, 'imagePath');
%   writetable(struct2table(csvAccu), savefileCsv);

end

function ds = setCosmoStructure(opt, ds, condLabelNb, condLabelName, modalityLabelNb, modalityLabelName)
  % sets up the target, chunk, labels by stimuli condition labels, runs,
  % number labels.

  % design info from opt
  nbRun = opt.mvpa.nbRun;
  betasPerCondition = opt.mvpa.nbTrialRepetition;

  % chunk (runs), target (condition), labels (condition names)
  conditionPerRun = length(condLabelNb);
  betasPerRun = betasPerCondition * conditionPerRun;

  chunks = repmat((1:nbRun)', 1, betasPerRun*conditionPerRun);
  chunks = chunks(:);

  targets = repmat(condLabelNb', 1, nbRun)';
  targets = targets(:);
  targets = repmat(targets, betasPerCondition, 2);
  targets = targets(:);%IQRA

  condLabelName = repmat(condLabelName', 1, nbRun)';
  condLabelName = condLabelName(:);
  condLabelName = repmat(condLabelName, betasPerCondition, 2);
  condLabelName = condLabelName(:);
  
  modalityLabelName = repmat(modalityLabelName', 1, nbRun*conditionPerRun)';
  modalityLabelName = modalityLabelName(:);
  modalityLabelName = repmat(modalityLabelName, betasPerCondition, 1);
  
  modalityLabelNb = repmat(modalityLabelNb', 1, nbRun*conditionPerRun)';
  modalityLabelNb = modalityLabelNb(:);
  modalityLabelNb = repmat(modalityLabelNb, betasPerCondition, 1);
  
%   modality = repmat(modalityLabelNb',1,nbRun*conditionPerRun)';
%   modality = modality(:) ;

  % assign our 4D image design into cosmo ds git
  ds.sa.targets = targets;
  ds.sa.chunks = chunks;
  ds.sa.labels = condLabelName;
  ds.sa.modality = modalityLabelNb;
 

  % figure; imagesc(ds.sa.chunks);

end

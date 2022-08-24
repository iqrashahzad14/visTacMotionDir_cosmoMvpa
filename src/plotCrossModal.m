
im='beta';%'t_maps', 'beta'
smooth='2';
voxNb='150';

subList={'004','005', '006', '007'};
roiList={'lhMT','rhMT','lS1','lPC', 'rPC', 'lMTt', 'rMTt'};
decodingConditionList = {'trainVisual_testTactile','trainTactile_testVision','both'};

visDir_lhMT=[]; visDir_rhMT=[]; visDir_lS1=[]; visDir_lPC=[]; visDir_rPC=[]; visDir_lMTt=[]; visDir_rMTt=[];
tacDir_lhMT=[]; tacDir_rhMT=[]; tacDir_lS1=[]; tacDir_lPC=[]; tacDir_rPC=[]; tacDir_lMTt=[]; tacDir_rMTt=[];
both_lhMT=[]; both_rhMT=[]; both_lS1=[]; both_lPC=[]; both_rPC=[]; both_lMTt=[]; both_rMTt=[];



for iAccu=1:length(accu)
    for iSub=1:length(subList)
        subID=subList(iSub);
        if strcmp(char({accu(iAccu).subID}.'),char(subID))==1

            if strcmp(char({accu(iAccu).image}.'), im)==1 && strcmp(num2str([accu(iAccu).ffxSmooth].'),smooth)==1 && strcmp(num2str([accu(iAccu).choosenVoxNb].'),voxNb)==1
                
                
                if strcmp(char({accu(iAccu).mask}.'), 'lhMT')==1
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'trainVisual_testTactile' )==1
                        visDir_lhMT = [visDir_lhMT;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'trainTactile_testVision' )==1
                        tacDir_lhMT = [tacDir_lhMT;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'both' )==1
                        both_lhMT = [both_lhMT;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'rhMT')==1
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'trainVisual_testTactile' )==1
                        visDir_rhMT = [visDir_rhMT;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'trainTactile_testVision' )==1
                        tacDir_rhMT = [tacDir_rhMT;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'both' )==1
                        both_rhMT = [both_rhMT;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'lS1')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'trainVisual_testTactile' )==1
                        visDir_lS1 = [visDir_lS1;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'trainTactile_testVision' )==1
                        tacDir_lS1 = [tacDir_lS1;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'both' )==1
                        both_lS1 = [both_lS1;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'lPC')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'trainVisual_testTactile' )==1
                        visDir_lPC = [visDir_lPC;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'trainTactile_testVision' )==1
                        tacDir_lPC = [tacDir_lPC;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'both' )==1
                        both_lPC = [both_lPC;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'rPC')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'trainVisual_testTactile' )==1
                        visDir_rPC = [visDir_rPC;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'trainTactile_testVision' )==1
                        tacDir_rPC = [tacDir_rPC;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'both' )==1
                        both_rPC = [both_rPC;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'lMTt')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'trainVisual_testTactile' )==1
                        visDir_lMTt = [visDir_lMTt;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'trainTactile_testVision' )==1
                        tacDir_lMTt = [tacDir_lMTt;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'both' )==1
                        both_lMTt = [both_lMTt;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'rMTt')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'trainVisual_testTactile' )==1
                        visDir_rMTt = [visDir_rMTt;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'trainTactile_testVision' )==1
                        tacDir_rMTt = [tacDir_rMTt;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'both' )==1
                        both_rMTt = [both_rMTt;[accu(iAccu).accuracy].'];
                    end
                end
            end    
       end
            

    end
end

decodAccu=[visDir_lhMT,tacDir_lhMT,both_lhMT, visDir_rhMT, tacDir_rhMT,both_rhMT,...
    visDir_lS1, tacDir_lS1, both_lS1, visDir_lPC, tacDir_lPC,both_lPC, visDir_rPC, tacDir_rPC,both_rPC,...
    visDir_lMTt, tacDir_lMTt, both_lMTt, visDir_rMTt, tacDir_rMTt, both_rMTt ];

meanDecodAccu=mean(decodAccu);
seDecodAccu=std(decodAccu)/sqrt(length(subList));

T = array2table(decodAccu,...
    'VariableNames',{'visDir_lhMT','tacDir_lhMT','both_lhMT', 'visDir_rhMT', 'tacDir_rhMT','both_rhMT',...
    'visDir_lS1', 'tacDir_lS1','both_lS1', 'visDir_lPC', 'tacDir_lPC','both_lPC', 'visDir_rPC', 'tacDir_rPC','both_rPC',...
    'visDir_lMTt', 'tacDir_lMTt', 'both_lMTt', 'visDir_rMTt', 'tacDir_rMTt','both_rMTt'});

figure(1)
% subNb=1;
% sub_series = [visDir_lhMT(subNb,:),tacDir_lhMT(subNb,:); visDir_rhMT(subNb,:), tacDir_rhMT(subNb,:);...
%     visDir_lS1(subNb,:), tacDir_lS1(subNb,:); visDir_lPC(subNb,:), tacDir_lPC(subNb,:); visDir_rPC(subNb,:), tacDir_rPC(subNb,:);...
%     visDir_lMTt(subNb,:), tacDir_lMTt(subNb,:); visDir_rMTt(subNb,:), tacDir_rMTt(subNb,:) ];

model_series = [mean(visDir_lhMT),mean(tacDir_lhMT), mean(both_lhMT); mean(visDir_rhMT), mean(tacDir_rhMT),mean(both_rhMT);...
    mean(visDir_lS1), mean(tacDir_lS1),mean(both_lS1); mean(visDir_lPC), mean(tacDir_lPC), mean(both_lPC); mean(visDir_rPC), mean(tacDir_rPC), mean(both_rPC);...
    mean(visDir_lMTt), mean(tacDir_lMTt), mean(both_lMTt); mean(visDir_rMTt), mean(tacDir_rMTt), mean(both_rMTt)];

model_error = [std(visDir_lhMT)/sqrt(length(subList)),std(tacDir_lhMT)/sqrt(length(subList)), std(both_lhMT)/sqrt(length(subList)); std(visDir_rhMT)/sqrt(length(subList)), std(tacDir_rhMT)/sqrt(length(subList)), std(both_rhMT)/sqrt(length(subList));...
    std(visDir_lS1)/sqrt(length(subList)), std(tacDir_lS1)/sqrt(length(subList)), std(both_lS1)/sqrt(length(subList)); std(visDir_lPC)/sqrt(length(subList)), std(tacDir_lPC)/sqrt(length(subList)), std(both_lPC)/sqrt(length(subList)); std(visDir_rPC)/sqrt(length(subList)), std(tacDir_rPC)/sqrt(length(subList)), std(both_rPC)/sqrt(length(subList));...
    std(visDir_lMTt)/sqrt(length(subList)), std(tacDir_lMTt)/sqrt(length(subList)), std(both_lMTt)/sqrt(length(subList)); std(visDir_rMTt)/sqrt(length(subList)), std(tacDir_rMTt)/sqrt(length(subList)), std(both_rMTt)/sqrt(length(subList))];

b = bar(model_series, 'grouped'); 

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
hold on
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');

scatter(repmat(x(1,1), length(visDir_lhMT), 1),visDir_lhMT,'MarkerEdgeColor','k');%,'MarkerFaceColor','k');
scatter(repmat(x(2,1), length(tacDir_lhMT), 1),tacDir_lhMT,'MarkerEdgeColor','k');
scatter(repmat(x(3,1), length(both_lhMT), 1),both_lhMT,'MarkerEdgeColor','k');

scatter(repmat(x(1,2), length(visDir_rhMT), 1),visDir_rhMT,'MarkerEdgeColor','k');
scatter(repmat(x(2,2), length(tacDir_rhMT), 1),tacDir_rhMT,'MarkerEdgeColor','k');
scatter(repmat(x(3,2), length(both_rhMT), 1),both_rhMT,'MarkerEdgeColor','k');

scatter(repmat(x(1,3), length(visDir_lS1), 1),visDir_lS1,'MarkerEdgeColor','k');
scatter(repmat(x(2,3), length(tacDir_lS1), 1),tacDir_lS1,'MarkerEdgeColor','k');
scatter(repmat(x(3,3), length(both_lS1), 1),both_lS1,'MarkerEdgeColor','k');

scatter(repmat(x(1,4), length(visDir_lPC), 1),visDir_lPC,'MarkerEdgeColor','k');
scatter(repmat(x(2,4), length(tacDir_lPC), 1),tacDir_lPC,'MarkerEdgeColor','k');
scatter(repmat(x(3,4), length(both_lPC), 1),both_lPC,'MarkerEdgeColor','k');

scatter(repmat(x(1,5), length(visDir_rPC), 1),visDir_rPC,'MarkerEdgeColor','k');
scatter(repmat(x(2,5), length(tacDir_rPC), 1),tacDir_rPC,'MarkerEdgeColor','k');
scatter(repmat(x(3,5), length(both_rPC), 1),both_rPC,'MarkerEdgeColor','k');

scatter(repmat(x(1,6), length(visDir_lMTt), 1),visDir_lMTt,'MarkerEdgeColor','k')
scatter(repmat(x(2,6), length(tacDir_lMTt), 1),tacDir_lMTt,'MarkerEdgeColor','k');
scatter(repmat(x(3,6), length(both_lMTt), 1),both_lMTt,'MarkerEdgeColor','k');

scatter(repmat(x(1,7), length(visDir_rMTt), 1),visDir_rMTt,'MarkerEdgeColor','k');%,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',1);
scatter(repmat(x(2,7), length(tacDir_rMTt), 1),tacDir_rMTt,'MarkerEdgeColor','k');
scatter(repmat(x(3,7), length(both_rMTt), 1),both_rMTt,'MarkerEdgeColor','k');

hold on
%plot chance level
yline(0.5, '--')

hold off
ylim([0 1])
legend({'trainVisual_testTactile','trainTactile_testVision','both' },'Location','northeast')
xlabel('ROI') 
ylabel('Decoding Accuracy')
xticklabels({'lhMT','rhMT','lS1','lPC', 'rPC', 'lMTt', 'rMTt'})
head1= 'CrossModal'; 
head2=strcat('image=', im,' smoothing=',smooth,' ', ' voxelNb=',voxNb);
title(head1, head2)
saveas(gcf,strcat('decoding_', 'CrossModal','_','image_',im,'_','smoothing',smooth,'_','voxelNb',voxNb,'.png'))
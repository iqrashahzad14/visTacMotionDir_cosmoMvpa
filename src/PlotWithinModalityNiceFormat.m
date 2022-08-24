
im='beta';%'t_maps', 'beta'
smooth='0';
voxNb='100';

subList={'004','005', '006', '007'};
roiList={'lhMT','rhMT','lS1','lPC', 'rPC', 'lMTt', 'rMTt'};
decodingConditionList = {'visual_vertical_vs_visual_horizontal', 'tactile_vertical_vs_tactile_horizontal'};

visDir_lhMT=[]; visDir_rhMT=[]; visDir_lS1=[]; visDir_lPC=[]; visDir_rPC=[]; visDir_lMTt=[]; visDir_rMTt=[];
tacDir_lhMT=[]; tacDir_rhMT=[]; tacDir_lS1=[]; tacDir_lPC=[]; tacDir_rPC=[]; tacDir_lMTt=[]; tacDir_rMTt=[];


for iAccu=1:length(accu)
    for iSub=1:length(subList)
        subID=subList(iSub);
        if strcmp(char({accu(iAccu).subID}.'),char(subID))==1

            if strcmp(char({accu(iAccu).image}.'), im)==1 && strcmp(num2str([accu(iAccu).ffxSmooth].'),smooth)==1 && strcmp(num2str([accu(iAccu).choosenVoxNb].'),voxNb)==1
                
                
                if strcmp(char({accu(iAccu).mask}.'), 'lhMT')==1
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'visual_vertical_vs_visual_horizontal' )==1
                        visDir_lhMT = [visDir_lhMT;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'tactile_vertical_vs_tactile_horizontal' )==1
                        tacDir_lhMT = [tacDir_lhMT;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'rhMT')==1
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'visual_vertical_vs_visual_horizontal' )==1
                        visDir_rhMT = [visDir_rhMT;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'tactile_vertical_vs_tactile_horizontal' )==1
                        tacDir_rhMT = [tacDir_rhMT;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'lS1')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'visual_vertical_vs_visual_horizontal' )==1
                        visDir_lS1 = [visDir_lS1;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'tactile_vertical_vs_tactile_horizontal' )==1
                        tacDir_lS1 = [tacDir_lS1;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'lPC')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'visual_vertical_vs_visual_horizontal' )==1
                        visDir_lPC = [visDir_lPC;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'tactile_vertical_vs_tactile_horizontal' )==1
                        tacDir_lPC = [tacDir_lPC;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'rPC')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'visual_vertical_vs_visual_horizontal' )==1
                        visDir_rPC = [visDir_rPC;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'tactile_vertical_vs_tactile_horizontal' )==1
                        tacDir_rPC = [tacDir_rPC;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'lMTt')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'visual_vertical_vs_visual_horizontal' )==1
                        visDir_lMTt = [visDir_lMTt;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'tactile_vertical_vs_tactile_horizontal' )==1
                        tacDir_lMTt = [tacDir_lMTt;[accu(iAccu).accuracy].'];
                    end
                end
                
                if strcmp(char({accu(iAccu).mask}.'), 'rMTt')==1 
                    varDecodCond={accu(iAccu).decodingCondition}.';
                    if strcmp(varDecodCond{1}{1},'visual_vertical_vs_visual_horizontal' )==1
                        visDir_rMTt = [visDir_rMTt;[accu(iAccu).accuracy].'];
                    elseif strcmp(varDecodCond{1}{1},'tactile_vertical_vs_tactile_horizontal' )==1
                        tacDir_rMTt = [tacDir_rMTt;[accu(iAccu).accuracy].'];
                    end
                end
            end    
       end
            

    end
end

decodAccu=[visDir_lhMT,tacDir_lhMT, visDir_rhMT, tacDir_rhMT,...
    visDir_lS1, tacDir_lS1, visDir_lPC, tacDir_lPC, visDir_rPC, tacDir_rPC,...
    visDir_lMTt, tacDir_lMTt, visDir_rMTt, tacDir_rMTt ];
meanDecodAccu=mean(decodAccu);
seDecodAccu=std(decodAccu)/sqrt(length(subList));

T = array2table(decodAccu,...
    'VariableNames',{'visDir_lhMT','tacDir_lhMT', 'visDir_rhMT', 'tacDir_rhMT',...
    'visDir_lS1', 'tacDir_lS1', 'visDir_lPC', 'tacDir_lPC', 'visDir_rPC', 'tacDir_rPC',...
    'visDir_lMTt', 'tacDir_lMTt', 'visDir_rMTt', 'tacDir_rMTt'});


%% PLOT

% set parameters for plots
shape='o';
%%set the width of the edge of the markers
LineWidthMarkers=1.5;
%%set the width of the edge of the mean line
LineWidthMean=4;
%%set the length of the mean line
LineLength=0.4; %the actual length will be the double of this value

%%%set the color for each condition in RGB (and divide them by 256 to be matlab compatible)
Col_A=[105 170 153]/256; % light green
Col_B=[24 96 88]/256;% dark green
Col_C=[255 158 74]/256; %light orangr
Col_D=[208 110 48]/256;% dark orange
Col_E=[198 131 239]/256;% light purple
Col_F=[121 57 195]/256; % dark purple
Colors=[Col_A;Col_B;Col_C;Col_D;Col_E;Col_F];

%%%set the transparency of the markers
Transparency=1;%0.7;
FontName='Avenir'; %set the style of the labels
FontSize=17; %%set the size of the labels
X_label = 'ROI';
Y_label='Decoding Accuracy'; %%
yLim=[0 1]; %put here your yLim

figure(1)
set(gcf,'color','w');
model_series = [mean(visDir_lhMT),mean(tacDir_lhMT); mean(visDir_rhMT), mean(tacDir_rhMT);...
    mean(visDir_lS1), mean(tacDir_lS1); mean(visDir_lPC), mean(tacDir_lPC); mean(visDir_rPC), mean(tacDir_rPC);...
    mean(visDir_lMTt), mean(tacDir_lMTt); mean(visDir_rMTt), mean(tacDir_rMTt)];

model_error = [std(visDir_lhMT)/sqrt(length(subList)),std(tacDir_lhMT)/sqrt(length(subList)); std(visDir_rhMT)/sqrt(length(subList)), std(tacDir_rhMT)/sqrt(length(subList));...
    std(visDir_lS1)/sqrt(length(subList)), std(tacDir_lS1)/sqrt(length(subList)); std(visDir_lPC)/sqrt(length(subList)), std(tacDir_lPC)/sqrt(length(subList)); std(visDir_rPC)/sqrt(length(subList)), std(tacDir_rPC)/sqrt(length(subList));...
    std(visDir_lMTt)/sqrt(length(subList)), std(tacDir_lMTt)/sqrt(length(subList)); std(visDir_rMTt)/sqrt(length(subList)), std(tacDir_rMTt)/sqrt(length(subList))];

% b = bar(model_series, 'grouped');
b = bar(model_series, 'grouped', 'FaceColor','flat', 'LineWidth',LineWidthMean );

b(1).EdgeColor = Colors(4,:); %group1
b(2).EdgeColor = Colors(2,:); %group2

b(1).CData(1,:) = 'w'; % group 1 1st bar
b(1).CData(2,:) = 'w'; % group 1 2nd bar
b(1).CData(3,:) = 'w'; % group 1 3rd bar
b(1).CData(4,:) = 'w'; % group 1 4th bar
b(1).CData(5,:) = 'w'; % group 1 5th bar
b(1).CData(6,:) = 'w'; % group 1 6th bar
b(1).CData(7,:) = 'w'; % group 1 7th bar

b(2).CData(1,:) = 'w'; % group 2 1st bar
b(2).CData(2,:) = 'w'; % group 2 2nd bar
b(2).CData(3,:) = 'w'; % group 2 3rd bar
b(2).CData(4,:) = 'w'; % group 2 4th bar
b(2).CData(5,:) = 'w'; % group 2 5th bar
b(2).CData(6,:) = 'w'; % group 2 6th bar
b(2).CData(7,:) = 'w'; % group 2 7th bar

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
hold on

sz=60;
scatter(repmat(x(1,1), length(visDir_lhMT), 1),visDir_lhMT,shape,'MarkerEdgeColor',Colors(4,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);%,'MarkerFaceColor','k');
scatter(repmat(x(2,1), length(tacDir_lhMT), 1),tacDir_lhMT,shape,'MarkerEdgeColor',Colors(2,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(1,2), length(visDir_rhMT), 1),visDir_rhMT,shape,'MarkerEdgeColor',Colors(4,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(2,2), length(tacDir_rhMT), 1),tacDir_rhMT,shape,'MarkerEdgeColor',Colors(2,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(1,3), length(visDir_lS1), 1),visDir_lS1,shape,'MarkerEdgeColor',Colors(4,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(2,3), length(tacDir_lS1), 1),tacDir_lS1,shape,'MarkerEdgeColor',Colors(2,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(1,4), length(visDir_lPC), 1),visDir_lPC,shape,'MarkerEdgeColor',Colors(4,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(2,4), length(tacDir_lPC), 1),tacDir_lPC,shape,'MarkerEdgeColor',Colors(2,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(1,5), length(visDir_rPC), 1),visDir_rPC,shape,'MarkerEdgeColor',Colors(4,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(2,5), length(tacDir_rPC), 1),tacDir_rPC,shape,'MarkerEdgeColor',Colors(2,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(1,6), length(visDir_lMTt), 1),visDir_lMTt,shape,'MarkerEdgeColor',Colors(4,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(2,6), length(tacDir_lMTt), 1),tacDir_lMTt,shape,'MarkerEdgeColor',Colors(2,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);
scatter(repmat(x(1,7), length(visDir_rMTt), 1),visDir_rMTt,shape,'MarkerEdgeColor',Colors(4,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);%,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',1);
scatter(repmat(x(2,7), length(tacDir_rMTt), 1),tacDir_rMTt,shape,'MarkerEdgeColor',Colors(2,:),'MarkerEdgeAlpha',Transparency,'LineWidth',LineWidthMarkers);

% scatter(repmat(x(1,1), length(visDir_lhMT), 1),visDir_lhMT,'MarkerEdgeColor','k');%,'MarkerFaceColor','k');
% scatter(repmat(x(2,1), length(tacDir_lhMT), 1),tacDir_lhMT,'MarkerEdgeColor','k');
% scatter(repmat(x(1,2), length(visDir_rhMT), 1),visDir_rhMT,'MarkerEdgeColor','k');
% scatter(repmat(x(2,2), length(tacDir_rhMT), 1),tacDir_rhMT,'MarkerEdgeColor','k');
% scatter(repmat(x(1,3), length(visDir_lS1), 1),visDir_lS1,'MarkerEdgeColor','k');
% scatter(repmat(x(2,3), length(tacDir_lS1), 1),tacDir_lS1,'MarkerEdgeColor','k');
% scatter(repmat(x(1,4), length(visDir_lPC), 1),visDir_lPC,'MarkerEdgeColor','k');
% scatter(repmat(x(2,4), length(tacDir_lPC), 1),tacDir_lPC,'MarkerEdgeColor','k');
% scatter(repmat(x(1,5), length(visDir_rPC), 1),visDir_rPC,'MarkerEdgeColor','k');
% scatter(repmat(x(2,5), length(tacDir_rPC), 1),tacDir_rPC,'MarkerEdgeColor','k');
% scatter(repmat(x(1,6), length(visDir_lMTt), 1),visDir_lMTt,'MarkerEdgeColor','k');
% scatter(repmat(x(2,6), length(tacDir_lMTt), 1),tacDir_lMTt,'MarkerEdgeColor','k');
% scatter(repmat(x(1,7), length(visDir_rMTt), 1),visDir_rMTt,'MarkerEdgeColor','k');%,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',1);
% scatter(repmat(x(2,7), length(tacDir_rMTt), 1),tacDir_rMTt,'MarkerEdgeColor','k');

hold on
% Plot the errorbars
eb=errorbar(x',model_series,model_error,'k','linestyle','none', 'LineWidth',LineWidthMean,'Color','k');
eb(1).Color=Colors(4,:); %group1
eb(2).Color=Colors(2,:); %group2

hold on
%plot chance level
yline(0.5, '--','LineWidth',LineWidthMean/2,'Color','k')

hold off
ylim([0 1])
legend({'vis-VerVsHor','tac-VerVsHor'},'Location','northeast')
legend boxoff 
xlabel(X_label) 
ylabel(Y_label)
XTickLabel={'lhMT','rhMT','lS1','lPC', 'rPC', 'lMTt', 'rMTt'};
xticklabels(XTickLabel)
head1= 'Within-Modality Direction Selectivity'; 
head2=strcat('image=', im,' smoothing=',smooth,' ', ' voxelNb=',voxNb);
title('Within-Modality Direction Selectivity')
% saveas(gcf,strcat('decoding_', 'Within-Modality','_','image_',im,'_','smoothing',smooth,'_','voxelNb',voxNb,'.png'))

ax=gca;
set(ax,'FontName',FontName,'FontSize',FontSize, 'FontWeight','bold',...
    'LineWidth',2.5,'TickDir','out', 'TickLength', [0,0],...
    'yLim',yLim,'XTickLabel', XTickLabel,'FontSize',FontSize);
xlabel(X_label,'FontSize',18,'FontAngle','italic');

ylabel(Y_label,'FontSize',FontSize,'FontAngle','italic');

box off
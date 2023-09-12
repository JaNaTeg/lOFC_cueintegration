% this script averages signal strength of the univariate contrast for compound versus signal from the lOFC and mOFC ROI
    % tests for compound effects and interaction

clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = load('fMRIdata_ROI_univar_compoundsingle.mat'); 
data= fname.data; % voxelwise values for each ROI

nsubs = length(data);
nroi = length(fname.roifile);

for r = 1:nroi
    for cond=1:4 %4 conditions: AX, B, CY, D 
        for s=1:nsubs
            res(s,cond,r) = mean(data{s,r}(cond,:));
        end   
    end
    figure, 
    ah = axes;
    y = res(:,:,r);
    x = [1:4];
    hold on
    boxplot(y);
    ylabel('fMRI response')
    title(strrep(fname.roifile{r}, '_', ' '));
    set(ah, 'xtick', x, 'xticklabel', {'AX', 'B', 'CY','D'})
    
end


%% lOFC
% AX vs B
[h,p,ci,stats] = ttest(res(:,1,1),res(:,2,1))
% %CY vs D
[h,p,ci,stats] = ttest(res(:,3,1),res(:,4,1))
% %Interaction
 [h,p,ci,stats] = ttest((res(:,1,1)-res(:,2,1)),(res(:,3,1)-res(:,4,1)))


%% mOFC
% AX vs B
[h,p,ci,stats] = ttest(res(:,1,2),res(:,2,2))
% %CY vs D
[h,p,ci,stats] = ttest(res(:,3,2),res(:,4,2))
% %Interaction
 [h,p,ci,stats] = ttest((res(:,1,2)-res(:,2,2)),(res(:,3,2)-res(:,4,2)))


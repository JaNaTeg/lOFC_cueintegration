% this script runs a leave-one-run-out SVM classfication analysis on
% outcome expectattions for sweet odor, savory odor and air across all 6
% experimental runs in the lOFC and mOFC ROIs
   % tests for above chance decoding accuracy
   % tests for specificity of sweet and savory odor representati

clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS 
% libsvm-3.25 needs to be downloaded and added to path (addpath())

ROIdata = load('data_ROI_6runs_Sw_Sa_No.mat');

nsubs = length(ROIdata.dataall);
nroi = length(ROIdata.roifile);
nruns = 6;

nconds = 3; %number of conditions in GLM (here = 3: Sw vs sa vs No)
targetlabel = repmat([1;2;3],6,1); %for confusion matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS
acc_decode = zeros(nsubs,nruns+1,length(nroi));
confmx = zeros(3,3,nroi,nsubs);
av_confmx = zeros(3,3,length(nroi));

for roi = 1:nroi
    for sub = 1:nsubs
   
    % GET DATA 
    data= ROIdata.dataall{sub,roi};
    
    
    % ROI decoding leave one out crossvalidation
    
    % define the labels for decoding for each run
      labels = repmat([1,2,3], nruns, 1);
      lPredict_sn=[];

        all=[1:nruns];
        mean_acc = [];
          
     % DO SOME MVPA ANALYSIS WITH CROSSVALIDATION HERE (this is up to you)
        for it = 1:nruns
            
            r_train = all(all~=it); %train on all runs -1
            r_test = it; %test on that run
            
            vectors_train = [];
            vectors_test = [];
            labels_train = [];
            labels_test = [];

            for cond = 1:3
               vectors_train = [vectors_train; squeeze(data(r_train,cond,:))];  % 1 row per condition & run (cond1, run1 & 2 followed by cond 2, run 1&2)
               labels_train = [labels_train; labels(r_train,cond)];   % 1 1 2 2

                vectors_test = [vectors_test, squeeze(data(r_test,cond,:))];
                labels_test = [labels_test; labels(r_test,cond)];
            end

            vectors_test = vectors_test';
            
            % train SVC
            model = svmtrain(labels_train,vectors_train,'-s 0 -t 0 -c 0.1 -q');  % function from LIBSVM toolbox -> needs to be in path!
            
            % test SVC
            %[predicted_label, accuracy, decision_values ] = svmpredict(labels_test, vectors_test, model, '-q');
            
            [labels_predict, accuracy, ~] = svmpredict(labels_test, vectors_test, model, '-q'); 

            acc_decode(sub,it,roi) = accuracy(1);
            lPredict_sn = [lPredict_sn;labels_predict];
            mean_acc = [mean_acc; accuracy(1)];
        end
        
        acc_decode(sub,it+1,roi) = mean(mean_acc); %average over all runs
        Cm_sn = confusionchart(targetlabel,lPredict_sn); %create confusion matrix
        Cm_sn.Normalization= 'row-normalized';          % normalize by row
        confmx(:,:,roi,sub) = Cm_sn.NormalizedValues;            % save normalized values

        %subtract chance
        acc_decode(sub,:,roi) = acc_decode(sub,:,roi) - (100/3);

        clear Cm_sn
    end % all rois
    
end

for ro = 1:nroi
    av_confmx(:,:,ro) = mean(confmx(:,:,ro,:),4)*100;
    figure ()
    % heatmap(av_confmx(:,:,ro),'Colormap',parula,'ColorScaling','scaledrows','ColorbarVisible','on')
    heatmap(av_confmx(:,:,ro),'Colormap',pink,'ColorLimits',[20 60],'ColorbarVisible','on') %'ColorMethod','none' % Fig 2 = lOFC, Fig 3 = mOFC
    %imagesc(av_confmx(:,:,ro))
end

% test whether average decoding accuracy is above chance
% lOFC
[~,p,~,stats] =  ttest(acc_decode(:,end,1))

% mOFC
[~,p,~,stats] =  ttest(acc_decode(:,end,2))

boxplot([acc_decode(:,end,2) acc_decode(:,end,1)])

% test whether representation of reward expectations are specific
R=confmx(1:2,1:2,:,:);
R_lofc_diag =  squeeze(mean( [R(1,1,1,:), R(2,2,1,:)]));
R_lofc_offdiag =  squeeze(mean( [R(1,2,1,:), R(2,1,1,:)]));
specific_acc_lOFC = R_lofc_diag - R_lofc_offdiag;

R_mofc_diag = squeeze( mean( [R(1,1,2,:), R(2,2,2,:)]));
R_mofc_offdiag =  squeeze(mean( [R(1,2,2,:), R(2,1,2,:)]));
specific_acc_mOFC = R_mofc_diag - R_mofc_offdiag;

[~,p,~,stats] = ttest(R_lofc_diag,R_lofc_offdiag,"Tail","right")
[~,p,~,stats] = ttest(R_mofc_diag,R_mofc_offdiag,"Tail","right")

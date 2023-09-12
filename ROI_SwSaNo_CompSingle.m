
% this script runs a SVM classification in the lOFC and mOFC ROIs to decode
% outcome expectations separately for compound and single cue conditions
    % resulting accuracies are compared
    % compound decoding is tested for specificity

clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS 
% libsvm-3.25 needs to be downloaded and added to path
% addpath('');

compdata = load('fMRIdata_ROI_SwSaNo_compound.mat');
singledata = load('fMRIdata_ROI_SwSaNo_single.mat');

nsubs = length(compdata.dataall_train);
nroi = length(compdata.roifile);

nruns_train = 6;
nruns_test = 4;
nconds = 3; %number of conditions in GLM (here = 2: Sw vs Sa, A1 vs A2)
targetlabel = repmat([1;2;3],4,1); %for confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% train sweet/savory/no, test compound

% results
  acc_decode_comp = zeros(nsubs,nruns_test+1,nroi);
   confmx_comp = zeros(3,3,nroi,nsubs);
   av_confmx_comp = zeros(3,3,nroi);

for sub = 1:nsubs    
   
  %% ROI
   for roi = 1:nroi
  
    % GET DATA 
    data_train= compdata.dataall_train{sub,roi};
    data_test= compdata.dataall_test{sub,roi};

    % ROI decoding leave one run out cross validation
      
    % define the labels for decoding for each run
      labels = repmat([1,2,3], nruns_train, 1);
       lPredict_sn=[];

      all=[1:nruns_train];
      mean_acc = [];
      
      % MVPA ANALYSIS WITH CROSSVALIDATION H
      
        for it = 1:nruns_test
            r_train = all(all~=(it+2)); %train on all runs -1
            r_test = it; %test on that run
            
            vectors_train = [];
            vectors_test = [];
            labels_train = [];
            labels_test = [];

            for cond = 1:3
               vectors_train = [vectors_train; squeeze(data_train(r_train,cond,:))];  % 1 row per condition & run (cond1, run1 & 2 followed by cond 2, run 1&2)
               labels_train = [labels_train; labels(r_train,cond)];   % 1 1 2 2

                vectors_test = [vectors_test, squeeze(data_test(r_test,cond,:))];
                labels_test = [labels_test; labels(r_test,cond)];
            end

             vectors_test = vectors_test';
             
            % train SVC
            model = svmtrain(labels_train,vectors_train,'-s 0 -t 0 -c 0.1 -q');  % function from LIBSVM toolbox -> needs to be in path!
            
            % test SVC
            [labels_predict, accuracy, ~] = svmpredict(labels_test, vectors_test, model, '-q'); 

            acc_decode(sub,it,roi) = accuracy(1);
            mean_acc = [mean_acc; accuracy(1)];   
            lPredict_sn = [lPredict_sn;labels_predict];
        end
%         
            acc_decode_comp(sub,it+1,roi) = mean(mean_acc); %average over all runs
            Cm_sn = confusionchart(targetlabel,lPredict_sn); %create confusion matrix
            Cm_sn.Normalization= 'row-normalized';          % normalize by row
            confmx_comp(:,:,roi,sub) = Cm_sn.NormalizedValues;            % save normalized values

        clear Cm_sn
   end
   
   %subtract chance
     acc_decode_comp(sub,:,:) = acc_decode_comp(sub,:,:) - (100/3);
     
end

for ro = 1:nroi
    av_confmx_comp(:,:,ro) = mean(confmx_comp(:,:,ro,:),4)*100;
    figure ()
    % heatmap(av_confmx(:,:,ro),'Colormap',parula,'ColorScaling','scaledrows','ColorbarVisible','on')
    %heatmap(av_confmx(:,:,ro),'Colormap',pink,'ColorbarVisible','off')
    heatmap(av_confmx_comp(:,:,ro),'Colormap',pink,'ColorLimits',[20 60],'ColorbarVisible','on') 
    %imagesc(av_confmx(:,:,ro))
end

R=confmx_comp(1:2,1:2,:,:);
R_lofc_diag =  squeeze(mean( [R(1,1,1,:), R(2,2,1,:)]));
R_lofc_offdiag =  squeeze(mean( [R(1,2,1,:), R(2,1,1,:)]));

R_mofc_diag = squeeze( mean( [R(1,1,2,:), R(2,2,2,:)]));
R_mofc_offdiag =  squeeze(mean( [R(1,2,2,:), R(2,1,2,:)]));

comp_diff_R_lofc = R_lofc_diag-R_lofc_offdiag;
comp_diff_R_mofc = R_mofc_diag-R_mofc_offdiag;


[~,p,~,stats] = ttest(R_lofc_diag,R_lofc_offdiag,"Tail","right") % specific outcome expectations for compounds in lOFC
[~,p,~,stats] = ttest(R_mofc_diag,R_mofc_offdiag, "Tail","right") % specific outcome expectations for compounds in mOFC


clear data_train data_test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% train sweet/savory/no, test single

% results
  acc_decode_single = zeros(nsubs,nruns_test+1,nroi);

for sub = 1:nsubs    
   
  %% ROI
   for roi = 1:nroi
  
    % GET DATA 
    data_train= singledata.dataall_train{sub,roi};
    data_test= singledata.dataall_test{sub,roi};

    % ROI decoding leave one run out cross validation
    % define the labels for decoding for each run
      labels = repmat([1,2,3], nruns_train, 1);
      lPredict_sn=[];

      all=[1:nruns_train];
      mean_acc = [];

         for it = 1:nruns_test
            
            r_train = all(all~=(it+2)); %train on all runs -1
            r_test = it; %test on that run
            
            vectors_train = [];
            vectors_test = [];
            labels_train = [];
            labels_test = [];

            for cond = 1:3
               vectors_train = [vectors_train; squeeze(data_train(r_train,cond,:))];  % 1 row per condition & run (cond1, run1 & 2 followed by cond 2, run 1&2)
               labels_train = [labels_train; labels(r_train,cond)];   % 1 1 2 2

                vectors_test = [vectors_test, squeeze(data_test(r_test,cond,:))];
                labels_test = [labels_test; labels(r_test,cond)];
            end

            vectors_test = vectors_test';
            
            % train SVC
            model = svmtrain(labels_train,vectors_train,'-s 0 -t 0 -c 0.1 -q');  % function from LIBSVM toolbox -> needs to be in path!
            
            % test SVC
            [labels_predict, accuracy, ~] = svmpredict(labels_test, vectors_test, model, '-q'); 

            acc_decode_single(sub,it,roi) = accuracy(1);
            mean_acc = [mean_acc; accuracy(1)];
           
        end
        
       acc_decode_single(sub,it+1,roi) = mean(mean_acc); %average over all runs
   end
   
   %subtract chance
     acc_decode_single(sub,:,:) = acc_decode_single(sub,:,:) - (100/3);
         
end

boxplot([acc_decode_comp(:,end,2),acc_decode_single(:,end,2),acc_decode_comp(:,end,1),acc_decode_single(:,end,2)]) 

% compound versus single in mOFC
[~,p,~,stats] = ttest(acc_decode_comp(:,end,2),acc_decode_single(:,end,2),"Tail","right")
% compound versus single in lOFC
[~,p,~,stats] = ttest(acc_decode_comp(:,end,1),acc_decode_single(:,end,1),"Tail","right")


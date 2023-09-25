library (Hmisc)
library (plotrix)
library (lme4)
library (afex)

# set directory if necessary
# setwd("/Users...")

######################################################## Experiment 1 - fMRI ##################################################################

# this section contains analyses of the behavioral data of the fMRI experiment: cond = last run of conditioning (before TMS); compound = compound phase after TMS (4 runs);

# data frame variables:
# rew: 1 = odor reward; 0 = no odor / air
# cue: 2 = compound; 1 = single 
# acc = % correct predictions
# acc_diff = prediction accuracy difference between compound and conditioning
# rt = average response time
# rt_diff = response time difference between compound and conditioning

df_fmri_pleas = read.table("fMRI_behavdata_pleasantness.txt", header=T )
df_fmri_cond = read.table("fMRI_behavdata_conditioning.txt", header=T )
df_fmri_comp = read.table("fMRI_behavdata_compound.txt", header=T )
df_fmri_specrewexp = read.table("fMRI_specific_rew_exp.txt", header=T )

# some variables
VP = unique(df_fmri_cond$Subj) 
N = length(VP)

### 1. odor pleasantness ratings ###
boxplot(df_fmri_pleas, ylab="pleasantness", horizontal=F) 
stripchart(data.frame(df_fmri_pleas),            # Data
           method = "jitter", # Random noise
           pch = 20,          # Pch symbols
           vertical = T,   # Vertical mode
           add = TRUE)        # Add it over

# odors versus air
t.test(rowMeans(df_fmri_pleas[1:2]),df_fmri_pleas$air,paired=T)

# sweet versus savory odor
t.test(df_fmri_pleas$sweet,df_fmri_pleas$savory,paired=T)


### 2. performance in last run of conditioning ###
df_fmri_cond$Subj = as.factor(df_fmri_cond$Subj)
df_fmri_cond$rew = as.factor(df_fmri_cond$rew)
df_fmri_cond$cue = as.factor(df_fmri_cond$cue)

bp=boxplot(cbind(df_fmri_cond$acc[df_fmri_cond$cond==1],df_fmri_cond$acc[df_fmri_cond$cond==2],df_fmri_cond$acc[df_fmri_cond$cond==3],df_fmri_cond$acc[df_fmri_cond$cond==4]), names=c('AX','B','CY','D'), ylim = c(0,100))
#points(seq(bp$n), colMeans(cbind(df_fmri_cond$acc[df_fmri_cond$cond==1],df_fmri_cond$acc[df_fmri_cond$cond==2],df_fmri_cond$acc[df_fmri_cond$cond==3],df_fmri_cond$acc[df_fmri_cond$cond==4])), col = "black", pch = 18)
stripchart(data.frame(cbind(df_fmri_cond$acc[df_fmri_cond$cond==1],df_fmri_cond$acc[df_fmri_cond$cond==2],df_fmri_cond$acc[df_fmri_cond$cond==3],df_fmri_cond$acc[df_fmri_cond$cond==4])),            # Data
           method = "jitter", # Random noise
           pch = 20,          # Pch symbols
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over


### 3. performance during compound 
# factors
df_fmri_comp$Subj = as.factor(df_fmri_comp$Subj)
df_fmri_comp$rew = as.factor(df_fmri_comp$rew)
df_fmri_comp$cue = as.factor(df_fmri_comp$cue)

## average change in accuracy per condition (across runs)
comp_acc_subj=tapply(df_fmri_comp$acc_diff,data.frame(df_fmri_comp$Subj,df_fmri_comp$cond),mean)

bp=boxplot(comp_acc_subj, names=c('AX','B','CY','D'))
# points(seq(bp$n), colMeans(comp_acc_subj), col = "black", pch = 18)
stripchart(data.frame(comp_acc_subj),            # Data
           method = "jitter", # Random noise
           pch = 20,          # Pch symbols
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over


## average change in prediction accuracy over time (4 runs of compound conditioning)
diff_comp_acc_mean = tapply(df_fmri_comp$acc_diff, data.frame(df_fmri_comp$cond,df_fmri_comp$run),mean )
diff_comp_acc_stderr = tapply(df_fmri_comp$acc_diff, data.frame(df_fmri_comp$cond,df_fmri_comp$run), std.error )

plot(diff_comp_acc_mean[1,],type='o',col='firebrick4',ylim=c(-5,25),lwd=3,ylab="accuracy (%)",xlab="run")
errbar(c(1:4),diff_comp_acc_mean[1,],diff_comp_acc_mean[1,]+diff_comp_acc_stderr[1,],diff_comp_acc_mean[1,]-diff_comp_acc_stderr[1,],add=T,cap=0,pch=1,errbar.col='firebrick4',col='firebrick4')
lines(diff_comp_acc_mean[2,],col="indianred2",lwd=3,type='o')
errbar(c(1:4),diff_comp_acc_mean[2,],diff_comp_acc_mean[2,]+diff_comp_acc_stderr[2,],diff_comp_acc_mean[2,]-diff_comp_acc_stderr[2,],add=T,cap=0,pch=1,errbar.col="indianred2",col="indianred2")
lines(diff_comp_acc_mean[3,],col='darkgreen',lwd=3,type='o')
errbar(c(1:4),diff_comp_acc_mean[3,],diff_comp_acc_mean[3,]+diff_comp_acc_stderr[3,],diff_comp_acc_mean[3,]-diff_comp_acc_stderr[3,],add=T,cap=0,pch=1,errbar.col='darkgreen',col='darkgreen')
lines(diff_comp_acc_mean[4,],col='darkseagreen3',lwd=3,type='o')
errbar(c(1:4),diff_comp_acc_mean[4,],diff_comp_acc_mean[4,]+diff_comp_acc_stderr[4,],diff_comp_acc_mean[4,]-diff_comp_acc_stderr[4,],add=T,cap=0,pch=1,errbar.col='darkseagreen3',col='darkseagreen3')
legend (1,25,c('AX+','B+','CY-','D-'),box.lty=0,ncol=4,
        fill=c('firebrick4','indianred2','darkgreen','darkseagreen3'))

### linear models ###
model.diff_acc_cond =lmer(acc_diff ~ 1 + rew*cue+run + (1|Subj), data = df_fmri_comp)
summary(model.diff_acc_cond)

# separate models for odor and no odor
model.diff_acc_condAXB =lmer(acc_diff ~ 1 + cue+run + (1|Subj), data = df_fmri_comp[df_fmri_comp$rew==1,])
summary(model.diff_acc_condAXB)

model.diff_acc_condCYD =lmer(acc_diff ~ 1 + cue+run + (1|Subj), data = df_fmri_comp[df_fmri_comp$rew==2,])
summary(model.diff_acc_condCYD)

# response times
model.diff_rt_cond =lmer(rt_diff ~ 1 + rew*cue+run + (1|Subj), data = df_fmri_comp)
summary(model.diff_rt_cond)

model.diff_rt_condAXB =lmer(rt_diff ~ 1 + cue+run + (1|Subj), data = df_fmri_comp[df_fmri_comp$rew==1,])
summary(model.diff_rt_condAXB)
model.diff_rt_condCYD =lmer(rt_diff ~ 1 + cue+run + (1|Subj), data = df_fmri_comp[df_fmri_comp$rew==2,])
summary(model.diff_rt_condCYD)

# specific reward expectations
diff_pred_acc_AX = (df_fmri_specrewexp$compAX_correct - df_fmri_specrewexp$condAX_correct) - (df_fmri_specrewexp$compAX_false - df_fmri_specrewexp$condAX_false)
diff_pred_acc_B = (df_fmri_specrewexp$compB_correct - df_fmri_specrewexp$condB_correct) - (df_fmri_specrewexp$compB_false - df_fmri_specrewexp$condB_false)

t.test(diff_pred_acc_AX ,diff_pred_acc_B, paired =T)

######################################################## Experiment 2 - TMS ##################################################################

# this section contains analyses of the behavioral data of the TMS experiment: cond = last run of conditioning (before TMS); compound = compound phase after TMS (4 runs);

# Stim = active stimulation group; Sham = control group 
# data frame variables:
  # rew: 1 = odor reward; 2 = no odor / air
  # cue: 2 = compound; 1 = single 
  # group: 1 = Stim; 0 = Sham
  # sex: 1 = male; 2 = female
  # acc = % correct predictions
  # acc_diff = prediction accuracy difference between compound and conditioning
  # rt = average response time
  # rt_diff = response time difference between compound and conditioning

df_tms_pleas = read.table("TMS_behavdata_pleasantness.txt", header=T )
df_tms_cond = read.table("TMS_behavdata_conditioning.txt", header=T )
df_tms_comp = read.table("TMS_behavdata_compound.txt", header=T )

# some variables
VPtms = unique(df_tms_comp$Subj)
Ntms = length(VPtms)
active_stim = c(2,4,5,6,7,12,13,16,21,24,27,33,38,40,41,44,46,61,63,65,66,72,74,77,78,79) #IDs of active Stim group

`%notin%` <- Negate(`%in%`)
sg = which(VPtms%in%active_stim)
cg = which(VPtms%notin%active_stim)

# 1. pleasantness

aov.pleas = data.frame(rep(df_tms_pleas$subj,3),rep(df_tms_pleas$stim,3),rep(1:3, each=Ntms),c(df_tms_pleas$sweet,df_tms_pleas$savory,df_tms_pleas$air))
colnames(aov.pleas) = c('Subj','group','cond','rating')

aov.pleas$Subj = as.factor(aov.pleas$Subj)
aov.pleas$group = as.factor(aov.pleas$group)
aov.pleas$cond = as.factor(aov.pleas$cond)

aov_groupxcond = aov(rating ~ (group*cond) + Error(Subj/cond), data = aov.pleas)
summary(aov_groupxcond)


### 2. performance in last run of conditioning 

# test for overall differences in accuracy after conditioning
av_acc = tapply(df_tms_cond$acc,df_tms_cond$Subj,mean)
t.test(av_acc[sg],av_acc[cg])
df_tms_cond$Subj = as.factor(df_tms_cond$Subj)
df_tms_cond$group = as.factor(df_tms_cond$group)

# factors
df_tms_cond$Subj = as.factor(df_tms_cond$Subj)
df_tms_cond$group = as.factor(df_tms_cond$group)
df_tms_cond$rew = as.factor(df_tms_cond$rew)
df_tms_cond$cue = as.factor(df_tms_cond$cue)

# ANOVAs to test for differences in performance at the end of conditioning
res.aov_accO= aov(acc ~ (group*cue) + Error(Subj/cue), data = df_tms_cond[df_tms_cond$rew==1,]) # odors only
summary(res.aov_accO)
res.aov_accNO= aov(acc ~ (group*cue) + Error(Subj/cue), data = df_tms_cond[df_tms_cond$rew==0,]) # no odors only
summary(res.aov_accNO)

res.aov_acc= aov(acc ~ (group*rew*cue) + Error(Subj/(rew*cue)), data = df_tms_cond) # prediction accuracy
summary(res.aov_acc)

res.aov_rt= aov(rt ~ (group*rew*cue) + Error(Subj/(rew*cue)), data = df_tms_cond) # response time
summary(res.aov_rt)


### 3. performance during compound
# factors
df_tms_comp$Subj = as.factor(df_tms_comp$Subj)
df_tms_comp$rew = as.factor(df_tms_comp$rew)
df_tms_comp$cue = as.factor(df_tms_comp$cue)
df_tms_comp$group = as.factor(df_tms_comp$group)

# average change in accuracy per condition 
comp_acc_subj_tms=tapply(df_tms_comp$acc_diff,data.frame(df_tms_comp$Subj,df_tms_comp$cue,df_tms_comp$rew),mean)

par(mfrow=c(1,2))
bp=boxplot(cbind(comp_acc_subj_tms[cg,2:1,2],comp_acc_subj_tms[cg,2:1,1]), names=c('AX','B','CY','D'), xlab = 'Sham')
points(seq(bp$n), colMeans(cbind(comp_acc_subj_tms[cg,2:1,2],comp_acc_subj_tms[cg,2:1,1])), col = "black", pch = 18)

bp=boxplot(cbind(comp_acc_subj_tms[sg,2:1,2],comp_acc_subj_tms[sg,2:1,1]), names=c('AX','B','CY','D'), xlab = 'Stim')
points(seq(bp$n), colMeans(cbind(comp_acc_subj_tms[sg,2:1,2],comp_acc_subj_tms[sg,2:1,1])), col = "black", pch = 18)


### linear models to test influence of stimulation on prediction accuracy 

model.diff_acc_cond_group =lmer(acc_diff ~ 1 + Group*rew*cue+run + (1|Subj), data = df_tms_comp)
summary(model.diff_acc_cond_group)

# separate models for odor predictions and air predictions
model.diff_acc_condAXB =lmer(acc_diff ~ 1 + group*cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$rew==1,]) # odors
summary(model.diff_acc_condAXB)

model.diff_acc_condCYD =lmer(acc_diff ~ 1 + group*cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$rew==0,]) # no odors
summary(model.diff_acc_condCYD)

# separate models for Stim and Sham
model.diff_acc_Stim =lmer(acc_diff ~ 1 + rew*cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$group==1,])
summary(model.diff_acc_Stim)
model.diff_acc_StimcondAXB =lmer(acc_diff ~ 1 + cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$group==1 & df_tms_comp$rew==1,])
summary(model.diff_acc_StimcondAXB)
model.diff_acc_StimcondCYD =lmer(acc_diff ~ 1 + cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$group==1 & df_tms_comp$rew==0,])
summary(model.diff_acc_StimcondCYD)

model.diff_acc_Sham =lmer(acc_diff ~ 1 + rew*cue+run + (1|Subj), data = df_tms_comp[df_tms_comp$group==0,])
summary(model.diff_acc_Sham)
model.diff_acc_ShamcondAXB =lmer(acc_diff ~ 1 + cue+run + (1|Subj), data = df_tms_comp[df_tms_comp$group==0 & df_tms_comp$rew==1,])
summary(model.diff_acc_ShamcondAXB)
model.diff_acc_ShamcondCYD =lmer(acc_diff ~ 1 + cue+run + (1|Subj), data = df_tms_comp[df_tms_comp$group==0 & df_tms_comp$rew==0,])
summary(model.diff_acc_ShamcondCYD)

### linear models to test influence of stimulation on prediction speed

model.diff_rt_cond_group =lmer(rt_diff ~ 1 + group*rew*cue+run + (1|Subj), data = df_tms_comp)
summary(model.diff_rt_cond_group)

# separate models for odor and no odor
model.diff_rt_condAXB =lmer(rt_diff ~ 1 + group*cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$rew==1,])
summary(model.diff_rt_condAXB)

model.diff_rt_condCYD =lmer(rt_diff ~ 1 + group*cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$rew==0,])
summary(model.diff_rt_condCYD)

# separate models for Stim and Sham
model.diff_rt_Stim =lmer(rt_diff ~ 1 + rew*cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$group==1,])
summary(model.diff_rt_Stim)
model.diff_rt_StimcondAXB =lmer(rt_diff ~ 1 + cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$group==1 & df_tms_comp$rew==1,])
summary(model.diff_rt_StimcondAXB)
model.diff_rt_StimcondCYD =lmer(rt_diff ~ 1 + cue+run + (1|Subj), data =df_tms_comp[df_tms_comp$group==1 & df_tms_comp$rew==0,])
summary(model.diff_rt_StimcondCYD)

model.diff_rt_Sham =lmer(rt_diff ~ 1 + rew*cue+run + (1|Subj), data = df_tms_comp[df_tms_comp$group==0,])
summary(model.diff_rt_Sham)
model.diff_rt_ShamcondAXB =lmer(rt_diff ~ 1 + cue+run + (1|Subj), data = df_tms_comp[df_tms_comp$group==0 & df_tms_comp$rew==1,])
summary(model.diff_rt_ShamcondAXB)
model.diff_rt_ShamcondCYD =lmer(rt_diff ~ 1 + cue+run + (1|Subj), data = df_tms_comp[df_tms_comp$group==0 & df_tms_comp$rew==0,])
summary(model.diff_rt_ShamcondCYD)



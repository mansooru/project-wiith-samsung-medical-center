%% Longitudinal analysis - Linear Mixed model
% LME modeling for cortical thinning for parkinson's disease
% based on freesurfer
% 
% created : Aug-28-2018
% Author: Mansukim, mansooru@skku.edu or mansooru.kim@gmail.com
% @ Sungkyunkwan University

clear all; clc;
hemi = 'rh';

input_dir = ['/store6/mskim/project/SMC/01_CT.analysis/p1_02_freesurfer/' ];
output_dir = '/store6/mskim/project/SMC/01_CT.analysis/p1_freesurfer/';

command1 = ['mris_preproc --qdec-long qdec.table.dat --target fsaverage --hemi ' hemi ' --meas thickness --out ' input_dir hemi '.thickness.mgh  --SUBJECTS_DIR ' input_dir];
command2 = ['mri_surf2surf --hemi ' hemi ' --s fsaverage --sval ' input_dir  hemi '.thickness.mgh --tval ' input_dir hemi '.thickness.fwhm10.fsaverage.mgh --fwhm-trg 10 --cortex  --noreshape'];
disp(command1); system(command1);
disp(command2); system(command2);

% lme modeling
clear all; clc;
hemi = 'rh';

input_dir = ['/store6/mskim/project/SMC/01_CT.analysis/p1_02_freesurfer/' ];
cd(input_dir)
[Y,mri] = fs_read_Y([hemi '.thickness.fwhm10.fsaverage.mgh']);

Qdec = fReadQdec('qdec.table.dat');
Qdec = rmQdecCol(Qdec,1);
sID = Qdec(2:end,1);
Qdec = rmQdecCol(Qdec,1);
M = Qdec2num(Qdec);
[M,Y,ni] = sortData(M,1,Y,sID);

lhsphere = fs_read_surf('/usr/local/freesurfer/subjects/fsaverage/surf/lh.sphere');
lhcortex = fs_read_label('/usr/local/freesurfer/subjects/fsaverage/label/lh.cortex.label');

% Design matrix X
% Intercept, time, group, time x group, sided, gender, age, MMSE, ICV
% We will consider, time x group and sided symptom information
X = [ones(length(M),1) M];

[lhTh0,lhRe] = lme_mass_fit_EMinit(X,[1 2],Y,ni,lhcortex,3,8);
[lhRgs,lhRgMeans] = lme_mass_RgGrow(lhsphere,lhRe,lhTh0,lhcortex,2,95);

% Check similarity between two esitimated reions
surf.faces =  lhsphere.tri;
surf.vertices =  lhsphere.coord';

figure; p1 = patch(surf);set(p1,'facecolor','interp','edgecolor','none','facevertexcdata',lhTh0(1,:)');
figure; p2 = patch(surf); set(p2,'facecolor','interp','edgecolor','none','facevertexcdata',lhRgMeans(1,:)');

% Spatiotemporal LME model fitting
lhstats = lme_mass_fit_Rgw(X,[1 2],Y,ni,lhTh0,lhRgs,lhsphere);

% Contrast matrix 
CM.C = [0 0 1 1 0 0 0 0];

F_lhstats = lme_mass_F(lhstats,CM,5);

nv=length(lhstats);
Beta2 = zeros(1,nv);
for i=1:nv
if ~isempty(lhstats(i).Bhat)
      Beta2(i) = lhstats(i).Bhat(4);
   end;
end;

mri1 = mri;
mri1.volsz(4) = 1;
fs_write_Y(Beta2,mri1,['long.beta.' hemi '.mgh' ]);

[detvtx,sided_pval,pth] = lme_mass_FDR2(F_lhstats.pval,F_lhstats.sgn,lhcortex,0.05,-1);

sided_pval = sided_pval*1000;
s_idx = (sided_pval > 1);
sided_pval(s_idx) = 1;

pcor = -log10(sided_pval);
fs_write_Y(pco,mri1,['long.fdr.pval.' hemi '.mgh' ] );

idx = (sided_pval <0.001);
nv=length(lhstats);
sig_beta = zeros(1,nv);
sig_beta = Beta2.*idx;

fs_write_Y(sig_beta,mri1,['long.sbeta.' hemi '.mgh' ]);


% m1 = [161:2:252];
% m2 = [162:2:252];
% 
% unit_diff = (Y(m2,:)-Y(m1,:))./X(m2,2);
% sig_unit_diff = unit_diff.*repmat(idx,45,1);
% fs_write_Y(mean(sig_unit_diff,1),mri1,'diff_ct_lh.mgh');
% 
% stat_unit_diff = reshape(unit_diff(:,idx),1,[]);
% tmp = mean(sig_unit_diff,1);
% tmp_idx = (tmp~=0);
% tmp(tmp_idx);
% hist(tmp(tmp_idx),100);
% 
% fs_write_Y(r,mri1,'r_rh.mgh');



%% Visulization
freeview -f /Applications/freesurfer/subjects/fsaverage/surf/lh.inflated:overlay=~/Documents/Project/SMC/01_result_lh/spval_sided.mgh:overlay_threshold=1.3,5:overlay=~/Documents/Project/SMC/01_result_lh/beta.mgh:overlay_threshold=0.0001,0.01


%%
clear all; clc;

[Y,mri] = fs_read_Y('rh.thickness_sm10.mgh');

Qdec = fReadQdec('qdec.table.dat');
Qdec = rmQdecCol(Qdec,1);
sID = Qdec(2:end,1);
Qdec = rmQdecCol(Qdec,1);
M = Qdec2num(Qdec);
[M,Y,ni] = sortData(M,1,Y,sID);

rhsphere = fs_read_surf('/usr/local/freesurfer/subjects/fsaverage/surf/rh.sphere');
rhcortex = fs_read_label('/usr/local/freesurfer/subjects/fsaverage/label/rh.cortex.label');

% Design matrix X
% Intercept, time, group, time x group, sided, gender, age, MMSE, ICV
% We will consider, time x group and sided symptom information
X = [ones(length(M),1) M];

[lhTh0,lhRe] = lme_mass_fit_EMinit(X,[1 2],Y,ni,rhcortex,3,30);
[lhRgs,lhRgMeans] = lme_mass_RgGrow(rhsphere,lhRe,lhTh0,rhcortex,2,95);

% Check similarity between two esitimated reions
surf.faces =  rhsphere.tri;
surf.vertices =  rhsphere.coord';

figure; p1 = patch(surf);set(p1,'facecolor','interp','edgecolor','none','facevertexcdata',lhTh0(1,:)');
figure; p2 = patch(surf); set(p2,'facecolor','interp','edgecolor','none','facevertexcdata',lhRgMeans(1,:)');

% Spatiotemporal LME model fitting
lhstats = lme_mass_fit_Rgw(X,[1 2],Y,ni,lhTh0,lhRgs,rhsphere);

% Contrast matrix 
CM.C = [0 0 1 1 0 0 0 0 0];

F_lhstats = lme_mass_F(lhstats,CM,30);

nv=length(lhstats);
Beta2 = zeros(1,nv);
for i=1:nv
if ~isempty(lhstats(i).Bhat)
      Beta2(i) = lhstats(i).Bhat(4);
   end;
end;

mri1 = mri;
mri1.volsz(4) = 1;
fs_write_Y(Beta2,mri1,'beta_rh.mgh');

[detvtx,sided_pval,pth] = lme_mass_FDR2(F_lhstats.pval,F_lhstats.sgn,rhcortex,0.05,0);

sided_pval = sided_pval.*50;
s_idx = (sided_pval > 1);
sided_pval(s_idx) = 1;

pcor = -log10(sided_pval);
fs_write_Y(pcor,mri1,'spval_rh.mgh');

idx = (sided_pval <0.001);
nv=length(lhstats);
sig_beta = zeros(1,nv);
sig_beta = Beta2.*idx;

fs_write_Y(sig_beta,mri1,'sig_beta_rh.mgh');



m1 = [161:2:252];
m2 = [162:2:252];

unit_diff = (Y(m2,:)-Y(m1,:))./X(m2,2);
sig_unit_diff = unit_diff.*repmat(idx,45,1);
fs_write_Y(mean(sig_unit_diff,1),mri1,'diff_ct_rh.mgh');

stat_unit_diff = reshape(unit_diff(:,idx),1,[]);
tmp = mean(sig_unit_diff,1);
tmp_idx = (tmp~=0);
tmp(tmp_idx);
hist(tmp(tmp_idx),100);




%%
[rhTh0,rhRe] = lme_mass_fit_EMinit(X,[1 2],Y,ni,rhcortex,3);
[rhRgs,rhRgMeans] = lme_mass_RgGrow(rhsphere,rhRe,rhTh0,rhcortex,2,95);

% Check similarity between two esitimated reions
surf.faces =  rhsphere.tri;
surf.vertices =  rhsphere.coord';

figure; p1 = patch(surf);set(p1,'facecolor','interp','edgecolor','none','facevertexcdata',rhTh0(1,:)');
figure; p2 = patch(surf); set(p2,'facecolor','interp','edgecolor','none','facevertexcdata',rhRgMeans(1,:)');

% Spatiotemporal LME model fitting

rhstats = lme_mass_fit_Rgw(X,[1 2],Y,ni,rhTh0,rhRgs,rhsphere);
CM.C = [0 0 1 1 0 0 0 0 0];

F_rhstats = lme_mass_F(rhstats,CM);

nv=length(rhstats);
Beta2 = zeros(1,nv);
for i=1:nv
if ~isempty(rhstats(i).Bhat)
      Beta2(i) = rhstats(i).Bhat(4);
   end;
end;

mri1 = mri;
mri1.volsz(4) = 1;
fs_write_Y(Beta2,mri1,'beta.mgh');

[detvtx,sided_pval,pth] = lme_mass_FDR2(F_rhstats.pval,F_rhstats.sgn,rhcortex,0.05,0);

pcor = -log10(sided_pval);
fs_write_Y(pcor,mri1,'spval_rh.mgh');

%%
clear all; clc;
load('avg_sig_jac.mat');
cov = csvread('Covariates_vols.csv',1,1);
avg_sig_jac(cov(:,2) ==0,:) = [];
cov(cov(:,2) ==0,:)=[];

Y = cov(:,1);
X = [10.^avg_sig_jac cov(:,3:8)];

for i = 1:31
    
    Xtrain = X;Ytrain = Y;
    Xtrain(i,:) = []; Ytrain(i,:) = [];
    
    Xtest = X(i,:); Ytest= Y(i,:);
    
    % 1st multiple linear regression model
    mdl=fitglm(Xtrain,Ytrain);
    predict_score(i) = predict(mdl,Xtest);
    
end
errors = predict_score - Y';

    
plot(predict_score,Y','o')
performance = rms(errors)
corr = corr(predict_score',Y)
xlabel('Predicted UdysRS'); ylabel('Actural UdysRS'); grid on; axis square; xlim([0 50])















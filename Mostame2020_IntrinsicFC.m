%% ----------------------------------------------- initialize
clear all; clc
% main directory
path_main='...'; cd(path_main);
% Personal functions directory
path_func='...\functions'; addpath(path_func);
% Raw data directory
path_root='...'; addpath(path_root);
% Fieldtrip directory
addpath( '...' );
ft_defaults;
% results path
path_results_rest=['...'];
path_results_task=['...'];
load Tasks

%% Call main function
Continuous_analysis = 0;
Update = 0;
for i=6
    Task=Tasks{i}; % for tasks with different variations, check if data/electrode files are loaded correctly
    for i_sub=[23]
        for freq=1 : 5
            clc
            i_sub
            freq
            if Update == 1
                Update_Files(i_sub, freq, 'Rest');
            else
                if Continuous_analysis
                    Main_rest(path_main, path_root, path_results_rest, freq, i_sub,Task);
                else
                    Main_tasks(path_main, path_root, path_results_task, freq, i_sub, Task);
                end
            end
        end
    end
end

%% ----------------------------------------------- Main function for rest (or any continous) data
function [] = Main_rest(path_main, path_root, path_results_rest, freq, i_sub, Task)
%% ----------------- load data & extract configurations
load Subjects
subject=Subjects{i_sub};
if i_sub<=10 % ------------------------------ Stanford data exclusion
    task='RestEyesOpen';
    if i_sub==10
        task='RestEyesClosed';
    end
    [~, edata, electrodes_coordinate, BadElecs]=load_clean_data(path_root,subject,task,1);
else
    cd([path_root '\' subject])
    load(sprintf('edata_Stanford_%s_%s.mat',Task, subject))
end % ------------------------------ Stanford data exclusion
% Extract Data specifics
numelec=numel(edata.label);
Fs=edata.fsample;

%% calculate electrode distance
if i_sub<=10 % ------------------------------ Stanford data exclusion
    load('BadElecs');
    dist_electrodes=nan(numelec,numelec);
    addpath('...');
    fil=xlsread('MNI_Coordinates', subject);
    if size(fil,2)>3
        fil(:,1:2)=[];
    end
    fil(BadElecs,:)=[]; fil(numelec+1:end,:)=[];
else
    cd([path_root '\' subject])
    load(sprintf('electrodes_%s_%s.mat', Task, subject)); fil=electrodes;
end % ------------------------------ Stanford data exclusion
cd(path_main)

for i=1:numelec
    for j=1:numelec
        dist_electrodes(i,j)=0.1*sqrt( (fil(i,1)-fil(j,1))^2 + (fil(i,2)-fil(j,2))^2 + (fil(i,3)-fil(j,3))^2 );
    end
    dist_electrodes(i,i)=nan;
end
dist_electrodes=round(dist_electrodes,2);

%% ----------------- calculate continuous PLV measure
h=waitbar(0,'Calculating phase coupling (PLV)...'); counter=0;
conn_PLV=[]; conn_PLV_Zscored=[];
for i=1:numelec
    for j=1:numelec
        if i<j
            counter=counter+1;
            waitbar(counter/(numelec*(numelec-1)/2));
            conn_PLV(i,j,:)=PLV_Sepideh(edata.trial{1}(i,:),edata.trial{1}(j,:), Fs, freq, 'ImC', 0, 1);
            conn_PLV_Zscored(i,j,:)=( conn_PLV(i,j,:)-nanmean(conn_PLV(i,j,:)) )/nanstd(conn_PLV(i,j,:));
            conn_PLV(j,i,:)=conn_PLV(i,j,:); conn_PLV_Zscored(j,i,:)=conn_PLV_Zscored(i,j,:);
        end
    end
    conn_PLV(i,i,:)=nan; conn_PLV_Zscored(i,i,:)=nan;
end
close(h)
% static FC
conn_PLV_static=nanmean(conn_PLV,3);
[conn_PLV_static_RegOut,~,~,~]=Dist_Reg_Out(conn_PLV_static,dist_electrodes);

%% ----------------- calculate Amp coupling measure
h=waitbar(0,'Calculating Amp coupling...'); counter=0;
conn_Amp=[]; conn_Amp_Zscored=[];
for i=1:numelec
    for j=1:numelec
        if i<j
            counter=counter+1;
            waitbar(counter/(numelec*(numelec-1)/2));
            conn_Amp(i,j,:)=abs(AmpCoupling(edata.trial{1}(i,:),edata.trial{1}(j,:), Fs, freq, 0, 1));
            conn_Amp_Zscored(i,j,:)=( conn_Amp(i,j,:)-nanmean(conn_Amp(i,j,:)) )/nanstd(conn_Amp(i,j,:));
            conn_Amp(j,i,:)=conn_Amp(i,j,:); conn_Amp_Zscored(j,i,:)=conn_Amp_Zscored(i,j,:);
        end
    end
    conn_Amp(i,i,:)=nan; conn_Amp_Zscored(i,i,:)=nan;
end
close(h)
% static FC
conn_Amp_static=nanmean(conn_Amp,3);
[conn_Amp_static_RegOut,~,~,~]=Dist_Reg_Out(conn_Amp_static,dist_electrodes);

%% ----------------- save data
clear oldfolder temp temp1 temp2
if i_sub<=10
    cd(path_results_rest);
    txt=sprintf('FC_ImC_Rest_%s_Freq%d.mat',subject,freq);
else
    if strcmp(Task,'fixation_pwrlaw') || strcmp(Task,'fixation_PAC')
        cd(path_results_rest);
        txt=sprintf('FC_ImC_%s_%s_freq%d.mat',Task,subject,freq);
    else
        cd('...');
        txt=sprintf('FC_ImC_%s_%s_freq%d.mat',Task,subject,freq);
    end
end
save(txt,'-v7.3')
cd(path_main);
end

%% ----------------------------------------------- Main function for task data
function []=Main_tasks(path_main, path_root, path_results_task, freq, i_sub, Task)
%% -------------------------------------------------
% DATA PREP SECTION
% -------------------------------------------------
load('Subjects');
subject=Subjects{i_sub};
addpath(strcat('...',subject));
%% Data loading
switch Task
    case 'Motor_Stanford'
        Tlim=[-2.5 2.5];
    case 'speech_basic'
        Tlim=[-1.5 1.5];
end
if i_sub<=10 % ------------------------------ Stanford data exclusion
    [edata,~,~,~]=load_clean_data(path_root,subject,Task,1);
else
    switch Task
        case 'Motor_Stanford'
            stim_code=11;
            load(sprintf('edata_Stanford_%s_Stim%d_%s.mat', Task, stim_code, subject));
        case 'speech_basic'
            load(sprintf('edata_Stanford_%s_%s.mat', Task, subject));
    end
end % ------------------------------ Stanford data exclusion

% Extract Data specifics
numelec=numel(edata.label);
numtrial=numel(edata.trial);
Electrodes=edata.label;
Fs=edata.fsample;

%% calculate electrode distance
if i_sub<=10 % ------------------------------ Stanford data exclusion
    load('BadElecs');
    dist_electrodes=nan(numelec,numelec);
    addpath('...');
    fil=xlsread('MNI_Coordinates', subject);
    if size(fil,2)>3
        fil(:,1:2)=[];
    end
    fil(BadElecs,:)=[]; fil(numelec+1:end,:)=[];
else
    cd([path_root '\' subject])
    load(sprintf('electrodes_%s_%s.mat',Task, subject)); fil=electrodes; clear electrodes_text
end % ------------------------------ Stanford data exclusion
cd(path_main)

for i=1:numelec
    for j=1:numelec
        dist_electrodes(i,j)=0.1*sqrt( (fil(i,1)-fil(j,1))^2 + (fil(i,2)-fil(j,2))^2 + (fil(i,3)-fil(j,3))^2 );
    end
    dist_electrodes(i,i)=nan;
end
dist_electrodes=round(dist_electrodes,2);

%% -------------------------------------------------
% DYNAMIC FC SECTION
% -------------------------------------------------
%% ----------------------------------------------- Time-freq analyses configurations
Time_step=0.02;
switch freq
    case 1
        Freqrange=[5 7]; Freq_step=1; Win_length=4;
    case 2
        Freqrange=[8 13]; Freq_step=2; Win_length=6;
    case 3
        Freqrange=[14 30]; Freq_step=3; Win_length=10;
    case 4
        Freqrange=[31 60]; Freq_step=4; Win_length=20;
    case 5
        Freqrange=[61 110]; Freq_step=5; Win_length=20;
end
% set cfg1
cfg1              = [];
cfg1.output       = 'fourier';
cfg1.method       = 'mtmconvol';
cfg1.pad         = 'nextpow2';
cfg1.padtype = 'zero';
% Time-Frequency range
cfg1.toi=edata.time{1}(1):Time_step:edata.time{1}(end);
cfg1.foi=Freqrange(1):Freq_step:Freqrange(2);
% Time window
cfg1.t_ftimwin(:,1)=Win_length./cfg1.foi;
if freq>=3
    cfg1.taper        = 'dpss';
    % frequency smoothing
    cfg1.tapsmofrq = 0.15*cfg1.foi;
else
    cfg1.taper        = 'hanning';
end
t=cfg1.toi;

%% ----------------------------------------------- Amp coupling
% estimate short time FFT
TFR=ft_freqanalysis(cfg1,edata);
% AmpAmp estimation
cfg2=[]; cfg2.removemean='yes'; cfg2.method='amplcorr'; cfg2.channelcmb= {Electrodes, Electrodes};
AmpAmp=ft_connectivityanalysis(cfg2,TFR);
FC_Amp=AmpAmp.amplcorrspctrm;
FC_Amp=squeeze(mean(FC_Amp,3)); clear temp;
FC_Amp=abs(FC_Amp);
FC_Amp_Zscored=FC_Amp;
for i=1:numelec
    for j=1:numelec
        if i<j
            FC_Amp_Zscored(i,j,:)=( FC_Amp(i,j,:)-nanmean(FC_Amp(i,j,cfg1.toi<0 & cfg1.toi>-0.5)) )/nanstd(FC_Amp(i,j,cfg1.toi<0 & cfg1.toi>-0.5));
            FC_Amp_Zscored(j,i,:)=FC_Amp_Zscored(i,j,:);
        end
    end
    FC_Amp(i,i,:)=nan; FC_Amp_Zscored(i,i,:)=nan;
end
% Baseline time-averaged static FC
FC_Amp_static_timeaveraged_baseline=nanmean(FC_Amp(:,:, t<0),3);
[FC_Amp_static_timeaveraged_baseline_regout,~,~,~]=Dist_Reg_Out(FC_Amp_static_timeaveraged_baseline,dist_electrodes);
% post-stimulus time-averaged static FC
FC_Amp_static_timeaveraged=nanmean(FC_Amp(:,:, t>0),3);
[FC_Amp_static_timeaveraged_regout,~,~,~]=Dist_Reg_Out(FC_Amp_static_timeaveraged,dist_electrodes);

%% ----------------------------------------------- PhC
% estimate PhC
cfg2=[]; cfg2.method= 'coh'; cfg2.complex = 'imag'; cfg2.channelcmb= {Electrodes, Electrodes};
PLV_TF=ft_connectivityanalysis(cfg2,TFR);
FC_PLV=squeeze(mean(PLV_TF.cohspctrm,2));
temp=zeros(numelec, numelec, size(FC_PLV,2));
for k=1:size(FC_PLV,1)
    i=find(strcmp(Electrodes, PLV_TF.labelcmb(k,1))>0);
    j=find(strcmp(Electrodes, PLV_TF.labelcmb(k,2))>0);
    temp(i,j,:)=FC_PLV(k,:);
end
FC_PLV=temp; clear temp; FC_PLV_Zscored=FC_PLV;
for i=1:numelec
    for j=1:numelec
        if i<j
            FC_PLV_Zscored(i,j,:)=( FC_PLV(i,j,:)-nanmean(FC_PLV(i,j,cfg1.toi<0 & cfg1.toi>-0.5)) )/nanstd(FC_PLV(i,j,cfg1.toi<0 & cfg1.toi>-0.5));
            FC_PLV_Zscored(j,i,:)=FC_PLV_Zscored(i,j,:);
        end
    end
    FC_PLV(i,i,:)=nan; FC_PLV_Zscored(i,i,:)=nan;
end
% Baseline time-averaged static FC
FC_PLV_static_timeaveraged_baseline=nanmean(FC_PLV(:,:, t<0),3);
[FC_PLV_static_timeaveragedbaseline_regout,~,~,~]=Dist_Reg_Out(FC_PLV_static_timeaveraged_baseline,dist_electrodes);
% post-stimulus time-averaged static FC
FC_PLV_static_timeaveraged=nanmean(FC_PLV(:,:, t>0),3);
[FC_PLV_static_timeaveraged_regout,~,~,~]=Dist_Reg_Out(FC_PLV_static_timeaveraged,dist_electrodes);

%% -------------------------------------------------
% STATIC FC SECTION
% -------------------------------------------------
%% ----------------------------------------------- 1s static connectivity Baseline
% set cfg0
cfg0              = [];
cfg0.output       = 'fourier';
cfg0.method       = 'mtmfft';
cfg0.pad         = 'nextpow2';
cfg0.padtype = 'zero';
cfg0.foi=Freqrange(1):Freq_step:Freqrange(2);
if freq>=3
    cfg0.taper        = 'dpss';
    % frequency smoothing
    cfg0.tapsmofrq = 5;
else
    cfg0.taper        = 'hanning';
end

% ------------ Baseline Freq analysis
cfg=[]; cfg.toilim=[Tlim(1) 0];

edata_trim=ft_redefinetrial(cfg, edata);
% estimate short time FFT
TFR=ft_freqanalysis(cfg0,edata_trim);

%% AmpC
cfg2=[]; cfg2.removemean='yes'; cfg2.method='amplcorr'; cfg2.channelcmb= {Electrodes, Electrodes};
AmpAmp=ft_connectivityanalysis(cfg2,TFR);
FC_Amp_static_baseline=AmpAmp.amplcorrspctrm;
FC_Amp_static_baseline=squeeze(mean(FC_Amp_static_baseline,3)); clear temp;
FC_Amp_static_baseline=abs(FC_Amp_static_baseline);
for elec=1:numelec
    FC_Amp_static_baseline(elec,elec)=nan;
end
[FC_Amp_static_baseline_regout,~,~,~]=Dist_Reg_Out(FC_Amp_static_baseline,dist_electrodes);

%% PLV
cfg2=[]; cfg2.method= 'coh'; cfg2.complex = 'imag';  cfg2.channelcmb= {Electrodes, Electrodes};
PLV_TF=ft_connectivityanalysis(cfg2,TFR);
FC_PLV_static_baseline=squeeze(mean(PLV_TF.cohspctrm,2));
temp=zeros(numelec, numelec);
for k=1:size(FC_PLV_static_baseline,1)
    i=find(strcmp(Electrodes, PLV_TF.labelcmb(k,1))>0);
    j=find(strcmp(Electrodes, PLV_TF.labelcmb(k,2))>0);
    temp(i,j)=FC_PLV_static_baseline(k);
end
FC_PLV_static_baseline=temp; clear temp;
for elec=1:numelec
    FC_PLV_static_baseline(elec,elec)=nan;
end
[FC_PLV_static_baseline_regout,~,~,~]=Dist_Reg_Out(FC_PLV_static_baseline,dist_electrodes);

%% ----------------------------------------------- 1s static FC post-stimulus
% post-stimlus
cfg=[]; cfg.toilim=[0 Tlim(2)];

edata_trim=ft_redefinetrial(cfg, edata);
% estimate short time FFT
TFR=ft_freqanalysis(cfg0,edata_trim);

%% AmpC
cfg2=[]; cfg2.removemean='yes'; cfg2.method='amplcorr'; cfg2.channelcmb= {Electrodes, Electrodes};
AmpAmp=ft_connectivityanalysis(cfg2,TFR);
FC_Amp_static=AmpAmp.amplcorrspctrm;
FC_Amp_static=squeeze(mean(FC_Amp_static,3)); clear temp;
FC_Amp_static=abs(FC_Amp_static);
for elec=1:numelec
    FC_Amp_static(elec,elec)=nan;
end
[FC_Amp_static_regout,~,~,~]=Dist_Reg_Out(FC_Amp_static,dist_electrodes);

%% PhC
cfg2=[]; cfg2.method= 'coh'; cfg2.complex = 'imag';  cfg2.channelcmb= {Electrodes, Electrodes};
PLV_TF=ft_connectivityanalysis(cfg2,TFR);
FC_PLV_static=squeeze(mean(PLV_TF.cohspctrm,2));
temp=zeros(numelec, numelec);
for k=1:size(FC_PLV_static,1)
    i=find(strcmp(Electrodes, PLV_TF.labelcmb(k,1))>0);
    j=find(strcmp(Electrodes, PLV_TF.labelcmb(k,2))>0);
    temp(i,j)=FC_PLV_static(k);
end
FC_PLV_static=temp; clear temp;
for elec=1:numelec
    FC_PLV_static(elec,elec)=nan;
end
[FC_PLV_static_regout,~,~,~]=Dist_Reg_Out(FC_PLV_static,dist_electrodes);

%% surr FC matrices
h=waitbar(0, 'Please wait...Generating surrogate correlation values...'); counter=0;
for repeat=1:1000
    counter=counter+1;
    waitbar( counter/ 1000 );
    FC_PLV_static_surr(:,:,repeat)=Phase_permute_2D(FC_PLV_static);
    FC_PLV_static_regout_surr(:,:,repeat)=Phase_permute_2D(FC_PLV_static_regout);
    FC_Amp_static_surr(:,:,repeat)=Phase_permute_2D(FC_Amp_static);
    FC_Amp_static_regout_surr(:,:,repeat)=Phase_permute_2D(FC_Amp_static_regout);
end
close(h); clear h counter
%% ----------------------------------------------- save data
clear temp
cd(path_results_task);
txt=sprintf('FC_ImC_%s_%s_freq%d.mat',Task, subject, freq);
save(txt)
cd(path_main);
end
% ----------------------------------------------- End of File
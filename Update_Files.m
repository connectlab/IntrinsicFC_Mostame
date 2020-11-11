function [] = Update_Files(subject, freq, Tasks)
switch freq
    case 1
        Freqrange=[5 7];
    case 2
        Freqrange=[8 13];
    case 3
        Freqrange=[14 30];
    case 4
        Freqrange=[31 60];
    case 5
        Freqrange=[61 110];
end
switch Tasks
    case 'Rest'
        cd('C:\Users\CONNE\Documents\MATLAB\UIUC\results\Rest\With_strips');
        txt=sprintf('FC_Rest_%s_Freq%d.mat',subject,freq);
    case 'CRM'
        cd('Y:\mostame2\UIUC\results\Task');
        txt=sprintf('AmpPhCorr_%s_Cond1_Freq%dto%d_smooth_abs.mat',subject,Freqrange(1),Freqrange(2));
    case 'Motor_Stanford'
        cd('Y:\mostame2\UIUC\results\Task');
        txt=sprintf('AmpPhCorr_%s_Cond1_Freq%dto%d_smooth_abs.mat',subject,Freqrange(1),Freqrange(2));
end
load(txt)

%% --------------------------------------------- do stuff
cd(['Y:\mostame2\ECOGimport' '\' subject])
stim_code=0;
load(sprintf('edata_Stanford_MotorBasic_Stim%d_%s.mat',stim_code,subject));
%% estimate AmpC for continuous data
conn_Amp_cont=[]; conn_Amp_Zscored_cont=[];
[conn_Amp_cont, conn_Amp_Zscored_cont, conn_Amp_static, conn_Amp_static_RegOut]=FC_estimate_AmpC(data_pre,numelec,Fs,freq,dist_electrodes);
%% estimate PhC for continuous data
conn_PLV_cont=[]; conn_PLV_Zscored_cont=[];
[conn_PLV_cont, conn_PLV_Zscored_cont, conn_PLV_static, conn_PLV_static_RegOut]=FC_estimate_PLV(data_pre,numelec,Fs,freq,dist_electrodes);
%% find out length of time window
% for PhC
A=[]; B=[];
for i=1:numelec
    for j=1:numelec
        if i<j
            temp=[];
            temp=squeeze(conn_PLV_Zscored_cont(i,j,:))'; temp=temp(~isnan(temp));
            [a w]=cpsd(temp,temp,[],[],1024); w=w/(2*pi);
            % find right side FWHM
            temp=[]; temp=abs((a-0.5*max(a))); temp(w<w(a==max(a)))=nan; temp=w( temp-nanmin(temp) == min( temp-nanmin(temp) ) );
            %             plot(w,a,'linewidth',2); line([temp temp],[0 a(w==temp)+2], 'color','k')
            A=[A w(a==max(a))];
            B=[B temp];
        end
    end
end
W1=median(A); W1_FWHM=median(B);
% for AmpC
A=[]; B=[];
for i=1:numelec
    for j=1:numelec
        if i<j
            temp=[];
            temp=squeeze(conn_Amp_Zscored_cont(i,j,:))'; temp=temp(~isnan(temp));
            [a w]=cpsd(temp,temp,[],[],1024); w=w/(2*pi);
            % find right side FWHM
            temp=[]; temp=abs((a-0.5*max(a))); temp(w<w(a==max(a)))=nan; temp=w( temp-nanmin(temp) == min( temp-nanmin(temp) ) );
            %             plot(w,a,'linewidth',2); line([temp temp],[0 a(w==temp)+2], 'color','k')
            A=[A w(a==max(a))];
            B=[B temp];
        end
    end
end
W2=median(A); W2_FWHM=median(B);
Corr_step=round(1/min(W1, W2),1,'significant');
if Corr_step>100
    fprintf('<<<<<<<<<<<<<<<WARNING! Corr_Step is longer than 100s = %d>>>>>>>>>>>>>>\n',Corr_step)
    Corr_step=100;
elseif Corr_step<20
    fprintf('<<<<<<<<<<<<<<<WARNING! Corr_Step is shorter than 20s = %d>>>>>>>>>>>>>>\n',Corr_step)
    Corr_step=20;
end
clear A B a w W1 W2
%% assess correlation between the two dynamics
h=waitbar(0,'Calculating dynamic correlations...'); counter=0;
Corr=nan(size(conn_PLV_static));
Corr_alltime=[];
for k=1:Corr_step/2:(size(conn_Amp_cont,3)-Corr_step+1)
    counter=counter+1;
    waitbar(counter/(numel(1:Corr_step/2:size(conn_Amp_cont,3)-Corr_step+1)));
    for i=1:numelec
        for j=1:numelec
            if i<j && Pair_significant(i,j)==1
                temp1=[]; temp1=conn_PLV_Zscored_cont(i,j,k:k+Corr_step-1); temp1=squeeze(temp1)'; temp1(isnan(temp1))=[];
                temp2=[]; temp2=conn_Amp_Zscored_cont(i,j,k:k+Corr_step-1); temp2=squeeze(temp2)'; temp2(isnan(temp2))=[];
                [acorr lag]=xcorr(temp1-nanmean(temp1),temp2-nanmean(temp2),'coeff');
                Corr(i,j)=acorr(lag==0); Corr(j,i)=Corr(i,j);
            end
        end
        Corr(i,i)=nan;
    end
    Corr_alltime=cat(3,Corr_alltime,Corr);
end
close(h); clear h;
%% permute FC dynamics for statistical test
R=500;
Corr_surr_repeat=nan([size(Corr_alltime), R]);
h=waitbar(0,'generating  correlations...'); counter=0;
for repeat=1:R
    counter=counter+1;
    waitbar(counter/R);
    Corr_surr_alltime=nan(size(Corr_alltime));
    for i=1:numelec
        for j=1:numelec
            if i<j && Pair_significant(i,j)==1
                temp=[]; temp=squeeze(conn_PLV_Zscored_cont(i,j,:))'; temp1=Phase_permute(temp);
                temp=[]; temp=squeeze(conn_Amp_Zscored_cont(i,j,:))'; temp2=Phase_permute(temp);
                count=0;
                for k=1:Corr_step/2:(size(conn_Amp_cont,3)-Corr_step+1)
                    count=count+1;
                    tmp1=temp1(k:k+Corr_step-1); tmp2=temp2(k:k+Corr_step-1);
                    tmp1(isnan(tmp1))=[]; tmp2(isnan(tmp2))=[];
                    [acorr lag]=xcorr(tmp1-nanmean(tmp1),tmp2-nanmean(tmp2),'coeff');
                    Corr_surr_alltime(i,j,count)=acorr(lag==0); Corr_surr_alltime(j,i,count)=Corr_surr_alltime(i,j,count);
                end
            end
        end
    end
    Corr_surr_repeat(:,:,:,repeat)=Corr_surr_alltime;
end
close(h); clear h PLV_perm Amp_perm

%% ----------------------------------------------- FDR control statistical test
% MCP Hochberg FDR for positive values
alpha=0.05;
pval=nan(size(Corr_alltime));
for k=1:size(Corr_alltime,3)
    for i=1:numelec
        for j=1:numelec
            if i<j && Pair_significant(i,j)==1
                clear temp; temp=sort( squeeze(Corr_surr_repeat(i,j,k,:)), 'ascend' )';
                pval(i,j,k)=1- (find(Corr_alltime(i,j,k)<[temp inf],1,'first')-1)/length(temp); pval(j,i,k)=pval(i,j,k);
            end
        end
    end
end
temp_pval=pval;
for i=1:numelec
    for j=1:numelec
        if i<=j
            temp_pval(i,j,:)=nan;
        end
    end
end
% positive tail
clear temp; temp=sort(temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
temp=temp./[1:length(temp)]*length(temp);
K_hochberg=find(temp<=alpha,1,'last');
if ~isempty(K_hochberg)
    clear temp; temp=sort(temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
    Corr_alltime_significant_hochberg=double(pval<=temp(K_hochberg));
else
    Corr_alltime_significant_hochberg=zeros(size(pval));
end
% negative tail
clear temp; temp=sort(1-temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
temp=temp./[1:length(temp)]*length(temp);
K_hochberg=find(temp<=alpha,1,'last');
if ~isempty(K_hochberg)
    clear temp; temp=sort(1-temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
    Corr_alltime_significant_hochberg=Corr_alltime_significant_hochberg -double((1-pval)<=temp(K_hochberg));
end
%exclude irrelevant electrodes
for i=1:size(Corr_alltime_significant_hochberg,3)
    temp=squeeze(Corr_alltime_significant_hochberg(:,:,i)); temp(Pair_significant==0)=nan;
    Corr_alltime_significant_hochberg(:,:,i)=temp;
end

%% ---------------------------------------------------
% save
clear temp
switch Tasks
    case 'Rest'
        temp=cd('C:\Users\CONNE\Documents\MATLAB\UIUC\results\Rest\With_strips');
        txt=sprintf('FC_Rest_%s_Freq%d.mat',subject,freq);
    case 'CRM'
        temp=cd('C:\Users\CONNE\Documents\MATLAB\UIUC\results\Task\With_strips');
        txt=sprintf('FC_%s_%s_freq%d.mat',Tasks,subject,freq);
    case 'Motor_Stanford'
        cd('Y:\mostame2\UIUC\results\Task');
        txt=sprintf('AmpPhCorr_%s_Cond1_Freq%dto%d_smooth_abs.mat',subject,Freqrange(1),Freqrange(2));
end
clear Tasks
save(txt)
cd('Y:\mostame2\UIUC');
end

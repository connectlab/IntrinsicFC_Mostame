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
        cd('...');
        txt=sprintf('FC_Rest_%s_Freq%d.mat',subject,freq);
    case 'Motor_Stanford'
        cd('...');
        txt=sprintf('AmpPhCorr_%s_Cond1_Freq%dto%d_smooth_abs.mat',subject,Freqrange(1),Freqrange(2));
end
load(txt)

%% --------------------------------------------- do stuff

%% ---------------------------------------------------
% save
clear temp
switch Tasks
    case 'Rest'
        temp=cd('...');
        txt=sprintf('FC_Rest_%s_Freq%d.mat',subject,freq);
    case 'Motor_Stanford'
        cd('...');
        txt=sprintf('AmpPhCorr_%s_Cond1_Freq%dto%d_smooth_abs.mat',subject,Freqrange(1),Freqrange(2));
end
clear Tasks
save(txt)
cd('...');
end

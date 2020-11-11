function [edata, data_pre, cfg, electrodes_coordinate] = load_clean_data(path_root,subject,Task,want_strips)
switch Task
    case 'RestEyesOpen'
        path_task='\01_RestEyesOpen';
    case 'CRM'
        path_task='\04_CRM';
    case 'ObjNaming_1'
        path_task='\02_ObjNaming';
    case 'ObjNaming_2'
        path_task='\02_ObjNaming';
    case 'RestEyesClosed'
        path_task='\01_RestEyesClosed';
end
% path_subject=strcat(path_root,subject); addpath(path_subject); load('BadElecs.mat');
% if want_strips==0
%     path_data=strcat(path_subject,path_task,'\Preprocessed');
% else
%     path_data=strcat(path_subject,path_task,'\Preprocessed_With_Strips');
% end
% --------------------------------------------------------------------------- test
path_data=path_root;
data_pre=[]; cfg=[];
% ---------------------------------------------------------------------------
addpath( path_data );
if strcmp(Task,'ObjNaming_1') || strcmp(Task,'ObjNaming_2')
    load(sprintf('edata_clean_%s_%s.mat','ObjNaming',subject))
else
    if ~want_strips
        load(sprintf('edata_clean_%s0.8_%s.mat',Task,subject))
        load(sprintf('edata_marked_%s_%s.mat',Task,subject))
        load(sprintf('edata_marked_%s_%s_cfg.mat',Task,subject))
    else
        load(sprintf('edata_clean_%s_%s.mat',Task,subject))
        load(sprintf('edata_marked_%s_%s.mat',Task,subject))
        load(sprintf('edata_marked_%s_%s_cfg.mat',Task,subject))
    end
end
% configure electrode coordinates
electrodes_coordinate=[];
if want_strips==0
    addpath([path_root '\MNI']); load(sprintf('electrodes_coordinate_Brain_%s.mat',subject));
    electrodes_coordinate(BadElecs(BadElecs<size(electrodes_coordinate,1)),:)=[];
end
end

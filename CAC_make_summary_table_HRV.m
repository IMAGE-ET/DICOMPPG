
%% this script will read all .mat files (Kubios preprocessed results) in a directory, 
% extract the relevant information and write them in a summary table, which
% will be saved in csv format for further processing (SPSS, R, etc.)
%
% CAC version
% michael.gaebler@gmail.com, 12.9.2014
%
% NB: assuming that file naming is CAC_12345_rest1*


datadir = 'F:\Directory1\test'; % ANPASSEN: add directory, where .mat files are - assuming, they're all in the same folder

bigdata = [];

    
    hrfiles = dir(fullfile(datadir,'*.mat'));
    
    for ifile = 1:length(hrfiles)
        
        
        clear Res
        
        load(fullfile(datadir,hrfiles(ifile).name));
        
        
        
        bigdata(ifile,1) = str2num(hrfiles(ifile).name(5:9)); % subject number !! ggf. ANPASSEN
        bigdata(ifile,2) = str2num(hrfiles(ifile).name(15)); % resting-state session number !! ggf. ANPASSEN
         
        bigdata(ifile,3) =  Res.HRV.Statistics.mean_RR;  % mean IBI
        bigdata(ifile,4) =  Res.HRV.Statistics.std_RR;   % SD IBI
        bigdata(ifile,5) =  Res.HRV.Statistics.mean_HRV; % mean HR
        bigdata(ifile,6) =  Res.HRV.Statistics.std_HRV;  % SD HR
        bigdata(ifile,7) =  Res.HRV.Statistics.RMSSD;
        bigdata(ifile,8) =  Res.HRV.Statistics.NN50;
        bigdata(ifile,9) =  Res.HRV.Statistics.pNN50;
        
        bigdata(ifile,10) =  Res.HRV.Frequency.Welch.LF_peak;
        bigdata(ifile,11) =  Res.HRV.Frequency.Welch.HF_peak;
        bigdata(ifile,12) =  Res.HRV.Frequency.Welch.LF_power;
        bigdata(ifile,13) =  Res.HRV.Frequency.Welch.HF_power * 1000000; % !!! HF-HRV power 
        
        bigdata(ifile,14) = Res.HRV.Frequency.Welch.HF_power_nu;  % HF normalized units
        bigdata(ifile,15) = Res.HRV.Frequency.Welch.HF_power_prc; % HF normalized units
        bigdata(ifile,16) = Res.HRV.Frequency.Welch.LF_HF_power;  % LF/HF ratio
        
        
    end
    
    bigdata =  [bigdata; bigdata];
    
    
    


csvwrite(fullfile(datadir,'CAC_HR_summary.csv'), bigdata);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Heart rate from Siemens .puls files
% Author: Michael Gaebler
% contact: michael.gaebler@gmail.com
% requires: peakdet.m from PhLEM toolbox: https://sites.google.com/site/phlemtoolbox/Home/downloads
% Version: 09/2014 (CAC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function CAC_HR_extract(subj)

clear all

bildflag  = input('Willst du Figures sehen? (1=ja; 0=nein): ', 's');

subjtmp  = input('Hier Nummer des Probanden eingeben: ', 's');
subj=['CAC_', subjtmp];

subjpath = fullfile('F:','Kathi\','pulse_data',subj); % ANPASSEN


scanlength_pts = 15750;% 5 min 15 sec = 315 sec = 315000/20 msec =15750 Datenpunkte


alldicoms = dir(fullfile(subjpath,['*rest*.dcm']));
pulsename = dir(fullfile(subjpath,['*rest*.puls']));

addpath('F:\Kathi\Matlab')


%% erstes rest

ima_hdr = spm_dicom_headers([fullfile(subjpath,alldicoms(1).name)]);
dicom_time_mdh = ima_hdr{1}.AcquisitionTime * 1000; % multiply by 1000 to get into same time scale


[cat,puls_time_mdh] = textread(fullfile(subjpath,pulsename(1).name), '%s%n%*[^\n]', 2, 'delimiter', ':', 'headerlines', 9);

for i = 1:2; % find 'LogStartMDHTime'
    if strcmp(cat(i),'LogStartMDHTime')
        % diff_mdh = int32(dicom_time_mdh) - puls_time_mdh(i); % onsets of DICOMS in ms
        diff_mdh = dicom_time_mdh - puls_time_mdh(i); % onset of DICOMS in ms
        
    end
end

diff_mdh_pts = floor(diff_mdh/20); % convert to data points (recorded with 50 Hz)


datdum = fopen(fullfile(subjpath,pulsename(1).name)); % öffne Datei und lies sie aus
d = textscan(datdum,'%4n');
fclose(datdum);

mdh_extract_start = diff_mdh_pts + 5; % 4 acquisition parameters at beginning of .puls file

ppg.values = d{1};

%ppg_val_clean = ppg_values(ppg_values < 5000); % werte ohne 5000er


scan_end = mdh_extract_start + scanlength_pts;% 5 min 15 sec = 315 sec = 315000/20 msec = 15750 Datenpunkte

ppg_data_tmp = ppg.values(mdh_extract_start:scan_end + 250); % approximation: 5 min * 50 bpm
ppg.data =  ppg.values(mdh_extract_start:scan_end+length(ppg_data_tmp(ppg_data_tmp==5000))); % add number of lost data points due to peaks
ppg.data_clean = ppg.data(ppg.data < 5000);

if str2num(bildflag)
    figure, plot(ppg.data), title('pulse data (rest1) with peak detection')
end

ppg.t_clean = 0:1/50:length(ppg.data_clean)/50-1/50;
ppg.t =  0:1/50:length(ppg.data)/50-1/50;

ppg.pks_clean = peakdet(ppg.data_clean,200);


if str2num(bildflag)
    
    figure;
    plot(ppg.t_clean,ppg.data_clean,'b',...
        ppg.t_clean(ppg.pks_clean(:,1)),ppg.pks_clean(:,2),...
        'rv','MarkerFaceColor','r') % plot PPG, clean and downsampled
    title(['Peakdet peak detection (rest1): ' num2str(length(ppg.pks_clean)) ' pks.'])
end

% extract peaks (coded as 5000 in Siemens files, if not or faulty, use
% peakdet.m)

%ppg_data = uint16(ppg_data);


ppg.siemens_pks = find(ppg.data==5000)*20; % peaks in ms
ppg.siemens_pks_pts = find(ppg.data==5000); % peaks in ms

if str2num(bildflag)
    
    figure;
    plot(ppg.t,ppg.data,'b',...
        ppg.t(ppg.siemens_pks_pts),ppg.data(ppg.siemens_pks_pts-1),...
        'rv','MarkerFaceColor','r') % plot PPG, clean and downsampled
    title(['Siemens peak detection (rest2): ' num2str(length(ppg.siemens_pks)) ' pks.'])
end

%[find_pks, find_locs] = findpeaks(ppg_data_clean);
ppg.ibis_peakdet = diff(ppg.pks_clean(:,1)); % IBIs (in ms)


ppg.ibis_siemens = diff(ppg.siemens_pks); % IBIs (in ms)

if str2num(bildflag)
    %     figure, hist(ibis), title(['IBIs in ms (rest1); avg. HR: ' num2str(length(siemens_pks)/5.15) ' bpm'])
    figure, hist(ppg.ibis_siemens), title(['IBIs in ms (rest1); avg. HR: ' num2str(60/(mean(ppg.ibis_siemens)/1000)) ' bpm'])
    figure, hist(ppg.ibis_peakdet), title(['IBIs in ms (rest1); avg. HR: ' num2str(60/(mean(ppg.ibis_peakdet)*.02)) ' bpm'])
    
end

%peaks = peakdet(ppg_data,2000);%2200); % delta just from experience - potentially adjust DELTA parameter to detect more/less peaks

ppg.oxyverlauf = ppg.data_clean;

if str2num(bildflag)
    figure, plot(ppg.oxyverlauf), title('pulse data (rest1) without peak detection')
end


dlmwrite(fullfile(subjpath, [subj '_rest1_oxyverlauf.txt']),ppg.oxyverlauf);
dlmwrite(fullfile(subjpath, [subj '_rest1_peaks_siemens.txt']),ppg.siemens_pks);
dlmwrite(fullfile(subjpath, [subj '_rest1_peaks_peakdet.txt']),ppg.pks_clean);
dlmwrite(fullfile(subjpath, [subj '_rest1_ibis_siemens.txt']),ppg.ibis_siemens);
dlmwrite(fullfile(subjpath, [subj '_rest1_ibis_peakdet.txt']),ppg.ibis_ibis);



%% zweites rest

clear diff* mdh_* ppg* oxyver* dicom_* puls_* ima_* *scan_* d* ibi*



ima_hdr = spm_dicom_headers([fullfile(subjpath,alldicoms(2).name)]);
dicom_time_mdh = ima_hdr{1}.AcquisitionTime * 1000; % multiply by 1000 to get into same time scale


[cat,puls_time_mdh] = textread(fullfile(subjpath,pulsename(1).name), '%s%n%*[^\n]', 2, 'delimiter', ':', 'headerlines', 9);

for i = 1:2; % find 'LogStartMDHTime'
    if strcmp(cat(i),'LogStartMDHTime')
        % diff_mdh = int32(dicom_time_mdh) - puls_time_mdh(i); % onsets of DICOMS in ms
        diff_mdh = dicom_time_mdh - puls_time_mdh(i); % onset of DICOMS in ms
        
    end
end

diff_mdh_pts = floor(diff_mdh/20); % convert to data points (recorded with 50 Hz)


datdum = fopen(fullfile(subjpath,pulsename(2).name)); % öffne Datei und lies sie aus
d = textscan(datdum,'%4n');
fclose(datdum);

mdh_extract_start = diff_mdh_pts + 5; % 4 acquisition parameters at beginning of .puls file

ppg.values = d{1};

%ppg_val_clean = ppg_values(ppg_values < 5000); % werte ohne 5000er


scan_end = mdh_extract_start + scanlength_pts;% 5 min 15 sec = 315 sec = 315000/20 msec = 15750 Datenpunkte

ppg_data_tmp = ppg.values(mdh_extract_start:scan_end + 250); % approximation: 5 min * 50 bpm
ppg.data =  ppg.values(mdh_extract_start:scan_end+length(ppg_data_tmp(ppg_data_tmp==5000))); % add number of lost data points due to peaks
ppg.data_clean = ppg.data(ppg.data < 5000);

if str2num(bildflag)
    figure, plot(ppg.data), title('pulse data (rest2) with peak detection')
end

ppg.t_clean = 0:1/50:length(ppg.data_clean)/50-1/50;
ppg.t =  0:1/50:length(ppg.data)/50-1/50;

ppg.pks_clean = peakdet(ppg.data_clean,200);


if str2num(bildflag)
    
    figure;
    plot(ppg.t_clean,ppg.data_clean,'b',...
        ppg.t_clean(ppg.pks_clean(:,1)),ppg.pks_clean(:,2),...
        'rv','MarkerFaceColor','r') % plot PPG, clean and downsampled
    title(['Peakdet peak detection (rest2): ' num2str(length(ppg.pks_clean)) ' pks.'])
end

% extract peaks (coded as 5000 in Siemens files, if not or faulty, use
% peakdet.m)

%ppg_data = uint16(ppg_data);


ppg.siemens_pks = find(ppg.data==5000)*20; % peaks in ms
ppg.siemens_pks_pts = find(ppg.data==5000); % peaks in ms

if str2num(bildflag)
    
    figure;
    plot(ppg.t,ppg.data,'b',...
        ppg.t(ppg.siemens_pks_pts),ppg.data(ppg.siemens_pks_pts-1),...
        'rv','MarkerFaceColor','r') % plot PPG, clean and downsampled
    title(['Siemens peak detection (rest2): ' num2str(length(ppg.siemens_pks)) ' pks.'])
end

%[find_pks, find_locs] = findpeaks(ppg_data_clean);
ppg.ibis_peakdet = diff(ppg.pks_clean(:,1)); % IBIs (in ms)


ppg.ibis_siemens = diff(ppg.siemens_pks); % IBIs (in ms)

if str2num(bildflag)
    %     figure, hist(ibis), title(['IBIs in ms (rest1); avg. HR: ' num2str(length(siemens_pks)/5.15) ' bpm'])
    figure, hist(ppg.ibis_siemens), title(['IBIs in ms (rest2); avg. HR: ' num2str(60/(mean(ppg.ibis_siemens)/1000)) ' bpm'])
    figure, hist(ppg.ibis_peakdet), title(['IBIs in ms (rest2); avg. HR: ' num2str(60/(mean(ppg.ibis_peakdet)*.02)) ' bpm'])
    
end

%peaks = peakdet(ppg_data,2000);%2200); % delta just from experience - potentially adjust DELTA parameter to detect more/less peaks

ppg.oxyverlauf = ppg.data_clean;

if str2num(bildflag)
    figure, plot(ppg.oxyverlauf), title('pulse data (rest2) without peak detection')
end


dlmwrite(fullfile(subjpath, [subj '_rest2_oxyverlauf.txt']),ppg.oxyverlauf);
dlmwrite(fullfile(subjpath, [subj '_rest2_peaks_siemens.txt']),ppg.siemens_pks);
dlmwrite(fullfile(subjpath, [subj '_rest2_peaks_peakdet.txt']),ppg.pks_clean);
dlmwrite(fullfile(subjpath, [subj '_rest2_ibis_siemens.txt']),ppg.ibis_siemens);
dlmwrite(fullfile(subjpath, [subj '_rest2_ibis_peakdet.txt']),ppg.ibis_ibis);






% Du kannst die Herzschläge (siemens_pks), die in ms angegeben sind, in TRs
% umrechnen, indem du durch 2020 (1 TR) dividierst.

% peaks_in_TR = siemens_pks/2020;

%% extract actual sampling frq. of pulse oxi and actual length of resting state scan
% michael.gaebler@gmail.com
% INPUT:
%       subjtmp - subject number as string
% OUTPUT: datout, including
%       sf - individual sampling frequency
%       name - resting-state session
%
% 11.9.2014

function  [datout] = CAC_HR_extract_length_and_oxyverlauf(subjtmp)

%[ppg.sf, scan_length_s] =
subj=['CAC_', subjtmp];


subjpath = fullfile('F:','Kathi\','pulse_data',subj); % ANPASSEN
%subjpath = fullfile('T:\Dokumente\','test1'); % ANPASSEN


start_dcm = dir(fullfile(subjpath,['*rest*start.dcm']));
end_dcm = dir(fullfile(subjpath,['*rest*end.dcm']));

pulsename = dir(fullfile(subjpath,['*rest*.puls']));

datout = [];


%% loop through resting state sessions

%for irest = 1:2
for irest = 1:length(start_dcm)
    try
        ima_hdr_start = spm_dicom_headers([fullfile(subjpath,start_dcm(irest).name)]);
        dicom_time_mdh_start = ima_hdr_start{1}.AcquisitionTime * 1000; % multiply by 1000 to get into same time scale (ms)
        
        
        ima_hdr_2 = spm_dicom_headers([fullfile(subjpath,end_dcm(irest).name)]);
        dicom_time_mdh_end = ima_hdr_2{1}.AcquisitionTime * 1000; % multiply by 1000 to get into same time scale (ms)
        
        scan_end = dicom_time_mdh_end + 2020;
        
        % calculate scan length
        scan_length_ms = scan_end - dicom_time_mdh_start; % scan length in ms
        %scan_length_pts = scan_length_ms/20; % scan length in datapoints
        scan_length_s = scan_length_ms/1000;
        scan_length_min = fix(scan_length_s/60);
        scan_length_sec = scan_length_s-60*scan_length_min;
        
        [cat,pulse_time_mdh] = textread(fullfile(subjpath,pulsename(irest).name), '%s%n%*[^\n]', 4, 'delimiter', ':', 'headerlines', 9);
        
        
        for i = 1:4; % find 'LogStartMDHTime'
            if strcmp(cat(i),'LogStartMDHTime')
                % diff_mdh = int32(dicom_time_mdh) - puls_time_mdh(i); % onsets of DICOMS in ms
                diff_mdh_start = dicom_time_mdh_start - pulse_time_mdh(i); % onset of DICOMS in ms
                pulse_start = pulse_time_mdh(i);
                
                % diff_mdh_end = dicom_time_mdh_end - puls_time_mdh(i); % onset LAST DICOM in ms
                
            elseif  strcmp(cat(i),'LogStopMDHTime')
                % diff_mdh = int32(dicom_time_mdh) - puls_time_mdh(i); % onsets of DICOMS in ms
                pulse_end = pulse_time_mdh(i);
                diff_mdh_end = pulse_end - dicom_time_mdh_end; % onset of DICOMS in ms
                
                % diff_mdh_end = dicom_time_mdh_end - puls_time_mdh(i); % onset LAST DICOM in ms
                %     else display('no onset/offset times provided !!')
            end
        end
        
      
        
        
        
        datdum = fopen(fullfile(subjpath,pulsename(irest).name)); % öffne Datei und lies sie aus
        d = textscan(datdum,'%4n');
        fclose(datdum);
        
        %mdh_extract_start = diff_mdh_start_pts + 5; % 4 acquisition parameters at beginning of .puls file
        
        ppg.values = d{1};
        ppg.pulse = ppg.values(5:end); % 4 acquisition parameters at beginning of .puls file
        
        ppg.pulse_clean = ppg.pulse(ppg.pulse < 5000);
        
        pulse_length_s = (pulse_end - pulse_start)/1000;
        
        ppg.sf = length(ppg.pulse_clean)/pulse_length_s;
        
        display([start_dcm(irest).name(1:15) ', pulse sampling frq.: ' num2str(ppg.sf)...
            ' Hz.; actual scan time: ' num2str(scan_length_min) ':' num2str(scan_length_sec) ' min.'])
        
        scan_length_pts = scan_length_s * ppg.sf;
        
          % pulse_start_pts = floor(pulse_start/20.1); % convert to data points (recorded with 50 Hz)
        diff_mdh_pts = floor(diff_mdh_start/((1/ppg.sf)*1000)); % convert to data points (recorded with individual sf)
        

        
        % control for case that pulse stopped before volume was written
        
        ctrl_end = min([uint32(diff_mdh_pts+scan_length_pts),length(ppg.pulse_clean)]);
        
        
        % ppg.pulse_clean_cropped = ppg.pulse_clean(uint32(diff_mdh_pts):uint32(diff_mdh_pts+scan_length_pts));
        ppg.pulse_clean_cropped = ppg.pulse_clean(uint32(diff_mdh_pts):ctrl_end);
        
        len_diff = uint32(diff_mdh_pts+scan_length_pts)-length(ppg.pulse_clean);
        
        if ctrl_end==length(ppg.pulse_clean)
            display(['ACHTUNG !!!: pulse file ' num2str(len_diff) ' pts. shorter than calculated length for ' start_dcm(irest).name(1:15)])
        else end
        
        % dlmwrite(fullfile(subjpath, [subj '_rest' num2str(irest) '_oxyverlauf.txt']),ppg.pulse_clean_cropped);
        dlmwrite(fullfile(subjpath, [start_dcm(irest).name(1:15) '_oxyverlauf_TEST.txt']),...
            ppg.pulse_clean_cropped);
       
%         transfer = ppg.pulse_clean_cropped;
%         save(fullfile(subjpath, [start_dcm(irest).name(1:15) '_oxyverlauf_TEST.txt']), ...
%             'transfer', '-ascii');
% 

% tmpfile = fopen(fullfile(subjpath, [start_dcm(irest).name(1:15) '_oxyverlauf_TEST.txt']),'a');
% 
% for iline = 1:length(ppg.pulse_clean_cropped)
% fprintf(tmpfile, '%u', ppg.pulse_clean_cropped(iline));
% fprintf(tmpfile, '\n' );
% end
% fclose(tmpfile)
       
% dlmwrite(fullfile(subjpath, [start_dcm(irest).name(1:15) '_oxyverlauf.txt']),ppg.pulse_clean_cropped);

        
        
        display(['@ Successfully wrote ' start_dcm(irest).name(1:15) '_oxyverlauf.txt']);
        
    catch     display(['!!! could NOT !!! write ' start_dcm(irest).name(1:15) '_oxyverlauf.txt']);
    end
    
    datout = [datout; str2num(start_dcm(irest).name(5:9)), str2num(start_dcm(irest).name(15)), ppg.sf];
    
    %%
    %     figure, plot(ppg.pulse_clean_cropped, 'LineWidth', 2)
    %     ppg.pulse_ibis = diff(find(ppg.pulse==5000))*20;
    %     ppg.mean_HR = 60/(mean(ppg.pulse_ibis)/1000);
    
    
end



end
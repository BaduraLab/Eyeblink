%% Code to analyze eyeblink traces per group
% Eyeblink_import needs to be run first
clearvars; 
% Specify the folder where the files live (matlab variables created with import). 
Folder = uigetdir('', 'Select Folder');
% Get a list of all files in the folder with the desired file name pattern and open them. 
filePattern = fullfile(Folder, '*.mat'); 
Files = dir(filePattern);
%Preallocation of the arrays that will be filled in with data from each mouse after every loop iteration.
% change this depending on how many sessions/trials per animal
% traces = zeros (1200,2000);
% sum_traces_US_trials = zeros (100,2000); 
% sum_traces_CSUS_trials = zeros (1000,2000); 
% sum_traces_CS_trials = zeros (100,2000);
% 
% sum_amplitudeCR_CSUS_sessions = double.empty; 
% sum_amplitudeCR_CS_sessions = double.empty;
% sum_amplitudeCR_CSUS_sessions_5 = double.empty;
% sum_amplitudeCR_CS_sessions_5 = double.empty;
% 
% sum_amplitudeCR_CSUS_sessions_std = double.empty; 
% sum_amplitudeCR_CS_sessions_std = double.empty; 
% sum_amplitudeCR_CSUS_sessions_std_5 = double.empty;
% sum_amplitudeCR_CS_sessions_std_5 = double.empty; 
% 
% sum_perc_CR_CSUS = double.empty;
% sum_perc_CR_CS = double.empty;
% sum_perc_CR_CSUS_std = double.empty;
% sum_perc_CR_CS_std = double.empty;
% 
% sum_traces_invalid_session = {};
% sum_traces_puff_session = {};
% sum_traces_paired_session = {};
% sum_traces_led_session = {};


for k = 1 : length(Files)
    baseFileName = Files(k).name;
    fullFileName = fullfile(Files(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', baseFileName);
    load (fullFileName); %name of the mouse variable file (MouseXXXallvariables.mat)
    
    %From here on we start calculating stuff
    traces_orig=traces;
    t=t-t(USstart(1)); % time 0 is when puff is triggered
    t_ms = t*100;
   
    %Baseline
    traces=bsxfun(@minus, traces, median(traces(:,t<-.3),2));
    sessions=unique(sesno); 
    
    % This section "removes" invalid traaces 
    traces_led = traces(ttype==3,:);
    traces(ttype==3,:)= nan; %replace CS only trials (led) for nans. We do this so that we do not apply the valid trial selection to CS trials.  

    %select valid and invalid trials only from trials that have puff (US, CSUS)
    invalid = any(median(traces(:,t>-.2 & t<.6),2)<0,2); %750-900 ms, response to puff (standard)
    valid = any(median(traces(:,t>-.2 & t<.6),2)>0,2); % this time windows can be changed if necessary
    traces_valid = traces(valid,:); 
    traces_invalid = traces(invalid,:);

    %invalid traces per session 
        for a=sessions'
            sel_=sesno==a&invalid==1; 
            id_sel = find(sel_); 
            inv_traces_session {a} = traces (id_sel,:);
        end
    nanss = isnan(traces); 
    traces(nanss) = traces_led; %add back the CS only trials
    traces (invalid,:) = nan; % "delete" the invalid trials for further analysis
      
    %Normalize data to puff trials (US, CSUS) within each session
    numSTD=3; %thic can be changed depending on how stringent you want to be
    tracesnorm=nan(size(traces));
        for a=sessions'
            sel=sesno==a&ttype<3;%sel are the US and CSUS trials in in the current session in the loop        
            idxses=find(sesno==a);
            USmean=nanmedian(traces(sel,:)); %USmean is the median of all traces of sel 
            URmean=median(USmean(t>0 & t<.3)); %mean response to puff in time window when puff happens 
            tracesnorm(idxses,:)= traces(idxses,:)/URmean; %normalization of all traces to UR mean
            %from here on is for determining treshold to exclude trials:
            dat=tracesnorm(idxses,t<-.25); %different time window, before US is presented 
            stds=nanstd(dat,[],2);%std of selection above, std of each trial
            mus=median(stds);%median value of all those stds
            iqrs=iqr(stds); %interquartile range (difference 75th - 25th percentiles of stds).   
            thresh=mus+numSTD*iqrs; %thresh = median of std + (3*iqrs-difference between interquartile range)
            tracesnorm(idxses(stds>thresh),:)=nan; 
            disp(['Excluded ' num2str(numel(find(stds>thresh))) ' Trials'])

        end

    %Select normalized traces for US only, CS-US and CS only trials (handy for the heatmaps later)
    traces_US_trials = tracesnorm(ttype==1,:);
    traces_CSUS_trials = tracesnorm(ttype==2,:); 
    traces_CS_trials = tracesnorm(ttype==3,:);
    
        for a=sessions'
            sel_puff=sesno==a&ttype==1; 
            id_sel_puff = find(sel_puff); 
            traces_puff_session {a} = tracesnorm (id_sel_puff,:);
         
            sel_paired=sesno==a&ttype==2; 
            id_sel_paired = find(sel_paired); 
            traces_paired_session {a} = tracesnorm (id_sel_paired,:);
                    
            sel_led=sesno==a&ttype==3; 
            id_sel_led = find(sel_led); 
            traces_led_session {a} = tracesnorm (id_sel_led,:);
                    
        end
  
    % calculate CR amplitude per session based on CS and CS-US trails
    amplitudesCR=nanmean(tracesnorm(:,t>-0.1 & t<-0.002),2); 
    
    % 5% of UR for CR detection
    meanUR = nanmean(tracesnorm(:,t>0 & t<0.6),2); %750-900
    amplitudesCR_orig = amplitudesCR;
    CRs_5 = amplitudesCR>(0.05*meanUR);  
    amplitudesCR_5 = amplitudesCR (CRs_5,:);
    ttype_5 = ttype(CRs_5,:);
    sesno_5 = sesno(CRs_5,:);
    
    %calculate CR amplitude per session based on CS and CS-US trails
    amplitudeCR_CSUS_sessions=grpstats(amplitudesCR(ttype>1),sesno(ttype>1),{'nanmean'});
    amplitudeCR_CSUS_sessions_std=grpstats(amplitudesCR(ttype>1),sesno(ttype>1),{'nanstd'});
    
    amplitudeCR_CSUS_sessions_5=grpstats(amplitudesCR_5(ttype_5>1),sesno_5(ttype_5>1),{'nanmean'});
    amplitudeCR_CSUS_sessions_std_5=grpstats(amplitudesCR_5(ttype_5>1),sesno_5(ttype_5>1),{'nanstd'});

    %calculate CR amplitude per session based on CS only trails
    amplitudeCR_CS_sessions=grpstats(amplitudesCR(ttype==3),sesno(ttype==3),{'nanmean'});
    amplitudeCR_CS_sessions_std=grpstats(amplitudesCR(ttype==3),sesno(ttype==3),{'nanstd'});
    
    amplitudeCR_CS_sessions_5=grpstats(amplitudesCR_5(ttype_5==3),sesno_5(ttype_5==3),{'nanmean'});
    amplitudeCR_CS_sessions_std_5=grpstats(amplitudesCR_5(ttype_5==3),sesno_5(ttype_5==3),{'nanstd'});
    
    %Calculate percentage of CRs (in 5%> UR mean)
    perc_CR5_CS=grpstats(CRs_5(ttype_5==3),sesno_5(ttype_5==3),{'mean'});
    perc_CR5_CS_std=grpstats(CRs_5(ttype_5==3),sesno_5(ttype_5==3),{'std'});

    perc_CR5_CSUS=grpstats(CRs_5(ttype_5<3),sesno_5(ttype_5<3),{'mean'});
    perc_CR5_CSUS_std=grpstats(CRs_5(ttype_5<3),sesno_5(ttype_5<3),{'std'});

    %save results
%     save([num2str(mouseno),'Results'],'traces','tracesnorm','traces_US_trials','traces_CSUS_trials'...
%         ,'traces_CS_trials', 'ttype','mouseno','sesno', 'CSstart','USstart','CSduration'...
%         ,'USduration','t','t_ms','amplitudesCR','amplitudeCR_CS_sessions','amplitudeCR_CSUS_sessions'...
%        ,'amplitudeCR_CSUS_sessions_std','amplitudeCR_CS_sessions_std','amplitudeCR_CS_sessions_5'...
%        ,'amplitudeCR_CSUS_sessions_5', 'amplitudeCR_CSUS_sessions_std_5','amplitudeCR_CS_sessions_std_5'...
%        ,'inv_traces_session','traces_valid','traces_invalid');

  %Fill the preallocated array with data from each loop (each mouse) 
  sum_traces_CS_trials = sum_traces_CS_trials + traces_CS_trials;
  sum_traces_CSUS_trials = sum_traces_CSUS_trials + traces_CSUS_trials;
  sum_traces_US_trials = sum_traces_US_trials + traces_US_trials;
  
  sum_amplitudeCR_CSUS_sessions = cat(3,sum_amplitudeCR_CSUS_sessions,amplitudeCR_CSUS_sessions);
  sum_amplitudeCR_CSUS_sessions_std = cat(3,sum_amplitudeCR_CSUS_sessions_std,amplitudeCR_CSUS_sessions_std);
  sum_amplitudeCR_CS_sessions = cat(3,sum_amplitudeCR_CS_sessions,amplitudeCR_CS_sessions);
  sum_amplitudeCR_CS_sessions_std = cat(3,sum_amplitudeCR_CS_sessions_std,amplitudeCR_CS_sessions_std);

  sum_amplitudeCR_CSUS_sessions_5 = cat(3,sum_amplitudeCR_CSUS_sessions_5,amplitudeCR_CSUS_sessions_5);
  sum_amplitudeCR_CSUS_sessions_std_5 = cat(3,sum_amplitudeCR_CSUS_sessions_std_5,amplitudeCR_CSUS_sessions_std_5);
  sum_amplitudeCR_CS_sessions_5 = cat(3,sum_amplitudeCR_CS_sessions_5,amplitudeCR_CS_sessions_5);
  sum_amplitudeCR_CS_sessions_std_5 = cat(3,sum_amplitudeCR_CS_sessions_std_5,amplitudeCR_CS_sessions_std_5);

  sum_perc_CR_CSUS = cat(3,sum_perc_CR_CSUS ,perc_CR5_CSUS);
  sum_perc_CR_CS = cat(3,sum_perc_CR_CS ,perc_CR5_CS);
  sum_perc_CR_CSUS_std = cat(3,sum_perc_CR_CSUS_std ,perc_CR5_CSUS_std);
  sum_perc_CR_CS_std = cat(3,sum_perc_CR_CS_std ,perc_CR5_CS_std);
    
  sum_traces_invalid_session = vertcat(sum_traces_invalid_session, inv_traces_session); 
  
  sum_traces_puff_session = vertcat(sum_traces_puff_session, traces_puff_session); 
  sum_traces_paired_session = vertcat(sum_traces_paired_session, traces_paired_session); 
  sum_traces_led_session = vertcat(sum_traces_led_session, traces_led_session); 

end

%% types of trials per session

for p = 1:size(sum_traces_puff_session,2) %p is number of sessions
    
    puff_traces_per_session{p}= vertcat(sum_traces_puff_session{:,p});
    paired_traces_per_session{p}= vertcat(sum_traces_paired_session{:,p});
    led_traces_per_session{p}= vertcat(sum_traces_led_session{:,p});
    
end


%% Calculate invalid trials per session 

for p = 1:size(sum_traces_invalid_session,2)
    
    invalid_traces_per_session{p} = vertcat(sum_traces_invalid_session{:,p});
    perc_invalid_session{p} = (length(invalid_traces_per_session{p}(:,1))/length(paired_traces_per_session{p}))*100; %change
    total_invalid{p} = sum(length(invalid_traces_per_session{p}(:,1)));
    
end

%percentage of invalid trials
total_traces_invalid = vertcat(invalid_traces_per_session{1,:});
total_invalid = cell2mat (total_invalid); 
total_invalid_group = sum (total_invalid); 
perc_invalid = (total_invalid_group/(240*p*k))*100; 
perc_invalid_session = cell2mat(perc_invalid_session);
perc_valid_session = 100 - perc_invalid_session; 

%% Calculate group average traces and group averages
average_traces_CS_trials = sum_traces_CS_trials/k;
average_traces_CSUS_trials = sum_traces_CSUS_trials/k;
average_traces_US_trials = sum_traces_US_trials/k;

average_amplitudeCR_CS_sessions = mean (sum_amplitudeCR_CS_sessions,3); 
std_amplitudeCR_CS_sessions = std(sum_amplitudeCR_CS_sessions,[],3);
sem_amplitudeCR_CS_sessions = std_amplitudeCR_CS_sessions/sqrt(length(average_amplitudeCR_CS_sessions)); 

average_amplitudeCR_CS_sessions_5 = mean (sum_amplitudeCR_CS_sessions_5,3); 
std_amplitudeCR_CS_sessions_5 = std(sum_amplitudeCR_CS_sessions_5,[],3);
sem_amplitudeCR_CS_sessions_5 = std_amplitudeCR_CS_sessions-5/sqrt(length(average_amplitudeCR_CS_sessions_5));

average_amplitudeCR_CSUS_sessions = mean (sum_amplitudeCR_CSUS_sessions,3); 
std_amplitudeCR_CSUS_sessions = std (sum_amplitudeCR_CSUS_sessions,[],3); 
sem_amplitudeCR_CSUS_sessions = std_amplitudeCR_CSUS_sessions/sqrt(length(average_amplitudeCR_CSUS_sessions)); 

average_amplitudeCR_CSUS_sessions_5 = mean (sum_amplitudeCR_CSUS_sessions_5,3); 
std_amplitudeCR_CSUS_sessions_5 = std (sum_amplitudeCR_CSUS_sessions_5,[],3); 
sem_amplitudeCR_CSUS_sessions_5 = std_amplitudeCR_CSUS_sessions_5/sqrt(length(average_amplitudeCR_CSUS_sessions_5)); 

average_percCR_CS = mean (sum_perc_CR_CS,3); 
average_percCR_CSUS = mean (sum_perc_CR_CSUS,3); 
average_percCR_CS_std = std (sum_perc_CR_CS,[],3);
average_percCR_CSUS_std = std (sum_perc_CR_CSUS,[],3);

%% Create structures to organize and save results
Traces_averages = struct('t_ms',t_ms,'average_traces_CS_trials',average_traces_CS_trials, 'average_traces_CSUS_trials'...
    ,average_traces_CSUS_trials,'average_traces_US_trials',average_traces_US_trials...
    ,'perc_invalid_session',perc_invalid_session,'perc_valid_session',perc_valid_session,'total_traces_invalid',total_traces_invalid...
    ,'perc_invalid',perc_invalid);

Traces_sessions = struct('paired_traces_per_session',paired_traces_per_session,'puff_traces_per_session',puff_traces_per_session...
    ,'led_traces_per_session',led_traces_per_session,'invalid_traces_per_session',invalid_traces_per_session);

Averages_CR_amplitude = struct('average_amplitudeCR_CS_sessions_5',average_amplitudeCR_CS_sessions_5,'std_amplitudeCR_CS_sessions_5'...
  ,std_amplitudeCR_CS_sessions_5,'std_amplitudeCR_CSUS_sessions_5',std_amplitudeCR_CSUS_sessions_5,'average_amplitudeCR_CSUS_sessions_5'...
  ,average_amplitudeCR_CSUS_sessions_5,'sem_amplitudeCR_CSUS_sessions_5',sem_amplitudeCR_CSUS_sessions_5,'sem_amplitudeCR_CS_sessions_5'...
  ,sem_amplitudeCR_CS_sessions_5);

Averages_CR_percentage = struct ('average_percCR_CS',average_percCR_CS,'average_percCR_CSUS',average_percCR_CSUS...
    ,'average_percCR_CS_std',average_percCR_CS_std,'average_percCR_CSUS_std',average_percCR_CSUS_std);

All_CR_amplitudes = struct('sum_amplitudeCR_CSUS_sessions_5',sum_amplitudeCR_CSUS_sessions_5,'sum_amplitudeCR_CSUS_sessions_std_5'...
    ,sum_amplitudeCR_CSUS_sessions_std_5,'sum_amplitudeCR_CS_sessions_5',sum_amplitudeCR_CS_sessions_5,'sum_amplitudeCR_CS_sessions_std_5'...
    ,sum_amplitudeCR_CS_sessions_std_5); 

All_CR_percentages = struct('sum_perc_CR_CSUS',sum_perc_CR_CSUS,'sum_perc_CR_CS',sum_perc_CR_CS,'sum_perc_CR_CSUS_std'...
    ,sum_perc_CR_CSUS_std,'sum_perc_CR_CS_std',sum_perc_CR_CS_std); 

save([Folder,'Results'],'Traces_averages','Traces_averages','Traces_sessions','Traces_sessions','Averages_CR_amplitude','Averages_CR_amplitude','Averages_CR_percentage','Averages_CR_percentage'...
    ,'All_CR_amplitudes','All_CR_amplitudes','All_CR_percentages','All_CR_percentages'); 


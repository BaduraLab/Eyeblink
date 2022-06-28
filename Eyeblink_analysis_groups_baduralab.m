%% Code to analyze eyeblink traces per group
% Eyeblink_import needs to be run first
clearvars; 
% Specify the folder where the files live (matlab variables created with import). 
Folder = 'Z:\Maria\Eyeblink Behavior Data\eyeblink_fabian\a'; 

% Get a list of all files in the folder with the desired file name pattern and open them. 
filePattern = fullfile(Folder, '*.mat'); 
Files = dir(filePattern);

%Preallocation of the arrays that will be filled in with
%data from each mouse after every loop iteration. 
sum_traces_US_trials = zeros (100,2000); 
sum_traces_CSUS_trials = zeros (1000,2000); 
sum_traces_CS_trials = zeros (100,2000);
sum_amplitudesCR = zeros (1200,1); 
sum_amplitudeCR_CSUS_sessions = zeros (5,1); 
sum_amplitudeCR_CS_sessions = zeros (5,1); 


for k = 1 : length(Files)
    baseFileName = Files(k).name;
    fullFileName = fullfile(Files(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    load (baseFileName); %baseFileName is the name of the mouse variable file (MouseXXXallvariables.mat)
    
    % From here on we start calculating stuff:
    traces_orig=traces;
    t=t-t(USstart(1));
    t_ms = t*1000;
    %Remove trace offset
    traces=bsxfun(@minus, traces, median(traces(:,t<-.3),2)); %substract the median between traces when t<-0.3 and 2 from the total traces. This is why we save traces_orig before.
    %traces=traces(:, t<1.24); % use only if something weird happens at the end of the traces

    %Normalize data to US only and US-CS trials within each session
    numSTD=3;
    sessions=unique(sesno); USpeakamp=[];
    tracesnorm=nan(size(traces));
    figure;
    hold all
    for a=sessions'
        sel=sesno==a&ttype<3;        
        idxses=find(sesno==a);
        USmean=nanmedian(traces(sel,:));
        plot(t,USmean); hold on
        USpeakamp=max(USmean(t>0.01 & t<.15));
        tracesnorm(idxses,:)=traces(idxses,:)/USpeakamp;
        dat=tracesnorm(idxses,t<-.25);    
        stds=nanstd(dat,[],2);
        mus=median(stds);
        iqrs=iqr(stds);    
        thresh=mus+numSTD*iqrs;
        tracesnorm(idxses(stds>thresh),:)=nan;
        disp(['Excluded ' num2str(numel(find(stds>thresh))) ' Trials'])
    end
    close
    %figure
    %plot(tracesnorm') %I plot this here just to have a quick view.

    %Select normalized traces for US only, CS-US and CS only trials (handy for the heatmaps later)
    traces_US_trials = tracesnorm(ttype==1,:);
    traces_CSUS_trials = tracesnorm(ttype==2,:); 
    traces_CS_trials = tracesnorm(ttype==3,:);

    %calculate CR amplitude per session based on CS and CS-US trails
    amplitudesCR=nanmean(tracesnorm(:,t>-0.04 & t<-0.002),2); %standard is t>-0.04 & t<-0.002
    amplitudeCR_CSUS_sessions=grpstats(amplitudesCR(ttype>1),sesno(ttype>1),{'mean'});

    %calculate and plot CR amplitude per session based on CS only trails
    amplitudesCR=nanmean(tracesnorm(:,t>-0.04 & t<-0.002),2); %standard is t>-0.04 & t<-0.002
    amplitudeCR_CS_sessions=grpstats(amplitudesCR(ttype==3),sesno(ttype==3),{'mean'});

    %save results
    save([num2str(mouseno),'summary'],'tracesnorm','traces_US_trials','traces_CSUS_trials', 'traces_CS_trials', 'ttype','mouseno','sesno'...
       ,'CSstart','USstart','CSduration','USduration','t','t_ms','amplitudesCR','amplitudeCR_CS_sessions','amplitudeCR_CSUS_sessions');

  %Fill the preallocated array with data from each loop (each mouse) 
  sum_traces_CS_trials = sum_traces_CS_trials + traces_CS_trials;
  sum_traces_CSUS_trials = sum_traces_CSUS_trials + traces_CSUS_trials;
  sum_traces_US_trials = sum_traces_US_trials + traces_US_trials;
  sum_amplitudesCR = sum_amplitudesCR + amplitudesCR; 
  sum_amplitudeCR_CSUS_sessions = sum_amplitudeCR_CSUS_sessions + amplitudeCR_CSUS_sessions;
  sum_amplitudeCR_CS_sessions = sum_amplitudeCR_CS_sessions + amplitudeCR_CS_sessions; 

end


%% Calculate group average traces and group average CR amp and save reults!

average_traces_CS_trials = sum_traces_CS_trials/k;
average_traces_CSUS_trials = sum_traces_CSUS_trials/k;
average_traces_US_trials = sum_traces_US_trials/k;
average_amplitudeCR_CSUS_sessions = sum_amplitudeCR_CSUS_sessions/k; 
average_amplitudeCR_CS_sessions = sum_amplitudeCR_CS_sessions/k; 
average_amplitudesCR = sum_amplitudesCR/k; 

save(('Results_groupX'),'average_amplitudeCR_CS_sessions', 'average_amplitudeCR_CSUS_sessions', 'average_amplitudesCR'...
    ,'average_traces_CS_trials', 'average_traces_CSUS_trials', 'average_traces_US_trials');

%% Plot group average traces.

figure
hold all
plot(t_ms,nanmedian(average_traces_US_trials));
plot(t_ms,nanmedian(average_traces_CSUS_trials)); 
plot(t_ms,nanmedian(average_traces_CS_trials)); 
legend('US','CS-US','CS'); 

%% Plot CR amplitude per session.

% Find a way to make boxplots prettier - or try some other
% graphs to represent CR amp per session. Although boxplots
% are handy to understand the spread of the data. 

figure
boxplot(average_amplitudesCR(ttype>1),sesno(ttype>1)); %boxplot shows median!
hold all
plot(sesno(ttype>1),average_amplitudesCR(ttype>1),'.')%individual dots are amplitude per trial
plot(sessions,average_amplitudeCR_CSUS_sessions,'.k','MarkerSize', 25)%mean is black
xlabel('Sessions')
ylabel('CR amplitude')
title ('CR amplitude CS-US trials')

figure
boxplot(average_amplitudesCR(ttype==3),sesno(ttype==3));
hold all
plot(sesno(ttype==3),average_amplitudesCR(ttype==3),'.')
plot(sessions,average_amplitudeCR_CS_sessions,'.k','MarkerSize', 25)%mean is black
xlabel('trials')
ylabel('CR amplitude')
title ('CR amplitude CS trials')

%% Plot group heatmaps - CS only trials
%Use the following line of code if you want to take a
%smaller time window for the heatmap: 
%tracesnorm_CSonly_trange = traces_CSonly_trials(:,t>=-0.25 & t<0.25); 

figure
h = heatmap (average_traces_CS_trials);
h.GridVisible = 'off';
h.MissingDataLabel = 'NaN'; 
h.MissingDataColor = [0.65,0.65,0.65];
h.Colormap = hot;
h.ColorLimits = [0,1];
h.XLabel = 'Time (ms)';
h.YLabel = 'CS only Trials'; 
title ('Timing of responses (CS trials)')

%This needs to change if the time window is changed. 
%Heatmap function is shitty bc you cannot customize the X and Yaxis. 
%This is a "dirty hack":
XLabels = 1:2000;
CustomXLabels = string(XLabels);
CustomXLabels(mod(XLabels,250) ~= 0) = " ";
h.XDisplayLabels = CustomXLabels;
 
YLabels = 1:100;
CustomYLabels = string(YLabels);
CustomYLabels(mod(YLabels,10) ~= 0) = " ";
h.YDisplayLabels = CustomYLabels;

%% Plot group heatmaps - CS-US only trials
%Use the following line of code if you want to take a
%smaller time window for the heatmap: 
%tracesnorm_CSonly_trange = traces_CSonly_trials(:,t>=-0.25 & t<0.25); 

figure
h = heatmap (average_traces_CSUS_trials);
h.GridVisible = 'off';
h.MissingDataLabel = 'NaN'; 
h.MissingDataColor = [0.65,0.65,0.65];
h.Colormap = hot;
h.ColorLimits = [0,1];
h.XLabel = 'Time (ms)';
h.YLabel = 'CS-US Trials'; 
title ('Timing of responses (CS-US trials)')

%This needs to change if the time window is changed. 
%Heatmap function is shitty bc you cannot customize the X and Yaxis. 
%This is a "dirty hack":
XLabels = 1:2000;
CustomXLabels = string(XLabels);
CustomXLabels(mod(XLabels,250) ~= 0) = " ";
h.XDisplayLabels = CustomXLabels;
 
YLabels = 1:1000;
CustomYLabels = string(YLabels);
CustomYLabels(mod(YLabels,100) ~= 0) = " ";
h.YDisplayLabels = CustomYLabels;



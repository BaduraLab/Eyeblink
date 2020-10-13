clearvars;

%% Code to analyze eyeblink traces
load('Mouse400allvariables.mat')
traces_orig=traces;
t=t-t(USstart(1)); 
figure
plot (traces_orig')

%% Remove trace offset
traces=bsxfun(@minus, traces, median(traces(:,t<-.3),2)); %substract the median between traces when t<-0.3 and 2 from the total traces. This is why we save traces_orig before.
%traces=traces(:, t<1.24); % use only if something weird happens at the end of the traces
figure
plot (traces')

%t=t(:, t<1.24); % use only if something weird happens at the end of the traces

%% Normalize data to US only and US-CS trials within each session
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
figure
plot(tracesnorm')

%% Select normalized traces for US only, CS-US and CS only trials (handy for the heatmaps later)
traces_US_trials = tracesnorm(ttype==1,:);
traces_CSUS_trials = tracesnorm(ttype==2,:); 
traces_CS_trials = tracesnorm(ttype==3,:);
t_ms = t*100;
figure
hold all
plot(t_ms,traces_CSUS_trials,'color',[0.85,0.33,0.10]); 
plot(t_ms,traces_US_trials,'color',[0.00,0.45,0.74]);
plot(t_ms,traces_CS_trials,'color',[0.93,0.69,0.13]);


%% Plot averages of CS, US and CS-US for the entire data set
figure
hold all

plot(t_ms,nanmedian(traces_US_trials));
plot(t_ms,nanmedian(traces_CSUS_trials));
plot(t_ms,nanmedian(traces_CS_trials));
legend('US','CS-US','CS'); 

%% save normalized traces
save([num2str(mouseno),'norm'],'tracesnorm','traces_US_trials','traces_CSUS_trials', 'traces_CS_trials', 'ttype','mouseno','sesno'...
   ,'CSstart','USstart','CSduration','USduration','t');

%% calculate and plot CR amplitude per session based on CS and CS-US trails
amplitudesCRs=nanmean(tracesnorm(:,t>-0.04 & t<-0.002),2); %standard is t>-0.04 & t<-0.002
amplitudeCRs_CSUS_sessions=grpstats(amplitudesCRs(ttype>1),sesno(ttype>1),{'mean'});
std_amplitudeCR_CSUS_sessions=grpstats(amplitudesCRs(ttype>1),sesno(ttype>1),{'std'});


figure
%boxplot(amplitudesCRs(ttype>1),sesno(ttype>1)); %when trial type is 2 or 3, CS and CS-US
plot(sessions,amplitudeCRs_CSUS_sessions); 
errorbar(amplitudeCRs_CSUS_sessions,std_amplitudeCR_CSUS_sessions); 

hold on
%plot(sesno(ttype>1),amplitudesCRs(ttype>1),'.')
xlim([.5 10.5]);
xticks([1:10]);
a = get(gca,'XTickLabel');  
set(gca,'TickDir','out');
xlabel('Sessions');
ylabel('CR amplitude')
box('off'); 


%% calculate and plot CR amplitude per session based on CS only trails
amplitudesCRs=nanmean(tracesnorm(:,t>-0.04 & t<-0.002),2); %standard is t>-0.04 & t<-0.002, when mice usually respond?

figure
boxplot(amplitudesCRs(ttype==3),sesno(ttype==3));
hold on
plot(sesno(ttype==3),amplitudesCRs(ttype==3),'.')
xlabel('trials')
ylabel('CR amplitude')
print ([num2str(mouseno),'CRamplitudeCSonlyTrails'],'-depsc')

amplitudeCRs_CSonly_sessions=grpstats(amplitudesCRs(ttype==3),sesno(ttype==3),{'mean'});

 %% save results
save([num2str(mouseno),'results'],'amplitudesCRs','amplitudeCRs_CSonly_sessions','amplitudeCRs_CSUS_sessions');


%% Heatmap plot for CS only trials
%Use the following line of code if you want to take a
%smaller time window for the heatmap: 
%tracesnorm_CSonly_trange = traces_CSonly_trials(:,t>=-0.25 & t<0.25); 

figure
h = heatmap (traces_CS_trials);
h.GridVisible = 'off';
h.MissingDataLabel = 'NaN'; 
h.MissingDataColor = [0.65,0.65,0.65];
h.Colormap = hot;
h.ColorLimits = [0,1];
h.XLabel = 'Time (ms)';
h.YLabel = 'CS only Trials'; 

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

%% Heatmap for CS-US trails
figure
h = heatmap (traces_CSUS_trials);
h.GridVisible = 'off';
h.MissingDataLabel = 'NaN'; 
h.MissingDataColor = [0.65,0.65,0.65];
h.Colormap = hot;
h.ColorLimits = [0,1];
h.XLabel = 'Time (ms)';
h.YLabel = 'CS-US Trials'; 

%This needs to change if the time window is changed. 
%Heatmap function is shitty bc you cannot customize the X and Yaxis. 
%This is a "dirty hack":XLabels = 1:2000;
CustomXLabels = string(XLabels);
CustomXLabels(mod(XLabels,250) ~= 0) = " ";
h.XDisplayLabels = CustomXLabels;
 
YLabels = 1:1000;
CustomYLabels = string(YLabels);
CustomYLabels(mod(YLabels,100) ~= 0) = " ";
h.YDisplayLabels = CustomYLabels;






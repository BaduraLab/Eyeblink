% clearvars;
% %% code to analyze eyeblink traces
% load('Mouse341allvariables.mat')
traces_orig=traces;
t=t-t(USstart(1));

%% remove pre tone trace offset
traces=bsxfun(@minus, traces,median(traces(:,t<-.3),2));
%traces=traces(:, t<1.24); % use only if something weird happens at the end of the traces
figure
plot (traces_orig')
figure
plot (traces')

%t=t(:, t<1.24); % use only if something weird happens at the end of the traces

%% loop though all sessions and compute US only average
%  us_alone_trace=nanmean(traces(ttype==1,:));
%  USpeakamp=max(us_alone_trace);
%  tracesnorm=traces/USpeakamp;
%  figure
%  plot(tracesnorm','color',[.8 .8 .8])
%  hold all
%  plot(nanmedian(tracesnorm(ttype==1,:)))
%  plot(nanmedian(tracesnorm(ttype==2,:)))
%  plot(nanmedian(tracesnorm(ttype==3,:)))
 
%% normalize data to US only and US-CS trials within each session
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
%% plot averages of CS, US and CS-US for the entire data set
close
hold all
figure
plot(t,nanmedian(tracesnorm(ttype==1,:)))
plot(t,nanmedian(tracesnorm(ttype==2,:)))
plot(t,nanmedian(tracesnorm(ttype==3,:)))
%% save normalized traces
% save([num2str(mouseno),'norm'],'tracesnorm','blockno','ttype','mouseno','sesno'...
%    ,'CSstart','USstart','CSduration','USduration','t');
%% calculate and plot CR amplitude per session based on CS and CS-US trails
amplitudesCRs=nanmean(tracesnorm(:,t>-0.04 & t<-0.002),2); %standard is t>-0.04 & t<-0.002
figure
boxplot(amplitudesCRs(ttype>1),sesno(ttype>1));
hold on
plot(sesno(ttype>1),amplitudesCRs(ttype>1),'.')

xlabel('trials')
ylabel('CR amplitude')
print ([num2str(mouseno),'CRamplitude_CS_CS-UStrails'],'-depsc')
%% calculate and plot CR amplitude per session based on CS only trails
amplitudesCRs=nanmean(tracesnorm(:,t>-0.04 & t<-0.002),2); %standard is t>-0.04 & t<-0.002
figure
boxplot(amplitudesCRs(ttype==3),sesno(ttype==3));
hold on
plot(sesno(ttype==3),amplitudesCRs(ttype==3),'.')

xlabel('trials')
ylabel('CR amplitude')
print ([num2str(mouseno),'CRamplitudeCSonlyTrails'],'-depsc')
%% EDITED TILL HERE



%% calculate and plot CR% per session Based on %of UR here 10% 
peakUR=nanmedian(tracesnorm(:,t>0.03 & t<0.07),2); 
%peakUR=nanmedian(tracesnorm(:,t>-0.0001 & t<0.25),2); %standard
CR10=amplitudesCRs>(0.1*peakUR);
perc_cr10=grpstats(CR10(ttype>1),sesno(ttype>1),{'mean'});

figure
plot(unique(sesno),perc_cr10)
xlabel('sessions')
ylabel('CR% 10%treshold')
ylim([0 100]);
print([num2str(mouseno),'CR%10%'],'-depsc')
%% calculate and plot CR% per session Based on %of UR here 15%
CR15=amplitudesCRs>(0.15*peakUR);
perc_cr15=grpstats(CR15(ttype>1),sesno(ttype>1),{'mean'});

figure
plot(unique(sesno),perc_cr15*100)
xlabel('sessions')
ylabel('CR%_15%treshold')
ylim([0 100]);
print ([num2str(mouseno),'CR%15%'],'-depsc')
%% calculate and plot CR% per session based on STD above basline
open_eye=nanmedian(tracesnorm(:,t<-.3),2);
open_eye_std=nanstd(tracesnorm(:,t<-.3),[],2);

sd=3; %change depending on what STD you want
thresh=nanmedian(open_eye + sd*open_eye_std); 
resp=(amplitudesCRs>thresh);
perc_cr_amp=grpstats(resp(ttype>1),sesno(ttype>1),{'mean'});

figure
plot(unique(sesno),perc_cr_amp)
xlabel('sessions')
ylabel('CR%')
ylim([0 100]);
print ([num2str(mouseno),'CRamplitudestd'],'-depsc')
%% ISSUE MAGNET
figure
plot(t,tracesnorm(amplitudesCRs<.1&sesno==10,:))

%% save statistics
save([num2str(mouseno),'stats'],'tracesnorm','blockno','ttype','mouseno','sesno',...
     't','perc_cr10','perc_cr15','perc_cr_amp','resp','CR15','CR10', 'peakUR','amplitudesCRs');
%% plot hist amplitude CRs
hold all
bins=[-.5:.05:1.5];
plot_hist(amplitudesCRs(muscimol==1 & ttype>1),bins, 'r');
plot_hist(amplitudesCRs(muscimol==2 & ttype>1),bins,'r','Linewidth',2);
plot_hist(amplitudesCRs(saline==1 & ttype>1),bins,'b');
plot_hist(amplitudesCRs(saline==2 & ttype>1),bins,'b','Linewidth',2);
histogram(amplitudesCRs(ttype==1),bins,'k');
histcounts(amplitudesCRs(ttype==1),'k');
legend('muscimol bs','muscimol','saline bs','saline','US alone')

% sess_musc=sesno(min(find(muscimol==2)));
% last_two=sessions(sess_musc-2:sess_musc-1);
% plot_hist(amplitudesCRs(sesno==last_two(1) & sesno==last_two(1) & ttype>1),bins);
% legend('muscimol bs','muscimol','saline bs','saline','US alone','2 sess bef musc')

print([num2str(mouseno),'CRamplitudeCumulative'],'-depsc')

%% for plotting and checking things (an example of checking one session)
plot(t,tracesnorm(sesno==7 & ttype==2,:))
%%plot(t,nanmedian(tracesnorm(sesno==16 & ttype==2,:)))

%% downsample taking care of nans
traces_ds=[];
hold all
num_nans=0;
for i=1:size(tracesnorm,1)    
    pp=tracesnorm(i,:);
    if numel(find(~isnan(pp)))>100
        pp(isnan(pp))=interp1(find(~isnan(pp)), pp(~isnan(pp)), find(isnan(pp)));
        traces_ds=[traces_ds; decimate(pp,10)];
    else
        num_nans=num_nans+1;
        disp(numel(find(~isnan(pp))))
        traces_ds=[traces_ds; downsample(pp,10)];
        disp(['downsampling without low-pass vector' num2str(i)]); 
        plot(pp)

    end
end
%% plot all lines with median downsampled
figure
plot(t_ds,traces_ds(sesno==2 & ttype==2,:),'Linewidth',0.1,'Color',[0.8,0.8,0.8])
hold on
plot(t_ds,nanmedian(traces_ds(sesno==2 &ttype==2,:)),'r','Linewidth',2);

print ([num2str(mouseno),'saline'],'-depsc')
%% plot all lines with median not downsampled
figure
plot(t,tracesnorm(sesno==2 & ttype==2,:),'Linewidth',0.1,'Color',[0.8,0.8,0.8])
hold on
plot(t,nanmedian(tracesnorm(sesno==16 &ttype>1,:)),'r','Linewidth',2);




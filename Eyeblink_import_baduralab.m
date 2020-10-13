%% clear workspace
clearvars;

%% Import eyeblink data from Excel file
datafolder='Z:\Maria\Eyeblink Behavior Data\eyeblink_fabian'; %for Maria's PC

mouseno=437; %mouse id
fsamp=250; %sampling frequency

%% Here we are saying which kind of file to open and how
xlsfname=sprintf('%s%s%03d%s',datafolder,[filesep 'Mouse'],mouseno,'.xlsx');
savefname=sprintf('%s%s%03d%s',datafolder,[filesep 'Mouse'],mouseno,'allvariables','.mat');

%% Defining variables
[A,B]=xlsread(xlsfname); 
coltitles=B{1,:}; %column titles are row 1 for all columns
B(1,:)=[]; %asign row 1 a null element
n=size(A,1); %define n

%% Remove empty trials
sel=true(n,1);
for a=1:n
    if isempty(B{a,7})
        sel(a)=false; 
    end
end
%Redifine variables after removing empty trials
A=A(sel,:); 
B=B(sel,:);
n=size(A,1);

%% select traces 
traces = A(:,14:end); %selecting only the traces from original matrix A 
m = size(traces,2); %m is the values per trail
c = 0:(m-1);
t = c/fsamp; %time depending on sampling frequency (time vector 1x2000)

%% Define session number and trial types

sesno=A(:,2); % session number

ttype=zeros(n,1); %trial type 
for a=1:n %do this for the total amount of trials (240)
    switch B{a,7} %for every trial type (puff, led and puff+led) asign 1, 2 or 3.
        case 'puff1'
            ttype(a)=1;
        case 'led1_puff2'
            ttype(a)=2;
        case 'led1'
            ttype(a)=3;
    end
end

%% Important variables to save
% Make sure that the selection is correct depending on how
% the excel file is organized.

blockno=A(:,6); 
CSstart=A(:,8);
CSduration=A(:,9);
USstart=A(:,12);
%USduration=A(:,13);
USduration=A(:,11); %change here bc excel files are different. US duration is in column 11 instead of 13.
clear A B;

save(savefname,'traces','blockno','ttype','mouseno','sesno'...
    ,'CSstart','USstart','CSduration','USduration','n','t','m');


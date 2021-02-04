% FFT_EEG_spectrum analysis (24hr,w/ normalization,REM>=40sec,pure R/W/NR epochs)
%
% more close to YH's method:
% (only count continuous 3 epochsm and delete its first & last epoch)
%
% My plan:
%   (M: the epochs with EEG artifacts)
%  [REM]NR,W,M_interruption->R, extract R>9 epochs, delete the 1st/last/NR_int 
%  [Wake,exclude W_int][NR,exclude NR_int]follow YH method

function[F1,F2]=FFTspec5(FileName,opts)
tic
% FileName=('C_181215Sat_181216Sun_A2AR.0002_Ch2_AE_cFFT.txt');
% opts=delimitedTextImportOptions('EmptyLineRule','read','VariableNamesLine',1,'DataLines',[2,Inf]);

T=readtable(FileName,opts); 
[r,c]=size(T); 
T=T(1:(r-1),1:(c-1)); %removing the extra one column and one row in T when read from .txt file
[r,c]=size(T); 
% r:row(stage), c:column(hz)
% stage starts from row-19; hz starts from column_4~65

stage1=zeros(r-18,1); % extract stages
stage1=table2array(T(19:r,2)); % table ->cell
[s,a]=size(stage1);% s:number of all stages/epochs, a=1

hz=str2double(table2array(T(19:r,4:c)));% extract hz, in double
[epoch_n,hz_n]=size(hz);
% hz_n=hz_n; % number of Hz calculating (ex. 0hz, 0.5hz, 1.0hz.....)

% %% calculating individual epoch mean from 0-30Hz (no matter W/NR/R)
% EPmean=zeros(s,1); % EPmean=mean of power(0-30Hz) in each epoch(EP:Epoch Power)
% for ii=1:s;
%     EPmean(ii,1)=mean(hz(ii,:));
% end
% EPmean_all=mean(EPmean); % mean of all epoch mean

%% calculating sum of each epoch from 0-30Hz (no matter W/NR/R)

% Find epochs with AE (M) to exclude from calculting normalization value(191127)
stage1_M=strrep(stage1,'M','1'); % cell
stage1_M=strrep(stage1_M,'W','0'); % cell
stage1_M=strrep(stage1_M,'NR','0'); % cell
stage1_M=strrep(stage1_M,'R','0'); % cell
stage1_M=str2double(stage1_M); % cell -> double

M_epoch=find(stage1_M==1);% epochs in hz should NOT be calculated
M_epoch_n=length(M_epoch); % number of M epochs
non_M_epoch=find(stage1_M==0);% epochs in hz should be calculated
n_calc_epoch=length(non_M_epoch); % number of non-M epochs should be used to calculate normalization value

EPsum=zeros(n_calc_epoch,1); %EPsum:sum of power (0-30Hz) in each epoch(EP:Epoch Power)
for ii=1:n_calc_epoch;
    EPsum(ii,1)=sum(hz((non_M_epoch(ii,1)),:));
end
EPsum_mean=mean(EPsum); % mean of epoch sum for normalization

%% replace stage with number
stage2=strrep(stage1,'W','0'); % cell
stage2=strrep(stage2,'NR','2'); % cell
stage2=strrep(stage2,'R','3'); % cell
stage2=strrep(stage2,'M','5'); % cell, M: artifacts epoch
stage2=str2double(stage2); % cell -> double

%% Calculating REM (only>=40sec)& delete 1st/last/NR_int
% replace 1 NR epoch interrupt -> R (result in stage6)
stage3_R=stage2; % stage2(before), stage3(after) NR,W,M replacement

NR_int=[]; % location of interrupted NR epoch
for ii=1:(s-1)
    if stage2(ii+1,1)-stage2(ii,1)==-1 & stage2(ii+1,1)-stage2(ii+2,1)==-1
        stage3_R(ii+1,1)=3;
        NR_int=[NR_int ii+1];
    end
end
NR_int=(NR_int)'+18;
n_NR_int=length(NR_int); %number of interrupted NR between REM

% replace W_int in REM episode
W_int=[];
for ii=1:(s-1)
    if stage2(ii+1,1)-stage2(ii,1)==-3 & stage2(ii+1,1)-stage2(ii+2,1)==-3
        stage3_R(ii+1,1)=3; 
        W_int=[W_int ii+1];
    end
end
W_int=(W_int)'+18; 
n_W_int=length(W_int);%number of interrupted W between REM

% replace M_int in REM episode
M_int=[];
for ii=1:(s-1)
    if stage2(ii+1,1)-stage2(ii,1)==5 & stage2(ii+1,1)-stage2(ii+2,1)==5
        stage3_R(ii+1,1)=3; 
        M_int=[M_int ii+1];
    end
end
M_int=(M_int)'+18;
n_M_int=length(M_int);%number of interrupted M between REM


% Distinguish REM from the others (all length)
R_epoch_all=stage3_R>2;
R_epoch_all=double(R_epoch_all); %logic -> double

[e,b]=find(R_epoch_all==1); %(e=no. of total R episode; b=1)
count=[e];
n_REMepoch=length(count); % total REM epoch number

% Find REM episode (extract REM episode, all length)
REMepisode=1;
Lcount(1)=count(1);
for jj=1:length(count)-1;
    if count(jj+1)-count(jj)>1
        REMepisode= REMepisode+1;
        
        Lcount(REMepisode)=count(jj+1);
        Ucount(REMepisode-1)=count(jj);
        REMepisode_length(REMepisode-1)=Ucount(REMepisode-1)-Lcount(REMepisode-1)+1;
        %episode_length=epoch numbers in each REM episode
    end
end
% define the last one;
Ucount(REMepisode)=count(length(count));
REMepisode_length=[REMepisode_length]';
n_REMepisode=REMepisode;

% Extract REM epoisode >= 40sec (which location information shows in long_REM_episode)
[r,d]=find(REMepisode_length>9); % r=no. of (long) REM episode
n_longREMepisode=length(r); % n_longREMepisode: total number of (long) REM episode
n_longREMepoch=sum(REMepisode_length(r)); % n_longREMepisode: total number of (long) REM epoch

% Extract long REM epochs
long_REM_epoch=[];
for hh=1:length(r);
    gg=r(hh);
    long_REM_epoch=[long_REM_epoch Lcount(gg):Ucount(gg)];
end
long_REM_epoch=(long_REM_epoch)'+18;

% delete from calculation: 1st/last R epoch & NR_int & W_int & M_int
R_epoch_del=([Lcount+18 Ucount+18 (NR_int)' (NR_int-1)' (NR_int+1)' (W_int)' (W_int-1)' (W_int+1)' (M_int)' (M_int-1)' (M_int+1)'])';% R_epoch_del: R epoch to be delete from calculation
R_epoch_del=sort(R_epoch_del,'ascend');
R_epoch_del=intersect(long_REM_epoch, R_epoch_del);% intersect: keep common elements without repetition 

R_epoch_calc=setdiff(long_REM_epoch, R_epoch_del);% R_epoch_cal: R epochs to calculate(returns the data in A that is not in B, with no repetitions)


% calculate longREM epoch mean in individual Hz range

Rm_hz=zeros(1,hz_n);
for hh=1:hz_n;
    C=hz((R_epoch_calc)-18,hh);
    Rm_hz(1,hh)=mean(C);
end

% Rm_hz=zeros(1,c-3);
% T_calc=str2double(table2array(T));
% for hh=4:c;
%     C=T_calc(R_epoch_calc,hh);
%     Rm_hz(1,hh-3)=mean(C);
% end

% % normalize with EPmean_all
% Rm_hz_nor=(Rm_hz)/(EPmean_all);

% normalize with EPsum_mean
Rm_hz_nor=(Rm_hz)/(EPsum_mean);
Rm_hz_nor=Rm_hz_nor*100;

% figure
x_axis=table2array(T(18,4:c)); 
X=strrep(x_axis,'Hz','0');
X=str2double(X);
X=[X]';

figure
subplot(4,1,1)
plot(X,Rm_hz_nor,'k')
ylabel('Power')
title('REM')
ax=gca;
ax.YLim=[0 20]
ax.XTick=[0:2:30]

%% Calculating Wake
% distinguish W from others (all length)
stage4_W=strrep(stage1,'W','1');
stage4_W=strrep(stage4_W,'NR','0');c
stage4_W=strrep(stage4_W,'R','0');
stage4_W=str2double(stage4_W); % cell->double

% extract the calculated W epochs
W_epoch_calc=[];
for ii=1:length(stage4_W)-2;
    if stage4_W(ii,1)* stage4_W(ii+1,1)* stage4_W(ii+2,1)==1
        W_epoch_calc=[W_epoch_calc (ii+1)];
    end
end
W_epoch_calc=(W_epoch_calc+18)'; % calculated W epochs

% calculate Wake epoch mean in individual Hz range

Wm_hz=zeros(1,hz_n);
for hh=1:hz_n;
    C=hz((W_epoch_calc)-18,hh);
    Wm_hz(1,hh)=mean(C);
end

% Wm_hz=zeros(1,c-3);
% for ii=4:c;
%     C=T_calc(W_epoch_calc,ii);
%     Wm_hz(1,ii-3)=mean(C);
% end

% % normalize with EPmean_all
% Wm_hz_nor=(Wm_hz)/(EPmean_all);

% normalize with EPsum_mean
Wm_hz_nor=(Wm_hz)/(EPsum_mean);
Wm_hz_nor=Wm_hz_nor*100;

subplot(4,1,2)
plot(X,Wm_hz_nor,'b')
ylabel('Power')
title('Wake')
ax=gca;
ax.YLim=[0 20]
ax.XTick=[0:2:30]


%% Calculating NREM
% distinguish NR from others (all length)
stage5_NR=strrep(stage1,'NR','1');
stage5_NR=strrep(stage5_NR,'W','0');
stage5_NR=strrep(stage5_NR,'R','0');
stage5_NR=str2double(stage5_NR); % cell->double

% extract the caclualted NR epochs
NR_epoch_calc=[];
for ii=1:length(stage5_NR)-2;
    if stage5_NR(ii,1)*stage5_NR(ii+1,1)*stage5_NR(ii+2,1)==1
        NR_epoch_calc=[NR_epoch_calc (ii+1)];
    end 
end
NR_epoch_calc=(NR_epoch_calc+18)'; % calculated NR epochs

% calculate NREM epoch mean in individual Hz range

NRm_hz=zeros(1,hz_n);
for hh=1:hz_n;
    C=hz((NR_epoch_calc)-18,hh);
    NRm_hz(1,hh)=mean(C);
end

% NRm_hz=zeros(1,c-3);
% for ii=4:c;
%     C=T_calc(NR_epoch_calc,ii);
%     NRm_hz(1,ii-3)=mean(C);
% end

% % normalize with EPmean_all
% NRm_hz_nor=(NRm_hz)/(EPmean_all);

% normalize with EPsum_mean
NRm_hz_nor=(NRm_hz)/(EPsum_mean);
NRm_hz_nor=NRm_hz_nor*100;

subplot(4,1,3)
plot(X,NRm_hz_nor,'r')
ylabel('Power')
title('NREM')
ax=gca;
ax.YLim=[0 20]
ax.XTick=[0:2:30]

%% plot final figure and save figure
subplot(4,1,4)
% plot(X,Rm_hz_nor,'k',X,Wm_hz_nor,'b',X,NRm_hz_nor,'r')
plot(X,Rm_hz_nor,'LineWidth',1.75,'Color','k')
hold on
plot(X,Wm_hz_nor,'LineWidth',1.75,'Color','b')
hold on
plot(X,NRm_hz_nor,'LineWidth',1.75,'Color','r')
xlabel('Freq (Hz)')
ylabel('Power')
title('All')
ax=gca;
ax.YLim=[0 20]
ax.XTick=[0:2:30]

FileName1=erase(FileName,'.xlsx'); %character:for file name
FileName2=string(erase(FileName,'.xlsx')); %to become string:for put into csv
saveas(gcf,strcat(FileName1,'.tif'))

%% create results file
FileName2(2:(c-3))=string(' ');
X=num2cell(X);
Rm_hz_nor=num2cell(Rm_hz_nor);
Wm_hz_nor=num2cell(Wm_hz_nor);
NRm_hz_nor=num2cell(NRm_hz_nor);

F1= [FileName2
    X'
    Rm_hz_nor
    Wm_hz_nor
    NRm_hz_nor]; % (F1)for File-1

% writetable(F1,strcat(FileName1,'.csv'))

% creat general information file
n_Lcount=length(Lcount);
n_Ucount=length(Ucount);
n_Repoch_del=length(R_epoch_del);
n_Repoch_calc=length(R_epoch_calc);
n_Wepoch_calc=length(W_epoch_calc);
n_NRepoch_calc=length(NR_epoch_calc);

toc
timeVal=toc;

% F2=[EPmean_all;
%     NaN; n_REMepisode; n_longREMepisode; 
%     NaN; n_REMepoch; n_longREMepoch; 
%     NaN; n_Lcount; n_Ucount;
%     NaN; n_NR_int; n_W_int;
%     NaN; n_Repoch_del; n_Repoch_calc;
%     NaN; n_Wepoch_calc;
%     NaN; n_NRepoch_calc;
%     NaN; timeVal];

F2=[EPsum_mean;
    NaN; n_REMepisode; n_longREMepisode; 
    NaN; n_REMepoch; n_longREMepoch; 
    NaN; n_Lcount; n_Ucount;
    NaN; n_NR_int; n_W_int; n_M_int;
    NaN; n_Repoch_del; n_Repoch_calc;
    NaN; n_Wepoch_calc;
    NaN; n_NRepoch_calc;
    NaN; timeVal];

F2=num2cell(F2);

end


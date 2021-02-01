% FFT_EEG_spectrum analysis (24hr,w/ normalization,REM>=40sec,pure R/W/NR epochs)
%
% YH's method:
% (only count continuous 3 epochs and delete its first & last epoch)
%
% FFTspec4: for analyzing .txt file!

function[F1,F2]=FFTspec6(FileName,opts)
tic
% FileName=('C_181215Sat_181216Sun_A2AR.0002_Ch3_AE_cFFT.txt');
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

%% Calculating REM (all the continuous REM epochs, all length)
stage2_R=strrep(stage1,'W','0'); % cell
stage2_R=strrep(stage2_R,'NR','0'); % cell
stage2_R=strrep(stage2_R,'R','1'); % cell
stage2_R=str2double(stage2_R); % cell -> double

%% extract the calculated R epochs
R_epoch_calc=[];
for ii=1:length(stage2_R)-2;
    if stage2_R(ii,1)* stage2_R(ii+1,1)* stage2_R(ii+2,1)==1
        R_epoch_calc=[R_epoch_calc (ii+1)];
    end
end
% R_epoch_calc=(R_epoch_calc+18)';% calculated R epochs
R_epoch_calc=(R_epoch_calc)';
R_epoch_calc_LP=R_epoch_calc(R_epoch_calc<10801);
R_epoch_calc_DP=R_epoch_calc(R_epoch_calc>10800);

% calculate REM epoch mean in individual Hz range
Rm_hz=zeros(1,hz_n);
Rm_hz_LP=zeros(1,hz_n);
Rm_hz_DP=zeros(1,hz_n);

for ii=1:hz_n;
    C=hz(R_epoch_calc,ii);
    C_LP=hz(R_epoch_calc_LP,ii);
    C_DP=hz(R_epoch_calc_DP,ii);
    
    Rm_hz(1,ii)=mean(C);
    Rm_hz_LP(1,ii)=mean(C_LP);
    Rm_hz_DP(1,ii)=mean(C_DP);
end
% % normalize with EPmean_all
% Rm_hz_nor=(Rm_hz)/(EPmean_all);

% normalize with EPsum_mean
Rm_hz_nor=((Rm_hz)/(EPsum_mean))*100;
Rm_hz_nor_LP=((Rm_hz_LP)/(EPsum_mean))*100;
Rm_hz_nor_DP=((Rm_hz_DP)/(EPsum_mean))*100;

% figure
x_axis=table2array(T(18,4:c)); 
X=strrep(x_axis,'Hz','0');
X=str2double(X);
X=[X]';

figure
subplot(4,1,1)
plot(X,Rm_hz_nor,'k', 'LineWidth',1.0)

hold on
plot(X,Rm_hz_nor_LP,'Color',[0.8500 0.3250 0.0980],'LineWidth',0.75,'LineStyle','-.')
plot(X,Rm_hz_nor_DP,'Color',[0 0.4470 0.7410],'LineWidth',0.75,'LineStyle','-.')

ylabel('Power')
title('REM')
ax=gca;
ax.YLim=[0 20]
ax.XTick=[0:2:30]
legend('All','LP','DP')

hold off

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
W_epoch_calc=(W_epoch_calc)';
W_epoch_calc_LP=W_epoch_calc(W_epoch_calc<10801);
W_epoch_calc_DP=W_epoch_calc(W_epoch_calc>10800);

% calculate Wake epoch mean in individual Hz range
Wm_hz=zeros(1,hz_n);
Wm_hz_LP=zeros(1,hz_n);
Wm_hz_DP=zeros(1,hz_n);

for ii=1:hz_n;
    C=hz(W_epoch_calc,ii);
    C_LP=hz(W_epoch_calc_LP,ii);
    C_DP=hz(W_epoch_calc_DP,ii);
    
    Wm_hz(1,ii)=mean(C);
    Wm_hz_LP(1,ii)=mean(C_LP);
    Wm_hz_DP(1,ii)=mean(C_DP);
end


% normalize with EPsum_mean
Wm_hz_nor=((Wm_hz)/(EPsum_mean))*100;
Wm_hz_nor_LP=((Wm_hz_LP)/(EPsum_mean))*100;
Wm_hz_nor_DP=((Wm_hz_DP)/(EPsum_mean))*100;

subplot(4,1,2)
plot(X,Wm_hz_nor,'b','LineWidth',1.0)

hold on
plot(X,Wm_hz_nor_LP,'Color',[0.8500 0.3250 0.0980],'LineWidth',0.75,'LineStyle','-.')
plot(X,Wm_hz_nor_DP,'Color',[0 0.4470 0.7410],'LineWidth',0.75,'LineStyle','-.')

ylabel('Power')
title('Wake')
ax=gca;
ax.YLim=[0 20]
ax.XTick=[0:2:30]
legend('All','LP','DP')

hold off

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
NR_epoch_calc=(NR_epoch_calc)';
NR_epoch_calc_LP=NR_epoch_calc(NR_epoch_calc<10801);
NR_epoch_calc_DP=NR_epoch_calc(NR_epoch_calc>10800);

% calculate NREM epoch mean in individual Hz range
NRm_hz=zeros(1,hz_n);
NRm_hz_LP=zeros(1,hz_n);
NRm_hz_DP=zeros(1,hz_n);

for ii=1:hz_n;
    C=hz(NR_epoch_calc,ii);
    C_LP=hz(NR_epoch_calc_LP,ii);
    C_DP=hz(NR_epoch_calc_DP,ii);
    
    NRm_hz(1,ii)=mean(C);
    NRm_hz_LP(1,ii)=mean(C_LP);
    NRm_hz_DP(1,ii)=mean(C_DP);
end
NRm_hz_nor=((NRm_hz)/(EPsum_mean))*100;
NRm_hz_nor_LP=((NRm_hz_LP)/(EPsum_mean))*100;
NRm_hz_nor_DP=((NRm_hz_DP)/(EPsum_mean))*100;

subplot(4,1,3)
plot(X,NRm_hz_nor,'r','LineWidth',1.75)

hold on
plot(X,NRm_hz_nor_LP,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.0,'LineStyle','-.')
plot(X,NRm_hz_nor_DP,'Color',[0 0.4470 0.7410],'LineWidth',1.0,'LineStyle','-.')

ylabel('Power')
title('NREM')
ax=gca;
ax.YLim=[0 20]
ax.XTick=[0:2:30]
legend('All','LP','DP')

hold off

%% plot final figure and save figure
subplot(4,1,4)
% plot(X,Rm_hz_nor,'k',X,Wm_hz_nor,'b',X,NRm_hz_nor,'r')
plot(X,Rm_hz_nor,'LineWidth',1.75,'Color','k')
hold on
plot(X,Wm_hz_nor,'LineWidth',1.75,'Color','b')
plot(X,NRm_hz_nor,'LineWidth',1.75,'Color','r')

xlabel('Freq (Hz)')
ylabel('Power')
title('All')
ax=gca;
ax.YLim=[0 20]
ax.XTick=[0:2:30]
legend('Wake','NREM','REM')

hold off

FileName1=erase(FileName,'.txt'); %character:for file name
FileName2=string(erase(FileName,'.txt')); %to become string:for put into csv
saveas(gcf,strcat(FileName1,'.tif'))

%% create results file
FileName2(2:hz_n)=string(' ');%the rest of columns were empty string
X=num2cell(X);
Rm_hz_nor=num2cell(Rm_hz_nor);
Wm_hz_nor=num2cell(Wm_hz_nor);
NRm_hz_nor=num2cell(NRm_hz_nor);

Rm_hz_nor_LP=num2cell(Rm_hz_nor_LP);
Wm_hz_nor_LP=num2cell(Wm_hz_nor_LP);
NRm_hz_nor_LP=num2cell(NRm_hz_nor_LP);

Rm_hz_nor_DP=num2cell(Rm_hz_nor_DP);
Wm_hz_nor_DP=num2cell(Wm_hz_nor_DP);
NRm_hz_nor_DP=num2cell(NRm_hz_nor_DP);

F1= [FileName2
    X'
    Rm_hz_nor
    Wm_hz_nor
    NRm_hz_nor
    Rm_hz_nor_LP
    Wm_hz_nor_LP
    NRm_hz_nor_LP
    Rm_hz_nor_DP
    Wm_hz_nor_DP
    NRm_hz_nor_DP]; % (F1)for File-1

% writetable(F1,strcat(FileName1,'.csv'))


%% creat general information file
n_Repoch_calc=length(R_epoch_calc);
n_Wepoch_calc=length(W_epoch_calc);
n_NRepoch_calc=length(NR_epoch_calc);

n_Repoch_calc_LP=length(R_epoch_calc_LP);
n_Wepoch_calc_LP=length(W_epoch_calc_LP);
n_NRepoch_calc_LP=length(NR_epoch_calc_LP);

n_Repoch_calc_DP=length(R_epoch_calc_DP);
n_Wepoch_calc_DP=length(W_epoch_calc_DP);
n_NRepoch_calc_DP=length(NR_epoch_calc_DP);

toc
timeVal=toc;

% F2=[EPmean_all;
%     NaN; n_Repoch_calc;
%     NaN; n_Wepoch_calc;
%     NaN; n_NRepoch_calc;
%     NaN; timeVal];

F2=[EPsum_mean;
    NaN; n_Repoch_calc;
    n_Repoch_calc_LP;
    n_Repoch_calc_DP;
    NaN; n_Wepoch_calc;
    n_Wepoch_calc_LP;
    n_Wepoch_calc_DP;
    NaN; n_NRepoch_calc;
    n_NRepoch_calc_LP;
    n_NRepoch_calc_DP;
    NaN; timeVal];

F2=num2cell(F2);

end


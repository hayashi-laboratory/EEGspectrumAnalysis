% For automatic Run all the files in a folder
File=dir('*.txt');
n=length(File);

%specify .txt reading porperties to meet the same pattern as reading.xlsx files
opts=delimitedTextImportOptions('EmptyLineRule','read','VariableNamesLine',1,'DataLines',[2,Inf]);

% information for output 
Results1=[]; %for output results(F1) into csv file
R1_column_title=({'File Name','Hz range','24h REM power','24h Wake power','24h NREM power','LP REM power','LP Wake power','LP NREM power','DP REM power','DP Wake power','DP NREM power'});

Results2=[]; %for output general info.(F2) into csv file 
R2_column_title=({'normalization value','','calculated 24h R epochs','calculated LP R epochs','calculated DP R epochs','','calculated W epochs','calculated LP W epochs','calculated DP W epochs','','calculated NR epochs','calculated LP NR epochs','calculated DP NR epochs','','analysis time'})';
R2_row_title={'',File.name};
for ii=1:n
    w=waitbar(ii/n,strcat('Please wait! Analysing file no. ',num2str(ii)))
    
    [F1,F2]=FFTspec6(File(ii).name,opts);
    
    F1=[R1_column_title;F1'];
    
    Results1=[Results1 F1];
    
    Results2=[Results2 F2];  
   
    close(w)
end


Results2=[R2_column_title Results2];
Results2=[R2_row_title
    Results2];

Results1=array2table(Results1);
Results2=cell2table(Results2);

writetable(Results1,strcat('01_FFT_Results_',date,'.csv'))
writetable(Results2,strcat('02_General Analysis Info_',date,'.csv'))




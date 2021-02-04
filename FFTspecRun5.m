% For automatic Run all the files in a folder
File=dir('*.txt');
n=length(File);

%specify .txt reading porperties to meet the same pattern as reading.xlsx files
opts=delimitedTextImportOptions('EmptyLineRule','read','VariableNamesLine',1,'DataLines',[2,Inf]);

Results1=[]; %for output results(F1) into csv file
R1_column_title=({'File Name','Hz range','long REM power','Wake power','NREM power'});

Results2=[]; %for output general info.(F2) into csv file 
R2_column_title=({'normalization value','','total REM episode','long REM episode','','total REM epochs','long REM epochs','','Lcount','Ucount','','NR interruption','W interruption','M interruption between REM','','deleted R epochs','calculated R epochs','','calculated W epochs','','calculated NR epochs','','analysis time'})';
R2_row_title={'',File.name};
for ii=1:n
    w=waitbar(ii/n,strcat('Please wait! Analysing file no. ',num2str(ii),', out of total:','',num2str(n)))
    
    [F1,F2]=FFTspec5(File(ii).name,opts);
    
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

writetable(Results1,strcat('01_FFT_Results_longREM_',date,'.csv'))
writetable(Results2,strcat('02_General Analysis Info_',date,'.csv'))
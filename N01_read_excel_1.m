%%%%%%%%%% ��excel����ÿһ������Ԫ���ݶ�ȡ��matlab��������mat�ļ���ʽ�洢
clear
clc

for nian=2008:2008
    
    file=['../WOD_',num2str(nian),'.xlsx']
    
    [DNA_series,txt,~]=xlsread(file);  %һ����Ҫ��cell���飬Ҫ��Ȼ���Կ���
    DNA_series=single(DNA_series);
    
    %%%�����ļ��� ���к�
    filename_info=char(txt(:,1));
    
    
    %%%�ֶ�����
    % variable_name=raw(1,2:end);
    variable_name={'accessin_number','dataset_id','lat','lon','year','month','day','probe_type','recorder','hour','minute','depth_number','maximum_depth','hasTemp','hasSal','hasOxygen','hasChlonophyll','country_id','GMT_time','WMO_id','dbase_orig','Project_name','plarform','vehicle','Institute','WOD_cruise_identifier','sum_temp','sum_salinity','sum_depth','std_depth','std_temp','std_salinity','corr_temp_depth','corr_sal_depth'};
    
    DNA_series(:,1)=[]; %ɾ��WOD_unique_id
    
    txt(:,1)=[];  %ɾ���ļ���
    txt(:,2)=[]; %ɾ��WOD_unique_id
    
    
    
    %%%%���ַ���ת��ΪASCILL�룬Ȼ�����
    istxt=all(isnan(DNA_series));
    istxt(26)=1;  %26. WOD_cruise_identifier���ַ���
    for i=1:length(istxt)
        i
        if(istxt(i))
            column=i;
            ascill_all=char(txt(:,i))+0;
            %�����пո�����Ϊnan
            ascill_all(ascill_all==' ')=NaN;
            %ÿһ���ַ������
            sum_ascill_all=sum(ascill_all,2,'omitnan');
            DNA_series(:,i)=sum_ascill_all;
        end
    end
    
    DNA_series(DNA_series==999)=NaN;  %%%%999 ȱ��ֵ����ΪNan
    
    %%%%%�������������ˣ�������Ԫ������Ϣ���������ֵ���ʽ�洢����һ����������У�����ÿһ�д���һ�����ߣ���Ӧ���ļ�����filename_info�洢�������Ϳ���һһ��Ӧ�ļ���
    eval(['save DNA_summary_',num2str(nian),'.mat DNA_series variable_name filename_info'])
end



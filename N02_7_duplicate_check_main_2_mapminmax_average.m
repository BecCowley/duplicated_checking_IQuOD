%%%%  ����ÿһ�й�һ��������Ȼ��������ƽ�������ž��������Ƚ��ĸ����ӽ�
%%%% ʹ��mapstd��Ԫ������Ϣ���б�׼��������Ȼ��������ƽ��
%%%  �������Ǳ�ڵġ����ܡ����ơ��ظ����ݶԵ��ļ�����һ��txt�ļ���
clear
clc

for nian=1995:1995
    
    eval(['load DNA_summary_',num2str(nian),'.mat'])
    
    %��׼������mapstd  �����ݱ�׼��Ϊ��ֵΪ0������Ϊ1������    mapminmax��һ��
    %��һ���ͱ�׼����������̽��
    DNA_series_copy=DNA_series;
    DNA_series_copy(:,20)=[]; %��20�� WMO_ID��Ϣ�����٣�ȥ��
    
    %���Զ�ÿһ�н��й�һ���������������ٵ�Ӱ��
    DNA_mapped_1=mapminmax(DNA_series_copy',0,1);
    DNA_mapped=DNA_mapped_1';
    DNA_mapped(:,5)=0;
    DNA_mapped(DNA_series_copy==0)=0;
    
    %��ÿһ�н��б�׼������
%     DNA_mapped_1=mapstd(DNA_series_copy',0,1);
%     DNA_mapped=DNA_mapped_1';
% %    DNA_mapped(:,5)=0;
%     DNA_mapped(DNA_series_copy==0)=0;
    
    %ÿһ��������ƽ��  -->������Լ�Ȩƽ��   ������Ծ���γ��Ȩ�ؼӴ�
    average_DNA=nanmean(DNA_mapped,2);
    
    %��average_DNA�����������򣬷�����������㷨�Ľ���
    [average_DNA,index]=sort(average_DNA);
    filename_info=filename_info(index,:);
    DNA_mapped=DNA_mapped(index,:);
    DNA_series=DNA_series(index,:);
    
    %%% ѭ������
    output_variables=['filename',variable_name];
    
    filename=['./potential_duplicates_output/',num2str(nian),'/potential_duplicat_mapminmax_',num2str(nian),'.txt'];
    if(exist(filename))
        delete(filename)
    end
    fid=fopen(filename,'w+');
    % for i=1:length(output_variables)
    %     fprintf(fid,'%s ',output_variables{i})
    % end
    % output_filename='potential_duplicates.xlsx';
    
    number_pairs=0;
    number_profiles=0;
    for i=1:length(average_DNA)
        i
        number1=average_DNA(i);
        difference=abs((number1-average_DNA)/number1*100);   %����ٷֱ�
        difference(1:i-1)=NaN;
        duplicate_number=sum(difference<0.0001);   %��ֵ0.01%      ��ֵ������֮���趨����
        if(duplicate_number>=2)
            %%%%�����ظ�
            %        pause
            %��difference�����С�����Ϊ0��
            difference(i)=NaN;
            id=[i;find(difference==nanmin(difference))];
            DNA_series_small=DNA_series(id,:);
            
            %%%%%%%%%%%%%%%%%%%����Ǹ���SUR/MRB/DRB���ݣ������������
            if(DNA_series(i,2)==5||DNA_series(i,2)==7||DNA_series(i,2)==9) %�ֱ����SUR MRB DRB����
                continue
            end
            
            %%%%%%�����Ȼ����¶Ȼ����ζ�ͬʱ���ϼ������ģ����߳�����С�Ŵ�һ�������ģ�ֱ����� ���ϵ���϶�����ȵ� %%%%��ʵ�������ȥ������һ��ģ��
            %����
            if(any(abs(DNA_series_small(1,[33,34])-DNA_series_small(2,[33,34]))<1e-4)) %���ϵ�����
                %             fprintf(fid,'%s\n','��Potential Duplicates pairs��:');
                %             for m=1:length(id)
                %                 fprintf(fid,'%s ',filename_info(id(m),:));
                %                 fprintf(fid,'%3d %.4f %.4f %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %.4f %.4f %3d %.4f %.4f %.4f %.4f %.4f\n',DNA_series(id(m),:));
                %             end
                %%%%%���ԭʼ�����ļ����ļ���
                for m=1:length(id)
                    fprintf(fid,'%s ',filename_info(id(m),:));
                end
                fprintf(fid,'\n');
                
                number_pairs=number_pairs+1;
                number_profiles=number_profiles+duplicate_number;
                continue
            end
            
            %%%%%%%%%����һ��һЩ�ų����ж� �����ж�����Ƭ���ж��ٸ���
            fragment_same_number=sum(abs(DNA_series_small(1,:)-DNA_series_small(2,:))<1e-5,'omitnan');
            if(fragment_same_number<25)  %û�����Ƶ�Ƭ�Σ���Ҫ������  %%%%%%%������Էּ�������׼ȷ�ظ����ϸ����32 ����31  �ĳ�27�����ҳ����󲿷�׼ȷ�ظ�
                continue   %25 26��
            end
            
            %%%%�����XBT CTD MBT BOT������һ�������������5�ȣ�����ͬһ��probe  �ų��ߺ������۲�
            %%%%type,platform,vehicle,����sum_temp,corr(temp,depth)��һ�������ж�Ϊͬһ������ ͬһ�����鴬/ƽ̨ �ϵĶ�ι۲�
            if((DNA_series_small(1,2)==4 && DNA_series_small(2,2)==4) || (DNA_series_small(1,2)==2 && DNA_series_small(2,2)==2) || (DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==3 && DNA_series_small(2,2)==3))
                index1=all(DNA_series_small(1,[5,6,8,23,24,26])==DNA_series_small(2,[5,6,8,23,24,26]));  %��Ҫһ��
                index2= abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.099; %sum_temp�����
                index3= abs(DNA_series_small(1,33)-DNA_series_small(2,33))>0.001;  %cor_temp_depth
                index4=any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<5) && any(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))>1e-5);
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            %%%%%�ų�����/�����㳤ʱ�������۲� ֻ��MRB Bottle SUR
            if((DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==7 && DNA_series_small(2,2)==7) || (DNA_series_small(1,2)==5 && DNA_series_small(2,2)==5))
                index1=all(DNA_series_small(1,[5,6,8,9,22,23,24])==DNA_series_small(2,[5,6,8,9,22,23,24]));  %��Ҫһ��
                index2=abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.05; %sum_temp�����
                index3=abs(DNA_series_small(1,29)-DNA_series_small(2,29))<1e-5;  %sum_depth���
                index4=all(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<0.01);  %���㣺��γ��С��0.01��
                if(index1 && index2 && index3 && index4)
                    continue
                end
            end
            
            
            %%%%%���ԭʼ�����ļ�
            %         fprintf(fid,'%s\n','��Potential Duplicates pairs��:');
            %         for m=1:length(id)
            %             fprintf(fid,'%s ',filename_info(id(m),:));
            %             fprintf(fid,'%3d %.4f %.4f %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %.4f %.4f %3d %.4f %.4f %.4f %.4f %.4f\n',DNA_series(id(m),:));
            %         end
            %         fprintf(fid,'\n');
            
            %%%%%���ԭʼ�����ļ�
            for m=1:length(id)
                fprintf(fid,'%s ',filename_info(id(m),:));
            end
            fprintf(fid,'\n');
            
            %                pause
            number_pairs=number_pairs+1;
            number_profiles=number_profiles+duplicate_number;
            
        end
    end
    
    number_pairs   %%%��ʱ��һ���㷨����ҳ��˶��ٸ��ظ���
    number_profiles
end

%%
% figure();
% plot(average_DNA2,'o');
% ylabel('Average DNA')
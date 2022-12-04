%%%%  ����ÿһ�й�һ��������Ȼ��������ƽ�������ž��������Ƚ��ĸ����ӽ�
%%%% �ص�������Ⱥͺ��¶Ⱥ���ȵ�
clear
clc

for nian=2008:2008
    
    eval(['load DNA_summary_',num2str(nian),'.mat'])
    
    
    DNA_series_temp_depth=DNA_series(:,[27,29]); %sum_temp,sum_depth
    DNA_series_temp_depth(abs(DNA_series_temp_depth)>1e6)=NaN; %%%��ȱֵmissing value
    % figure()
    % plot(DNA_series_temp_depth(:,2),'o')
    
    %%%%%%%%%%  ����Ⱥͺ��¶Ⱥ͵�ƽ��
    average_DNA=nanmean(DNA_series_temp_depth,2);
    
    %��average_DNA�����������򣬷�����������㷨�Ľ���
    [average_DNA,index]=sort(average_DNA);
    filename_info=filename_info(index,:);
    DNA_series=DNA_series(index,:);
    
    %%% ѭ������
    output_variables=['filename',variable_name];
    
    filename=['./potential_duplicates_output/',num2str(nian),'/potential_duplicate_',num2str(nian),'_depth_temp.txt'];
    if(exist(filename))
        delete(filename)
    end
    fid=fopen(filename,'w+');

    
    number_pairs=0;
    number_profiles=0;
    for i=1:length(average_DNA)
        i
        number1=average_DNA(i);
        difference=abs((number1-average_DNA)/number1*100);   %����ٷֱ�
        difference(1:i-1)=NaN;
        duplicate_number=sum(difference<0.0001);   %��ֵ0.001%      ��ֵ������֮���趨����
        if(duplicate_number>=2)
            %%%%�����ظ�
            %        pause
            %��difference�����С�����Ϊ0��
            difference(i)=NaN;
            id=[i;find(difference==nanmin(difference))];
            DNA_series_small=DNA_series(id,:);
            %%%%%%%%%%%%%%%%%%%����Ǹ���MRB���ݣ������� ���xxxxxxx
            if(DNA_series(i,2)==7)
                continue
            end
            
            if(DNA_series(i,12)<=3)  %%%%%%%%����۲����:��ȸ���С��3
                continue
            end
            
            
            %%%%%%%%%����һ��һЩ�ų����ж� �����ж�����Ƭ���ж��ٸ���
            fragment_same_number=sum(abs(DNA_series_small(1,:)-DNA_series_small(2,:))<1e-5,'omitnan');
            if(fragment_same_number<26)  %û�����Ƶ�Ƭ�Σ���Ҫ������  %%%%%%%������Էּ�������׼ȷ�ظ����ϸ����32 ����31  �ĳ�27�����ҳ����󲿷�׼ȷ�ظ�
                continue
            end
            
            %%%%�ų���Ⱥ����ܴ��
            for m=2:length(id)
                if(abs(DNA_series_small(1,29)-DNA_series_small(m,29))>1)
                    id(m)=NaN;
                end
            end
            id(isnan(id))=[];
            if(length(id)<=1)
                continue
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
            
            %%%%%%%%�ų��ζȺͲ���ȵģ��Ҳ�Ϊ0
            if((abs(DNA_series_small(1,28)-DNA_series_small(2,28))>1e-3) && (DNA_series_small(1,28)>1e-6))
                continue
            end
            
            
            %         %%%%%���ԭʼ�����ļ�
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
    
    number_pairs
end
%%
% figure();
% plot(average_DNA2,'o');
% ylabel('Average DNA')
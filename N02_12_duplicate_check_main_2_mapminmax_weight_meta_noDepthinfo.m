%%%%%% �˳�������ʶ��һЩ��ȸ�������ȼ�¼��һ�����ظ����� ������depth number, maximum depth��
%%%%  ����ÿһ�й�һ��������Ȼ��������ƽ�������ž��������Ƚ��ĸ����ӽ�


clear
clc

load('DNA_summary_1995.mat')

DNA_series_meta=DNA_series(:,[1:26]);
DNA_series_meta(:,[12,13,20])=[]; %��12 13��depth_number, maximum depth; ��20�� WMO_ID��Ϣ������
%��׼������mapstd  �����ݱ�׼��Ϊ��ֵΪ0������Ϊ1������    mapminmax��һ��
%��һ���ͱ�׼����������̽��
% DNA_mapped=mapstd(DNA_series,0,1);
%���Զ�ÿһ�н��й�һ���������������ٵ�Ӱ��
DNA_mapped_1=mapminmax(DNA_series_meta',0,1);
DNA_mapped=DNA_mapped_1';
DNA_mapped(:,5)=0;
DNA_mapped(DNA_series_meta==0)=0;

%ÿһ�н��б�׼��
% DNA_mapped_1=mapstd(DNA_series_meta',0,1);
% DNA_mapped=DNA_mapped_1';
% %DNA_mapped(:,5)=0;
% DNA_mapped(DNA_series_meta==0)=0;

%%%%%%%%%%��Ȩ��
[weight]=entropy_weight(DNA_series_meta);
figure();bar(weight)

%ÿһ��������ƽ��  -->������Լ�Ȩƽ��   ������Ծ���γ��Ȩ�ؼӴ�
% average_DNA=nanmean(DNA_mapped,2);
%%%%%%%%%%  ��Ȩƽ��
average_DNA_single=NaN(size(DNA_mapped));
for i=1:length(weight)
    average_DNA_single(:,i)=DNA_mapped(:,i)*weight(i);
end
average_DNA=sum(average_DNA_single,2,'omitnan');


%��average_DNA�����������򣬷�����������㷨�Ľ���
[average_DNA,index]=sort(average_DNA);
filename_info=filename_info(index,:);
DNA_mapped=DNA_mapped(index,:);
DNA_series=DNA_series(index,:);
DNA_series_meta=DNA_series_meta(index,:);

% figure();plot(average_DNA,'o')
%%% ѭ������
output_variables=['filename',variable_name];
filename='./potential_duplicates_output/1995/potential_duplicate_1995_mapminmax_weight_meta_noDepthInfo.txt';
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
    duplicate_number=sum(difference<0.0001);   %��ֵ0.001%      ��ֵ������֮���趨����
    if(duplicate_number>=2)
        %%%%�����ظ�
        %        pause
        %��difference�����С�����Ϊ0��
        difference(i)=NaN;
        id=[i;find(difference==nanmin(difference))];
        DNA_series_small=DNA_series(id,:);
        DNA_series_small_meta=DNA_series_meta(id,:);
        %%%%%%%%%%%%%%%%%%%����Ǹ���MRB���ݣ������� ���xxxxxxx
%         if(DNA_series(i,1)==7)
%             continue
%         end
         
        %%%%%%%%%����һ��һЩ�ų����ж� �����ж�����Ƭ���ж��ٸ���
        %û�����Ƶ�Ƭ�Σ���Ҫ������
        fragment_same_number=sum(abs(DNA_series_small_meta(1,:)-DNA_series_small_meta(2,:))<1e-5,'omitnan');
        if(DNA_series(i,2)==7 || DNA_series(i,2)==9 || DNA_series(i,2)==5)  %��������խһЩ DRB MRB SUR
            if(fragment_same_number<23)
                continue
            end
        else
            if(fragment_same_number<20)  %%%����ֻ��20��Ƭ�������Ƶģ�������23��Ƭ��
                continue
            end
        end
       
        %%%%%�ų�����/�����㳤ʱ�������۲� ֻ��MRB Bottle SUR
        if((DNA_series_small(1,2)==1 && DNA_series_small(2,2)==1) || (DNA_series_small(1,2)==7 && DNA_series_small(2,2)==7) || (DNA_series_small(1,2)==5 && DNA_series_small(2,2)==5))
            index1=all(DNA_series_small(1,[5,6,8,9,22,23,24])==DNA_series_small(2,[5,6,8,9,22,23,24]));  %��Ҫһ��
            index2=   abs(DNA_series_small(1,27)-DNA_series_small(2,27))>0.05; %sum_temp�����
            index3= abs(DNA_series_small(1,29)-DNA_series_small(2,29))<1e-5;  %sun_depth���
            index4=all(abs(DNA_series_small(1,[3,4])-DNA_series_small(2,[3,4]))<0.01);  %���㣺��γ��С��0.01��
            if(index1 && index2 && index3 && index4)
                continue
            end
        end
        

        %%%%%���ԭʼ�����ļ����ļ���
        for m=1:length(id)
            fprintf(fid,'%s ',filename_info(id(m),:));
        end
        fprintf(fid,'\n');
        
        %%%%%���ԭʼ�����ļ�
%         fprintf(fid,'%s\n','��Potential Duplicates pairs��:');
%         for m=1:length(id)
%             fprintf(fid,'%s ',filename_info(id(m),:));
%             fprintf(fid,'%3d %.4f %.4f %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %.4f %.4f %3d %.4f %.4f %.4f %.4f %.4f\n',DNA_series(id(m),:));
%         end
%         fprintf(fid,'\n');
%                pause
        number_pairs=number_pairs+1;
        number_profiles=number_profiles+duplicate_number;
        
    end
end

number_pairs
number_profiles


%%
% figure();
% plot(average_DNA2,'o');
% ylabel('Average DNA')


% for i=1:length(average_DNA)
%     if(contains(filename_info(i,:),'CASv1_T_S_19950724_00859_CTD.nc'))
%         i
%         pause
%     end
% end


clear
clc

total_profiles=142537;  %%%1975���������������
duplicate_profile=433;
duplicate_precentage=duplicate_profile/total_profiles*100


filename='./potential_duplicates_output/1975/DuplicateList_potential_duplicate_ALL_1975_unique_1117.xlsx';
[num,txt,~] = xlsread(filename,'A2:BM433');




% 4. ͳ�ƹ��Ҽ�������� ��41��
% 5. ͳ��platform code��������� ��49��
% 6. ͳ�� WOD_cruise_identifier1 ���ͼ�����/ռ�� ��53��
%%%%%%%%%% 7. ͳ��ÿ���ظ����͵����� ����ռ���ظ������ı���
% ��5�У� �ߺ��۲���ͬʱ�䲻ͬλ��_1  ��6�У��۲ⱻ����_2 ��7�У��ߺ��۲�ͬһʱ��ͬһλ�ö����ͬ�۲�_3
% ��8�У��۲ⱻ��������_4  ��9�У��������λ����Ϣ_5 ��10�У�����������Ϣ6 ��11�У�����ʱ����Ϣ7 ��12�У�������Ҵ���8
% ��13�У�������������9  ��14�У�������Ϣ��������ȫ���10



%%% 1.ͳ��accession number ÿ������ж��ٸ��������Ƿ�������������⣨�����ࣩ��access_number
access_num_all=num(:,16);
[unique(access_num_all');histc(access_num_all',unique(access_num_all'));histc(access_num_all',unique(access_num_all'))./length(lat_all)*100]

access_num_copy=access_num_all;

%%%% 2. ������λ�õ�ͼ 
lat_all=num(:,18);
lon_all=num(:,20);
plot_geo_location(lon_all,lat_all,'Location of potential duplicated pairs')
saveas(gcf,'./pics/locaiton_1975.png')

% 3. ͳ��ÿ�����������ظ����ݵĸ���  ��16�У������ǲ�����һЩ�������ر���
instrument_all=txt(:,16);
instrument_type_unique=unique(instrument_all);
for i=1:length(lat_all)
    
    for ii=1:length(instrument_type_unique)
        if(instrument_all{i}==instrument_type_unique{ii})
            instrument_all_order(i)=ii;
        end
    end
end
stat_instrument=[unique(instrument_all_order);histc(instrument_all_order,unique(instrument_all_order));histc(instrument_all_order,unique(instrument_all_order))./length(lat_all)*100]
instrument_type_unique

yan'z

country_all=num(:,39)';
stat_country=[unique(country_all);histc(country_all,unique(country_all));histc(country_all,unique(country_all))./length(lat_all)*100]


platform_code_all=txt(:,49);
platform_unique=unique(platform_code_all)
% for i=1:length(lat_all)
%     
%     for ii=1:length(platform_unique)
%         if(platform_code_all{i}==platform_unique{ii})
%             platform_unique_order(i)=ii;
%         end
%     end
% end
% unique(platform_code_all)
% stat_platform=[unique(platform_unique_order);histc(platform_unique_order,unique(platform_unique_order));histc(platform_unique_order,unique(platform_unique_order))./length(lat_all)*100]
% 


WOD_cruise_identifier_all=txt(:,53);
WOD_cruise_identifier_unique=unique(WOD_cruise_identifier_all)
% for i=1:length(lat_all)
%     
%     for ii=1:length(WOD_cruise_identifier_unique)
%         if(WOD_cruise_identifier_all{i}==WOD_cruise_identifier_unique{ii})
%             WOD_cruise_identifier_order(i)=ii;
%         end
%     end
% end
% unique(WOD_cruise_identifier_all)
% stat_cruise_id=[unique(WOD_cruise_identifier_order);histc(WOD_cruise_identifier_order,unique(WOD_cruise_identifier_order));histc(WOD_cruise_identifier_order,unique(WOD_cruise_identifier_order))./length(lat_all)*100]
% 

% ��5�У� �ߺ��۲���ͬʱ�䲻ͬλ��_1  ��6�У��۲ⱻ����_2 ��7�У��ߺ��۲�ͬһʱ��ͬһλ�ö����ͬ�۲�_3
% ��8�У��۲ⱻ��������_4  ��9�У��������λ����Ϣ_5 ��10�У�����������Ϣ6 ��11�У�����ʱ����Ϣ7 ��12�У�������Ҵ���8
% ��13�У�������������9  ��14�У�������Ϣ��������ȫ���10
duplicate1=num(:,3);
duplicate2=num(:,4);
duplicate3=num(:,5);
duplicate4=num(:,6);
duplicate5=num(:,7);
duplicate6=num(:,8);
duplicate7=num(:,9);
duplicate8=num(:,10);
duplicate9=num(:,11);
duplicate10=num(:,12);

sum(duplicate1~=0),sum(duplicate1~=0)./length(duplicate1)
sum(duplicate2~=0),sum(duplicate2~=0)./length(duplicate2)
sum(duplicate3~=0),sum(duplicate3~=0)./length(duplicate3)
sum(duplicate4~=0),sum(duplicate4~=0)./length(duplicate4)
sum(duplicate5~=0),sum(duplicate5~=0)./length(duplicate5)
sum(duplicate6~=0),sum(duplicate6~=0)./length(duplicate6)
sum(duplicate7~=0),sum(duplicate7~=0)./length(duplicate7)
sum(duplicate8~=0),sum(duplicate8~=0)./length(duplicate8)
sum(duplicate9~=0),sum(duplicate9~=0)./length(duplicate9)
sum(duplicate10~=0),sum(duplicate10~=0)./length(duplicate10)



%% �洢wod_unique_id  �޳�ȫ��id�����ѡȡһ��id���޳� �Է������OHC
clear
clc

filename='./potential_duplicates_output/DuplicateList_potential_duplicate_ALL_1995_unique.xlsx';
[num,txt,~] = xlsread(filename,'A2:BL312');

all_unique_id=[num(:,1),num(:,2)];
all_unique_id=all_unique_id(:);
remove_unique_id=all_unique_id;

save ./wod_unique_id_forOHC/all_duplicated_unique_ID_1995.mat remove_unique_id
clear remove_unique_id

half_unique_id=num(:,1);
remove_unique_id=half_unique_id;
save ./wod_unique_id_forOHC/half_duplicated_unique_ID_1995.mat remove_unique_id



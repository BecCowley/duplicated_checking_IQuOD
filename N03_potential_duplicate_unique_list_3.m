%%%%%%%% duplicate list ȥ��ȡΨһ
%%%%%%%% N02�����г��������˺ܶ��Ǳ�ڵ��ظ����ļ���������Ҫ����Щ�ļ�������ļ��Խ��кϲ������txt�ļ��ϲ���һ�����txt�ļ���
%%%% �Ա�����һ��duplicate ����н��д������Ǻϲ���ʱ����ܻ����ظ����ļ��ԣ����ʱ����Ҫ��ȥ�ء�����֤�ļ�����Ψһ��
clear
clc

files_txt=dir('./potential_duplicates_output/2008/potential_dup*.txt');

filenames_column1={};
filenames_column2={};
filenames_combines={};
m=1;
for i=1:length(files_txt)
    
    file1=[files_txt(i).folder,'/',files_txt(i).name];
    fid=fopen(file1,'r');
    
    while ~feof(fid)
        str=fgetl(fid);
        str=strtrim(str);
        s=regexp(str,'\s+','split');
        for k=2:length(s)
            filename1=s{1};
            filename2=s{k};
            
            filename_combine=[filename1,' ',filename2];
            
            filenames_column1{m}=filename1;
            filenames_column2{m}=filename2;
            filenames_combines{m}=filename_combine;
            m=m+1;
        end
    end
end

%%%%ȥ��
filename_unique_pairs=unique(filenames_combines);


%%%%%% ��������µ�potential_list
fid=fopen('./potential_duplicate_ALL_2008_unique_1119.txt','w+');
%%%%%���ԭʼ�����ļ�
for m=1:length(filename_unique_pairs)
    s=regexp(filename_unique_pairs{m},'\s+','split');
    output_filename1=s{1};
    output_filename2=s{2};
    
    fprintf(fid,'%s   %s',output_filename1,output_filename2);
    fprintf(fid,'\n');
    clear s
end

length(filename_unique_pairs)

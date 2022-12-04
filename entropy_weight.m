
function [weight]=entropy_weight(DNA_series)

%��׼������mapstd  �����ݱ�׼��Ϊ��ֵΪ0������Ϊ1������    mapminmax��һ��
%��һ���ͱ�׼����������̽��
% DNA_mapped=mapstd(DNA_series,0,1);
DNA_series(1,4)=1994;
DNA_mapped=mapminmax(DNA_series',0,1);
DNA_mapped=DNA_mapped';
DNA_mapped(DNA_series==0)=0;


B=DNA_mapped;
%B�����е�ÿһ�о��Ѿ���һ����0-1
B(B==0)=0.00001;
B(B==1)=0.99999;
 
[n,m]=size(B); % n������, m��ָ��
%%�����j��ָ���£���i������ռ��ָ��ı���p(i,j)
p=NaN(n,m);
dd=sum(B,1,'omitnan');
for j=1:m
    p(:,j)=B(:,j)./dd(j);
end
% for i=1:n
%     i
%     for j=1:m
%         p(i,j)=B(i,j)/sum(B(:,j),'omitnan');
%     end
% end
%%�����j��ָ�����ֵe(j)
k=1/log(n);
% e=NaN();
for j=1:m
%     j
    e(j)=-k*sum(p(:,j).*log(p(:,j)),'omitnan');
end
d=ones(1,m)-e; %������Ϣ�������
weight=d./sum(d,'omitnan'); %��Ȩֵw ��Ҫ��Ҫ�Ľ�������
s=100*weight*B'; %���ۺϵ÷�

end
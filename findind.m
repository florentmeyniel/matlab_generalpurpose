function [r]=findind(a,b)
% [r]=findind(a,b) r is the logical indexing, r(i) is 1 if a(i) is a value 
% of b, and 0 if no element of b is equal to r(i)
% r has the same dimension as a
% a and b can be row or column vector
% 
% Florent Meyniel

if min(size(a))>1 || min(size(b))>1
    error('a and b should both be one-dimensional')
end

if size(a,2)<size(a,1)
    a_iscol=1;
    a=a';
else
    a_iscol=0;
end

if size(b,2)<size(b,1)
    b=b';
end

temp=[];
for i=b
    temp(i,:)=[a==i];
end

r=logical(sum(temp));
if a_iscol
    r=r';
end

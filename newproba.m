function newconn = newproba(conn)

N=length(conn);

for kkk=1:N

    Cu=ones(N,N); 
    inds = find(conn~=0);
           
    Cu(inds) = Cu(inds) ./ conn(inds);  %1/x
 
    newconn = log10(Cu); %log(1/x)
     
end
 
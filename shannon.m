
% Author: Pardhu M
% Contact: pardhu.madipalli@outlook.com
% The code corresponds to the method implemented in the following paper:
% Blunt, Shannon D., and Karl Gerlach. "Adaptive pulse compression via MMSE estimation." IEEE Transactions on Aerospace and Electronic Systems 42.2 (2006): 572-584.
% Author is not reponsible for errors in the code



function xcap3 = shannon(tran, recv, L)
%Output: Filter coefficients
%received samples should go from -2(N-1) to 0 to L-1 to L-1+3(N-1)
%trans signal length is N
%window length is L


y=recv;
s=tran;
N=length(tran);
c=zeros(N,N);
%finding initial stage MMSE filter coefficients W
for i=0:N-1
    for j=0:N-1
        if i>=j
            p=[zeros(1,i-j) s(1:N-(i-j))];
        else
            p=[s(j-i+1:N) zeros(1, j-i)];
        end
        
        c(i+1,j+1)=s*ctranspose(p);
    end
end


w=inv(c)*s'; %first MMSE filter coeffs
size(w);


%finding initial x coefficients
xcap=zeros(1, L+4*(N-1));
pho=zeros(1, L+4*(N-1));

for i=1:L+4*(N-1)
    xcap(i)=sum(y(i:N+i-1).*ctranspose(w)); %output for 1st stage
    pho(i)=abs(xcap(i))^2;
end

c=zeros(N,N, L+(2*N)-2);

%finding 2nd stage MMSE filter coefficinets
for l=1:L+(2*N)-2
    
    for n=-N+1:N-1
        if n>=0
            p=[zeros(1,n) s(1:N-n)]';
        else
            p=[s(1-n:N) zeros(1,-n)]';
        end
        
        c(:,:,l)=c(:,:,l)+pho(N+l+n-1)*(p*ctranspose(p));
    end
end

w=zeros(N, L+(2*N)-2);

for l=1:L+(2*N)-2
    w(:,l)=pho(N+1-1)*inv(c(:,:,l))*s'; %second stage filter coefficents
end


y=y(N:length(y)-N+1);

xcap2=zeros(1, L+2*(N-1)); %output for 2nd stage
pho2=zeros(1, L+2*(N-1));

for l=1:L+2*(N-1)
    xcap2(l)=sum(y(l:N+l-1).*ctranspose(w(:,l)));
    pho2(l)=abs(xcap2(l))^2;
end


%for finding 3rd stage w3 and xcap3 
c3=zeros(N,N, L);

for l=1:L
    
    for n=-N+1:N-1
        if n>=0
            p=[zeros(1,n) s(1:N-n)]';
        else
            p=[s(1-n:N) zeros(1,-n)]';
        end
        
        c3(:,:,l)=c3(:,:,l)+pho2(N+l+n-1)*(p*ctranspose(p));
    end
end

w3=zeros(N, L);

for l=1:L
    w3(:,l)=pho2(N+1-1)*inv(c3(:,:,l))*s';
end

y=y(N:length(y)-N+1);

xcap3=zeros(1, L);

for l=1:L
    xcap3(l)=sum(y(l:N+l-1).*ctranspose(w3(:,l))); %final output
end
        

function [G,Gidx]=freqmat(a,b,c,d,w);
% function [G,Gidx]=freqmat(a,b,c,d,w);
% CALCULATE MODEL FREQ RESPONSES
nw=length(w);
[ncp,nia]=size(d);
j=sqrt(-1);
G=zeros(ncp,nia*nw);
Gidx=(1:nia:nw*nia)-1;


for i=1:nia,
%  [re,im]=nyquist(a,b,c,d,i,w);
%  G(:,Gidx+i)=re'+j*im';
 G(:,Gidx+i)=freqresp(a,b,c,d,i,j*w).';
end

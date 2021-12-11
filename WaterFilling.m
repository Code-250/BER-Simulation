function [gamat]=WaterFilling(Lamda,SNR,nT)
r=length(Lamda);
gama = zeros(1,r);
 index=[1:r];
indext=index;
p=1;

while p<r
    ir=[1:r-p+1].'; 
    temp= sum(1./Lamda(indext(ir)));
    mu = nT/(r-p+1.)*(1+1/SNR*temp);
    gama(indext(ir))=mu-nT./(SNR*Lamda(indext(ir)));
    if min(gama(indext))<0
      i=find(gama==min(gama));  ii=find(indext==i);
      indext1=[indext([1:ii-1]) indext([ii+1:end])];
           indext=indext1;
      p=p+1;
      clear gama;
     else
  p=r;
    end
  end
gamat=zeros(1,length(Lamda));
gamat(indext)=gama(indext);

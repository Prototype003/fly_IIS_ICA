function [ SugiCorr , SugiY , SugiX , origY , origX ] = Sugi( X , Y , tau , E  )

L=length(X);
T=1+(E-1)*tau;
Xm=zeros((L-T+1),E);
Ym=zeros((L-T+1),E);
SugiN=E+1;
N = L-T+1;

%% RECONTRUCTIONS OF ORIGINAL SYSTEMS

for t=1:(L-T+1)
    Xm(t,:)=X((T+t-1):-tau:(T+t-1-(E-1)*tau));
    Ym(t,:)=Y((T+t-1):-tau:(T+t-1-(E-1)*tau));
end
%%

SugiX=zeros(N,1);
SugiY=zeros(N,1);

origY=Y(T:end);
origX=X(T:end);

parfor j=1:N
%% neighborhood search 

[n1,d1]=knnsearch(Xm,Xm(j,:),'k',E+2);
[n2,d2]=knnsearch(Ym,Ym(j,:),'k',E+2);
susY=origY(n1(2:end)); 
susX=origX(n2(2:end));

%% CMM

SugsusY=susY(1:SugiN);
SugsusX=susX(1:SugiN);
Sugid1=d1(:,2:SugiN+1);
Sugid2=d2(:,2:SugiN+1);
u1=exp(-Sugid1./(Sugid1(:,1)*ones(1,SugiN)));
u2=exp(-Sugid2./(Sugid2(:,1)*ones(1,SugiN)));
w1=u1./(sum(u1,2)*ones(1,SugiN));
w2=u2./(sum(u2,2)*ones(1,SugiN));
SugiY(j)= w1*SugsusY;
SugiX(j)= w2*SugsusX;

end

SugiCorr1=corrcoef(origY,SugiY);
SugiCorr(2,1)=SugiCorr1(1,2);

SugiCorr2=corrcoef(origX,SugiX);
SugiCorr(1,1)=SugiCorr2(1,2);

end
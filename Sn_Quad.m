function [out]= Sn_Quad(fun,n)
%% 
Octants=[1 1 1; -1 1 1; -1 -1 1; 1 -1 1;1 1 -1; -1 1 -1; -1 -1 -1; 1 -1 -1];
types=[4 6 8 12 16];
M1=[0.3500212 0.2666355 0.2182179 0.1672126 0.1389568];
Pos=zeros(1,n/2);
Pos(1)=M1(types==n);
for i=2:n/2
    Pos(i)=sqrt(Pos(1)^2+2*(i-1)*(1-3*Pos(1)^2)/(n-2));
end
%weights
W=zeros(5,36);
W(1,1:3)=0.3333333;

W(2,[1,3,6])=0.1761263;
W(2,[2,4,5])=0.1572071;

W(3,[1 4 10])=0.1209877;
W(3,[2 3 5 7 8 9])=0.0907407;
W(3,6)=0.0925926;

W(4,[1 6 21])=0.0707626;
W(4,[2 5 7 11 19 20])=0.0558811;
W(4,[3 4 12 15 16 18])=0.0373377;
W(4,[8 10 17])=0.0502819;
W(4,[9 13 14])=0.0258513;

W(5,[1 8 36])=0.0489872;
W(5,[2 7 9 15 34 35])=0.0413296;
W(5,[3 6 16 21 31 33])=0.0212326;
W(5,[4 5 22 26 27 30])=0.0256207;
W(5,[10 14 32])=0.0360486;
W(5,[11 13 17 20 28 29])=0.0144589;
W(5,[12 23 25])=0.0344958;
W(5,[18 19 24])=0.0085179;


W=W*pi/2; %normalization
Wspecific=W(types==n,:);
%%%%%%%%

%Points
Triplets=unique(nchoosek([Pos Pos Pos],3),'rows');
j=1;
%%
clear AbsPoints
clear TripletsCorrected;
for i=1:size(Triplets,1)
    if sum(Triplets(i,:).^2)<1.00001 && sum(Triplets(i,:).^2)>0.99999
        TripletsCorrected(j,:)=Triplets(i,:);
        j=j+1;
    end
end
AbsPoints=TripletsCorrected;

% m=1;
% for i=1:size(TripletsCorrected,1)
% AbsPoints(m:m+2,:)=perms(TripletsCorrected(i,:));
% m=m+3;
% end
%% Integration
value=0;
for oct=1:size(Octants,1)
    for i=1:size(AbsPoints,1)
        xyz=AbsPoints(i,:).*Octants(oct,:);
        value=value+fun(xyz(1), xyz(2), xyz(3))*Wspecific(i);
    end
end
   out=value;
end

clear all
clc
%% Introduction of variables.
Coordinates= [0 5; 5 10];
L1= Coordinates(1,2)-Coordinates(1,1);
L2= Coordinates(2,2)-Coordinates(2,1);
Sigma_s=[45; 8];
Sigma_g=[4; 50];
Sigma_t= Sigma_s + Sigma_g;
N=100000;
Tally_capture=0;
Tally_FWS=0;

% HP: void boundary condition

%% Histories
for hist=1:N
    hist
    % Generate the  particle
    x= rand*L1 + Coordinates(1,1);
   
    % particles "in the game"
          % Sample direction
    if rand>0.5
        dir=1; %right
    else
        dir=0; %left
    end
    flag_bound=0;
    absorbed=0;
    lost=0;
    while (absorbed==0 && lost==0)
       % where is the particle? 
       if ne(x,Coordinates(1,2))
       region=(x>=Coordinates(1,2))+1;
       else
           region=dir+1;
       end
       % determine closest boundary
       % get number of mfp until boundary
       sb=abs(x-Coordinates(region,(dir+1)));
       nb=sb*Sigma_t(region); 
       if flag_bound==0
       %sample the number of mfp until next collision
       nc=-log(1-rand);
       % get the distance
       sc=nc/Sigma_t(region);
       end

           % If collision
           if nc<nb
               %update location x
               x=x+sign(dir-0.5)*sc;
               %sample reaction type
               rt=rand;
               %if absorbed
               if rt<Sigma_g(region)/Sigma_t(region)
                   absorbed=1;
                   if region==2 
                       Tally_capture=Tally_capture +1; %(weight equal to one, so just increase counter)
                   end
               else
                   flag_bound=0;
                   % sample new direction
                   if rand<0.5 % change direction
                       dir = 1-dir;
                   else %preserve direction and score if in region 1 (weight equal to one, so just increase counter)
                       if region==1
                       Tally_FWS=Tally_FWS+1;
                       end
                   end
               end
           % If not collision
           else 
               x=x+sign(dir-0.5)*sb;
               if (x<=Coordinates(1,1) || x>= Coordinates(2,2))
                   lost=1;
               end  
               nc=nc-nb;
               flag_bound=1;
           end
    end
end

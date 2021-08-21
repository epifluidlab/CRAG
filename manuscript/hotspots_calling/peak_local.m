function p=peak_local(i,density,flag_bin,left_bin,right_bin,window_width)
%%%%%%%%calculate the significance of local 5kb for the current window
%%i:the location of the current window
%%density: IFS for the windows of the  whole chromosome
%%flag_bin, vector: whether the windows were a candidate: 1:ture; 0: false
%%left_bin,right_bin: the first and the last window of the cyrrent chromosme
%%window_width: window size


shift_left=fix((window_width-20)/20);
shift_right=fix(window_width/20);

left1=max(left_bin,i-125); %% the furst window of local 5kb
left2=max(left_bin,i-shift_left);

right1=min(right_bin,i+shift_right);
right2=min(right_bin,i+126);


if (left1 == i) %%the
    temp1=density(right1:right2,1);
    flag_temp=flag_bin(right1:right2,1);
else
    if (right1 == i) %%the
        temp1=density(left1:left2,1);
        flag_temp=flag_bin(left1:left2,1);
    else
        temp1=[density(left1:left2,1);density(right1:right2,1)];
        flag_temp=[flag_bin(left1:left2,1);flag_bin(right1:right2,1)];
    end
end

temp1=temp1(flag_temp(:,1)==1,1); %%the IFS score of the local 5kb windows
if length(temp1)>=30
    if (mean(ceil(temp1))>=std(ceil(temp1))^2)
        lamda1=poissfit(ceil(temp1));
        p=poisscdf(density(i,1),lamda1); %% poisson test
    else
        lamda1=fitdist(ceil(temp1),'NegativeBinomial');
        p=cdf(lamda1,density(i,1)); %% poisson test
    end
    
else
    p=1;
end

end

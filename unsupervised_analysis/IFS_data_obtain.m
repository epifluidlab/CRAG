function [feature_data,peak_origin]=IFS_data_obtain(path_name,peak,peak_type)
%%%%%%Obtain the z-score transformed IFS in four hotspot sets for the sample 'path_name'

file_name=strcat(path_name,'/data_n/');

feature=zeros(length(peak),1);
peak_origin=zeros(length(peak),3);
cou=0;

if (peak_type==1)
    
    para=zeros(22,2);
    for i=1:22
        t_file=strcat('chr',num2str(i));
        t_file=strcat(file_name,t_file);
        
        load(t_file);
        null=zeros(100000,1);
        region_len=length(loc);
        for j=1:100000
            while (1)
                %%%%%%%%%Obtain 100,000 random regions and get 1000,000 IFS as
                %%%%%%%%%barkground
                a=randi(region_len-200,1);
                b=a+199;
                if (sum(dark_flag(a:b,1))==200)
                    break;
                end
            end
            null(j,1)=sum(loc(a:b,1));
        end
        [mu,sigma] = normfit(null); %%Get the mu and sigma for z-socre tranform
        para(i,1)=mu;
        para(i,2)=sigma;
    end
    
    %%%%%%%%Get the IFS for the peak and z-score transformed
    for i=1:22
        t_file=strcat('chr',num2str(i));
        t_file=strcat(file_name,t_file);
        load(t_file);
        t_peak=peak(peak(:,1)==i,2:3);
        m=length(t_peak);
        for j=1:m
            base1=fix((t_peak(j,1)+t_peak(j,2))/2);
            left1=base1-100;
            right1=base1+99;
            val1=sum(loc(left1:right1,1));
            cou=cou+1;
            feature(cou,1)=(val1-para(i,1))/para(i,2);
            peak_origin(cou,1)=i;
            peak_origin(cou,2:3)=t_peak(j,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear t_peak;
    end
else
    for i=1:22
        chrm=strcat('chr',num2str(i));
        t_file=strcat(file_name,chrm);
        
        G_name='./Basic_info/GC/';
        G_name=strcat(G_name,chrm);
        load (G_name);
        loc_g=loc;
        
        %%%%%%%load mappability
        path(path,'./Basic_info/mappability/');
        m_name='./Basic_info/mappability/';
        m_name=strcat(m_name,chrm);
        m_name=strcat(m_name,'_m.mat');  %%get the region with no mappbility
        load (m_name);
        load(t_file);
        
        region_len=length(loc);
        step=200;
        loop_end=region_len-200;
        region=zeros(ceil(region_len/step),2);
        flag_bin=zeros(ceil(region_len/step),1);
        cou_candidate=0;
        for j=1:step:loop_end
            a=j;
            b=j+199;
            cou_candidate=cou_candidate+1;
            region(cou_candidate,1)=sum(loc(a:b,1));
            region(cou_candidate,2)=ceil(mean(loc_g(a:b,1)));
            if (region(cou_candidate,2)==0)
                region(cou_candidate,2)=1;   %%%For unexpection data, set it as 1
            end
            te=mean(map_b(a:b,1));
            if (sum(dark_flag(a:b,1)) == 200) && (te >= 0.90) %%filter the region with dark regoion or no mappbility
                flag_bin(cou_candidate,1)=1;
            end
        end
        
        t_val=region(flag_bin==1,:);
        data = smooth(t_val(:,2),t_val(:,1),0.75,'loess');
        
        res=t_val(:,1)-data;
        mean_ave=mean(t_val(:,1));
        val=res+mean_ave;   %%%The loess corrected IFS (coverage)
        [mu,sigma] = normfit(val);
        
        %%%%%%%%%%%%%%%%%GC information and fitted value for each GC context
        pre_GC=zeros(100,1)+mean_ave; %%GC from 1 - 100.Set the default value as mean_val
        %%If the GC content is not exiested in the
        %%sliding window, then the correct value for
        %%the GC contect is set as the original value of
        %%that window (don't corrected that window).
        
        [fit,index]=unique(data); %%Find the fit value for each GC content
        temp=[t_val(index,2) fit]; %% The fitted value for each GC content
        num=length(temp(:,1));
        for j=1:num
            pre_GC(temp(j,1),1)=temp(j,2);  %%The fitted value for each GC content
        end
        
        t_peak=peak(peak(:,1)==i,2:3);
        m=length(t_peak);
        for j=1:m
            base1=fix((t_peak(j,1)+t_peak(j,2))/2);
            left1=base1-100;
            right1=base1+99;
            val1=sum(loc(left1:right1,1));
            GC_val=ceil(mean(loc_g(left1:right1,1)));
            if (GC_val==0)
                GC_val=1;
            end
            est_val1=val1-pre_GC(GC_val,1); %%residual for this window
            est_val1=est_val1+mean_ave; %%GC xorrected IFS (coverage) for the current window
            cou=cou+1;
            feature(cou,1)=(est_val1-mu)/sigma;
            peak_origin(cou,1)=i;
            peak_origin(cou,2:3)=t_peak(j,:);
        end
    end
end

feature_data=feature;


end

function Hotspot_call(path_name,peak_type,loop_id,threshold_p,local_p,threshold_q)
%%%%%%path_name:the file name which you want to call hotspots: i.e.BH01
%%%%%peak_type:call the hotspots based on original fragmentation (1) or GC-bias
%%%%%corrected fragmentation (2).
%%%%%loop_id: chromsome id
%%%%threshold_p: global p-value cut-off
%%%%local_p:local p-value cut-off
%%%%threshold_q: fdr cut-off
if (nargin==5)
    threshold_q=0.01;   %%FDR cut-off
end
if (nargin==4)
    local_p=0.00001; %%P-value cut-off
    threshold_q=0.01;
end
if (nargin==3)
    threshold_p=0.00001;
    local_p=0.00001;
    threshold_q=0.01;
end

path_name=strcat(path_name,'/');

produce_file=strcat(path_name,'data_n/'); %%input file of the IFS

out_path=strcat(path_name,'outpeak_n/'); %%output peak files (before merging process)
if ~exist(out_path,'dir')
system(['mkdir ' out_path]);
end

path(path,produce_file);

listing = dir(produce_file);


cou=0;
for i=1:length(listing)   %%obtain all the names of the usefule files
    lo=strfind(listing(i).name,'ch');
    if ~isempty(lo)
        cou=cou+1;
        file_name{cou,1}=listing(i).name;
    end
end
%%%%%window_width:the sliding window size
window_width=200;

if (peak_type == 1) %%%%Call peaks using IFS
    ii=loop_id;
    %%%%%%%%load IFS data of the current chromosome
    chr=strcat('chr',num2str(ii));
    t_chr=strcat(chr,'.mat');
    
    t_name=strcat(produce_file,chr);
    
    lo=strcmpi(file_name,t_chr);
    if (sum(lo)==0)
        disp('The chromosome does not exist!');
        return;
    end
    load (t_name);
    
    %%%%%%%load mappability score of the current window
    m_name='Basic_info/mappability/';
    
    m_name=strcat(m_name,chr);
    m_name=strcat(m_name,'_m.mat');  %%get the region with no mappbility
    load (m_name);
    region_len=length(loc); %%chromosome length
    
    step=20; %% the step of sliding window
    loop_end=region_len-window_width;
    region=zeros(ceil(region_len/step),4);
    flag_bin=zeros(ceil(region_len/step),1);
    cou_candidate=0;
    for i=1:step:loop_end
        %% scan the whole chromose to get the bais information
        %%for all the sliding windows
        a=i;
        b=i+window_width-1;
        cou_candidate=cou_candidate+1;
        region(cou_candidate,1)=a; %%start location
        region(cou_candidate,2)=b; %%end location
        region(cou_candidate,3)=sum(loc(a:b,1));%%the total IFS of the current window
        if (sum(dark_flag(a:b,1)) == window_width && mean(map_b(a:b,1)) >= 0.9 && region(cou_candidate,3) > 0)
            %%filter the region with dark regoion or low mappbility (avergae mappbiliaty less than 0.9)
            %%Only the windows with at least one fragment were considered.
            flag_bin(cou_candidate,1)=1;
        end
    end
    
    region=region(1:cou_candidate,:);
    flag_bin=flag_bin(1:cou_candidate,1);
    temp_density=region(flag_bin==1,3); %%IFS scores for all the candidate windows
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pd=fitdist(temp_density,'Poisson');
    %fit ladma for (global) possion distribution ;
    
    density=region(:,3);
    
    %%%%%%%%%%%%%detect peak region%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_local_cut=local_p;  %threshold of p-value in the local region
    n=length(region);
    final_region=zeros(sum(flag_bin),6);
    cou=0;
    for i=1:n
        if flag_bin(i,1)==1 %%For the candidate windows
            probability = cdf(pd,region(i,3));%%p-value for global possion test
            if probability <= threshold_p %%Only the windows pass the cut-off were considered further
                cou=cou+1;
                final_region(cou,1:3)=region(i,1:3);
                final_region(cou,4)=probability; %%p-value for global
                final_region(cou,5)=peak_local(i,density,flag_bin,1,n,p_local_cut,window_width);
                %%Using function 'peak_local' to calculate the local
                %%p-value
            end
        end
    end
    
    n=length(final_region);
    if cou<n
        final_region((cou+1):n,:)=[];
    end
    %%%%%%%%%%%%%%%%%%%FDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p=ones(sum(flag_bin),1);
    cou=length(final_region(:,1));
    p(1:cou,1)=final_region(:,4);
    fdr=mafdr(p,'BHFDR',true);
    %%Using Benjamini and Hochberg correction
    %%to obtain the FDR for global p-value
    final_region(:,6)=fdr(1:cou,1);
    peak=final_region(final_region(:,6) <= threshold_q & final_region(:,5) <= p_local_cut,:);
    %%The sliding windows with fdr, local p-value and global p-value
    %%pass the cut-off were set as hotspots.
    
    peak_file=strcat(out_path,chr);
    peak_file=strcat(peak_file,'_peak.mat');
    save((peak_file),'peak','region_len','-v7.3');
    %%%Saving the peak file for the current window (./outpeak/chrid_peak.mat)
    clear region;
else
    if (peak_type == 2)   %%%%Call peaks using fragment length (GC bias corrected IFS)
        ii=loop_id; %%chromosome id
        %%%%%%%%%%%chromosome id%%%%%
        chr=strcat('chr',num2str(ii));
        %%%%%%%%%%%load GC content%%%%%%%%%%%
        G_name='Basic_info/GC/';
        G_name=strcat(G_name,chr);
        load (G_name);
        loc_g=loc;
        
        %%%%%%%%%%load IFS
        t_name=strcat(produce_file,chr);
        t_chr=strcat(chr,'.mat');
        lo=strcmpi(file_name,t_chr);
        if (sum(lo)==0)
            disp('The chromosome does not exist!');
            return;
        end
        load (t_name);
        
        %%%%%%%load mappability
        m_name='Basic_info/mappability/';
        m_name=strcat(m_name,chr);
        m_name=strcat(m_name,'_m.mat');
        load (m_name);
        
        region_len=length(loc);
        step=20;
        loop_end=region_len-window_width;
        region=zeros(ceil(region_len/step),6);
        flag_bin=zeros(ceil(region_len/step),1);
        cou_candidate=0;
        for i=1:step:loop_end
            a=i;
            b=i+window_width-1;
            cou_candidate=cou_candidate+1;
            region(cou_candidate,1)=a;  %%start location
            region(cou_candidate,2)=b;  %%end location
            region(cou_candidate,3)=sum(loc(a:b,1)); %%IFS score in the current window
            region(cou_candidate,4)=ceil(mean(loc_g(a:b,1))); %%average GC content of the current window
            if (region(cou_candidate,4)==0)
                region(cou_candidate,4)=1;   %%%For the window without GC, set it as 1
            end
            
            region(cou_candidate,5)=mean(map_b(a:b,1));%%average mappability score
            if (sum(dark_flag(a:b,1)) == window_width) && region(cou_candidate,3)>0 && region(cou_candidate,5)>=0.90 %%filter the region with dark regoion or no mappbility
                %%filter the region with dark regoion or low mappbility (avergae mappbiliaty less than 0.9)
                %%Only the windows with at least one fragment were considered.
                flag_bin(cou_candidate,1)=1;
            end
        end
        
        region=region(1:cou_candidate,:);
        flag_bin=flag_bin(1:cou_candidate,1);
        
        index=(1:length(flag_bin))';
        index=index(flag_bin==1,1); %% The IFS for the candidate widows
        data = smooth(region(index,4),region(index,3),0.75,'loess');
        %%Using loess fit to fit the GC contents and IFS socres for all the
        %%candidate windows
        
        res=region(index,3)-data; %%Residual value for all the candidate windows
        mean_ave=mean(region(index,3)); %%Averge IFS
        
        val=res+mean_ave; %%For all the candidate windows, using Residual value
        %%and average IFS to add back the values.
        %%val is the corrected IFS for all the candidate windows
        
        val(val<0,1)=0;
        
        region(index,3)=val; %%GC bias corrected IFS
        temp_density=region(flag_bin==1,3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pd=fitdist(temp_density,'Poisson');
        density=region(:,3);
        
        %%%%%%%%%%%%%detect peak region%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_local_cut=local_p;  %threshold of p-value in the local region
        n=length(region);
        final_region=zeros(sum(flag_bin),6);
        cou=0;
        for i=1:n
            if flag_bin(i,1)==1
                probability = cdf(pd,region(i,3));
                if probability <= threshold_p
                    cou=cou+1;
                    final_region(cou,1:2)=region(i,1:2);
                    final_region(cou,3)=region(i,3);
                    final_region(cou,4)=probability;
                    final_region(cou,5)=peak_local(i,density,flag_bin,1,n,p_local_cut,window_width);
                    %%Using function 'peak_local' to calculate the local
                    %%p-value
                end
            end
        end
        
        n=length(final_region);
        if cou<n
            final_region((cou+1):n,:)=[];
        end
        %%%%%%%%%%%%%%%%%%%FDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p=ones(sum(flag_bin),1);
        cou=length(final_region(:,1));
        p(1:cou,1)=final_region(:,4);
        fdr=mafdr(p,'BHFDR',true);
        final_region(:,6)=fdr(1:cou,1);
        peak=final_region(final_region(:,6) <= threshold_q & final_region(:,5) <= p_local_cut ,:);
        %%The sliding windows with fdr, local p-value and global p-value
        %%pass the cut-off were set as hotspots.
        peak_file=strcat(out_path,chr);
        peak_file=strcat(peak_file,'_peak.mat');
        save((peak_file),'peak','region_len','-v7.3');
        clear region;
    end
end
end


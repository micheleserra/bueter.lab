clear all
close all

%% give initial thresholds - here you can modify the parameters
Fs=1000; %Hz; Sampling frequency of the signal
FOI=0.5; %Hz; cut off frequency for filtering (in Hz)
artefact_speed_threshold=20; %ml/sec; Speed signal is considered as an artefact above the threshold (absolute value). The detected interval is deleted and interpolated linearly.
min_amplitude_of_suck=0.5; %ml/sec; speed is set to zero below this threshold
min_peak_distance=1; %sec; the minimal distance between peaks sucks (prevent to detect noine as suck) 
min_suck_amount=1; %ml; below this threshold the sucks considered as noise and not included in the analysis.
min_burst_interval=3.5; %sec; set IBI
refill_pause=min_burst_interval; % the refill intervall is replace by the IBI
name_of_the_refill_table='refill_table.xlsx';
refill_table=0; %1=yes, 0=no
want_figures=1; %1=yes, 0=no



%% read the files from folder

dir_name = uigetdir;
Path = dir_name;
dir_file = dir(dir_name);
num_file = size(dir_file);

%% read refill data
if refill_table
    [interval_refill,filenames_refill] = xlsread(name_of_the_refill_table);
    filenames_refill=filenames_refill(2:end,1);
end

Suck_all_time=[];
addpath(dir_name)
cntr_subj=0;
 cntr_no_subj=0;
for i=1:num_file(1)
    File = dir_file(i).name;
    [~,name,ext] = fileparts([dir_name,'/',File]);
    
    xls_e = strcmp(ext, '.tdms');
    
    if (xls_e > 0)
        
        
        
        
        cntr_subj=cntr_subj+1;
        
        FileName{cntr_subj}=File;
        TDMS = TDMS_readTDMSFile(fullfile([dir_name,'/',File]));
        
        
            nofSamples=size(TDMS.data{3},2);
            Data1=zeros(size(TDMS.data,2)-3,nofSamples);
            for idx=4:size(TDMS.data,2)
            Data1(idx-3,:)=TDMS.data{idx};
            end

            t=(1:nofSamples)./Fs;

            fidInfo=fopen([dir_name,'/', name,'.info']);
            frewind(fidInfo)
            C=textscan(fidInfo,'%s %f' , 'CommentStyle',{'%'},'Delimiter',{':', '\n'});
            fclose(fidInfo);

            weight=zeros(4,nofSamples);
            height=zeros(4,nofSamples);
            temp=zeros(4,nofSamples);

            for IDX_unit=1:4
                %apply the parameters form the output linearization
                weight(IDX_unit,:) = (Data1(1+(IDX_unit-1)*3,:)).*C{1,2}...
                (26+(IDX_unit-1)*2)+C{1,2}(27+(IDX_unit-1)*2)-C{1,2}...
                (46+(IDX_unit-1)*2);
                height(IDX_unit,:) = (Data1(2+(IDX_unit-1)*3,:)).*C{1,2}...
                (18+(IDX_unit-1)*2)+C{1,2}(19+(IDX_unit-1)*2)-C{1,2}...
                (47+(IDX_unit-1)*2);
                temp(IDX_unit,:) = (Data1(3+(IDX_unit-1)*3,:)).*C{1,2}...
                (34+(IDX_unit-1)*2)+C{1,2}(35+(IDX_unit-1)*2);
            end

            rate=1;
            weight_ds=zeros(4,ceil(nofSamples/rate));
            height_ds=zeros(4,ceil(nofSamples/rate));
            temp_ds=zeros(4,ceil(nofSamples/rate));
            for IDX_unit=1:4
            weight_ds(IDX_unit,:) = weight(IDX_unit,:);
            height_ds(IDX_unit,:) = height(IDX_unit,:);
            temp_ds(IDX_unit,:) = temp(IDX_unit,:);
            end

            weight_comp=zeros(4,ceil(nofSamples/rate));
            A_tube=2*pi*(C{1,2}(10)/2)^2-2*pi*(C{1,2}(11)/2)^2;
            A_hose=2*pi*(C{1,2}(12)/2)^2-2*pi*(C{1,2}(13)/2)^2;
            A=A_tube+A_hose;
            for IDX_unit=1:4
            Density(IDX_unit)=C{1,2}(42+(IDX_unit-1));
            weight_comp(IDX_unit,:)=A.*height_ds(IDX_unit,:).*Density(IDX_unit);
            end
            % Appliy correction for weight difference by buonancy
            Volume=weight_ds-weight_comp;
            Volume=ones(size(Volume)).*max(Volume)-Volume;
            t=rate/Fs:rate/Fs:nofSamples/(Fs/rate)/rate;
            Volume=Volume(1,:);
            V=Volume(1,:);

%% refill
if refill_table
            L=zeros(1,length(filenames_refill));
            for b=1:length(filenames_refill)
                
                tf = regexpi(File,filenames_refill{b});
                if ~isempty(tf)
                    L(b)=1;
                end
                
            end
            t_r=[];
            if sum(L)>0
                hol=find(L);
                for b=1:sum(L)
                    t_r(b)=interval_refill(hol(b),1);
                    
                    % deal with the refill
                    Volume((interval_refill(hol(b),2)*Fs):end) = Volume((interval_refill(hol(b),2)*Fs):end) + ...
                        (Volume(interval_refill(hol(b),1)*Fs)-Volume(interval_refill(hol(b),2)*Fs)); %
                    
                    % in small interval around the refill make the signal
                    % equal to prevent spurious sucks
                     h1=(interval_refill(hol(b),1)*Fs)-500:(interval_refill(hol(b),1)*Fs);
                     h2=(interval_refill(hol(b),2)*Fs):(interval_refill(hol(b),2)*Fs)+500;
                    Volume([h1,h2])=mean(Volume([h1,h2]));
                        Volume((interval_refill(hol(b),1)*Fs):(interval_refill(hol(b),2)*Fs)) = mean(Volume([h1,h2]));
                end
                
            end
end
%% end refill

%% filter
       
            data_filt=bst_bandpass_fft(Volume(1,:)*1000, 1000, 0,FOI, 1, 1);
            speed=gradient(data_filt)' .*(Fs/rate);
            amount_consumed=data_filt;
%% interpolation
            [ind,val]=find(abs(speed)>artefact_speed_threshold);
            vege=ind;
            interp_speed=speed;
            interp_amount_consumed=amount_consumed;
            while ~isempty(vege)
                %% find close artefacts together
                ind_interval=vege(1)-2000:vege(1)+2000; %mekkora intervallum legyen, ahol a közöen elõforduló mûtermékeket interpolálja?
                indi=0;
                ia_e=1;
                while indi==0
                [C,ia,ib] = intersect(vege,ind_interval); %ia=vege; ib=ind_interval
                if (length(ia)-length(ia_e))>0
                    ia_e=ia;
                    ind_interval=(ind_interval(1)) : (ind_interval(end)+2000);
                    
                else
                    indi=1;
                    
                end
                end
                vege(ia)=[];

                time_interp=t(ind_interval);
                time_ok=setdiff(t,time_interp);

                %speed
                good_speed=speed;
                good_speed(ind_interval)=[];
                signal_interp = interp1(time_ok,good_speed,time_interp,'linear');
                interp_speed(ind_interval)=signal_interp;

                %amount_consumed
                good_amount_consumed=amount_consumed;
                good_amount_consumed(ind_interval)=[];
                signal_interp = interp1(time_ok,good_amount_consumed,time_interp,'linear');
                interp_amount_consumed(ind_interval)=signal_interp;

            end
            interp_speed(interp_speed<min_amplitude_of_suck)=0;
            interp_speed(end)=0;
            
            
            %% find zeros, make zeros
%             error('here')
if refill_table          
           if sum(L)>0
                hol=find(L);
                for b=1:sum(L) %here the accidently arising sucks (due to noise) were deleted
                     h1=(interval_refill(hol(b),1)*Fs)-500:(interval_refill(hol(b),1)*Fs);
                     h2=(interval_refill(hol(b),2)*Fs):(interval_refill(hol(b),2)*Fs)+500;
                     interp_speed(h1(1):h2(end))=0;
                     ind=0;
                     cntr=1;
                     while ind==0
                         if interp_speed(h1(1)-cntr)==0
                             ind=1;
                         else
                             interp_speed(h1(1)-cntr)=0;
                             cntr=cntr+1;
                         end
                     end
                     
                     ind=0;
                     cntr=1;
                     while ind==0 
                         if interp_speed(h2(end)+cntr)==0 | h2(end)+cntr==length(interp_speed)-1
                             ind=1;
                         else
                             interp_speed(h2(end)+cntr)=0;
                             cntr=cntr+1;
                         end
                     end
                     
                end
                %%
               interp_speed=interp_speed';
                hh=0;
                hol2=[];
                for b=1:sum(L) %find closest nonzero speeds around refills
                    
                     h1=(interval_refill(hol(b),1)*Fs-hh)-500:(interval_refill(hol(b),1)*Fs-hh);
                     h2=(interval_refill(hol(b),2)*Fs-hh):(interval_refill(hol(b),2)*Fs-hh)+500;
                     ind=0;
                     cntr=1;
                     while ind==0
                         if interp_speed(h1(1)-cntr)~=0 | h1(1)-cntr==1
                             ind=1;
                             hol2(b,1)=h1(1)-cntr;
                         else
                             cntr=cntr+1;
                         end
                     end
                     
                     ind=0;
                     cntr=1;
                     while ind==0
                         if interp_speed(h2(end)+cntr)~=0 | h2(end)+cntr==length(interp_speed)-1
                             ind=1;
                             hol2(b,2)=h2(end)+cntr;
                         else
                             cntr=cntr+1;
                         end
                     end
                     
                      % delete the pauses around refill sucks and add 3000
                      % datapoints
                     interp_speed=[(interp_speed(1:hol2(b,1))),(zeros(1,refill_pause*Fs)),(interp_speed(hol2(b,2):end))];
                     interp_amount_consumed=[(interp_amount_consumed(1:hol2(b,1))),(ones(1,refill_pause*Fs)*interp_amount_consumed(hol2(b,1))),(interp_amount_consumed(hol2(b,2):end))];
                     hh=hh+(hol2(b,2)-hol2(b,1))-refill_pause*Fs;
                end
                
           end
end
            t=(1:length(interp_amount_consumed))/Fs;
            
            
%%

        [pks,locs,w,proms] = findpeaks(interp_speed,t,'MinPeakProminence',...
             min_amplitude_of_suck,'MinPeakDistance',min_peak_distance,'Annotate','peaks'); %,'MinPeakWidth',min_peak_width
        
        [C,ai,ib] = intersect(t,locs);
        
        
        
        Basis_ind=zeros(2,length(ai));
        
        for a=1:length(ai)
            
            ind=0;
            cntr=1;
            while ind==0
                
                if interp_speed(ai(a)+cntr)>0
                    cntr=cntr+1;
                else
                    Basis_ind(1,a)=ai(a)+cntr;
                    ind=1;
                end
                
                if ismember(ai(a)+cntr,ai)
                    
                    [M,I]=min(interp_speed(ai(a):ai(a)+cntr));
                    Basis_ind(1,a)=ai(a)+I-1;
                    ind=1;
                    
                end
                
            end
            
            
            ind=0;
            cntr=1;
            while ind==0
                
                if interp_speed(ai(a)-cntr)>0
                    cntr=cntr+1;
                else
                    Basis_ind(2,a)=ai(a)-cntr;
                    ind=1;
                end
                
                if ismember(ai(a)-cntr,ai)

                    [M,I]=min(interp_speed(ai(a)-cntr:ai(a)));
                    Basis_ind(2,a)=ai(a)-length(interp_speed(ai(a)-cntr:ai(a)))+I;
                    ind=1;
                    
                end
                
            end
        end
        
              %% consumed for sucks
        suck_consumed=(interp_amount_consumed(Basis_ind(1,:)-1)-interp_amount_consumed(Basis_ind(2,:)+1))';
        

        
        %%
        T_suck=(t(Basis_ind(1,:)-1)-t(Basis_ind(2,:)+1))';
        small_suck=find(suck_consumed<min_suck_amount);
        suck_consumed(small_suck)=[];
        T_suck(small_suck)=[];
        Vmax=pks;
        Vmax(small_suck)=[];
        TimeVamx=locs;
        locs(small_suck)=[];
        TimeVamx(small_suck)=[];
        Basis_ind(:,small_suck)=[];
        
        start_end_point(1)=min(Basis_ind(2,:));
        start_end_point(2)=max(Basis_ind(1,:));
        
        %% intersuck interval
        suck_interval=zeros(1,size(Basis_ind,2)-1);
        
        for a=1:size(Basis_ind,2)-1
            
            suck_interval(a)=t(Basis_ind(2,a+1))-t(Basis_ind(1,a));
            
        end
        O(cntr_subj,1)=sum(suck_consumed);
        O(cntr_subj,2)=interp_amount_consumed(end);
        O(cntr_subj,3)=sum(suck_interval)./max(t);


timestamp=t;
%%
if want_figures==1   
        
        f=figure;
        findpeaks(interp_speed,timestamp,'MinPeakProminence',...
        min_amplitude_of_suck,'MinPeakDistance',min_peak_distance,'Annotate','peaks'); %'MinPeakWidth',min_peak_width,
        
        hold on
        
        plot(timestamp(Basis_ind(1,:)),interp_speed(Basis_ind(1,:)),'r*')
        
        hold on
        
        plot(timestamp(Basis_ind(2,:)),interp_speed(Basis_ind(2,:)),'g*')
        
        xlabel('time (s)')
        ylabel('speed (ml/s)')
        title (File)
        
        saveas(f,[File,'.png'])
        close(f)
        
end





%%  all table

     
% inter_sucking interval
OUT_all.inter_suck_interval(cntr_subj,:)= suck_interval;
      
    end
    
    end
%% write all
Sum_out_all=cat(2,OUT_all.inter_suck_interval);

Label={'inter_sucking_interval', cntr_subj};

W=cat(1,Label,num2cell(Sum_out_all));
W=cat(2,(['FileNames', FileName])',W);

xlswrite([erase(['Full_Consumption','_',num2str(min_amplitude_of_suck),'_',num2str(min_peak_distance),'_',num2str(min_suck_amount),'_',...
    num2str(min_burst_interval),'_',date,],'.'),'.xlsx'],W) 

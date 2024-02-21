%%%%%%% read metadata and secondary data from WOD18 netCDF files
clear
clc

%%%%%% path with the WOD18 netCDf files
filepath='/Users/cow074/code/CARS/CODA/';  %%% This path doesn't include all the year 1995 netCDF files. It is justed a demo
%%%%%%% You need to download the netCDF files from WOD https://www.ncei.noaa.gov/access/world-ocean-database-select/dbsearch.html

filenames=dir(filepath);
filenames([filenames.isdir])=[];
n_files=length(filenames);

DNA_series = []; files = [];

meta_name={'wod_unique_id','accession_number','dataset_name','lat','lon',...
    '','','','probe_type','recorder','','',...
    'depth_number','maximum_depth','hasTemp','hasSalinity','hasOxygen',...
    'hasChlorophyll','country_id','GMT_time','WMO_ID','dbase_orig',...
    'project_name','Platform','ocean_vehicle',...
    'Institute','WOD_cruise_identifier','sum_temp',...
    'sum_salinity','sum_depth','std_depth','std_temp','std_salinity',...
    'corr_temp_depth','corr_sal_depth'};

var_name = {'WOD_id','Access_no','dataset','latitude','longitude',...
    '','','','Temperature_Instrument','Recorder','','','','',...
    'Temperature','Salinity','Oxygen',...
    'Chlorophyll','country','GMT_time','','dbase_orig',...
    'Project','Platform','Ocean_Vehicle',...
    'Institute','WOD_cruise_identifier','',...
    '','','','','',...
    '',''};
variable_name={'wod_unique_id','accession_number','dataset_id','lat','lon',...
    'yyyy','mm','dd','probe_type','recorder','HH','MM',...
    'depth_number','maximum_depth','hasTemp','hasSalinity','hasOxygen',...
    'hasChlorophyll','country_id','GMT_time','WMO_ID','dbase_orig',...
    'project_name','Platform','ocean_vehicle',...
    'Institute','WOD_cruise_identifier','sum_temp',...
    'sum_salinity','sum_depth','std_depth','std_temp','std_salinity',...
    'corr_temp_depth','corr_sal_depth'};

for nf = n_files
    file=[filepath,filenames(nf).name];
    n_prof = length(ncread(file,'cast'));
    ds=single(NaN(n_prof,35));

    % Reading data from the NetCDF file

    for ivar = 1:length(meta_name)
        if ~isempty(meta_name{ivar})
            switch ivar
                case [9, 10, 19, 22:27] %text format
                    dat = repmat('', n_prof,1);
                otherwise
                    dat = NaN*ones(n_prof,1);
            end
            try
                eval(['dat = ncread(file,''' var_name{ivar} ''');'])
                if ischar(dat)
                    dat = dat';
                end
            catch

            end
            eval([meta_name{ivar} ' = dat;'])
        end
    end

    % Sort out the time into each column
    time_units = ncreadatt(file, 'time','units');
    time_calendar = ncreadatt(file,'time','calendar');
    % get the full time information:
    [time,timezone]=cdfdate2num(time_units,time_calendar,ncread(file,'time'));
    yyyy = year(time);
    mm = month(time);
    dd = day(time);
    HH = hour(time);
    MM = minute(time);

    % depth information
    depth = ncread(file, 'Temperature_z');
    depth(depth > 12000 | depth < -10) = NaN;
    % sums and stdev:
    sum_depth = round(sum(depth,'omitnan'),4)';
    std_depth = round(std(depth,'omitnan'),4)';
    depth_number = sum(~isnan(depth))';
    maximum_depth = max(depth)';

    %temperature information
    try
        temp = ncread(file, 'Temperature');
        temp(temp > 40 | temp < -2.5) = NaN;
        sum_temp = round(sum(temp,'omitnan'),4)';
        std_temp = round(std(temp,'omitnan'),4)';

        % set zero sums to NaN (ie, no temperature data)
        sum_temp(sum_temp == 0.0) = NaN;
        std_temp(std_temp == 0.0) = NaN;
        hasTemp = ~isnan(sum_temp)';
    catch
        temp = NaN*zeros(size(depth));
        [sum_temp,std_temp,cor_temp_depth,hasTemp] = deal(zeros(n_depth,n_prof));
    end

    % Read salinity and apply conditions
    try
        sal = ncread(file, 'Salinity');
        sal(sal > 43 | sal < -1) = NaN;

        sum_sal = round(sum(sal,'omitnan'), 4);
        std_sal = round(std(sal,'omitnan'), 4);
        hasSalinity = ~isnan(sum_sal)';

        % set zero sums to NaN (ie, no temperature data)
        sum_temp(sum_sal == 0.0) = NaN;
        std_temp(std_sal == 0.0) = NaN;

    catch
        sal = NaN*zeros(size(depth));
        [sum_sal,std_sal,cor_sal_depth,hasSalinity] = deal(zeros(1,n_prof));
    end

    try
        oxygen = ncread(file, 'Oxygen');
        sum_oxy = sum(oxygen,'omitnan');
        hasOxygen = ~isnan(sum_oxy)';
    catch
        hasOxygen = zeros(n_prof,1);
    end
    try
        chlorophyll = ncread(file, 'Chlorophyll');
        sum_chl = sum(chlorophyll,'omitnan');
        hasChlorophyll = ~isnan(sum_chl)';
    catch
        hasChlorophyll = zeros(n_prof,1);
    end

    dataset_id = NaN*ones(n_prof,1);
    dataset_name = string(deblank(lower(dataset_name)));
    % Assign a number for each dataset type
    pat = ["bod", "bottle", "ocean station", "osd", "low-resolution", "low resolution"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 1;
    pat = ["towed","uor","undulating"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 10;
    pat = ["ctd","xctd"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 2;
    pat = ["mbt","mechanica","mb"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 3;
    pat = ["xbt","xb","expendable"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 4;
    pat = ["sur"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 5;
    pat = ["apb","autonomous","animal"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 6;
    pat = ["mrb","moored","tao"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 7;
    pat = ["pfl","argo","profiling"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 8;
    pat = ["drb","drifting"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 9;
    pat = ["gld","glider"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 11;
    pat = ["dbt"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 12;
    pat = ["std"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 13;
    pat = ["microbt"];
    ii = contains(dataset_name, pat);
    dataset_id(ii) = 14;

    % Get filename from the path
    [save_path, filename, ext] = fileparts(file);

    %%%%%%% put into the array and cell
    for a = 1:length(var_name)
        dat = eval(variable_name{a});
        if isstring(dat) || ischar(dat)
            ascill_all=char(dat)+0;
            %%% Set all NULL to nan
            ascill_all(ascill_all==' ')=NaN;
            %%% Sum per line
            sum_ascill_all=sum(ascill_all,2,'omitnan');
            ds(:,a)=sum_ascill_all;
        else
            ds(:,a) = dat;
        end
    end
    % add to the master arrays
    DNA_series = [DNA_series; ds];
    files = [files,file];
end

yr = datestr(now,'yyyymmdd');
save([filepath 'Metadata_summary_' yr '.mat'], 'DNA_series', 'variable_name','files')

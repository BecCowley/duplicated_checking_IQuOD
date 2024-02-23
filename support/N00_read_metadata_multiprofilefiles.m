%%%%%%% read metadata and secondary data from CODA netCDF files
% CODA format files have each variable in a separate file and multiple
% profiles in each file
% could be adapted to work on ragged array netcdf files from WOD
clear
clc

%%%%%% path with the CODA format netCDf files
filepath='/oa-decadal-climate/work/observations/CARSv2_ancillary/CODA/CODA_test/2010/'; 

filenames=dir(filepath);
filenames([filenames.isdir])=[];
n_files=length(filenames);

DNA_series = []; files = {}; fileid = [];

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

for nf = 1:n_files
    file=[filepath,filenames(nf).name];
    % Get filename from the path
    [save_path, filename, ext] = fileparts(file);

    disp(filename)
    n_prof = length(ncread(file,'cast'));
    disp(n_prof)
    ds=double(NaN(n_prof,35));

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
    if ~isdouble(time)
        disp('time not double')
        time = double(time);
    end
    yyyy = year(time);
    mm = month(time);
    dd = day(time);
    HH = hour(time);
    MM = minute(time);

    % depth information
    var = strsplit(filename,'_');
    depth = ncread(file, [var{5} '_z']);
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
        hasTemp = zeros(n_prof,1);
        [sum_temp,std_temp,corr_temp_depth] = deal(NaN*zeros(n_prof,1));
    end

    % Read salinity and apply conditions
    try
        sal = ncread(file, 'Salinity');
        sal(sal > 43 | sal < -1) = NaN;

        sum_salinity = round(sum(sal,'omitnan'), 4);
        std_salinity = round(std(sal,'omitnan'), 4);
        hasSalinity = ~isnan(sum_salinity)';

        % set zero sums to NaN (ie, no temperature data)
        sum_salinity(sum_salinity == 0.0) = NaN;
        std_salinity(std_salinity == 0.0) = NaN;
        

    catch
        sal = NaN*zeros(size(depth));
        hasSalinity = zeros(n_prof,1);
        [sum_salinity,std_salinity,corr_sal_depth] = deal(NaN*zeros(n_prof,1));
    end
    % get the correlation coefficients for the t/z and s/z relationships
    for ip = 1:n_prof
        inan = isnan(temp(:,ip) .* depth(:,ip));
        if sum(~inan) > 0
            temp2 = temp(~inan,ip);depth2 = depth(~inan,ip);
            corr_temp_depth(ip) = round(corr(temp2, depth2), 5);
        end
        inan = isnan(sal(:,ip) .* depth(:,ip));
        if sum(~inan) > 0
            sal2 = sal(~inan,ip); depth2 = depth(~inan,ip);
            corr_sal_depth(ip) = round(corr(sal2, depth2), 5);
        end
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
    files = [files;[filename ext]];
    fileid = [fileid;repmat(nf,n_prof,1)];

end
% combine metadata from the same wod_unique_ids since we are reading
% separate files for each variable
[C, ia, ic] = unique(DNA_series(:,1));
new_dna = DNA_series(ia,:);
for a = 1:length(C)
    ii = find(DNA_series(:,1) == C(a));
    if length(ii) >1
        % only for multiple records
        % cols 1:14 should be identical
        dd = diff(DNA_series(ii,[1:14,19:27,30:31]));
        if any(dd) > 0
            disp('Different metadata!!')
            disp(DNA_series(ii,1))
            disp(files(fileid(ii)))
            continue
        end
        % assign ~NaN values to first version of the record
        cols = [15:18,28:29,32:35];

        % now assign to new_dna
        for b = 2:length(ii)
            dat = (DNA_series(ii(b),cols));
            ij = ~isnan(dat) & dat > 0;
            new_dna(a,cols(ij)) = DNA_series(ii(b),cols(ij));
        end
    end
end
fileid = fileid(ia);
DNA_series = new_dna;

yr = datestr(now,'yyyymmdd');
save([filepath 'Metadata_summary_' yr '.mat'], 'DNA_series', 'variable_name','files','fileid')

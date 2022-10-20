function [dtstruct, data, cell_id, cell_num, filename] = h5_file_readout(file_path)
%This function loops through a folder and read out the header and the
%sweeps from the h5 data file and organize them by cell IDs (indicated in
%the original file name.

% Input:
%   file_path- path of the data folder

% Output:
%   dtstruct- includes both the header file and all the sweeps
%   data- all sweeps of a h5 file saved by its cell ID
%   cell_id- cell numbers and trace numbers within each cell
%   cell_num- cell IDs in numerical format
%   filename- original names of all files within the specified
%   folder


cd(file_path)
all_files = dir;
file_num = numel(all_files)-2;

filename = cell(1,file_num);
cell_id = cell(1,file_num);
cell_num = NaN(1,file_num);
data = cell(1,file_num);
dtstruct = cell(1,file_num);

if file_num == 0
    dtstruct = NaN;
    data = NaN;
    cell_id = NaN;
    cell_num = NaN;
    filename = NaN;
else

    for f = 1:file_num
        filename{f} = all_files(f+2).name;
        num_filename = numel(filename{f});

        if isnan(str2double(filename{f}(6)))
            cellID = str2double(filename{f}(5));
        else
            cellID = str2double(filename{f}(5:6));
        end

        cell_num(f) = cellID;

        if isempty(strfind(filename{f}, '-'))
            range(1) = str2double(filename{f}(num_filename-6 : num_filename-3));
            range(2) = range(1);
        else
            range(1) = str2double(filename{f}(num_filename - 11:num_filename - 8)); %first sweep #
            range(2) = str2double(filename{f}(num_filename - 6:num_filename - 3)); %last sweep #
        end

        %format file names
        for tii = range(1):range(2)

            trace_id = tii;
            cell_id{1,cellID}(trace_id,1) = trace_id;

            if tii < 10
                sname = strcat('000',num2str(tii));
            elseif tii >= 10 && tii < 100
                sname = strcat('00',num2str(tii));
            elseif tii >= 100 && tii < 1000
                sname = strcat('0',num2str(tii));
            else
                sname = num2str(tii);
            end


            dtstruct_temp = ws.loadDataFile(filename{f});
            dtstruct{1,cellID} = dtstruct_temp;
            prefix = strcat('dtstruct_temp.sweep_',sname,'.analogScans');
            data_temp = eval(prefix); 
            data{1,cellID}(1:numel(data_temp(:,1)),trace_id) = ...,
                data_temp(1:numel(data_temp(:,1)),1);
        end

    end

    %%% Replace zeroes in each cell with NaN
    for gi = 1:max(cell_num)
        cell_id{1,gi} = nonzeros(cell_id{1,gi});
        if isempty(cell_id{1,gi}) == 1
            continue
        else
            for ri = 1:cell_id{1,gi}(end,1)
                if ismember(ri,cell_id{1,gi}) == 0
                    data{1,gi}(:,ri) = NaN;
                end
            end
        end
    end


end
end
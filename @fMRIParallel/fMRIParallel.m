classdef fMRIParallel <handle
    properties
        subject_array;
        global_mask;
    end
    
    methods (Access = public)
        function obj  = fMRIParallel(sa)
            obj.subject_array = sa;
        end
        
        function global_mask = createMask(obj)
            % creates a mask for all runs in the subject_array.
            runs = obj.subject_array.get_runs;
            n = size(runs,1);
            scan_header = spm_vol(runs{1}.scans{1});
            global_mask = ones(1,prod(scan_header.dim));
            parfor i = 1:n
                mask = obj.calcMask(runs{i});
                global_mask = global_mask & mask;
                disp('done mask');
            end
            obj.global_mask = global_mask;
        end
        
        function createZ(obj, mask)
            if nargin<2
                mask = obj.global_mask;
            end
            runs = obj.subject_array.get_runs;
            n = size(runs,1);
            for i = 1:n
                obj.calculateZ(runs{i}, mask);
                disp('done Z');
            end
        end
        
        % doesn't really fit into this class
        function createG(obj, onsets, num_bins, cond_names, num_scans)
            runs = obj.subject_array.get_runs;
            n = size(runs,1);
            G = obj.calcG(onsets, num_bins, cond_names, num_scans);
            [G, ~, ~] = zscore(G, 1);
            for i = 1:n
                folder = obj.getFolder(runs{i});
                save([folder filesep 'G'], 'G');
                runs{i}.add_associated_matrix('G', [folder filesep 'G']);
            end
        end
        
        function createC(obj)
            runs = obj.subject_array.get_runs;
            n = size(runs,1);
            for i = 1:n
                try
                    obj.calcC(runs{i});
                    disp('done C');
                catch e
                    disp(e);
                    error(['Error at: ' num2str(i)]);
                end
            end
            
        end
        
        function XX = createCC(obj)
            runs = obj.subject_array.get_runs;
            n = size(runs,1);
            Cmat = load(runs{1}.get_associated_matrix('C'));
            XX = zeros(n*size(Cmat.C,1));
            numProcs = 4;
            parfor i = 1:numProcs
                matfile_objs = cellfun(@(x)matfile(x.get_associated_matrix('C')), runs, 'UniformOutput', 0);
                range = obj.calcBlock(i, size(matfile_objs{1}, 'C', 2), numProcs);
                arr = cellfun(@(x)(x.C(:, range)), matfile_objs, 'UniformOutput', 0);
                X = cell2mat(arr);
                XX = XX+X*X';
            end
        end
        
        function [U, D, V, P] = extractComponents(obj, CC)
            runs = obj.subject_array.get_runs;
            [U, D, V] = svds(CC, 3);
            G = cellfun(@(x)(getfield(load(x.get_associated_matrix('G')), 'G')), runs, 'UniformOutput', 0);
            P = sqrt(size(G,1))*sqrtm(inv(G'*G))*U;
            U = G*P;
        end
        
    end
    
    methods (Access = protected)
        function mask = calcMask(~, run)
            % takes each scan in the run and reads into scan_matrix.
            % then takes the grand mean of the scan matrix
            % voxels with values greater than the grand mean are set to
            % true. If over the time period 95% of the voxels are true,
            % this will make it into the final mask
            scans = run.get_scans;
            scan_matrix = cell2mat(cellfun(@(x)reshape(spm_read_vols(spm_vol(x)), 1, []), scans, 'UniformOutput', 0));
            grand_mean = mean(mean(scan_matrix));
            mask = scan_matrix>grand_mean;
            mask = floor(sum(mask)/(size(scan_matrix, 1)*0.95));
        end
        
        function Z = calculateZ(obj, run, mask)
            scans = run.get_scans;
            scan_matrix = cell2mat(cellfun(@(x)reshape(spm_read_vols(spm_vol(x)), 1, []), scans, 'UniformOutput', 0));
            Z = scan_matrix(:, mask);
            [Z, ~, ~] = zscore(Z, 1);
            folder = obj.getFolder(run);
            save([folder filesep 'Z'], 'Z');
            run.add_associated_matrix('Z', [folder filesep 'Z'])
        end
        
        function G = calcG(~, onsets, num_bins, cond_names, num_scans)
            designMatrix = Design(onsets, num_bins, cond_names, num_scans);
            G = designMatrix.createG;
        end
        
        function calcC(obj, run)
            folder = obj.getFolder(run);
            Zpath = run.get_associated_matrix('Z');
            Gpath = run.get_associated_matrix('G');
            Zmatrix = load(Zpath);
            Gmatrix = load(Gpath);
            C = Gmatrix.G\Zmatrix.Z;
            C = sqrtm(Gmatrix.G'*Gmatrix.G)*C;
            save([folder filesep 'C'], 'C', '-v7.3');
            run.add_associated_matrix('C', [folder filesep 'C']);
        end
        
        function range = calcBlock(~, index, numN, numP)
            blockSize = ceil(numN/(numP));
            startInd = ((index-1)*blockSize)+1;
            endInd = min(numN, (index*blockSize));
            range = startInd:endInd;
        end
        
        function folder = getFolder(~, run)
            scans = run.get_scans;
            folder = fileparts(scans{1});
        end
    end
    
end
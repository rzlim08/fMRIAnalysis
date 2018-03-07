 
% fs = NiftiFS(pwd);
% set_functional_dirstruct(fs, '{tOp_leVel}/{subject}/{scans}');
% set_structural_dirstruct(fs);
% set_subject_strmatch(fs, 's*');
% set_run_strmatch(fs, 'r*');
% set_scan_strmatch(fs, 'fs*');
% set_subjects(fs);
% %set_runs(fs);


% fs = NiftiFS(pwd);
% set_functional_dirstruct(fs, '{tOp_leVel}/{subject}/Functional/{runs}/{scans}');
% set_structural_dirstruct(fs);
% set_subject_strmatch(fs, 'SESOCD*');
% set_run_strmatch(fs, 'fMRI_ToL*');
% set_scan_strmatch(fs, 'g*');
% set_subjects(fs);
% set_runs(fs);
% set_is_nii(fs, 1);
% 
% remove_runs(fs, 'abort');
% rm(fs, {'clin';'abort';'SESOCD051';'nofMRI'})
% get_functional_scans(fs);


fs = NiftiFS(pwd);
set_functional_dirstruct(fs, '{top_level}/*/{groups}/{subjects}/{runs}/{scans}');
set_scan_strmatch(fs, 'f*');
set_subject_strmatch(fs, 's*');
set_subjects(fs);
set_runs(fs);
set_groups(fs);
tic;
run_on_cell_list(fs, @spm_get_space, get_functional_scans(fs));
toc;
tic;
run_on_cell_list_serial(fs, @spm_get_space, get_functional_scans(fs));
toc;
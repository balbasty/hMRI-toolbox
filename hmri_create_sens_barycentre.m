function jobsub = hmri_create_sens_barycentre(jobsub)
% Relative RF sensitivity calculation using a barycentric appproach.
% 
% Yael Balbastre, Wellcome Centre or Human Neuroimaging, 2020.

param = parse_param(jobsub);
% . struct{contrast}{echoe}         - MPM raw data
% . sens{contrast}{1=array|2=body}  - SMaps raw data
% . reg(1=array|2=body)             - Regularisation factors
% . dircalc                         - RF-specific directory
% . dirsuppl                        - Supplementary directory
% . json                            - ?
% . tag{contrast}                   - Contrast name (T1w, PDw, MTw)

% - Temporarily add multibias stuff to the path
if ~isdeployed
    dirscript  = fileparts(which('hmri_create_sens_barycentre'));
    oldpath    = path;
    cleanupObj = onCleanup(@() path(oldpath));
    path(fullfile(dirscript, 'multi-bias'), path);
    path(fullfile(dirscript, 'multi-bias', 'sub'), path);
end

% - Estimate bias fields
opt         = struct;
opt.coreg   = true;
opt.verbose = 3;
opt.output  = {'z' 'r' 'r0'};
opt.lambda  = [];
dat         = {};
sub         = logical([]);
for c=1:numel(param.sens)
    for m=1:numel(param.sens{c})
        opt.lambda(end+1) = param.reg(m);
        dat{end+1}        = param.sens{c}{m};
        if m == 1, sub(end+1) = true;
        else,      sub(end+1) = false; end
    end
end
[z,R,R0] = multibias(dat(:), opt);
z = z(:,:,:,sub);   % log-bias fields in mean space
R = R(:,:,sub);     % world smap-to-mean rigid matrix

% - Co-register structurals to smaps
M  = zeros(4,4,numel(param.sens));
dm = zeros(1,3,numel(param.sens));
for c=1:numel(param.sens)
    ref         = spm_vol(param.sens{c}{1});
    mov         = spm_vol(param.struct{c}{1});
    M(:,:,c)    = mov.mat;
    dm(1,:,c)   = mov.dim(1:3);
    q           = spm_coreg(ref, mov, struct('graphics', false));
    R(:,:,c)    = spm_matrix(q) * R(:,:,c);
end

% - Warp bias fields and apply
for c=1:numel(param.sens)
    
    nii0 = nifti(param.struct{c}{1});
    
    % - Warp bias field and save
    b    = exp(pull(z(:,:,:,c), R0\(R(:,:,c)\M(:,:,c)), dm(1,:,c)));
    nii  = nii0;
    sensname      = spm_file(spm_file(param.sens{c}{1}, 'basename'), ...
                        'suffix', '_sens', 'ext', '.nii');
    sensname      = fullfile(param.dircalc, sensname);
    nii.dat.fname = sensname;
    nii.dat.dtype = 'float32';
    nii.descrip   = 'Relative RF sensitivity';
    nii.dat.dim   = size(b);
    create(nii);
    nii.dat(:,:,:) = b*100;
    
    % - Metadata
    proc          = struct;
    proc.descrip  = ['hMRI toolbox - ' mfilename '.m - RF sensitivity correction'];
    proc.version  = hmri_get_version;
    proc.params   = param;
    output.imtype = 'sensitivity map';
    output.units  = 'p.u.';
    header        = init_output_metadata_structure(param.sens{c}{1}, proc, output);
    header.history.output.imtype = sprintf('Relative RF sensitivity map (HC/mean) for %sw images',param.tag{c});
    header.history.output.units  = 'p.u.';
    set_metadata(sensname,header,param.json);
        
    % - Copy file to supplementary folder
    copyfile(sensname,fullfile(param.dirsuppl,spm_file(sensname,'filename')));
    try copyfile([spm_str_manip(sensname,'r') '.json'],fullfile(dirsuppl,[spm_file(sensname,'basename'), '.json'])); end %#ok<*TRYNC>
    
    
    % - Apply bias field to each echo and save
    structname = cell(1,numel(param.struct{c}));
    for e=1:numel(param.struct{c})
        
        nii0 = nifti(param.struct{c}{e});
        nii  = nii0;
        structname{e} = spm_file(spm_file(nii.dat.fname, 'basename'), ...
                            'suffix', '_RFSC', 'ext', '.nii');
        structname{e} = fullfile(param.dircalc, structname{e});
        nii.dat.fname = structname{e};
        create(nii);
        nii.dat(:,:,:) = nii0.dat()./b;
        
        % - Metadata
        proc          = struct;
        proc.descrip  = ['hMRI toolbox - ' mfilename '.m - RF sensitivity correction'];
        proc.version  = hmri_get_version;
        proc.params   = param;
        output.imtype = 'sensitivity map';
        output.units  = 'p.u.';
        header        = init_output_metadata_structure(param.struct{c}{e}, proc, output);
        header.history.output.imtype = sprintf('RF sensitivity corrected %s-weighted echo',param.tag{c});
        header.history.output.units  = 'a.u.';
        header.acqpar = struct(...
            'RepetitionTime',   get_metadata_val(param.struct{c}{e},'RepetitionTime'), ...
            'EchoTime',         get_metadata_val(param.struct{c}{e},'EchoTime'), ...
            'FlipAngle',        get_metadata_val(param.struct{c}{e},'FlipAngle'));
        set_metadata(structname{e},header,param.json);
    end
    jobsub.raw_mpm.(param.tag{c}) = char(structname);
end

% -------------------------------------------------------------------------

function param = parse_param(job)

senstype = fieldnames(job.sensitivity);
senstype = senstype{1};
if job.sensitivity.(senstype).mode ~= 2
    error('We should be in barycentre mode.');
end

param.json     = hmri_get_defaults('json');
param.dircalc  = job.path.rfsenspath;
param.dirsuppl = job.path.supplpath;
param.reg      = [hmri_get_defaults('RFsens.reg.array')
                  hmri_get_defaults('RFsens.reg.body')];
if isempty(param.reg)
    warning(['Could not find default regularisation parameters. ' ...
             'Using [%g %g] instead.'], 1E7, 1E9);
end
              
param.struct = {};
param.sens   = {};
param.tag    = {};
possible_contrasts = {'MT' 'PD' 'T1'};
for c=1:numel(possible_contrasts)
    tag    = possible_contrasts{c};
    struct = job.raw_mpm.(tag);
    struct = spm_file(char(struct), 'number', '');
    if isempty(struct), continue; end
    sens   = job.sensitivity.(senstype).(sprintf('raw_sens_%s',tag));
    if isempty(sens), error('Missing sensitivity for %s', tag); end
    
    param.struct{end+1} = cellstr(struct);
    param.sens{end+1}   = cellstr(sens);
    param.tag{end+1}    = tag;
end
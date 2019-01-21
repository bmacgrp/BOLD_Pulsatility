% Requires: Tools for NIFTI and ANALYZE image to run
% Available at: http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% Also requires FSeries: https://www.mathworks.com/matlabcentral/fileexchange/31013-simple-real-fourier-series-approximation/content/Fseries.m
% and rsquare: https://www.mathworks.com/matlabcentral/fileexchange/34492-r-square--the-coefficient-of-determination
% Required Inputs: Vector of participant IDs, file of cardiac positions
% and fMRI file in .nifti format. Include entire directory of prefix if not
% in current path.
% Optional Inputs: flag_prefix: prefix of flag file for which volumes  to
%                  remove due to poor cardiac (1 - keep, 0 - discard). If
%                  not provided, program assumes all volumes are to be kept
%                  number of volumes to normalize analysis to. Defaults to
%                  smallest number of volumes in flag file. Value must be
%                  lower than smallest value in totalvolumes from creating
%                  cardiac files - cardiacfiles-general.r
%                  number of repetitions: how many permutations to conduct.
%                  Defaults to 45000.
%                  Threshold: Minimum intensity of voxels to analyze,
%                  adjust based on image intensity. Set to 0 if you do not wish to
%                  threshold. Defaults to 2000.
%                  Significance threshold: Minimum number of standard deviations
%                  from random mean where fit is considered significant.
%                  Defaults to 5.
%                  Output prefix: Name prefix of output files. Defaults to
%                  input file prefix.
% Outputs: .nifti file with Fourier series values and R^2 value, .nifti file with SD values, number 
% of significant voxels, total number of voxels and mean voxel value (output as vectors)
%Created by Athena Theyers, modified by Sarah Atwi, January 8th, 2019.

function [sigvox, numvox, meanval] = mapcardiacpulsatility_fourservalues(ID,input_prefix,cardiac_prefix,flag_prefix,volumes,reps,thresh,sig,output_prefix)

	if nargin > 9
		error('Too many inputs, requires at most 9')
	end

	if nargin < 3
		error('Too few inputs, requires at minimum ID vector, 4D fMRI scan, cardiac file')
	end

	switch nargin
        case 3
			flag_prefix = -1; % no flag file included, assume max volumes to be used. Set to -1 so can change at later steps.
            volumes = -1; % if no volumes entered, default to maximum. Set to -1 so can change at later steps.
            reps = 45000; % number of permutations to calulate
			thresh = 2000; % threshold for input scan
			sig = 5; % significance threshold
			output_prefix = input_prefix;
		case 4
			volumes = -1; % if no volumes entered, default to maximum. Set to -1 so can change at later steps.
            reps = 45000; % number of permutations to calulate
			thresh = 2000; % threshold for input scan
			sig = 5; % significance threshold
			output_prefix = input_prefix;
		case 5
			reps = 45000; % number of permutations to calulate
			thresh = 2000; % threshold for input scan
			sig = 5; % significance threshold
			output_prefix = input_prefix;
		case 6
			thresh = 2000; % threshold for input scan
			sig = 5; % significance threshold
			output_prefix = input_prefix;
		case 7
			sig = 5; % significance threshold
			output_prefix = input_prefix;
		case 8
			output_prefix = input_prefix;
		case 9
		
	end

	r2 = zeros(1,reps); 

	sigvox=zeros(numel(ID),1);
	numvox=zeros(numel(ID),1);
    meanval=zeros(numel(ID),1);
     
	name = strcat(input_prefix,num2str(ID(1)),'.nii.gz');
	nii = load_nii(name,[],[],[],[],[],0.5);
    if flag_prefix ~= -1
        read = repmat('%f ',1,numel(ID));
        fileID = fopen(strcat(flag_prefix,'.txt'),'r'); % file of volumes to remove (1-keep, 0-remove)
        flags = textscan(fileID, read);
        fclose(fileID);
    else
        flags = ones(nii.hdr.dime.dim(5),numel(ID));
        flags = num2cell(flags,1);
        volumes = nii.hdr.dime.dim(5);
    end
    
    testflags = flags{1};
    maxvol = numel(testflags(testflags==1));
    
    for l = 2:numel(ID)
        testflags = flags{l};
        if numel(testflags(testflags==1)) < maxvol
            maxvol =  numel(testflags(testflags==1));
        end
    end
    
    if volumes == -1
        volumes = maxvol;
    end

	if volumes > maxvol
        volumes = maxvol;
        warning('The number of volumes you have chosen, exceeds maximum number allowed with current flags. Reducing to maximum available.');
	end

	name4 = strcat(cardiac_prefix,num2str(ID(1)),'.txt');
    slices = nii.hdr.dime.dim(4);
	read = repmat('%f ',1,slices);
	fileID = fopen(name4,'r');
	cardiac = textscan(fileID, read);
	fclose(fileID);
	flag = flags{1};
	xres = nii.hdr.dime.dim(2);
	yres = nii.hdr.dime.dim(3);
	i = round(xres/2);
	j = round(yres/2);
	k = round(slices/2);
	card = round(cardiac{k}*100)/100;
	card = card(flag~=0);
	y = randsample(numel(card),volumes);
	card = card(y);
	point = squeeze(nii.img(i,j,k,:));
	point = point(flag~=0);
	point = point(y);
	for n=1:reps
		card1 = datasample(card,numel(card),'Replace',false); 
		[~,~,yfit] = Fseries(card1,point,3);
		[r2(n), ~] = rsquare(point,yfit); 
	end

	for m=1:numel(ID)
		name = strcat(input_prefix,num2str(ID(m)),'.nii.gz'); % input fMRI file
		name2 = strcat(output_prefix,num2str(ID(m)),'R2.nii.gz'); % output fMRI file containing R^2 value and Fourier series values
		name3 = strcat(output_prefix,num2str(ID(m)),'Std.nii.gz'); % output fMRI file containing SD value
		name4 = strcat(cardiac_prefix,num2str(ID(m)),'.txt'); % cardiac position file
		
		nii = load_nii(name,[],[],[],[],[],0.8);
		xres = nii.hdr.dime.dim(2);
		yres = nii.hdr.dime.dim(3);
		nii2 = nii;
		nii2.img = zeros(xres,yres,slices,8);
		nii2.hdr.dime.dim(5) = 8;
		nii3 = nii2;
		nii3.img = zeros(xres,yres,slices,1);
		nii3.hdr.dime.dim(5) = 1;
		read = repmat('%f ',1,slices);
		fileID = fopen(name4,'r');
		cardiac = textscan(fileID, read);
		fclose(fileID);
		flag = flags{m};
		for k=1:slices
			cardiac{k} = round(cardiac{k}*100)/100;
			card = cardiac{k};
			card = card(flag~=0);
			y = randsample(numel(card),volumes);
			card = card(y);
			for i=1:xres
				for j = 1:yres
					if min(squeeze(nii.img(i,j,k,:))) > thresh
						point = squeeze(nii.img(i,j,k,:));
						point = point(flag~=0);
						point = point(y);
						[a,b,yfit] = Fseries(card,point,3);
						[nii2.img(i,j,k,1), ~] = rsquare(point,yfit); 
						nii2.img(i,j,k,2:2:8) = a;
						nii2.img(i,j,k,3:2:7) = b;
					end
				end
			end
		end
		save_nii(nii2,name2);

		for k=1:slices
			for i=1:xres
				for j = 1:yres
					if nii2.img(i,j,k,1) ~= 0
						nii3.img(i,j,k,1) = (nii2.img(i,j,k,1)-mean(r2))/std(r2);
					end
				end
			end
		end
		save_nii(nii3,name3);

		sigvox(m)=numel(nii3.img(nii3.img>sig)); %number of voxels considered to have significant pulsatility
		numvox(m)=numel(nii3.img(nii3.img~=0)); %number of voxels with signal
		meanval(m) = mean(nii3.img(nii3.img~=0)); %mean number of standard deviations across entire scan
	end

end
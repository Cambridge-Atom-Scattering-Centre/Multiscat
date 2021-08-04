classdef Multiscat
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function PreparePotentialFiles(potStructArray)
            
            if ~isstruct(potStructArray)
                V=potStructArray; clear potStructArray;
                potStructArray(1).V=V;
            end
            
            for i=1:length(potStructArray)
                
                [Nx,Ny,Nz] = size(potStructArray(i).V);
                V_FT = fft2(potStructArray(i).V)/(Nx*Ny); % Assuming V(:,:,i) is the potential on the XY plane at Z(:,:,i)
                shifted_V_FT = fftshift(fftshift(V_FT,1),2);
                
                filePot = fopen(['./pot' num2str(10000+i) '.in'],'w');
                % Multiscat is written to treat first 5 lines as description lines
                fprintf(filePot,'Dummy line1\nDummy line2\nDummy line3\nDummy line4\nDummy line5\n');
                
                for n_i=1:Nx
                    for m_i=1:Ny
                        singleFC = squeeze(shifted_V_FT(m_i,n_i,:));
                        for l_i=1:Nz
                            potStr=sprintf('(%+0.6e, %+0.6e)\n', real(singleFC(l_i)), imag(singleFC(l_i)));
                            fprintf(filePot,'%s', potStr);
                        end
                    end
                end
                fclose(filePot);                
            end            
        end
        
        function prepareFourierLabels(V)
            
            [Nx,Ny,Nz] = size(V); % Assuming V(:,:,i) is the potential on the XY plane at Z(:,:,i)
            fileFourierLabels = fopen('FourierLabels.in','w');
            for n_i=1:Nx
                n=-Nx/2+n_i-1;
                for m_i=1:Ny
                    m=-Ny/2+m_i-1;
                    fprintf(fileFourierLabels,'%d %d\n', m, n);
                end
            end
            fclose(fileFourierLabels);
        end
        
        function confStruct = createConfigStruct(potStructArray)

            
            % realSpaceStructArray - fields: X, Y, Z (which describe the spatial space for V            
            a1_tmp=potStructArray(1).a1; a2_tmp=potStructArray(1).a2;
            angl = acos(a1_tmp*a2_tmp'/norm(a1_tmp)/norm(a2_tmp))*180/pi;
            a1=norm(a1_tmp);
            a2= cosd(angl)*norm(a2_tmp);
            b2= sind(angl)*norm(a2_tmp);
            

            [Nx,Ny,Nz] = size(potStructArray(1).V); % Assuming V(:,:,i) is the potential on the XY plane at Z(:,:,i)
            for i=1:length(potStructArray)
                tmp_minV(i) = min(min(min(potStructArray(i).V)));
            end
            minV = min(tmp_minV);
            
            % multiscat parameters
            confStruct.itestComment = 'itest=1 enables output of each diffraction intensity; itest=0 outputs specular only';
            confStruct.itest = 1;
            confStruct.gmresPreconditionerFlagComment = 'gmres preconditioner flag (ipc)';
            confStruct.gmresPreconditionerFlag = 0;
            confStruct.numSigFigConvergenceComment = 'number of significant figures convergence (nsf)';
            confStruct.numSigFigConvergence = 5;	
            confStruct.numOfFCComment = 'total number of fc';
            confStruct.numOfFC = Nx*Ny;
            confStruct.zIntegrationRangeComment = 'integration range (zmin,zmax)';
            confStruct.zIntegrationRange = [potStructArray(1).zmin potStructArray(1).zmax];
            confStruct.vminComment = 'potential well depth (vmin)';
            confStruct.vmin = abs(minV)+2;
            confStruct.dmaxComment = 'max -ve energy of closed channels (dmax)';
            confStruct.dmax = 120;
            confStruct.imaxComment = 'max index of channels (imax)';
            confStruct.imax = 120;
            confStruct.energiesToCalcComment = 'energies to calculate at (initial energy, step, number of steps)';
            confStruct.scatsub_a1Comment = 'a1 (see subroutine basis in "scatsub.f")';
            confStruct.scatsub_a1 = a1;
            confStruct.scatsub_a2Comment = 'a2';
            confStruct.scatsub_a2 = a2;
            confStruct.scatsub_b2Comment = 'b2';
            confStruct.scatsub_b2 = b2;
            confStruct.nzfixedComment = 'number of fixed z points,nzfixed';
            confStruct.nzfixed = potStructArray(1).zPoints;
            confStruct.stepzminComment = 'stepzmin  (max and min z values of fixed z points)';
            confStruct.stepzmin = potStructArray(1).zmin;
            confStruct.stepzmaxComment = 'stepzmax';
            confStruct.stepzmax = potStructArray(1).zmax;
            confStruct.startIndexComment = 'startindex';
            confStruct.startIndex = 10001;
            confStruct.endIndexComment = 'endindex';
            confStruct.endIndex = 10000 + length(potStructArray);
            confStruct.incidentParticleMassComment = 'helium mass';
            confStruct.incidentParticleMass = 4;
            disp('Mass set to 4')
        end
       

        function prepareConfigFile(confStruct)
            
            fileConf = fopen('Multiscat.conf','w');
            
            % Print conf file
            fprintf(fileConf,'% s \n', 'FourierLabels.in 	!The fourier labels input file');
            fprintf(fileConf,'% s \n', 'scatCond.in	! The scattering conditions input file');
 % the above two lines are written by Boyao on 17 Nov 2020 to make it
 % compatible with the latest version of multiscat
            fprintf(fileConf,'%d       !%s\n', confStruct.itest, confStruct.itestComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.gmresPreconditionerFlag, confStruct.gmresPreconditionerFlagComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.numSigFigConvergence, confStruct.numSigFigConvergenceComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.numOfFC, confStruct.numOfFCComment);
            fprintf(fileConf,'%d,%d       !%s\n', confStruct.zIntegrationRange(1), confStruct.zIntegrationRange(2), confStruct.zIntegrationRangeComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.vmin, confStruct.vminComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.dmax, confStruct.dmaxComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.imax, confStruct.imaxComment);
%             fprintf(fileConf,'%s       !%s\n', confStruct.energiesToCalc, confStruct.energiesToCalcComment);
%             fprintf(fileConf,'%s       !%s\n', confStruct.energiesArrayToCalc, confStruct.energiesArrayToCalcComment);
%             fprintf(fileConf,'%s       !%s\n', confStruct.thetaValuesToCalc, confStruct.thetaValuesToCalcComment);
%             fprintf(fileConf,'%s       !%s\n', confStruct.thetaArrayToCalc, confStruct.thetaArayToCalcComment);
%             fprintf(fileConf,'%s       !%s\n', confStruct.phiValuesToCalc, confStruct.phiValuesToCalcComment);
%             fprintf(fileConf,'%s       !%s\n', confStruct.phiArrayToCalc, confStruct.phiArrayToCalcComment);
            % the above five line are commented by Boyao on 17 Nov 2020 because they
            % seem incompatible with the latest version of multiscat
            fprintf(fileConf,'%d       !%s\n', confStruct.scatsub_a1, confStruct.scatsub_a1Comment);
            fprintf(fileConf,'%d       !%s\n', confStruct.scatsub_a2, confStruct.scatsub_a2Comment);
            fprintf(fileConf,'%d       !%s\n', confStruct.scatsub_b2, confStruct.scatsub_b2Comment);
            fprintf(fileConf,'%d       !%s\n', confStruct.nzfixed, confStruct.nzfixedComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.stepzmin, confStruct.stepzminComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.stepzmax, confStruct.stepzmaxComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.startIndex, confStruct.startIndexComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.endIndex, confStruct.endIndexComment);
            fprintf(fileConf,'%d       !%s\n', confStruct.incidentParticleMass, confStruct.incidentParticleMassComment);
            
            fclose(fileConf);
        end
        
        function [dpars]=readDiffracOut(filename)
            
            if exist('filename','var')==0
                error('MATLAB:DAN:NoFile','Please provide a filename when calling read_sr')
            end
            
            line_counter=0;
            fid = fopen(filename,'r');
            
            tline = fgetl(fid);
            line_counter=line_counter+1;
            
            ind1=0;
            while (1)
                if (size(tline,2)>46) && strcmp(tline(1:47),' Required number of z grid points, m =         ') %start of new block
                    ind1=ind1+1;
                    numstr = tline(48:size(tline,2));
                    dpars.m(ind1)=str2double(numstr);
                elseif (size(tline,2)>44) && strcmp(tline(1:44),' Number of diffraction channels, n =        ')
                    numstr = tline(45:size(tline,2));
                    dpars.n(ind1)=str2double(numstr);
                elseif (size(tline,2)>24) && strcmp(tline(1:25),' # Beam energy           ')
                    % Identify the number of points
                    numstr = tline(30:size(tline,2));
                    dpars.ei(ind1) = str2double(numstr);
                elseif (size(tline,2)>31) && strcmp(tline(1:32),' # Polar angle        theta =   ')
                    % Identify Imax
                    numstr = tline(33:44);
                    dpars.theta(ind1) = str2double(numstr);
                    numstr = tline(46:end);
                    dpars.phi(ind1) = str2double(numstr);

                    tline = fgetl(fid); ind2=1;
                    while strcmp(tline(1),'#')
                        numstr=sscanf(tline(2:end),'%f%f%f');
                        
                        if ind1>1 && length(dpars.preB2(ind1-1,:)) < ind2,
                            dpars.preB2(1:ind1-1,ind2)=NaN;
                            dpars.preB1(1:ind1-1,ind2)=NaN;
                            dpars.Idiff(1:ind1-1,ind2)=NaN;
                        end
                        
                        dpars.preB2(ind1,ind2)=numstr(1);                        
                        dpars.preB1(ind1,ind2)=numstr(2);                        
                        dpars.Idiff(ind1,ind2)=numstr(3);
                        ind2=ind2+1;
                        tline = fgetl(fid);
                        line_counter=line_counter+1;
                        dpars.openChannels(ind1) = ind2 - 1;
                        if (size(tline,2)>46) && strcmp(tline(1:47),' Required number of z grid points, m =         ') %start of new block

                            ind1=ind1+1;
                            numstr = tline(48:size(tline,2));
                            dpars.m(ind1)=str2double(numstr);
                        end
                    end
                end

                tline = fgetl(fid);
                line_counter=line_counter+1;
                if tline==-1
                    break
                end


            end

            A = unique([reshape(dpars.preB2',numel(dpars.preB2),1) reshape(dpars.preB1',numel(dpars.preB2),1)],'rows');

            [noOfChannels, ~] = size(A); 

            eis = unique(dpars.ei);
            noOfEis = numel(eis);
            for channelInd = 1:noOfChannels
                
                % m and n for this row
                m = A(channelInd, 1);
                n = A(channelInd, 2);

                eiInd = 1;
                while eiInd <= noOfEis
                    % energy is eis(eiInd)

                    % find value of I_mn at this energy
                    ind2s = find(dpars.preB2(eiInd,1:dpars.openChannels(eiInd)) == m);
                    ind2 = find(dpars.preB1(eiInd,ind2s) == n);
                    ind2 = ind2s(ind2);

                    if isempty(ind2) 
                        A(channelInd,eiInd+2) = NaN;
                    else
                        A(channelInd,eiInd+2) = dpars.Idiff(eiInd,ind2);
                    end

                    eiInd = eiInd + 1;

                end

            end

            dpars.channelTable = A;
            
        end
        
        function fortObj = readOutputs(filename)
            %% Import data from text file filename
            if exist(filename,'file') == 2
                
                fprintf([filename ', reading file...  ']);

                % Initialize variables.
                delimiter = ' ';
                startRow = 37;

                % Format string for each line of text:
                %   column3: double (%f)
                %       column4: text (%s)
                %   column5: double (%f)
                %       column6: text (%s)
                %   column7: double (%f)
                %       column8: text (%s)
                %   column9: double (%f)
                %       column10: text (%s)
                %   column11: double (%f)
                % For more information, see the TEXTSCAN documentation.
                formatSpec = '%*s%*s%f%s%f%s%f%s%f%s%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

                % Open the text file.
                fileID = fopen(filename,'r');

                % Read columns of data according to format string.
                % This call is based on the structure of the file used to generate this
                % code. If an error occurs for a different file, try regenerating the code
                % from the Import Tool.
                dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);

                % Close the text file.
                fclose(fileID);

                % Post processing for unimportable data.
                % No unimportable data rules were applied during the import, so no post
                % processing code is included. To generate code which works for
                % unimportable data, select unimportable cells in a file and regenerate the
                % script.

                % Allocate imported array to column variable names
                fortObj.E = dataArray{:, 1};
                fortObj.theta = dataArray{:, 3};
                fortObj.phi = dataArray{:, 5};
                fortObj.I00 = dataArray{:, 7};
                fortObj.sum = dataArray{:, 9};

                fprintf('Done\n');

                % Clear temporary variables
                clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;
                
            else
                % If fort.6 does not exist
                fprintf('Error: fort.6 does not exist\n');
                fortObj = 0;
            end
        end
    end
    
end


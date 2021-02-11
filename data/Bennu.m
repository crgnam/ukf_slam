classdef Bennu < handle
    properties (Access = public)
        % Gravity information:
        G = 6.67430*10^-11; %(m3⋅kg–1⋅s–2) Gravitational Constant
        M = (7.329*10^10); %(kg) Mass of Bennu
        mu % Standard gravitational constant of Bennu
        
        % landmark locations
        lmks %3XN of landmark locations
        lmk_norms %3XN of landmark surface normals
        
        % Visualization stuff:
        ax % handle for the plotting axis
        p  % handle for patch object
        l  % handle for light object
        lmk_vis % handle for visible landmarks
        lmk_inv % handle for unseen landmarks
        plotted_lmks = false
        
        % Components of the parsed .obj file
        fname = 'Bennu-Radar.obj'
        v
        vt
        vn
        f
    end
   
    %% Constructor:
    methods (Access = public)
        function [self] = Bennu(num_lmks)
            self.mu = self.G*self.M;
            self.readObj(self.fname);
            
            % Rescale from km to m:
            self.v = 1000*self.v;
            
            % Randomly sample shape for lmks:
            idx = randperm(size(self.v,1));
            self.lmks = self.v(idx(1:num_lmks),:)';
            self.lmk_norms = normc(self.lmks);
        end
    end
    
    %% Public Methods:
    methods (Access = public)
        function [] = drawBody(self)
            self.p = patch('Faces',self.f.v,'Vertices',self.v,...
                           'FaceColor',[0.5 0.5 0.5],'EdgeColor','None',...
                           'FaceLighting','gouraud','AmbientStrength',0.5,...
                           'SpecularStrength',0);
            self.ax = gca;
            self.l = light(self.ax);
            hold on; axis equal
            set(gca,'Color','k')
        end
        
        function [] = drawBodyStandalone(self,axs)
            patch(axs,'Faces',self.f.v,'Vertices',self.v,...
                           'FaceColor',[0.5 0.5 0.5],'EdgeColor','None',...
                           'FaceLighting','gouraud','AmbientStrength',0.5,...
                           'SpecularStrength',0);
            self.l = light(axs);
            hold on; axis equal
            set(axs,'Color','k')
        end
        
        function [] = drawLmks(self,visible,varargin)
            if numel(visible) == 1
                inds = true(1,size(self.lmks,2));
            elseif length(visible) == size(self.lmks,2)
                inds = visible;
            else
                error('Invalid input for VISIBLE')
            end
            if ~self.plotted_lmks
                self.lmk_vis = plot3(self.lmks(1,inds),...
                                     self.lmks(2,inds),...
                                     self.lmks(3,inds),'.g',varargin{:});
                self.lmk_inv = plot3(self.lmks(1,~inds),...
                                     self.lmks(2,~inds),...
                                     self.lmks(3,~inds),'.b',varargin{:});
                self.plotted_lmks = true;
            else
                set(self.lmk_vis,'XData',self.lmks(1,inds),...
                                 'YData',self.lmks(2,inds),...
                                 'ZData',self.lmks(3,inds));
                set(self.lmk_inv,'XData',self.lmks(1,~inds),...
                                 'YData',self.lmks(2,~inds),...
                                 'ZData',self.lmks(3,~inds));
            end
        end
    end
    
    %% Private properties and methods  
    methods (Access = private)
        function [self] = readObj(self,fname)
            self.v = []; self.vt = []; self.vn = []; self.f.v = []; self.f.vt = []; self.f.vn = [];
            
            % Parse .obj file: 
            fid = fopen(fname);
            while 1    
                tline = fgetl(fid);
                if ~ischar(tline),   break,   end  % exit at end of file 
                 ln = sscanf(tline,'%s',1); % line type 
                 %disp(ln)
                switch ln
                    case 'v'   % mesh vertexs
                        self.v = [self.v; sscanf(tline(2:end),'%f')'];
                    case 'vt'  % texture coordinate
                        self.vt = [self.vt; sscanf(tline(3:end),'%f')'];
                    case 'vn'  % normal coordinate
                        self.vn = [self.vn; sscanf(tline(3:end),'%f')'];
                    case 'f'   % face definition
                        fv = []; fvt = []; fvn = [];
                        str = textscan(tline(2:end),'%s'); str = str{1};

                       nf = length(findstr(str{1},'/')); % number of fields with this face vertices
                       [tok str] = strtok(str,'//');     % vertex only
                        for k = 1:length(tok) fv = [fv str2num(tok{k})]; end

                        if (nf > 0) 
                        [tok str] = strtok(str,'//');   % add texture coordinates
                            for k = 1:length(tok) fvt = [fvt str2num(tok{k})]; end
                        end
                        if (nf > 1) 
                        [tok str] = strtok(str,'//');   % add normal coordinates
                            for k = 1:length(tok) fvn = [fvn str2num(tok{k})]; end
                        end
                         self.f.v = [self.f.v; fv]; self.f.vt = [self.f.vt; fvt]; self.f.vn = [self.f.vn; fvn];
                end
            end
            fclose(fid);
        end
    end
end
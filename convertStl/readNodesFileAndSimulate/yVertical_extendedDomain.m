clear;
filename = 'nodes.out';

TX = 764; TY = 29; TZ = 82;

PX=10;
PY=1;
PZ=1;

fid = fopen(filename,'r');
v = fscanf(fid, '%f');
fclose(fid);

assert(numel(v) == TX*TY*TZ, '数据量不匹配：期望 %d，实际 %d', TX*TY*TZ, numel(v));
tmp = reshape(v, [TZ, TY, TX]);   % tmp(k, j, i)

A = permute(tmp, [1 2 3]);       % A(i, j, k) 尺寸: (nx, ny, nz)
TZ = 42
disp(size(A));
disp([min(A(:)) max(A(:))]);

for i=1:TX
    for j=1:TY
        for k=1:TZ
            G(TX-i+1,j,k)=A(k,j,i);
        end
    end
end

%for i = 765:800
%  for j = 1:TY
%    for k = 1:TZ
%        if j < 5
%          G(i,j,k) = 2;
%        endif
%      endfor
%    endfor
%  endfor

for i=1:TX
    for j=1:TY
        for k=1:TZ
%            if G(i,j,k)==2
%                G(i,j,k)=3;  %solid with contact angle 90
%            end
            if G(i,j,k)==1
                G(i,j,k)=2;  %solid with contact angle 30
            end
        end
    end
end

for i=1:TX   %liquid initialization
    for j=1:TY
        for k=1:TZ
%            if i>TX/2-40 && i<TX/2+40 && j<TY/2+20 && G(i,j,k)==0
%                G(i,j,k)=1;
%            end

            if G(i,j,k)==-1
                G(i,j,k)=1;
            end
        end
    end
end

Nx_extra = 100
G_ext = ones(TX + Nx_extra, TY, TZ);

G_ext(Nx_extra:Nx_extra-1+TX,:,:) = G;

TX_old = TX;
TX = TX + Nx_extra;

G = G_ext;

for i = 1:100
    for j = 1:10
        for k = 1:TZ
            G(i,j,k) = 2;
        end
    end
end

isosurface(G);
%image(60*G(:,:,82)');
axis equal;

mkdir('data');
for a=1:PX*PY*PZ
        rankz=ceil(a/(PX*PY));
        ranky=ceil((a-(rankz-1)*PX*PY)/PX);
        rankx=a-(rankz-1)*PX*PY-(ranky-1)*PX;
        rank(a,1)=rankx;
        rank(a,2)=ranky;
        rank(a,3)=rankz;

        if rankx<=mod(TX,PX)
            lengthx=ceil(TX/PX);
            startx=(rankx-1)*lengthx;
        else
            lengthx=floor(TX/PX);
            startx=(rankx-1)*lengthx+mod(TX,PX);
        end

        if ranky<=mod(TY,PY)
            lengthy=ceil(TY/PY);
            starty=(ranky-1)*lengthy;
        else
            lengthy=floor(TY/PY);
            starty=(ranky-1)*lengthy+mod(TY,PY);
        end

        if rankz<=mod(TZ,PZ)
            lengthz=ceil(TZ/PZ);
            startz=(rankz-1)*lengthz;
        else
            lengthz=floor(TZ/PZ);
            startz=(rankz-1)*lengthz+mod(TZ,PZ);
        end

        start(a,1)=startx;
        start(a,2)=starty;
        start(a,3)=startz;
        length(a,1)=lengthx;
        length(a,2)=lengthy;
        length(a,3)=lengthz;

        NX=lengthx+2;
        NY=lengthy+2;
        NZ=lengthz+2;

        data=zeros(NX,NY,NZ);

        for i=1:NX
            for j=1:NY
                for k=1:NZ
                    id=startx+i-1;
                    if id==0
                        id=TX;
                    end
                    if id==TX+1
                        id=1;
                    end

                    jd=starty+j-1;
                    if jd==0
                        jd=TY;
                    end
                    if jd==TY+1
                        jd=1;
                    end

                    kd=startz+k-1;
                    if kd==0
                        kd=TZ;
                    end
                    if kd==TZ+1
                        kd=1;
                    end

                    data(i,j,k)=G(id,jd,kd);
                end
            end
        end

        mpirank=num2str(a-1,'%04d');
        filename=['./data/data',mpirank,'.dat'];
        fid=fopen(filename,'w');
        for i=1:NX
            for j=1:NY
                for k=1:NZ
                    fprintf(fid,'%i\n',data(i,j,k));
                end
            end
        end
        fclose(fid);
end

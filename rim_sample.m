function h=rim_sample(mi, so, ro, be, Np, Nr, Tw, Fc)
% mi (microphone), so (source) and ro (room) are three-dimensional column vectors.
% Np: samples of the RIR.
% Nr: no. of random samples (Nr=0 for original IM).
% Tw: samples of low-pass filter, Fc: cut-off freq.
% All quantities above are in sample periods.
% be: matrix of refl. coeff. [x1,y1,z1;x2,y2,z2]
h=zeros(Np,1); ps=perm([0,1],[0,1],[0,1]);
Rps=repmat(so,[1,8])+(2.*ps-1).*repmat(mi,[1,8]);
or=floor(Np./(ro.*2))+1;
rs=perm(-or(1):or(1),-or(2):or(2),-or(3):or(3));
for i=1:size(rs'); r=rs(:, i);
 for j=1:8; p=ps(:,j); Rp=Rps(:,j);
  d=norm(2*ro.*r+Rp)+1+Nr*(2*rand-1);
  if round(d)>Np || round(d)<1; continue; end
  am=be(1,:)'.^abs(r+p).*be(2,:)'.^abs(r);
  if Tw==0; n=round(d); else
    n=(max(ceil(d-Tw/2),1):min(floor(d+Tw/2),Np))';
  end
  s=(1+cos(2*pi*(n-d)/Tw)).*sinc(Fc*(n-d))/2;
  s(isnan(s))=1; h(n)=h(n)+s*prod(am)/(4*pi*(d-1));
end; end;
function res=perm(varargin)
[res{1:nargin}]=ndgrid(varargin{1:nargin});
res=reshape(cat(nargin+1,res{:}),[],nargin)';
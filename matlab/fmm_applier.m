function out=fmm_applier(kh,source_in,charge,target_in)
%those are the source points
%they will be the points in the circle
%source should be 2Xnsource
[msource,nsource]=size(source_in);
source=source_in;
if (msource~=2) 
    nsource=msource;
    source=source_in';
end

srcinfo = [];
srcinfo.sources = source;
srcinfo.charges = charge(:).';

pg = 0;
pgt = 1;

%charge are the total charges in each source point

%those are the target points
%here the points can be in any configuration matrix or line
[mtarget,ntarget]=size(target_in);
target=target_in;
if (mtarget~=2) 
    ntarget=mtarget;
    target=target_in';
end
eps = 1e-7;
U = hfmm2d(eps,kh,srcinfo,pg,target,pgt);

%setting the solution
out=U.pottarg;

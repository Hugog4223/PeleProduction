

c ::: SCCS stuff "@(#)FABUTIL_2D.F        3.1\t6/25/93"
c ::: --------------------------------------------------------------
      subroutine cinterp2d (fine,floi1,floi2,fhii1,fhii2, fblo,fbhi,nvar
     &,lratio, crse,clo,chi,cblo,cbhi,fslo,fshi, cslope,clen,fslope,fd
     &at,flen,voff)

      integer floi1,floi2, fhii1,fhii2

      integer fblo(2), fbhi(2)
      integer cblo(2), cbhi(2)
      integer fslo(2), fshi(2)
      integer lratio, nvar, clen, flen, clo, chi
      DOUBLE PRECISION fine(floi1 :fhii1 ,floi2 :fhii2, nvar)
      DOUBLE PRECISION crse(clo:chi, nvar)
      DOUBLE PRECISION cslope(clo:chi, 2)
      DOUBLE PRECISION fslope(flen, 2)
      DOUBLE PRECISION fdat(flen)
      DOUBLE PRECISION voff(flen)

c ::: NOTE: data must be sent in so that
c ::: cslope(1,*) and crse(1,*) are associated with
c ::: the same cell
c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ::: cinterp:   conservative interpolation from coarse grid to
c ::: subregion of fine grid defined by (fblo,fbhi)
c :::
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio       =>  (const)  refinement ratio between levels
c :::
c ::: crse         =>  (const)  coarse grid data widended by 1 zone
c ::: and unrolled
c ::: clo,chi      =>  (const)  1.0D0 dimensional limits of crse grid
c ::: cslo,cshi    =>  (const)  coarse grid index limits where
c ::: slopes are to be defined. This is
c ::: the projection of (fblo,fbhi) down
c ::: to the coarse level
c ::: fslo,fshi    =>  (const)  fine grid index limits where
c ::: slopes are needed.  This is the
c ::: refinement of (cslo,cshi) and
c ::: contains but may not be identical
c ::: to (fblo,fbhi).
c ::: cslope       =>  (modify) temp array coarse grid slopes
c ::: clen         =>  (const)  length of coarse gtid slopes
c ::: fslope       =>  (modify) temp array for fine grid slope
c ::: flen         =>  (const)  length of fine grid slope array
c ::: fdat         =>  (const)  temp array for fine grid data
c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ::: local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      DOUBLE PRECISION hafrat, volratio
      DOUBLE PRECISION cen, forw, back, slp
      DOUBLE PRECISION xoff, yoff
      integer ncbx, ncby
      integer ncsx, ncsy
      integer islo, jslo
      integer icc, istart, iend
      integer lenx, leny, maxlen

      hafrat = 0.5D0*float(lratio-1)
      volratio = 1.0D0/float(lratio)

      ncbx = cbhi(1)-cblo(1)+1
      ncby = cbhi(2)-cblo(2)+1
      ncsx = ncbx+2
      ncsy = ncby+2
      ist = 1
      jst = ncsx
      islo = cblo(1)-1
      jslo = cblo(2)-1
      lenx = fbhi(1)-fblo(1)+1
      leny = fbhi(2)-fblo(2)+1
      maxlen = max(lenx,leny)
      if (maxlen .eq. lenx) then
          do 100 i = fblo(1), fbhi(1)
              fn = i-fslo(1)+1
              ioff = mod(fn-1,lratio)
              voff(fn) = float(ioff)-hafrat
100       continue
      else
          do 110 j = fblo(2), fbhi(2)
              fn = j-fslo(2)+1
              joff = mod(fn-1,lratio)
              voff(fn) = float(joff)-hafrat
110       continue
      end if
      do 120 n = 1, nvar

c ::: ::::: compute slopes in x direction
          do 130 i = 1, clen
              cen = 0.5D0*(crse(i+ist,n)-crse(i-ist,n))
              forw = crse(i+ist,n)-crse(i,n)
              back = crse(i,n)-crse(i-ist,n)
              slp = sign(1.0D0,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,1)=merge(slp,0.0D0,forw*back>=0.0d0)
130       continue
c ::: ::::: compute slopes in y direction
          do 140 i = 1, clen
              cen = 0.5D0*(crse(i+jst,n)-crse(i-jst,n))
              forw = crse(i+jst,n)-crse(i,n)
              back = crse(i,n)-crse(i-jst,n)
              slp = sign(1.0D0,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,2)=merge(slp,0.0D0,forw*back>=0.0d0)
140       continue
          if (maxlen .eq. lenx) then
              do 150 jc = cblo(2), cbhi(2)

c ::: ..,.......::::: strip out a fine grid slope vector
                  do 160 ioff = 1, lratio
                      icc = clo + ist + jst*(jc-jslo)
                      istart = ioff
                      iend = ioff + (ncbx-1)*lratio
                      do 170 fn = istart, iend, lratio
                          fslope(fn,1) = cslope(icc,1)
                          fslope(fn,2) = cslope(icc,2)
                          fdat(fn) = crse(icc,n)
                          icc = icc + ist
170                   continue
160               continue

                  do 180 joff = 0, lratio-1
                      j = lratio*jc + joff
                      if (j .lt. fblo(2)) then
                          goto 180
c                         --- next ---
                      end if
                      if (j .gt. fbhi(2)) then
                          goto 181
c                         --- break ---
                      end if
                      yoff = float(joff)-hafrat

                      do 190 i = fblo(1), fbhi(1)
                          fn = i-fslo(1)+1
                          fine(i,j,n) = fdat(fn) + volratio* (voff(fn)*f
     &slope(fn,1)+yoff*fslope(fn,2))
190                   continue
180               continue
181               continue
150           continue
          else
              do 200 ic = cblo(1), cbhi(1)

c ::: ..,.......::::: strip out a fine grid slope vector
                  do 210 joff = 1, lratio
                      icc = clo + ist*(ic-islo) + jst
                      istart = joff
                      iend = joff + (ncby-1)*lratio
                      do 220 fn = istart, iend, lratio
                          fslope(fn,1) = cslope(icc,1)
                          fslope(fn,2) = cslope(icc,2)
                          fdat(fn) = crse(icc,n)
                          icc = icc + jst
220                   continue
210               continue

                  do 230 ioff = 0, lratio-1
                      i = lratio*ic + ioff
                      if (i .lt. fblo(1)) then
                          goto 230
c                         --- next ---
                      end if
                      if (i .gt. fbhi(1)) then
                          goto 231
c                         --- break ---
                      end if
                      xoff = float(ioff)-hafrat

                      do 240 j = fblo(2), fbhi(2)
                          fn = j-fslo(2)+1
                          fine(i,j,n) = fdat(fn) + volratio* (xoff*fslop
     &e(fn,1)+voff(fn)*fslope(fn,2))
240                   continue
230               continue
231               continue
200           continue
          end if
120   continue

      return
      end

c ::: --------------------------------------------------------------
      subroutine pcinterp2d (fine,floi1,floi2,fhii1,fhii2,fblo,fbhi,lrat
     &,nvar,crse,cloi1,cloi2,chii1,chii2,cblo,cbhi,temp,tloi,thii)

      integer floi1,floi2
      integer fhii1,fhii2
      integer cloi1,cloi2
      integer chii1,chii2

      integer fblo(2), fbhi(2)
      integer cblo(2), cbhi(2)
      integer lrat, nvar, tloi, thii
      DOUBLE PRECISION fine(floi1 :fhii1 ,floi2 :fhii2, nvar)
      DOUBLE PRECISION crse(cloi1 :chii1 ,cloi2 :chii2, nvar)
      DOUBLE PRECISION temp(tloi:thii + 1)
c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ::: pcinterp:  use piecewise constant interpolation to define
c ::: values on the subregion of the fine FAB defined
c ::: by (fblo,fbhi).
c :::
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fab to get interpolated values
c ::: flo,fhi      =>  (const)  index limits of fine
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: crse         =>  (const)  fab holding coarse grid values
c ::: clo,chi      =>  (const)  index limits of src
c ::: cblo,cbhi    =>  (const)  subregion of coarse grid holding values
c ::: temp         =>  (modify) temporary space for vectorization
c ::: tlo,thi      =>  (const)  index limits of temp space
c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ::: local var
      integer i,j,k,ic,jc,kc,ioff,n
      integer ixproj,ii,ll
      ixproj(ii,ll) = (ii + ll*iabs(ii))/ll - iabs(ii)

      do 120 j = fblo(2), fbhi(2)
          jc = ixproj(j,lrat)
          do 130 n = 1, nvar
              do 140 ioff = 0, lrat-1
                  do 150 ic = cblo(1),cbhi(1)
                      i = lrat*ic + ioff
                      temp(i) = crse(ic,jc,n)
150               continue
140           continue
              do 160 i = fblo(1), fbhi(1)
                  fine(i,j,n) = temp(i)
160           continue
130       continue
120   continue

      return
      end

c ::: --------------------------------------------------------------
      subroutine cartgridminmax2d (data, lo1, lo2, hi1, hi2,vfracdata, v
     &feps, dmin, dmax)
      implicit none

      integer lo1, lo2, hi1, hi2
      DOUBLE PRECISION data(lo1:hi1 ,lo2:hi2)
      DOUBLE PRECISION vfracdata(lo1:hi1 ,lo2:hi2)
      DOUBLE PRECISION vfeps, dmin, dmax

      integer i, j

      dmax = -1.0D30
      dmin = 1.0D30
      do 420 j = lo2, hi2
          do 430 i = lo1, hi1
c      print *, "i j vfracdata(i,j) = ",i,j,vfracdata(i,j)
              if ( .not. (vfracdata(i,j).lt.vfeps)) then
                dmax = max(dmax,data(i,j))
                dmin = min(dmin,data(i,j))
              endif
430       continue
420   continue

      return
      end


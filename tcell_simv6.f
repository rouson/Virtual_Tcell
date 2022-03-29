      program tcsim
      implicit none

      logical actual_frames
      
      integer nside
      integer num_movies,ncellst
      integer i,j,k,jj,iii,jsave
      integer num_volume_levels
      integer nintervals_vol_v,jchoose
      integer nintervals_vol,inbox,maxnumframes
      integer num_tcells_allmovies,print_distributions
      integer num_volume_bin,num_total_frames,nframes_vol
      integer ii,ij,nit,total_frames,nminutes
      integer kk,numbin,numbinb,ncube,maxnf,ncube_reduce
      integer ijka,ivoltype,ncroot
      integer iig,jjg,kkg,ig,jg,kg,iangle,ki,kij,ijt
      integer sumframes_total,ntimes,idss
      integer id,nlines,ijtaa,maxnfa
      integer ivolt,irr,icell,itemp,num_left,nbinsize
      integer numf(100000)
      
      double precision lengthsim,maxrange,facab,time,time_end
      double precision max_speed_v,phi_new,pi,conv
      double precision xrefa,yrefa,zrefa,average_ang
      double precision volumetcells,xpos,ypos,zpos
      double precision vx,vy,vz,rmag,timee,vell
      double precision vxn,vyn,vzn,rmagn
      double precision phi_exact,vel_exact,width
      double precision velb,vele,sumcm
      double precision speed_max,speed_min
      double precision speed_sum,speed2_sum
      double precision xdisp,ydisp,zdisp,tdisp
      double precision disp,disp2
      double precision xsims,ysims,zsims,tsims
      double precision dispvprev,dispvafter
      double precision xdispnb,ydispnb,zdispnb
      double precision xdispna,ydispna,zdispna,dot,angle
      double precision minvolume,offset
      double precision xsimjj,ysimjj,zsimjj,tsimjj
      double precision rnumbin,rnumbinb,dspeed_sec
      double precision standard_dev
      double precision ang_min,ang_max,dang
      double precision vrx,vry,vrz,vnmag2
      double precision vnx,vny,vnz,rr,theta
      double precision xx,yy,zz
      double precision rrw,vel_new
      double precision angb,ange
      double precision dirx,diry,dirz,dist,dist2,dista
      double precision distance2,distsum2,dxg,dxg2
      double precision radius_tcell_squared,diameter_tcell
      double precision rr1,rr2,rr3
      double precision speed_avg,speedb,speede,timea
      double precision timedif,timen,timeobv
      double precision timesum2
      double precision vec1x,vec1y,vec1z
      double precision vec2x,vec2y,vec2z
      double precision velinitial,vg
      double precision vpx,vpy,vpz
      double precision xref,yref,zref
      double precision xc,yc,zc
      double precision xmin,ymin,zmin
      double precision xmax,ymax,zmax
      double precision xmina,ymina,zmina
      double precision xmaxa,ymaxa,zmaxa            
      double precision xminb,yminb,zminb
      double precision xmaxb,ymaxb,zmaxb      
      double precision xminc,yminc,zminc
      double precision xmaxc,ymaxc,zmaxc
      double precision vxp,vyp,vzp,rmagp      
      double precision rangex,rangey,rangez
      double precision time_step
      double precision dxs,dxs2,eps
      double precision speedaa,speedbb
      integer nside2,nside3

      integer nstpsum
      double precision timestpsum,timstp
      double precision speedfreqvoltot,angfreqvoltot
      double precision total_time,total_disp,vel_avg
      
      parameter (ncube = 168)
      
      logical grid(ncube,ncube,ncube)
c23456789012345678901234567890123456789012345678901234567890123456789012      
c     Argument is movie id      
c      integer, allocatable :: number_tcellsa(:),number_framesa(:)
c     Argument is total cell id
      integer, allocatable :: istore(:),ivoltypefull(:)
      integer, allocatable :: numframesa(:),numframesa_gen(:)
      integer, allocatable :: map_sim(:),numframes_v(:)

      integer, allocatable :: moved(:),iassign(:)
      integer, allocatable :: bincount(:),mark(:),mark_gen(:)
      integer, allocatable :: list(:,:),volume_type(:)
      integer, allocatable :: list_gen(:,:)
      integer, allocatable :: sorted(:),isum(:),sumframes(:)
      integer, allocatable :: indexsort(:),minframe(:),maxframe(:)
      integer, allocatable :: num_cells_level(:),numvf(:)
      integer, allocatable :: num_cells_level_gen(:)
      integer, allocatable :: num_volume_bin_gen(:)
c     23456789012345678901234567890123456789012345678901234567890123456789012      

      double precision, allocatable :: tcellmaxtime(:),tcellmintime(:)
      double precision, allocatable :: turnang(:),velangle(:)
      double precision, allocatable :: turnangc(:),velanglec(:)      
      double precision, allocatable :: xsim(:,:),ysim(:,:)
      double precision, allocatable :: zsim(:,:),tsim(:,:)
      double precision, allocatable :: xsimc(:,:),ysimc(:,:)
      double precision, allocatable :: zsimc(:,:),tsimc(:,:)
      double precision, allocatable :: xsimo(:,:),ysimo(:,:)
      double precision, allocatable :: zsimo(:,:),tsimo(:,:)      
c      double precision, allocatable :: dwidth(:,:)
      double precision, allocatable :: dspeed_seck(:)
      double precision, allocatable :: volumeall(:)
      double precision, allocatable :: rmark(:)
      double precision, allocatable :: aa(:),bb(:)
      double precision, allocatable :: speedbefore(:)
      double precision, allocatable :: speedafter(:)
      double precision, allocatable :: speedavg(:)
      double precision, allocatable :: angle_a(:)
      double precision, allocatable :: speedfreq_vol(:,:,:)
      double precision, allocatable :: angfreq_vol(:,:,:)
      double precision, allocatable :: speedbinfull_vol(:,:,:)
      double precision, allocatable :: cum_speedfreq_vol(:,:,:)
      double precision, allocatable :: angbinfull_vol(:,:,:)
      double precision, allocatable :: cum_angfreq_vol(:,:,:)
      double precision, allocatable :: speed_bin_after(:)
      double precision, allocatable :: ang_bin_after(:)
      double precision, allocatable :: volumetcellsim(:)
      double precision, allocatable :: frac_numframes_class(:)
      double precision, allocatable :: xp(:),yp(:),zp(:)
      double precision, allocatable :: xpp(:),ypp(:),zpp(:)      
      
      open(unit=97,file='flu_volume')
      open(unit=319,file='positions_new')
      open(unit=303,file='positions_sim')
      open(unit=220,file='angle_sim')
      open(unit=223,file='angle_empirical')      
      open(unit=197,file='volume_sim')
      open(unit=198,file='volume_empirical')      
      open(unit=876,file='volume_patrolled')

      actual_frames = .true.
      conv = 180.d0/acos(-1.d0)
      
c     Run simulation for nminutes
      nminutes = 90

      print_distributions = 1

      maxrange = 420.d0
      
c     Number of classes of T cells      
      num_volume_levels = 4
c      num_volume_levels = 1

c     Number of subgroups for different angle and velocity distributions      
      nintervals_vol_v = 1

c     Number of bins in distributions      
      nintervals_vol = 15

c     Set jchoose = 1 to print distributions      
      jchoose = 1

c      read(319,*) num_movies
c      allocate(number_tcellsa(num_movies))
c      allocate(number_framesa(num_movies))


      pi = acos(-1.d0)
      eps = 1.d-12
c     Diameter of a tcell in microns
      diameter_tcell = 10.d0
      radius_tcell_squared = (diameter_tcell/2.d0)**2


      nlines = 0
      DO
        READ (319,*, END=10)
        nlines = nlines + 1
      END DO
 10   REWIND (319)

      print*,'Number of lines in input positions file = ',nlines

      numf = 0
      maxnf = 0
      ncellst = 0
      do i = 1,nlines
         read(319,*) id,xsims,ysims,zsims,tsims
c         print*,id,xsims,ysims,zsims,tsims
         numf(id) = numf(id) + 1 
         ncellst = max(ncellst,id)
         maxnf = max(maxnf,numf(id))
      end do

      num_tcells_allmovies = ncellst
c      maxnfa = maxnf
      maxnfa = 2*maxnf
      print*,'Total number tcells = ',ncellst
      print*,'Maximum number of frames over all T cells = ',maxnf
      
      REWIND(319)

      allocate(numvf(num_volume_levels))
      allocate(volume_type(ncellst))
      allocate(xsim(ncellst,maxnfa),ysim(ncellst,maxnfa),
     &         zsim(ncellst,maxnfa),tsim(ncellst,maxnfa))
      allocate(xsimc(ncellst,maxnfa),ysimc(ncellst,maxnfa),
     &         zsimc(ncellst,maxnfa),tsimc(ncellst,maxnfa))
      allocate(xsimo(ncellst,maxnfa),ysimo(ncellst,maxnfa),
     &         zsimo(ncellst,maxnfa),tsimo(ncellst,maxnfa))      
      allocate(numframesa(ncellst))
      allocate(numframes_v(ncellst))      
      allocate(numframesa_gen(ncellst))      
      allocate(volumetcellsim(ncellst))
      allocate(volumeall(ncellst))
      allocate(tcellmaxtime(ncellst))
      allocate(tcellmintime(ncellst))

      numframesa = 0
      do i = 1,nlines
         read(319,*) id,xsims,ysims,zsims,tsims         
         numframesa(id) = numframesa(id) + 1
         xsim(id,numframesa(id)) = xsims
         ysim(id,numframesa(id)) = ysims
         zsim(id,numframesa(id)) = zsims
         tsim(id,numframesa(id)) = tsims
      end do
      
      xmin = 1.d+20
      ymin = 1.d+20
      zmin = 1.d+20
      maxnumframes = -1
      do i = 1,num_tcells_allmovies
         maxnumframes = max(maxnumframes,numframesa(i))
         do j = 1,numframesa(i)
            xmin = min(xsim(i,j),xmin)
            ymin = min(ysim(i,j),ymin)
            zmin = min(zsim(i,j),zmin)
         end do
      end do
      print*,'maxnumframes maxnfa = ',
     &        maxnumframes,maxnfa


      xminb = 1.d+20
      yminb = 1.d+20
      zminb = 1.d+20      
      xmaxb = -1.d+20
      ymaxb = -1.d+20
      zmaxb = -1.d+20
      do i = 1,num_tcells_allmovies
         do j = 1,numframesa(i)
            xsim(i,j) = xsim(i,j) - xmin
            ysim(i,j) = ysim(i,j) - ymin
            zsim(i,j) = zsim(i,j) - zmin

            xminb = min(xsim(i,j),xminb)
            yminb = min(ysim(i,j),yminb)
            zminb = min(zsim(i,j),zminb)                        
            xmaxb = max(xsim(i,j),xmaxb)
            ymaxb = max(ysim(i,j),ymaxb)
            zmaxb = max(zsim(i,j),zmaxb)
            
         end do
      end do
      rangex = xmaxb - xminb
      rangey = ymaxb - yminb
      rangez = zmaxb - zminb

      print*,'Minimum and Maximum position values of T cells'
      print*,'Minimum values: xmin ymin zmin = ',xminb,yminb,zminb
      print*,'Maximum values: xmax ymax zmax = ',xmaxb,ymaxb,zmaxb
      print*,'Range: rangex rangey rangez = ',
     &       rangex,rangey,rangez

      if (rangex .ge. maxrange .or.
     &    rangey .ge. maxrange .or.
     &    rangez .ge. maxrange) then
         print*,'rangex rangey rangez = ',
     &        rangex,rangey,rangez
         print*,'maxrange = ',maxrange
         stop
      end if


c     Calculate the volume traversed by each T cell
      print*,'Calculate volume patrolled by empirical cells'
c     We assume our grid ranges between
c     (-10:410,-10:410,-10:410)

      dxg = maxrange/dble(ncube)
c     dxg should be 2.5 microns
      dxg2 = dxg/2.d0
c     Volume of 2.5um x 2.5um x 2.5um cube
      vg = dxg**3

      xref = -10.d0
      yref = -10.d0
      zref = -10.d0

      do i = 1,num_tcells_allmovies
         grid = .false.
         tcellmaxtime(i) = -1.d0
         tcellmintime(i) = 1.d+10         
         do j = 1,numframesa(i)

c           (iig,jjg,kkg) are the logical coordinates of the T cell

            if (tsim(i,j) .gt. -1.d-6) then
               tcellmaxtime(i) = max(tcellmaxtime(i),tsim(i,j))
               tcellmintime(i) = min(tcellmintime(i),tsim(i,j))
               iig = (xsim(i,j) - xref)/dxg
               jjg = (ysim(i,j) - yref)/dxg
               kkg = (zsim(i,j) - zref)/dxg

               do ig = iig-2,iig+2
                  do jg = jjg-2,jjg+2
                     do kg = kkg-2,kkg+2
                        xc = dxg*dble(ig) + xref + dxg2
                        yc = dxg*dble(jg) + yref + dxg2
                        zc = dxg*dble(kg) + zref + dxg2
                        distance2 = (xsim(i,j) - xc)**2 +
     &                              (ysim(i,j) - yc)**2 +
     &                              (zsim(i,j) - zc)**2
                        if (ig .gt. ncube .or.
     &                      jg .gt. ncube .or.
     &                      kg .gt. ncube .or.
     &                      ig .lt. 0 .or.
     &                      jg .lt. 0 .or.
     &                      kg .lt. 0) then                           
                           print*,'ig jg kg = ',ig,jg,kg
                           stop
                        end if
                        if (distance2 .lt. 
     &                      radius_tcell_squared) then
                           grid(ig,jg,kg) = .true.
                        end if
                     end do
                  end do
               end do
            else
               print*,'tsim(',i,j,') = ',tsim(i,j)
               stop
            end if

         end do
                  
         volumetcellsim(i) = 0.d0
         do kkg = 1,ncube
            do jjg = 1,ncube
               do iig = 1,ncube
                  if (grid(iig,jjg,kkg)) then
                     volumetcellsim(i) = 
     &               volumetcellsim(i) + vg
                  end if
               end do
            end do
         end do

c         timeobv = tsim(i,numframesa(i))-tsim(i,1)
         timeobv = tcellmaxtime(i) - tcellmintime(i)         
         if (timeobv .gt. 0.d0) then
            volumetcellsim(i) = volumetcellsim(i)/timeobv
            volumeall(i) = volumetcellsim(i)
            write(198,*) volumetcellsim(i)
         end if

      end do

      print*,'Finished calculating volume patrolled by empirical cells'

      
      xsimc = xsim
      ysimc = ysim
      zsimc = zsim
      tsimc = tsim

      xsimo = xsim
      ysimo = ysim
      zsimo = zsim
      tsimo = tsim
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccVolume treatment ccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      allocate(minframe(num_volume_levels))
      allocate(maxframe(num_volume_levels))
      allocate(num_volume_bin_gen(num_volume_levels))
      allocate(isum(num_volume_levels))
      allocate(sumframes(num_volume_levels))      
      allocate(frac_numframes_class(num_volume_levels))
      allocate(num_cells_level(num_volume_levels))
      allocate(num_cells_level_gen(num_volume_levels))      
      allocate(bincount(num_volume_levels))
      allocate(mark(num_volume_levels+1))
      allocate(mark_gen(num_volume_levels+1))      
      allocate(rmark(num_volume_levels+1))      
      allocate(list(num_volume_levels,num_tcells_allmovies))
      allocate(list_gen(num_volume_levels,num_tcells_allmovies))       
      allocate(sorted(num_tcells_allmovies))  
      allocate(indexsort(num_tcells_allmovies))  

      allocate(aa(nintervals_vol),bb(nintervals_vol))
      aa = 0.d0
      bb = 0.d0
c     23456789012345678901234567890123456789012345678901234567890123456789012
      allocate(dspeed_seck(num_volume_levels))
c      allocate(dwidth(num_volume_levels,nintervals_vol_v))
      allocate(speedfreq_vol(num_volume_levels,
     &                       nintervals_vol_v,nintervals_vol))
      allocate(angfreq_vol(num_volume_levels,
     &                     nintervals_vol_v,nintervals_vol))
      allocate(speedbinfull_vol(num_volume_levels,
     &         nintervals_vol_v,nintervals_vol))
      allocate(cum_speedfreq_vol(num_volume_levels,
     &         nintervals_vol_v,nintervals_vol+1))         
      allocate(angbinfull_vol(num_volume_levels,
     &         nintervals_vol_v,nintervals_vol))
      allocate(cum_angfreq_vol(num_volume_levels,
     &         nintervals_vol_v,nintervals_vol+1))

c23456789012345678901234567890123456789012345678901234567890123456789012               

      print*,'Sorting volumes'      
      do i = 1,num_tcells_allmovies
         sorted(i) = 0
      end do

c     Sort volume patrolled by each T cell from least to greatest
      do i = 1,num_tcells_allmovies
         minvolume = 1.d+10
         do j = 1,num_tcells_allmovies
            if (volumeall(j) .lt. minvolume .and.
     &          sorted(j) .eq. 0) then
               minvolume = volumeall(j)
               jsave = j
            end if
         end do
         indexsort(i) = jsave
         sorted(jsave) = 1
      end do

c     Divide number of T cells equally into bins      
      num_volume_bin = int(dble(num_tcells_allmovies)/
     &                     dble(num_volume_levels))

      num_total_frames = 0
      do j = 1,num_tcells_allmovies
         do k = 1,numframesa(j)
            num_total_frames = num_total_frames + 1
         end do
      end do

c      print*,'num_total_frames = ',num_total_frames

      allocate(speed_bin_after(num_total_frames))
      allocate(ang_bin_after(num_total_frames))
      allocate(turnang(num_total_frames))
      allocate(velangle(num_total_frames))
      allocate(turnangc(num_total_frames))
      allocate(velanglec(num_total_frames))      
      allocate(speedbefore(num_total_frames))
      allocate(speedafter(num_total_frames))
      allocate(speedavg(num_total_frames))      
      allocate(angle_a(num_total_frames))
      
c     Construct lists of cells in each volume class
      do ii = 1,num_volume_levels+1
         mark(ii) = num_volume_bin*(ii-1) + 1
      end do
      mark(num_volume_levels+1) = num_tcells_allmovies+1

      do i = 1,num_volume_levels
         ij = 0
         do j = mark(i),mark(i+1)-1
            ij = ij + 1
            list(i,ij) = indexsort(j)
         end do
         num_cells_level(i) = ij
      end do

c     Units are 20 um/min
      max_speed_v = 20.d0
      
      do k = 1,num_volume_levels
         print*,'Creating bins of cells for volume per time class ',k
c        Create lists of speedbefore, speedafter, angle_a in class k
         speed_max = -1.d+20
         speed_min = 1.d+20
         speed_sum = 0.d0
         speed2_sum = 0.d0
         bincount(k) = 0
         minframe(k) = 1.d+8
         maxframe(k) = -1.d+8
         do ii = 1,num_cells_level(k)
            i = list(k,ii)
            minframe(k) = min(minframe(k),numframesa(i))
            maxframe(k) = max(maxframe(k),numframesa(i))
            do j = 2,numframesa(i)-1
               if (tsimc(i,j-1) .gt. -1.d-6 .and.
     &             tsimc(i,j) .gt. -1.d-6 .and.
     &             tsimc(i,j+1) .gt. -1.d-6) then

                  xdisp = xsimc(i,j) - xsimc(i,j-1)
                  ydisp = ysimc(i,j) - ysimc(i,j-1)
                  zdisp = zsimc(i,j) - zsimc(i,j-1)
                  tdisp = tsimc(i,j) - tsimc(i,j-1)
                  disp2 = xdisp**2 + ydisp**2 + zdisp**2 
                  disp = sqrt(disp2)+1.d-12
c                 Displacement velocity in um/s
                  dispvprev = disp/tdisp
                  total_disp = disp
                  total_time = tdisp

c                 Normalized vector directions from j-1 to j                  
                  xdispnb = xdisp/disp
                  ydispnb = ydisp/disp
                  zdispnb = zdisp/disp

                  xdisp = xsimc(i,j+1) - xsimc(i,j)
                  ydisp = ysimc(i,j+1) - ysimc(i,j)
                  zdisp = zsimc(i,j+1) - zsimc(i,j)
                  tdisp = tsimc(i,j+1) - tsimc(i,j)
                  disp2 = xdisp**2 + ydisp**2 + zdisp**2
                  disp = sqrt(disp2)+1.d-12
                  total_disp = total_disp + disp
                  total_time = total_time + tdisp
                  dispvafter = disp/tdisp

c                 Normalized vector directions from j to j+1                                    
                  xdispna = xdisp/disp
                  ydispna = ydisp/disp
                  zdispna = zdisp/disp

                  dot = xdispnb*xdispna + 
     &                  ydispnb*ydispna +
     &                  zdispnb*zdispna

                  angle = acos(dot)

                  bincount(k) = bincount(k) + 1
c                 Speeds are in um/s
                  speedavg(bincount(k)) = total_disp/(total_time+1.d-12)
                  speedbefore(bincount(k)) = dispvprev
                  speedafter(bincount(k)) = dispvafter
                  speed_max = max(speed_max,dispvprev)
                  speed_min = min(speed_min,dispvprev)
                  speed_sum = speed_sum + dispvprev
                  speed2_sum = speed2_sum + dispvprev**2
c                 Angles are in radians
                  angle_a(bincount(k)) = angle

               end if
            end do
         end do

         rnumbin = 1.d0/dble(bincount(k))

c        Dimensions of dspeed_sec are um/s      
c         dspeed_sec = (max_speed_v/60.d0 - 0.d0/60.d0)/
c     &        dble(nintervals_vol_v)


         speed_avg = speed_sum*rnumbin
         speed_min = 0.d0
         standard_dev = speed2_sum + rnumbin*(speed_sum**2)
         standard_dev = dsqrt(standard_dev/dble(bincount(k)-1))
         speed_max = min(speed_avg + 2.d0*standard_dev,40.d0/60.d0)
c        changed
         speed_max = max_speed_v/60.d0
c        Units of speed_min and speed_max are um/second         

         dspeed_sec = (speed_max - speed_min)/dble(nintervals_vol_v)

         dspeed_seck(k) = dspeed_sec

c        Create distribution in each class k         
         do j = 1,nintervals_vol_v

c           Creating velocity bin
            velb = speed_min + dble(j-1)*dspeed_sec
            vele = speed_min + dble(j)*dspeed_sec

            kk = 0
            do i = 1,bincount(k)
               if (speedbefore(i) .ge. velb .and.
     &             speedbefore(i) .lt. vele) then
                  kk = kk + 1
c                 units of um/min
                  speed_bin_after(kk) = speedafter(i)
c                  ang_bin_after(kk) = angle_a(i)
               end if
            end do
             
c           Construct a distribution with speeds in bin
c           numbin: number of items in bin j
            numbin = kk
            rnumbin = 1.d0/dble(numbin)

            kk = 0
            do i = 1,bincount(k)
c               if (speedafter(i) .ge. velb .and.
c     &              speedafter(i) .lt. vele) then
               if (speedavg(i) .ge. velb .and.
     &             speedavg(i) .lt. vele) then                  
                  kk = kk + 1
c                 units of um/min
c                  speed_bin_after(kk) = speedafter(i)
                  ang_bin_after(kk) = angle_a(i)
               end if
            end do
             
c           Construct a distribution with speeds in bin
c           numbin: number of items in bin j
            numbinb = kk
            rnumbinb = 1.d0/dble(numbinb)            
            
            speed_min = 1.d+20
            speed_max = -1.d+20
            speed_sum = 0.d0
            speed2_sum = 0.d0
            do i = 1,numbin
               speed_sum = speed_sum + speed_bin_after(i)
               speed2_sum = speed2_sum + speed_bin_after(i)**2
               speed_min = min(speed_min,speed_bin_after(i))
               speed_max = max(speed_max,speed_bin_after(i))
            end do
            standard_dev = speed2_sum + rnumbin*(speed_sum**2)
            standard_dev = dsqrt(standard_dev/dble(numbin-1))
            speed_avg = speed_sum*rnumbin

            speed_min = 0.d0
            speed_max = min(speed_avg + 2.d0*standard_dev,40.d0/60.d0)
changed            
            speed_max = max_speed_v/60.d0
            
c            dwidth(k,j) = (speed_max - speed_min)/dble(nintervals_vol)
            width = (speed_max - speed_min)/dble(nintervals_vol)

            speedfreqvoltot = 0.d0
            do i = 1,nintervals_vol
c               speedb = speed_min + dble(i-1)*dwidth(k,j)
c               speede = speed_min + dble(i)*dwidth(k,j)

               speedb = speed_min + dble(i-1)*width
               speede = speed_min + dble(i)*width
               
c              units of speedbinfullavg is um/min 
               speedbinfull_vol(k,j,i) = (speedb+speede)/2.d0
               speedfreq_vol(k,j,i) = 0.d0
               do ki = 1,numbin
                  if (speed_bin_after(ki) .ge. speedb .and.
     &                speed_bin_after(ki) .lt. speede) then
                     speedfreq_vol(k,j,i) =
     &               speedfreq_vol(k,j,i) + 1.d0
                     speedfreqvoltot = speedfreqvoltot + 1.d0
                  end if
               end do
            end do
            do i = 1,nintervals_vol
               speedfreq_vol(k,j,i) =
     &         speedfreq_vol(k,j,i)/speedfreqvoltot
            end do


            ang_min = 1.d+20
            ang_max = -1.d+20
            do i = 1,numbinb
               ang_min = min(ang_min,ang_bin_after(i))
               ang_max = max(ang_max,ang_bin_after(i))
            end do
            ang_min = 0.d0
            ang_max = acos(-1.d0)

            dang = (ang_max - ang_min)/dble(nintervals_vol)

            angfreqvoltot = 0.d0
            do i = 1,nintervals_vol
               angb = ang_min + dble(i-1)*dang
               ange = ang_min + dble(i)*dang
               angbinfull_vol(k,j,i) = (angb+ange)/2.d0
               angfreq_vol(k,j,i) = 0.d0
               do ki = 1,numbinb
                  if (ang_bin_after(ki) .ge. angb .and.
     &                ang_bin_after(ki) .lt. ange) then
                     angfreq_vol(k,j,i) = angfreq_vol(k,j,i) + 1.d0
                     angfreqvoltot = angfreqvoltot + 1.d0
                  end if
               end do
            end do
            do i = 1,nintervals_vol
               angfreq_vol(k,j,i) = angfreq_vol(k,j,i)/angfreqvoltot
            end do

         end do


c     k loop         
      end do


      print*,'Creating cumulative distributions'
c     Create cumulative distributions      
      do k = 1,num_volume_levels
                        
         do j = 1,nintervals_vol_v

            do i = 1,nintervals_vol
               aa(i) = aa(i) + speedbinfull_vol(k,j,i)*60.d0
               bb(i) = bb(i) + speedfreq_vol(k,j,i)

c               print*,'j jchoose nintervals_vol_v = ',
c     &              j,jchoose,nintervals_vol_v
c               stop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               
c              Comment if one does not wish to print distributions
c              if (k .eq. 1 .and. print_distributions .eq. 1) then
               if (j .eq. jchoose) then
               write(600+k,*) speedbinfull_vol(k,j,i)*60.d0,
     &                 speedfreq_vol(k,j,i)
               end if
c               end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
            end do
            
c           Dimensions of cumulative distribution are um/s
            cum_speedfreq_vol(k,j,1) = 0.d0
            sumcm = 0.d0
            do i = 1,nintervals_vol
               cum_speedfreq_vol(k,j,i+1) = cum_speedfreq_vol(k,j,i) +
     &              speedfreq_vol(k,j,i)
               sumcm = sumcm + speedfreq_vol(k,j,i)
            end do
            print*,'k = ',k,sumcm

c           Dimensions of culumlative distribution are radians
            cum_angfreq_vol(k,j,1) = 0.d0

            average_ang = 0.0
            do i = 1,nintervals_vol


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               
c              Comment if one does not wish to print distributions               
c               if (k .eq. 1 .and. print_distributions .eq. 1) then
                  write(700+k,*) angbinfull_vol(k,j,i)*conv,
     &                    angfreq_vol(k,j,i)
c               end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               
               average_ang = average_ang +
     &         angbinfull_vol(k,j,i)*conv*angfreq_vol(k,j,i)
               
               cum_angfreq_vol(k,j,i+1) = cum_angfreq_vol(k,j,i) +
     &         angfreq_vol(k,j,i)
            end do
            write(872,*) average_ang,60.*dspeed_seck(k)*dble(j)
            
         end do
         
      end do


c     Distance between cell centers
      dxs = 40.d0
      dxs2 = dxs/2.d0
      nside = 30
c     Large box is of size  [0:1200,0:1200,0:1200] microns
      nside2 = nside**2
      nside3 = nside**3
      
c     Positions of virtual cells

      allocate(map_sim(0:nside-1))
      allocate(xp(nside3),yp(nside3),zp(nside3))
      allocate(xpp(nside3),ypp(nside3),zpp(nside3))      
      
      do i = 1,nside-1
         map_sim(i) = i
      end do
      map_sim(0) = nside

      nbinsize = dble(nside3)/dble(num_volume_levels)

      allocate(istore(nside3),iassign(nside3))
      allocate(ivoltypefull(nside3))

      print*,'Assign volume classes to cells'
      num_left = nside3
      do i = 1,nside3
         istore(i) = i
         iassign(i) = 0
      end do

      offset = 400.d0
c     Adjust T cells to reside in box of size  [-400:800,-400:800,-400:800]
      do i = 1,nside3
         ii = map_sim(mod(i,nside))
         jj = map_sim(mod((i-1)/nside + 1,nside))
         kk = map_sim(mod((i-1)/nside2 + 1,nside))

         xpp(i) = dble(ii)*dxs - dxs2 - offset
         ypp(i) = dble(jj)*dxs - dxs2 - offset
         zpp(i) = dble(kk)*dxs - dxs2 - offset
         
         xp(i) = dble(ii)*dxs - dxs2  - offset
         yp(i) = dble(jj)*dxs - dxs2  - offset
         zp(i) = dble(kk)*dxs - dxs2  - offset
         
         iii = (dble(i)-1.d-8)/dble(nbinsize)
         ivolt = min(iii+1,num_volume_levels)
         
         call random_number(rr)
         irr = rr*dble(num_left)
         irr = irr + 1
         icell = istore(irr)
         
         itemp = istore(num_left)
         istore(num_left) = istore(irr)
         istore(irr) = itemp
         num_left = num_left - 1

         iassign(icell) = 1

c        Assign each T cell a volume class         
         ivoltypefull(icell) = ivolt
         
      end do

      do i = 1,nside3
         if (iassign(i) .eq. 0) then
            print*,'iassign(',i,') = ',iassign(i)
         end if
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
ccccccccccccccccccccccccc Begin virtual cell simulations cccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      print*,'Performing virtual cell simulations using'
      print*,'initial empirical positions'
      print*,'Creating speed, angle, and volume'
      print*,'distributions for comparison'
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
ccccccccccccccccccccccccc Move virtual cells ccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ijka = 0

      allocate(moved(num_tcells_allmovies))
      do j = 1,num_tcells_allmovies
         moved(j) = 0
      end do

      isum = 0
      sumframes = 0
      sumframes_total = 0
      
      timestpsum = 0.d0
      nstpsum = 0
      do kij = 1,num_volume_levels
         ivoltype = kij

         do ii = 1,num_cells_level(kij)
            j = list(kij,ii)
            moved(j) = 1

            do k = 2,numframesa(j)
               
               if (tsimo(j,k) .gt. -1.d-6) then

                  ijka = ijka + 1
                  do iii = 1,num_volume_levels
                     if (ivoltype .eq. iii) then
                        isum(iii) = isum(iii) + 1
                     end if
                  end do
                  
               end if

               if (tsimo(j,k) .gt. -1.d-6 .and.
     &             tsimo(j,k-1) .gt. -1.d-6) then
                  timstp = tsimo(j,k)-tsimo(j,k-1)
                  timestpsum = timestpsum + timstp
                  nstpsum = nstpsum + 1
               end if
               
            end do

         end do

      end do

      timstp = timestpsum/dble(nstpsum)
      print*,'timstp = ',timstp
      
      do iii = 1,num_volume_levels
c         print*,'isum(',iii,') = ',isum(iii),
c     &   isum(iii)/num_cells_level(iii)
         numvf(iii) = isum(iii)/num_cells_level(iii)
         print*,'numvf(',iii,') = ',numvf(iii)
      end do
      
      do kij = 1,num_volume_levels
         do ii = 1,num_cells_level(kij)
            j = list(kij,ii)
            do k = 2,numvf(kij)
               tsimo(j,1) = 0.d0
               tsimo(j,k) = timstp*dble(k-1)
            end do
         end do
      end do
      
      print*,'Moving virtual cells using exact number of cells'
      print*,'and the exact number of frames for each cell'
      do kij = 1,num_volume_levels
         ivoltype = kij

         do ii = 1,num_cells_level(kij)
            j = list(kij,ii)
c            print*,'list(',kij,ii,') = ',list(kij,ii)
            moved(j) = 1

            volume_type(j) = kij
c23456789012345678901234567890123456789012345678901234567890123456789012                              
            sumframes(kij) = sumframes(kij) + numframesa(j)
            sumframes_total = sumframes_total + numframesa(j)
            numframes_v(j) = numvf(kij)
c            write(303,*) j,xsimo(j,1),ysimo(j,1),zsimo(j,1),
c     &      idss,tsimo(j,1)
            do k = 2,numvf(kij)
c            do k = 2,numframesa(j)
c            do k = 2,2*numframesa(j)
c               print*,'tsimo(',j,k,') = ',tsimo(j,k)
               
               if (tsimo(j,k) .gt. -1.d-6) then
c     &             .or. k .gt. numframesa(j)) then

                  ijka = ijka + 1
                  do iii = 1,num_volume_levels
                     if (ivoltype .eq. iii) then
                        isum(iii) = isum(iii) + 1
                     end if
                  end do                  
                  
                  jj = 0
                  if (k .eq. 2) then
                     call random_number(rr2)                  
                     jj=int(max(rr2*dble(nintervals_vol_v)-eps,0.d0))+1
                  else
                     timee = tsimo(j,k-1) - tsimo(j,k-2)
                     if (tsimo(j,k-2) .lt. -1.d-6 .or.
     &                   tsimo(j,k-1) .lt. -1.d-6) then
                        call random_number(rr2)                  
                        jj =
     &                  int(max(rr2*dble(nintervals_vol_v)-eps,0.d0))+1
                     end if
                  end if

                  if (jj .eq. 0) then
c                    Compute previous velocity

                     call random_number(rr1)
                     call random_number(rr2)
                     call random_number(rr3)
                     
                     vx = xsimo(j,k-1) - xsimo(j,k-2) + eps*rr1
                     vy = ysimo(j,k-1) - ysimo(j,k-2) + eps*rr2
                     vz = zsimo(j,k-1) - zsimo(j,k-2) + eps*rr3
 
                     rmag = dsqrt(vx**2 + vy**2 + vz**2) + 1.d-12
                     if (rmag .lt. 1.d-8) then
                        print*,'Cell did not move first'
                        print*,'jj = ',jj
                        print*,'vx vy vz = ',vx,vy,vz
                        print*,'j k = ',j,k-1,k-2
                        print*,xsimo(j,k-1),xsimo(j,k-2)
                        print*,ysimo(j,k-1),ysimo(j,k-2)
                        print*,zsimo(j,k-1),zsimo(j,k-2)                        
c                       stop
                     end if
                     if (timee .lt. 1.d-8) then
                        print*,'timee = ',timee
c                        stop
                     end if
                     vell = rmag/timee
                     rmag = 1.d0/rmag
c                    Unit vector (vx,vy,vz)
                     vx = vx*rmag
                     vy = vy*rmag
                     vz = vz*rmag
                  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                    Use exact positions xsimc,ysimc,zsimc
c                     vxp = xsimc(j,k-1) - xsimc(j,k-2)
c                     vyp = ysimc(j,k-1) - ysimc(j,k-2)
c                     vzp = zsimc(j,k-1) - zsimc(j,k-2)
c                     rmagp = dsqrt(vxp**2 + vyp**2 + vzp**2) + 1.d-12
c                     rmagp = p1.d0/rmagp
c                     vxp = vxp*rmagp
c                     vyp = vyp*rmagp
c                     vzp = vzp*rmagp
c                     
c                     vxn = xsimc(j,k) - xsimc(j,k-1)
c                     vyn = ysimc(j,k) - ysimc(j,k-1)
c                     vzn = zsimc(j,k) - zsimc(j,k-1)
c                     rmagn = dsqrt(vxn**2 + vyn**2 + vzn**2) + 1.d-12;
c                     vxn = vxn/rmagn
c                     vyn = vyn/rmagn
c                     vzn = vzn/rmagn
c                        
c                     dot = vxp*vxn + vyp*vyn + vzp*vzn
c
cc                    Exact angle calculated from empirical data
c                     phi_exact = acos(dot)
c                     timen = tsimc(j,k) - tsimc(j,k-1)
c
cc                    Exact speed calculated from empirical data
c                     vel_exact = rmagn/timen
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               

c                    vell is the previous velocity            
                     jj=
     &               min(int(vell/dspeed_seck(kij))+1,nintervals_vol_v)

c                     print*,'jj = ',jj,k,dspeed_seck(kij)
c                     stop
                     
                  else
                     
                     call random_number(rr1)
                     call random_number(rr2)
                     call random_number(rr3)
                     vx = rr1
                     vy = rr2
                     vz = rr3

                     rmag = dsqrt(vx**2 + vy**2 + vz**2) + 1.d-12
                     rmag = 1.d0/rmag
                     vx = vx*rmag
                     vy = vy*rmag
                     vz = vz*rmag
                  end if

                  call random_numbere(rr)
                  width = speedbinfull_vol(ivoltype,jj,2) -
     &                    speedbinfull_vol(ivoltype,jj,1)

                  do i = 1,nintervals_vol

c23456789012345678901234567890123456789012345678901234567890123456789012                              
                     if (rr .ge. cum_speedfreq_vol(ivoltype,jj,i) .and. 
     &                   rr.lt. cum_speedfreq_vol(ivoltype,jj,i+1)) then
                  
cc                      speedbinfullavg has units of um/min
            
                        call random_numbere(rrw)
                        rrw = rrw - 0.5d0
                        vel_new =
     &                       speedbinfull_vol(ivoltype,jj,i)+width*rrw
                     end if
                  end do

                  vel_avg = (vell + vel_new)/2.d0

c                  if (jj .eq. jchoose) then
c                     write(420+kij,*) 0.d0,vel_new*60.d0
c                  end if

c                  jj =
c     &            min(int(vel_new/dspeed_seck(kij))+1,nintervals_vol_v)

                  jj =
     &            min(int(vel_avg/dspeed_seck(kij))+1,nintervals_vol_v)                  

                  
                  width = angbinfull_vol(ivoltype,jj,2) -
     &                    angbinfull_vol(ivoltype,jj,1)            
                  
                  call random_numbere(rr)
                  do i = 1,nintervals_vol               

                     if (rr .ge. cum_angfreq_vol(ivoltype,jj,i) .and. 
     &                   rr .lt. cum_angfreq_vol(ivoltype,jj,i+1)) then               
               
c                       phi is in degrees
                        call random_numbere(rrw)
                        rrw = rrw - 0.5d0
                  
                        phi_new = angbinfull_vol(ivoltype,jj,i) +
     &                  rrw*width 
                     end if
                  end do

                  nit  = 0
                  vnmag2 = 0.0
                  do while (vnmag2 .lt. 1.d-2 .and. nit .le. 100)
                     nit = nit + 1
                     call random_numbere(rr1)
                     call random_numbere(rr2)
                     call random_numbere(rr3)
                     vrx = rr1
                     vry = rr2
                     vrz = rr3
                     rmag = dsqrt(vrx**2 + vry**2 + vrz**2) + 1.d-12
                     rmag = 1.d0/rmag
                     vrx = vrx*rmag
                     vry = vry*rmag
                     vrz = vrz*rmag
                     call cross(vx,vy,vz,vrx,vry,vrz,vnx,vny,vnz)
                     vnmag2 = vnx**2 + vny**2 + vnz**2
                  end do
                  if (nit .gt. 100) then
                     print*,'Too many iterations in coordinate system'
                     stop
                  end if

                  rmagn = 1.0/dsqrt(vnmag2)
                  vnx = vnx*rmagn
                  vny = vny*rmagn
                  vnz = vnz*rmagn

                  
                  call cross(vx,vy,vz,vnx,vny,vnz,vpx,vpy,vpz)

                  call random_numbere(rr)
                  theta = rr*2.d0*acos(-1.d0)

c                 These coordinates are in the frame of the previous velocity vector
c                 (vx,vy,vz,vnx,vny,vnz,vpx,vpy,vpz)
                  zz = cos(phi_new)
                  xx = sin(phi_new)*cos(theta)
                  yy = sin(phi_new)*sin(theta)

c                 Map the coordinates to regular cartesian coordinates               
                  dirx = xx*vnx + yy*vpx + zz*vx
                  diry = xx*vny + yy*vpy + zz*vy
                  dirz = xx*vnz + yy*vpz + zz*vz
c                  print*,'dirx diry dirz = ',dirx,diry,dirz
                  
                  lengthsim = vel_new*(tsimo(j,k) - tsimo(j,k-1))
                  lengthsim = max(lengthsim,1.d-10);

                  if (lengthsim .lt. 1.d-14) then
                     print*,'lengthsim = ',lengthsim
                     print*,'vel_new = ',vel_new
                     print*,'t diff = ',
     &               tsimo(j,k),tsimo(j,k-1)
                     print*,'j k = ',j,k
c                     stop
                  end if

       
                  xsimo(j,k) = xsimo(j,k-1) + dirx*lengthsim
                  ysimo(j,k) = ysimo(j,k-1) + diry*lengthsim
                  zsimo(j,k) = zsimo(j,k-1) + dirz*lengthsim
c                  write(303,*) j,
c     &            xsimo(j,k),ysimo(j,k),zsimo(j,k),idss,tsimo(j,k)
                  
               end if

            end do

         end do
c         stop
      end do

      
      do j = 1,num_tcells_allmovies         
         if (moved(j) .eq. 0) then
            print*,'j not moved = ',j,moved(1)
            stop
         end if
      end do

c     Analyze velocities and turning angles of empirical cells
      iangle = 0
      do i = 1,num_tcells_allmovies
         do j = 1,numframesa(i)-2

            if (tsimc(i,j) .gt. -1.d-6 .and.
     &          tsimc(i,j+1) .gt. -1.d-6 .and.
     &          tsimc(i,j+2) .gt. -1.d-6) then

               vec1x = xsimc(i,j+1) - xsimc(i,j)
               vec1y = ysimc(i,j+1) - ysimc(i,j)
               vec1z = zsimc(i,j+1) - zsimc(i,j)
               dist2 = vec1x**2 + vec1y**2 + vec1z**2 
               dist = sqrt(dist2)+1.d-16
               dista = dist
               timedif = tsimc(i,j+1) - tsimc(i,j)
               if (timedif .lt. 0.d0) then
                  print*,'timedif = ',timedif
                  print*,tsimc(i,j+1),tsimc(i,j)
                  print*,i,j
                  stop
               end if
               timea = timedif
               velinitial = dist/timedif
               vec1x = vec1x/dist
               vec1y = vec1y/dist
               vec1z = vec1z/dist

               vec2x = xsimc(i,j+2) - xsimc(i,j+1)
               vec2y = ysimc(i,j+2) - ysimc(i,j+1)
               vec2z = zsimc(i,j+2) - zsimc(i,j+1)
               dist2 = vec2x**2 + vec2y**2 + vec2z**2 
               dist = sqrt(dist2)+1.d-16
               timedif = tsimc(i,j+2)-tsimc(i,j+1)
               timesum2 = timea + timedif
               distsum2 = dista + dist
               vec2x = vec2x/dist
               vec2y = vec2y/dist
               vec2z = vec2z/dist
                  
               dot = vec1x*vec2x+vec1y*vec2y+vec1z*vec2z
               angle = acos(dot)*180.d0/pi
               
               if (dist .gt. 1.d-12 .and.
     &             dista .gt. 1.d-12) then
                  iangle = iangle + 1
                  turnangc(iangle) = angle
                  velanglec(iangle) = distsum2/timesum2
                  velanglec(iangle) = dista/timea
                  velanglec(iangle) = dist/timedif                  
                  write(223,*) turnangc(iangle),velanglec(iangle)*60.
               end if

            end if
            
         end do
      end do

c      numframesa = numframesa_gen
c      do i = 1,num_tcells_allmovies
c         do j = 1,numframesa(i)
c            tsimc(i,j) = time_step*dble(j)
c         end do
c     end do

      if (actual_frames) then
         xsim = xsimo
         ysim = ysimo
         zsim = zsimo
         tsim = tsimo
      end if

c     Analyze velocities and turning angles of virtual cells      
      iangle = 0
      do i = 1,num_tcells_allmovies
c     do j = 1,numframesa(i)-2
        do j = 1,numframes_v(i)
           write(303,*) i,
     &     xsim(i,j),ysim(i,j),zsim(i,j),j,tsim(i,j)
        end do
         
         do j = 1,numframes_v(i)-2            
            
            vec1x = xsim(i,j+1) - xsim(i,j)
            vec1y = ysim(i,j+1) - ysim(i,j)
            vec1z = zsim(i,j+1) - zsim(i,j)
            timedif = tsim(i,j+1) - tsim(i,j)
            dist2 = vec1x**2 + vec1y**2 + vec1z**2 
            dist = sqrt(dist2)+1.d-16
            dista = dist
            
            if (tsim(i,j) .gt. -1.d-6 .and.
     &          tsim(i,j+1) .gt. -1.d-6 .and.
     &          tsim(i,j+2) .gt. -1.d-6) then
            
               timea = timedif
               velinitial = dist/timedif
               vec1x = vec1x/dist
               vec1y = vec1y/dist
               vec1z = vec1z/dist

               vec2x = xsim(i,j+2) - xsim(i,j+1)
               vec2y = ysim(i,j+2) - ysim(i,j+1)
               vec2z = zsim(i,j+2) - zsim(i,j+1)
               dist2 = vec2x**2 + vec2y**2 + vec2z**2 
               dist = sqrt(dist2)+1.d-16
               timedif = tsim(i,j+2)-tsim(i,j+1)
               timesum2 = timea + timedif
               distsum2 = dista + dist
               vec2x = vec2x/dist
               vec2y = vec2y/dist
               vec2z = vec2z/dist
                  
               dot = vec1x*vec2x+vec1y*vec2y+vec1z*vec2z
               angle = acos(dot)*180.d0/pi

               if (dist .gt. 1.d-12 .and.
     &             dista .gt. 1.d-12) then
                  iangle = iangle + 1
                  turnang(iangle) = angle
                  velangle(iangle) = distsum2/timesum2
                  velangle(iangle) = dista/timea
                  velanglec(iangle) = dist/timedif                                    
                  write(220,*) turnang(iangle),velangle(iangle)*60.

                  speedaa = speed_min + dble(jchoose-1)*dspeed_sec
                  speedbb = speed_min + dble(jchoose)*dspeed_sec

c                  print*,'speedaa speedbb = ',speedaa,speedbb
c                  print*,'dspeed_sec = ',dspeed_sec
c     stop
                  
                  do k = 1,num_volume_levels
                     if (volume_type(i) .eq. k .and.
     &                   (velinitial .gt. speedaa .and.
     &                    velinitial .le. speedbb) ) then
c                        write(420+k,*) turnang(iangle),velangle(iangle)*60.
                     write(420+k,*) turnang(iangle),dist*60.d0/timedif
                     end if
                  end do
                        
               end if

            end if
            
         end do
      end do
      
c      xmin = 1.d+20
c      ymin = 1.d+20
c      zmin = 1.d+20
c      do i = 1,num_tcells_allmovies
c         do j = 1,numframesa(i)
c            xmin = min(xsim(i,j),xmin)
c            ymin = min(ysim(i,j),ymin)
c            zmin = min(zsim(i,j),zmin)
c         end do
c      end do
cc      print*,'xmin ymin zmin = ',xmin,ymin,zmin
c
c      do i = 1,num_tcells_allmovies
c         do j = 1,numframesa(i)
c            xsim(i,j) = xsim(i,j) - xmin
c            ysim(i,j) = ysim(i,j) - ymin
c            zsim(i,j) = zsim(i,j) - zmin            
c         end do
c      end do


      xmin = 1.d+20
      ymin = 1.d+20
      zmin = 1.d+20
      do i = 1,num_tcells_allmovies
         do j = 1,numframes_v(i)            
            xmin = min(xsim(i,j),xmin)
            ymin = min(ysim(i,j),ymin)
            zmin = min(zsim(i,j),zmin)
         end do
      end do
c      print*,'xmin ymin zmin = ',xmin,ymin,zmin

      xmina = 1.d+20
      ymina = 1.d+20
      zmina = 1.d+20      
      xmaxa = -1.d+20
      ymaxa = -1.d+20
      zmaxa = -1.d+20
      do i = 1,num_tcells_allmovies
         do j = 1,numframes_v(i)                     

            xsim(i,j) = xsim(i,j) - xmin
            ysim(i,j) = ysim(i,j) - ymin
            zsim(i,j) = zsim(i,j) - zmin

            xmina = min(xsim(i,j),xmina)
            ymina = min(ysim(i,j),ymina)
            zmina = min(zsim(i,j),zmina)                        
            xmaxa = max(xsim(i,j),xmaxa)
            ymaxa = max(ysim(i,j),ymaxa)
            zmaxa = max(zsim(i,j),zmaxa)
            
         end do
      end do
      print*,'Minimum and Maximum of virtual comparision cells'
      print*,'xmina ymina zmina = ',xmina,ymina,zmina
      print*,'xmaxa ymaxa zmaxa = ',xmaxa,ymaxa,zmaxa
      print*,'rangex rangey rangez = ',
     &       xmaxa-xmina,ymaxa-ymina,zmaxa-zmina
      rangex = xmaxa - xmina
      rangey = ymaxa - ymina
      rangez = zmaxa - zmina
      if (rangex .ge. maxrange .or.
     &    rangey .ge. maxrange .or.
     &    rangez .ge. maxrange) then
         print*,'rangex rangey rangez = ',
     &        rangex,rangey,rangez
         stop
      end if
      
c      stop

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calculate the volume traversed by each T cell
      print*,'Calculate volume patrolled by virtual cells'

      facab = 1.d0
      do i = 1,num_tcells_allmovies
         grid = .false.
         tcellmaxtime(i) = -1.d0
         tcellmintime(i) = 1.d+10
         inbox = 1
         
c         do j = 1,numframesa(i)
         do j = 1,numframes_v(i)            

c           (iig,jjg,kkg) are the logical coordinates of the T cell

c           Using initial domain of empirical cells for volume calculation
            if (xsim(i,j) .ge. xminb/facab .and.
     &          xsim(i,j) .le. xmaxb/facab .and.
     &          ysim(i,j) .ge. yminb/facab .and.
     &          ysim(i,j) .le. ymaxb/facab .and.
     &          zsim(i,j) .ge. zminb/facab .and.
     &          zsim(i,j) .le. zmaxb/facab .and.
     &          inbox .eq. 1) then
               inbox = 1
            else
               inbox = 0
c               print*,'xsim(',i,j,') = ',xsim(i,j),ysim(i,j),zsim(i,j)
            end if
            inbox = 1
            
c            print*,'i j = ',i,j

            if (tsim(i,j) .gt. -1.d-6 .and. inbox .eq. 1) then

               tcellmaxtime(i) = max(tcellmaxtime(i),tsim(i,j))
               tcellmintime(i) = min(tcellmintime(i),tsim(i,j))
               
               iig = (xsim(i,j) - xref)/dxg
               jjg = (ysim(i,j) - yref)/dxg
               kkg = (zsim(i,j) - zref)/dxg
c               print*,'iig jjg kkg = ',iig,jjg,kkg,ncube
               
               do ig = iig-2,iig+2
                  do jg = jjg-2,jjg+2
                     do kg = kkg-2,kkg+2
c                        print*,'ig jg kg = ',ig,jg,kg
                        xc = dxg*dble(ig) + xref + dxg2
                        yc = dxg*dble(jg) + yref + dxg2
                        zc = dxg*dble(kg) + zref + dxg2
                        distance2 = (xsim(i,j) - xc)**2 +
     &                              (ysim(i,j) - yc)**2 +
     &                              (zsim(i,j) - zc)**2
                        if (ig .gt. ncube .or.
     &                      jg .gt. ncube .or.
     &                      kg .gt. ncube .or.
     &                      ig .le. 0 .or.
     &                      jg .le. 0 .or.
     &                      kg .le. 0) then
c                           print*,'ig jg kg = ',ig,jg,kg
c                           stop
                        else
                           if (distance2 .lt. 
     &                          radius_tcell_squared) then
c                              print*,'ig jg kg = ',ig,jg,kg
                              grid(ig,jg,kg) = .true.
                           end if
                        end if
                     end do
                  end do
               end do
c               print*,'End triple loop'
            end if

         end do

c         print*,'Volume'
         
         volumetcellsim(i) = 0.d0
         do kkg = 1,ncube
            do jjg = 1,ncube
               do iig = 1,ncube
                  if (grid(iig,jjg,kkg)) then
                     volumetcellsim(i) = 
     &               volumetcellsim(i) + vg
                  end if
               end do
            end do
         end do

c        timeobv = tsim(i,numframesa(i))-tsim(i,1)
         timeobv = tcellmaxtime(i) - tcellmintime(i)
         if (timeobv .gt. 0.d0) then
            volumetcellsim(i) = volumetcellsim(i)/timeobv
            write(197,*) volumetcellsim(i)
         end if

      end do
      print*,'Finished virtual simulation'


      return
      end


      subroutine cross(ax,ay,az,bx,by,bz,cx,cy,cz)
      implicit none
      double precision ax,ay,az,bx,by,bz,cx,cy,cz
      cx = ay*bz - by*az
      cy = bx*az - ax*bz
      cz = ax*by - ay*bx
      end

      subroutine random_numbere(rr)
      implicit none
      double precision rr
      call random_number(rr)
      end 

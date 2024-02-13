program Goes_ReBroadcast_nc2bufr
!
! Purpose: Convert GOES ReBroadcast netCDF files to BUFR format.
!          Currently only processes bands 7-16.
!
! Author: Jamie Bresch NCAR/MMM
!
!ifort -o GOES16_nc2bufr.exe GOES_GRB_nc2bufr.f90 \
!   -L${NETCDF}/lib -lnetcdf -lnetcdff -lm -I${NETCDF}/include \
!   -L${GSI_COREDIR}/lib -lbufr_i4r8
!
! gfortran -o GOES16_nc2bufr.exe -ffree-line-length-none  GOES_GRB_nc2bufr.f90 \
!    -L${NETCDF}/lib -lnetcdf -lnetcdff -lm -I${NETCDF}/include  \
!    -L${GSI_COREDIR}/lib -lbufr_i4r8
!
! input files:
!    (1) flist.txt: contains a list of nc files (exclude path) to be processed
!                     GoesReBroadcast file
!                     (optional) Clear Sky Mask output of cspp-geo-aitf package
!    (2) namelist.goes_nc2bufr
!        &data_nml
!          nc_list_file = 'flist.txt'
!          data_dir = '/parc/hclin/data/goes', ! path of the GRB nc files
!          data_id = 'OR_ABI-L1b-RadC-M3'      ! prefix of the downloaded GRB nc files
!          sat_id = 'G16'
!          bufr_tbl_file = 'bufrtab_NC021040.txt'
!          n_subsample = 1
!        /
!    (3) bufrtab_NC021040.txt: (temporary local) BUFR table
!
! Modified by Paola Corrales to work with GSI 3.7

   implicit none
   include 'netcdf.inc'

   integer, parameter  :: r_single = selected_real_kind(6)  ! single precision
   integer, parameter  :: r_double = selected_real_kind(15) ! double precision
   integer, parameter  :: i_byte   = selected_int_kind(1)   ! byte integer
   integer, parameter  :: i_short  = selected_int_kind(4)   ! short integer
   integer, parameter  :: i_long   = selected_int_kind(8)   ! long integer
   integer, parameter  :: i_kind   = i_long                 ! default integer
   integer, parameter  :: r_kind   = r_double               ! default real

   ! prefix of Clear Sky Mask (Binary Cloud Mask) output of cspp-geo-aitf package
   character(len=14), parameter :: BCM_id = 'OR_ABI-L2-ACMF'

   integer(i_kind) :: nband      = 10  ! IR bands 7-16
   integer(i_kind) :: band_start = 7
   integer(i_kind) :: band_end   = 16

   real(r_kind) :: pi, deg2rad, rad2deg

   real(r_kind), allocatable :: glat(:,:)    ! grid latitude (nx,ny)
   real(r_kind), allocatable :: glon(:,:)    ! grid longitude (nx,ny)
   real(r_kind), allocatable :: gzen(:,:)    ! satellite zenith angle (nx,ny)
   real(r_kind), allocatable :: solzen(:,:)  ! solar zenith angle (nx,ny)

   real(r_kind),    allocatable :: rad_2d(:,:)  ! radiance(nx,ny)
   real(r_kind),    allocatable :: bt_2d(:,:)   ! brightness temperature(nx,ny)
   real(r_kind),    allocatable :: sdtb_2d(:,:) ! std_dev (nx,ny) !P: add sdtb 2d
   integer(i_kind), allocatable :: qf_2d(:,:)   ! quality flag(nx,ny)
   integer(i_kind), allocatable :: cm_2d(:,:)   ! cloud_mask(nx,ny)

   type rad_type
      real(r_kind),    allocatable :: rad(:,:,:)  ! radiance(nband,nx,ny)
      real(r_kind),    allocatable :: bt(:,:,:)   ! brightness temperature(nband,nx,ny)
      integer(i_kind), allocatable :: qf(:,:,:)   ! quality flag(nband,nx,ny)
      real(r_kind),    allocatable :: sd(:,:,:)   ! std_dev(nband,nx,ny) P: change to 3D
      integer(i_kind), allocatable :: cm(:,:)     ! cloud mask(nx,ny)
   end type rad_type
   type(rad_type), allocatable  :: rdata(:)  ! (ntime)

   character(len=22), allocatable :: time_start(:)  ! (ntime) 2017-10-01T18:02:19.6Z

   integer(i_kind) :: ncid, nf_status
   integer(i_kind) :: nx, ny
   integer(i_kind) :: it, ib, ii, i, j
   integer(i_kind) :: ntime
   integer(i_kind) :: t_index
   integer(i_kind) :: band_id

   namelist /data_nml/ nc_list_file, data_dir, data_id, sat_id, bufr_tbl_file, n_subsample

   integer(i_kind)      :: nml_unit = 81
   integer(i_kind)      :: tbl_unit = 86
   integer(i_kind)      :: iunit    = 87

   character(len=256)              :: nc_list_file  ! the text file that contains a list of netcdf files to process
   character(len=256)              :: data_dir
   character(len=18)               :: data_id
   character(len=3)                :: sat_id
   character(len=256)              :: bufr_tbl_file
   integer(i_kind)                 :: n_subsample

!   real(r_kind)                    :: sdtb ! to be done  !P: replaced by sdtb_2b
   integer(i_kind)                 :: istat
   integer(i_kind)                 :: nfile, ifile, nlen
   logical                         :: isfile
   logical                         :: found_time
   logical                         :: got_grid_info
   logical, allocatable            :: valid(:), is_BCM(:)
   character(len=256), allocatable :: nc_fnames(:)
   character(len=256)              :: fname
   character(len=256)              :: bufr_out_fname
   character(len=256)              :: txtbuf
   character(len=18)               :: finfo
   character(len=2)                :: mode_id, scan_mode
   character(len=22), allocatable  :: scan_time(:) ! 2017-10-01T18:02:19.6Z
   integer(i_kind),   allocatable  :: fband_id(:)
   integer(i_kind),   allocatable  :: ftime_id(:)
   integer(i_kind),   allocatable  :: julianday(:)

   continue

   pi = acos(-1.0)
   deg2rad = pi/180.0
   rad2deg = 1.0/deg2rad
   !
   ! initialize namelist variables
   !
   nc_list_file   = 'flist.txt'
   data_dir       = '.'
   data_id        = 'OR_ABI-L1b-RadF-M3'
   sat_id         = 'G16'
   bufr_tbl_file  = 'bufrtab_NC021040.txt'
   n_subsample    = 1
   !
   ! read namelist
   !
   open(unit=nml_unit, file='namelist.goes_nc2bufr', status='old', form='formatted')
   read(unit=nml_unit, nml=data_nml, iostat=istat)
   write(0,nml=data_nml)
   if ( istat /= 0 ) then
      write(0,*) 'Error reading namelist data_nml'
      stop
   end if

   inquire(file=trim(bufr_tbl_file), exist=isfile)
   if ( .not. isfile ) then
      write(0,*) 'Error: File BUFR Table '//trim(bufr_tbl_file)//' not found'
      stop
   end if

   ! get file names from nc_list_file
   nfile  = 0  ! initialize the number of netcdf files to read
   inquire(file=trim(nc_list_file), exist=isfile)
   if ( .not. isfile ) then
      write(0,*) 'File not found: nc_list_file '//trim(nc_list_file)
      stop
   else
      open(unit=iunit, file=trim(nc_list_file), status='old', form='formatted')
      !first find out the number of netcdf files to read
      istat = 0
      do while ( istat == 0 )
         read(unit=iunit, fmt='(a)', iostat=istat) txtbuf
         if ( istat /= 0 ) then
            exit
         else
            nfile = nfile + 1
         end if
      end do
      if ( nfile > 0 ) then
         allocate (nc_fnames(nfile))
         !read the nc_list_file again to get the netcdf file names
         rewind(iunit)
         do ifile = 1, nfile
            read(unit=iunit, fmt='(a)', iostat=istat) nc_fnames(ifile)
         end do
      else
         write(0,*) 'File not found from nc_list_file '//trim(nc_list_file)
         stop
      end if
      close(iunit)
   end if !nc_list_file

   allocate (ftime_id(nfile))
   allocate (scan_time(nfile))
   allocate (julianday(nfile))
   allocate (fband_id(nfile))
   allocate (valid(nfile))
   allocate (is_BCM(nfile))
   valid( :) = .false.
   is_BCM(:) = .false.

   nlen = len_trim(data_id)
   mode_id = data_id(nlen-1:nlen)

   ! parse the file list
   t_index = 0
   file_loop1: do ifile = 1, nfile

      fname = trim(data_dir)//'/'//trim(nc_fnames(ifile))
      inquire(file=trim(fname), exist=isfile)
      if ( .not. isfile ) then
         write(0,*) 'File not found: '//trim(fname)
         cycle file_loop1
      end if

      ! OR_ABI-L1b-RadC-M3C16_G16_s20172741802196_e20172741804580_c20172741805015.nc
      ! retrieve some basic info from the netcdf filename itself
      call decode_nc_fname(trim(nc_fnames(ifile)),finfo, scan_mode, is_BCM(ifile), fband_id(ifile), scan_time(ifile), julianday(ifile))

      ! all files must be the same mode
      if ( scan_mode /= mode_id ) then
         cycle file_loop1
      end if

      if ( .not. is_BCM(ifile) ) then
         ! id of the file name must match specified data_id
         if ( finfo /= data_id ) then
            cycle file_loop1
         else
            ! only process band 7-16
            if ( fband_id(ifile) < band_start .or. fband_id(ifile) > band_end ) then
               cycle file_loop1
            end if
         end if
      end if

      valid(ifile) = .true.

      ! group files of the same scan time
      if ( t_index == 0 ) then
         t_index = t_index + 1
         ftime_id(ifile) = t_index
      else
         found_time = .false.
         find_time_loop: do ii = ifile-1, 1, -1
            if ( valid(ii) ) then
               if ( scan_time(ifile) == scan_time(ii) ) then
                  ftime_id(ifile) = ftime_id(ii)
                  found_time = .true.
                  exit find_time_loop
               end if
            end if
         end do find_time_loop
         if ( .not. found_time ) then
            t_index = t_index + 1
            ftime_id(ifile) = t_index
         end if
      end if

      ntime = t_index

   end do file_loop1

   if ( ntime <= 0 ) then
      write(0,*) 'ntime = ', ntime
      write(0,*) 'No valid files found from nc_list_file '//trim(nc_list_file)
      stop
   end if

   allocate (time_start(ntime))
   allocate (rdata(ntime))

   got_grid_info = .false.
   file_loop2: do ifile = 1, nfile

      if ( valid(ifile) ) then

         fname = trim(data_dir)//'/'//trim(nc_fnames(ifile))
         nf_status = nf_OPEN(trim(fname), nf_NOWRITE, ncid)
         if ( nf_status == 0 ) then
            write(0,*) 'Reading '//trim(fname)
         else
            write(0,*) 'ERROR reading '//trim(fname)
            cycle file_loop2
         end if

         if ( .not. got_grid_info ) then
            call read_GRB_dims(ncid, nx, ny)
            allocate (glat(nx, ny))
            allocate (glon(nx, ny))
            allocate (gzen(nx, ny))
            allocate (solzen(nx, ny))
            write(0,*) 'Calculating lat/lon from fixed grid x/y...'
            call read_GRB_grid(ncid, nx, ny, glat, glon, gzen)
            call calc_solar_zenith_angle(nx, ny, glat, glon, scan_time(ifile), julianday(ifile), solzen)
            got_grid_info = .true.
            allocate (rad_2d(nx, ny))
            allocate (bt_2d(nx, ny))
            allocate (qf_2d(nx, ny))
            allocate (cm_2d(nx, ny))
            allocate (sdtb_2d(nx, ny))
         end if

         it = ftime_id(ifile)
         ib = fband_id(ifile)

         if ( .not. is_BCM(ifile) ) then

            call read_GRB(ncid, nx, ny, rad_2d, bt_2d, qf_2d, sdtb_2d, band_id, time_start(it))

            if ( band_id /= ib ) then
               write(0,*) 'ERROR: band_id from the file name and the file content do not match.'
               cycle file_loop2
            end if

            if ( time_start(it) /= scan_time(ifile) ) then
               write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
               cycle file_loop2
            end if

            if ( .not. allocated(rdata(it)%rad) ) allocate (rdata(it)%rad(nband,nx,ny))
            if ( .not. allocated(rdata(it)%bt) )  allocate (rdata(it)%bt(nband,nx,ny))
            if ( .not. allocated(rdata(it)%qf) )  allocate (rdata(it)%qf(nband,nx,ny))
            if ( .not. allocated(rdata(it)%sd) )  allocate (rdata(it)%sd(nband,nx,ny)) !P: change to 3D

            do j = 1, ny
               do i = 1, nx
                  ! convert band id 7-16 to array index 1-10
                  rdata(it)%rad(ib-band_start+1,i,j) = rad_2d(i,j)
                  rdata(it)%bt(ib-band_start+1,i,j)  = bt_2d(i,j)
                  rdata(it)%qf(ib-band_start+1,i,j)  = qf_2d(i,j)
                  rdata(it)%sd(ib-band_start+1,i,j)  = sdtb_2d(i,j) !P: change sdtb to sdtb_2d
               end do
            end do

         else

            call read_L2_BCM(ncid, nx, ny, cm_2d, time_start(it))

            if ( time_start(it) /= scan_time(ifile) ) then
               write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
               cycle file_loop2
            end if

            if ( .not. allocated(rdata(it)%cm) )  allocate (rdata(it)%cm(nx,ny))
            rdata(it)%cm(:,:) = cm_2d(:,:)

         end if

         nf_status = nf_CLOSE(ncid)

      end if

   end do file_loop2

   if ( allocated(rad_2d) ) deallocate(rad_2d)
   if ( allocated(bt_2d) )  deallocate(bt_2d)
   if ( allocated(qf_2d) )  deallocate(qf_2d)
   if ( allocated(cm_2d) )  deallocate(cm_2d)
   if ( allocated(sdtb_2d) )  deallocate(sdtb_2d) !P: add allocatable sdtb_2d

   do it = 1, ntime
      bufr_out_fname = trim(data_id)//'_'//sat_id//'_'//time_start(it)//'.bufr'
      write(0,*) 'Writing ', trim(bufr_out_fname)
      if ( allocated(rdata(it)%cm) ) then
         call write_bufr(trim(bufr_out_fname), time_start(it), nx, ny, nband, &
            glat, glon, gzen, solzen, rdata(it)%bt, rdata(it)%qf, rdata(it)%sd, rdata(it)%cm)
      else
         call write_bufr(trim(bufr_out_fname), time_start(it), nx, ny, nband, &
            glat, glon, gzen, solzen, rdata(it)%bt, rdata(it)%qf, rdata(it)%sd)
      end if
   end do

   if ( allocated(glat) )   deallocate(glat)
   if ( allocated(glon) )   deallocate(glon)
   if ( allocated(gzen) )   deallocate(gzen)
   if ( allocated(solzen) ) deallocate(solzen)

   do it = 1, ntime
      if ( allocated(rdata(it)%rad) ) deallocate (rdata(it)%rad)
      if ( allocated(rdata(it)%bt)  ) deallocate (rdata(it)%bt)
      if ( allocated(rdata(it)%qf)  ) deallocate (rdata(it)%qf)
      if ( allocated(rdata(it)%cm)  ) deallocate (rdata(it)%cm)
      if ( allocated(rdata(it)%sd)  ) deallocate (rdata(it)%sd)
   end do
   deallocate(rdata)
   deallocate(time_start)

   deallocate(nc_fnames)
   deallocate(ftime_id)
   deallocate(scan_time)
   deallocate(julianday)
   deallocate(fband_id)
   deallocate(valid)
   deallocate(is_BCM)

contains

subroutine read_GRB_dims(ncid, nx, ny)
   implicit none
   integer(i_kind), intent(in)  :: ncid
   integer(i_kind), intent(out) :: nx, ny
   integer(i_kind)              :: dimid
   integer(i_kind)              :: nf_status(4)
   continue
   nf_status(1) = nf_INQ_DIMID(ncid, 'x', dimid)
   nf_status(2) = nf_INQ_DIMLEN(ncid, dimid, nx)
   nf_status(3) = nf_INQ_DIMID(ncid, 'y', dimid)
   nf_status(4) = nf_INQ_DIMLEN(ncid, dimid, ny)
   if ( any(nf_status /= 0) ) then
      write(0,*) 'Error reading dimensions'
      stop
   end if
   return
end subroutine read_GRB_dims

!NC_BYTE 8-bit signed integer
!NC_SHORT 16-bit signed integer
!NC_INT (or NC_LONG) 32-bit signed integer
!NC_FLOAT 32-bit floating point
!NC_DOUBLE 64-bit floating point

subroutine read_GRB_grid(ncid, nx, ny, glat, glon, gzen)
   implicit none
   integer(i_kind), intent(in)    :: ncid
   integer(i_kind), intent(in)    :: nx, ny
   real(r_kind),    intent(inout) :: glat(nx,ny)
   real(r_kind),    intent(inout) :: glon(nx,ny)
   real(r_kind),    intent(inout) :: gzen(nx,ny)
   integer(i_kind)                :: varid, i, j
   integer(i_kind)                :: nf_status
   integer(i_kind)                :: istart(1), icount(1)
   integer(i_short), allocatable  :: itmp_short_1d(:)
   real(r_kind),     allocatable  :: x(:)
   real(r_kind),     allocatable  :: y(:)
   real(r_single) :: scalef, offset
   real(r_double) :: dtmp
   real(r_double) :: r_eq    ! GRS80 semi-major axis of earth
   real(r_double) :: r_pol   ! GRS80 semi-minor axis of earth = (1-f)*r_eq
   real(r_double) :: lon_sat ! satellite longitude, longitude_of_projection_origin
   real(r_double) :: h_sat   ! satellite height
   real(r_double) :: a, b, c, rs, sx, sy, sz
   real(r_kind)   :: rlat, rlon, lon_diff, tmp1, theta1, theta2
   continue

!int goes_imager_projection ;
!  goes_imager_projection:long_name = "GOES-R ABI fixed grid projection" ;
!  goes_imager_projection:grid_mapping_name = "geostationary" ;
!  goes_imager_projection:perspective_point_height = 35786023. ;
!  goes_imager_projection:semi_major_axis = 6378137. ;
!  goes_imager_projection:semi_minor_axis = 6356752.31414 ;
!  goes_imager_projection:inverse_flattening = 298.2572221 ;
!  goes_imager_projection:latitude_of_projection_origin = 0. ;
!  goes_imager_projection:longitude_of_projection_origin = -89.5 ;
!  goes_imager_projection:sweep_angle_axis = "x" ;

   nf_status = nf_INQ_VARID(ncid, 'goes_imager_projection', varid)
   nf_status = nf_GET_ATT_DOUBLE(ncid, varid, 'semi_major_axis',  dtmp)
   r_eq = dtmp
   nf_status = nf_GET_ATT_DOUBLE(ncid, varid, 'semi_minor_axis',  dtmp)
   r_pol = dtmp
   nf_status = nf_GET_ATT_DOUBLE(ncid, varid, 'perspective_point_height',  dtmp)
   h_sat = dtmp + r_eq  ! perspective_point_height + semi_major_axis
   nf_status = nf_GET_ATT_DOUBLE(ncid, varid, 'longitude_of_projection_origin',  dtmp)
   lon_sat = dtmp * deg2rad

!short x(x) ;
!  x:scale_factor = 5.6e-05f ;
!  x:add_offset = -0.075012f ;
!  x:units = "rad" ;
!  x:axis = "X" ;
!  x:long_name = "GOES fixed grid projection x-coordinate" ;

   istart(1) = 1
   icount(1) = nx
   allocate(itmp_short_1d(nx))
   nf_status = nf_INQ_VARID(ncid, 'x', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:1), icount(1:1), itmp_short_1d(:))
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   allocate(x(nx))
   do i = 1, nx
      x(i) = offset + itmp_short_1d(i) * scalef
   end do
   deallocate(itmp_short_1d)

!short y(y) ;
!  y:scale_factor = -5.6e-05f ;
!  y:add_offset = 0.126532f ;
!  y:units = "rad" ;
!  y:axis = "Y" ;
!  y:long_name = "GOES fixed grid projection y-coordinate" ;
!  y:standard_name = "projection_y_coordinate" ;

   istart(1) = 1
   icount(1) = ny
   allocate(itmp_short_1d(ny))
   nf_status = nf_INQ_VARID(ncid, 'y', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:1), icount(1:1), itmp_short_1d(:))
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   allocate(y(ny))
   do i = 1, ny
      y(i) = offset + itmp_short_1d(i) * scalef
   end do
   deallocate(itmp_short_1d)

   ! Product Definition and User's Guide (PUG) Volume 3, pp. 19-21
   ! from fixed grid x/y to geodetic lat/lon
   do j = 1, ny
      do i = 1, nx
         a = sin(x(i))*sin(x(i)) + cos(x(i))*cos(x(i)) * &
             (cos(y(j))*cos(y(j))+(r_eq/r_pol)*(r_eq/r_pol)*sin(y(j))*sin(y(j)))
         b = -2.0*h_sat*cos(x(i))*cos(y(j))
         c = h_sat*h_sat - r_eq*r_eq
         rs = (-1.0*b - sqrt(b*b-4.0*a*c)) / (2.0*a)
         sx = rs * cos(x(i)) * cos(y(j))
         sy = -1.0 * rs * sin(x(i))
         sz = rs * cos(x(i)) * sin(y(j))
         !glat(i,j) = (atan((r_eq/r_pol)*(r_eq/r_pol)*(sz/sqrt((h_sat-sx)*(h_sat-sx)+sy*sy)))) * rad2deg
         !glon(i,j) = (lon_sat - atan(sy/(h_sat-sx))) * rad2deg
         glat(i,j) = atan((r_eq/r_pol)*(r_eq/r_pol)*(sz/sqrt((h_sat-sx)*(h_sat-sx)+sy*sy)))
         glon(i,j) = lon_sat - atan(sy/(h_sat-sx))
      end do
   end do

   deallocate(x)
   deallocate(y)

   ! calculate geostationary satellite zenith angle
   do j = 1, ny
      do i = 1, nx
         rlat = glat(i,j) ! in radian
         rlon = glon(i,j) ! in radian
         lon_diff = abs(rlon-lon_sat)
         tmp1 = sqrt((2.0*r_eq*sin(lon_diff/2.)-r_eq*(1.0-cos(rlat))*sin(lon_diff/2.))**2 &
           +(2.0*r_eq*sin(rlat/2.))**2-(r_eq*(1.0-cos(rlat))*sin(lon_diff/2.))**2)
         theta1 = 2.0*asin(tmp1/r_eq/2.)
         theta2 = atan(r_eq*sin(theta1)/((h_sat-r_eq)+r_eq*(1.0-sin(theta1))))
         gzen(i,j) = (theta1+theta2) * rad2deg
         !gzen(i,j) = 90.0 - atan((cos(lon_diff)*cos(rlat)-0.1512)/(sqrt(1.0-cos(lon_diff)*cos(lon_diff)*cos(rlat)*cos(rlat)))) * rad2deg
         glat(i,j) = glat(i,j) * rad2deg
         glon(i,j) = glon(i,j) * rad2deg
      end do
   end do

   return
end subroutine read_GRB_grid

subroutine read_GRB(ncid, nx, ny, rad, bt, qf, sd, band_id, time_start)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   integer(i_kind),   intent(out)   :: band_id
   real(r_kind),      intent(out)   :: sd(nx,ny) !P: change to 2D
   real(r_kind),      intent(inout) :: rad(nx,ny)
   real(r_kind),      intent(inout) :: bt(nx,ny)
   integer(i_kind),   intent(inout) :: qf(nx,ny)
   character(len=22), intent(out)   :: time_start  ! 2017-10-01T18:02:19.6Z
   integer(i_byte),  allocatable    :: itmp_byte_1d(:)
   integer(i_byte),  allocatable    :: itmp_byte_2d(:,:)
   integer(i_short), allocatable    :: itmp_short_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_short)                 :: ifill
   real(r_single)                   :: rfill
   real(r_single)                   :: rtmp
   real(r_single)                   :: planck_fk1, planck_fk2
   real(r_single)                   :: planck_bc1, planck_bc2
   real(r_single)                   :: scalef, offset
   real(r_kind)                     :: rmiss = -999.0
   integer(i_kind)                  :: imiss = -999
   real(r_kind)                     :: temp(9)!for sd calculations
   integer(i_kind)                  :: w = 3  !window size for sd
   real(r_kind)                     :: m,s    !for sd calculations
   continue

   ! time_start is the same for all bands, but time_end is not
   nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_start', time_start)
   !nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_end',   time_end)

   istart(1) = 1
   icount(1) = 1
   allocate(itmp_byte_1d(1))
   nf_status = nf_INQ_VARID(ncid, 'band_id', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:1), icount(1:1), itmp_byte_1d(:))
   band_id = itmp_byte_1d(1)
   deallocate(itmp_byte_1d)

   !nf_status = nf_INQ_VARID(ncid, 'std_dev_radiance_value_of_valid_pixels', varid)
   !nf_status = nf_GET_VAR_REAL(ncid, varid, rtmp)
   !sd = rtmp  !P: comment this to calculate sd later

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'DQF', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:2), icount(1:2), itmp_byte_2d(:,:))
   qf(:,:) = imiss
   do j = 1, ny
      do i = 1, nx
         qf(i,j) = itmp_byte_2d(i,j)
      end do
   end do
   deallocate(itmp_byte_2d)

   nf_status = nf_INQ_VARID(ncid, 'planck_fk1', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, planck_fk1)
   nf_status = nf_GET_ATT_REAL(ncid, varid, '_FillValue',  rfill)

   nf_status = nf_INQ_VARID(ncid, 'planck_fk2', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, planck_fk2)

   nf_status = nf_INQ_VARID(ncid, 'planck_bc1', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, planck_bc1)

   nf_status = nf_INQ_VARID(ncid, 'planck_bc2', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, planck_bc2)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_short_2d(nx, ny))
   nf_status = nf_INQ_VARID(ncid, 'Rad', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:2), icount(1:2), itmp_short_2d(:,:))
   nf_status = nf_GET_ATT_INT2(ncid, varid, '_FillValue',  ifill)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   rad(:,:) = rmiss
   bt(:,:)  = rmiss
   do j = 1, ny
      do i = 1, nx
         if ( itmp_short_2d(i,j) /= ifill ) then
            rad(i,j) = offset + itmp_short_2d(i,j) * scalef
            if ( planck_fk1 /= rfill .and. planck_fk2 /= rfill .and. &
                 planck_bc1 /= rfill .and. planck_bc2 /= rfill ) then
               bt(i,j) = (planck_fk2/(log((planck_fk1/rad(i,j))+1.0))-planck_bc1)/planck_bc2
            end if
         end if
      end do
   end do
   deallocate(itmp_short_2d)

! Calculates sd for 3x3 pixels 
   
   sd(:,:) = rmiss
   w       = 3 !window size 
   do j = 2, ny-1 !skip border
      do i = 2, nx-1 ! skip border
         temp = reshape(bt(i:i+w-1, j:j+w-1), (/9/))  ! bt on a 3x3 pixels window
         m    = sum(temp)/size(temp)
         s    = sqrt(sum((temp - m)**2)/size(temp))
         sd(i,j) = s
      end do
   end do 

   return
end subroutine read_GRB

subroutine read_L2_BCM(ncid, nx, ny, cm, time_start)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   integer(i_kind),   intent(inout) :: cm(nx,ny)
   character(len=22), intent(out)   :: time_start  ! 2017-10-01T18:02:19.6Z
   integer(i_byte),  allocatable    :: itmp_byte_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_kind)                  :: imiss = -999
   integer(i_kind)                  :: qf(nx,ny)
   continue

   ! time_start is the same for all bands, but time_end is not
   nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_start', time_start)
   !nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_end',   time_end)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'DQF', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:2), icount(1:2), itmp_byte_2d(:,:))
   qf(:,:) = imiss
   do j = 1, ny
      do i = 1, nx
         qf(i,j) = itmp_byte_2d(i,j)
      end do
   end do
   deallocate(itmp_byte_2d)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'BCM', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:2), icount(1:2), itmp_byte_2d(:,:))
   cm(:,:) = imiss
   do j = 1, ny
      do i = 1, nx
         if ( qf(i,j) == 0 ) then ! good quality
            cm(i,j) = itmp_byte_2d(i,j)
         end if
      end do
   end do
   deallocate(itmp_byte_2d)

   return
end subroutine read_L2_BCM

subroutine decode_nc_fname(fname, finfo, scan_mode, is_BCM, band_id, start_time, jday)
   implicit none
   character(len=*),  intent(in)  :: fname
   character(len=18), intent(out) :: finfo
   character(len=2),  intent(out) :: scan_mode
   logical,           intent(out) :: is_BCM
   integer(i_kind),   intent(out) :: band_id
   character(len=22), intent(out) :: start_time
   integer(i_kind),   intent(out) :: jday
   integer(i_kind) :: year, month, day, hour, minute, sec1, sec2

   if ( fname( 1:14) == BCM_id ) then
      is_BCM = .true.
      band_id = -99
   else
      is_BCM = .false.
   end if
   !CG_ABI-L2-ACMC-M3_G16_s20180351202275_e20180351205060_c20180351205106.nc
   !OR_ABI-L1b-RadC-M3C16_G16_s20172741802196_e20172741804580_c20172741805015.nc
   !1234567890123456789012345678901234567890123456789012345678901234567890123456
   if ( .not. is_BCM ) then
      read(fname( 1:18), '(a18)') finfo
      read(fname(17:18), '(a2)')  scan_mode
      read(fname(20:21), '(i2)')  band_id
      read(fname(28:31), '(i4)')  year
      read(fname(32:34), '(i3)')  jday
      read(fname(35:36), '(i2)')  hour
      read(fname(37:38), '(i2)')  minute
      read(fname(39:40), '(i2)')  sec1   ! integer part of second
      read(fname(41:41), '(i1)')  sec2   ! decimal part of second
      ! get month and day from julian day
      call get_date(year, jday, month, day)
      ! 2017-10-01T18:02:19.6Z
      write(start_time,'(i4.4,4(a,i2.2),a,i2.2,a,i1,a)') &
            year, '-', month, '-', day, 'T', hour, ':',  minute, ':', sec1, '.', sec2, 'Z'
   else
      read(fname( 1:17), '(a17)') finfo
      read(fname(16:17), '(a2)')  scan_mode
      read(fname(24:27), '(i4)')  year
      read(fname(28:30), '(i3)')  jday
      read(fname(31:32), '(i2)')  hour
      read(fname(33:34), '(i2)')  minute
      read(fname(35:36), '(i2)')  sec1   ! integer part of second
      read(fname(37:37), '(i1)')  sec2   ! decimal part of second
      ! get month and day from julian day
      call get_date(year, jday, month, day)
      ! 2017-10-01T18:02:19.6Z
      write(start_time,'(i4.4,4(a,i2.2),a,i2.2,a,i1,a)') &
            year, '-', month, '-', day, 'T', hour, ':',  minute, ':', sec1, '.', sec2, 'Z'
   end if
   return
end subroutine decode_nc_fname

subroutine get_date(ccyy, jday, month, day)
   implicit none
   integer(i_kind), intent(in)  :: ccyy, jday
   integer(i_kind), intent(out) :: month, day
   integer(i_kind) :: mmday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer(i_kind) :: i, jdtmp
   continue

   if ( MOD(ccyy,4) == 0 ) then
      mmday(2) = 29
      if ( MOD(ccyy,100) == 0 ) then
         mmday(2) = 28
      end if
      if ( MOD(ccyy,400) == 0 ) then
         mmday(2) = 29
      end if
   end if

   jdtmp = 0
   do i = 1, 12
      jdtmp = jdtmp + mmday(i)
      if ( jday <= jdtmp ) then
         month = i
         day = jday - ( jdtmp - mmday(i) )
         exit
      end if
   end do

   return
end subroutine get_date

subroutine write_bufr(bufr_fname, time_start, nx, ny, nband, lat, lon, sat_zen, sun_zen, bt, qf, sdtb, cloudmask)

   implicit none

   character(len=*),   intent(in) :: bufr_fname
   character(len=22),  intent(in) :: time_start
   integer(i_kind),    intent(in) :: nx, ny, nband
   real(r_kind),       intent(in) :: lat(nx,ny)
   real(r_kind),       intent(in) :: lon(nx,ny)
   real(r_kind),       intent(in) :: sat_zen(nx,ny)
   real(r_kind),       intent(in) :: sun_zen(nx,ny)
   real(r_kind),       intent(in) :: bt(nband,nx,ny)
   integer(i_kind),    intent(in) :: qf(nband,nx,ny)
   real(r_kind),       intent(in) :: sdtb(nband,nx,ny) !P: change to 3D
   integer(i_kind),    intent(in), optional :: cloudmask(nx,ny)

   integer(i_kind)  :: cm(nx,ny)
   character(len=8) :: subset
   integer(i_kind)  :: iunit_bufr, iunit_dx, iret, itmp
   integer(i_kind)  :: iline, isample, iband, ldate
   integer(i_kind)  :: year, month, day, hour, min, sec
   real(r_double)   :: timeinfo(6)
   real(r_double)   :: latlon(2)
   real(r_double)   :: satinfo(10)
   real(r_double)   :: cminfo(10)
   !hcl real(r_double)   :: radout(13,nband)
   real(r_double)   :: radout(2,nband)
   real(r_double)   :: sdtbout(nband)  !P: add sdtbout 
   real(r_double), parameter :: bmiss = 10.0e10
   character(80)    :: hdrstr
   real(r_double)   :: hdrabi(13)

   iunit_bufr = 50
   iunit_dx = 51

   subset = 'NC021046' ! clear sky id
   timeinfo = bmiss
   latlon = bmiss
   satinfo = bmiss
   radout = bmiss
   sdtbout = bmiss
   cminfo = bmiss

   open(iunit_bufr,file=trim(bufr_fname),iostat=iret,form='unformatted',status='unknown')
   if ( iret /= 0 ) then
      write(0,*) 'error opening file ', trim(bufr_fname)
      stop
   end if
   open(iunit_dx,file=trim(bufr_tbl_file),iostat=iret,form='formatted',status='old')
   if ( iret /= 0 ) then
      write(0,*) 'error opening bufr table'
      stop
   end if

   if ( present(cloudmask) ) then
      cm(:,:) = cloudmask(:,:)
   else
      cm(:,:) = -99
   end if

   call openbf(iunit_bufr,'OUT',iunit_dx)

   read(time_start( 1: 4), '(i4)') year
   read(time_start( 6: 7), '(i2)') month
   read(time_start( 9:10), '(i2)') day
   read(time_start(12:13), '(i2)') hour
   read(time_start(15:16), '(i2)') min
   read(time_start(18:19), '(i2)') sec

   timeinfo(1) = year
   timeinfo(2) = month
   timeinfo(3) = day
   timeinfo(4) = hour
   timeinfo(5) = min
   timeinfo(6) = sec
   ldate = year*1000000+month*10000+day*100+hour
   satinfo(1) = 270 ! GOES-16

   hdrstr='SAID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH SAZA BEARAZ SOZA SOLAZI'

   do iline = 1, ny, n_subsample
      do isample = 1, nx, n_subsample
         if ( all(bt(:,isample,iline)<0.0) ) cycle
         latlon(1) = lat(isample,iline)
         latlon(2) = lon(isample,iline)
         satinfo(8) = sat_zen(isample,iline)
         satinfo(9) = sun_zen(isample,iline)
         select case ( cm(isample,iline) )
            ! set Cloud amount in segment based on cloud mask
            case ( 0 )
               ! clear or probably clear
               cminfo = 100.0  ! It was 0.0 but the cloud mask is different
            case ( 1 )
               ! cloudy or probably cloudy
               cminfo = 0.0
         end select
         do iband = 1, nband
            radout(1,iband) = iband + 6
            ! qf (DQF, Data Quality Flag)
            ! 0:good, 1:conditionally_usable, 2:out_of_range, 3:no_value
            ! keep only qf=0,1 pixels
            if ( qf(iband,isample,iline) > 1 ) cycle
            if ( bt(iband,isample,iline) > 0.0 ) then
               !hcl radout(8,iband) = bt(iband,isample,iline)
               radout(2,iband) = bt(iband,isample,iline)
               sdtbout(iband) = sdtb(iband,isample,iline) !P: add sdtb
            end if
            !if ( sdtb(iband) > 0.0 ) then
               !hcl radout(12,iband) = sdtb(iband)
            !end if
            select case ( qf(iband,isample,iline) )
             case ( 0 )
                !hcl radout(13,iband) = 100.0
             case ( 1 )
                !hcl radout(13,iband) = 60.0
            end select
         end do
 
         hdrabi(1)=satinfo(1)
         hdrabi(2:7)=timeinfo
         hdrabi(8:9)=latlon
         hdrabi(10)=satinfo(8)
         hdrabi(11)=bmiss
         hdrabi(12)=satinfo(9)
         hdrabi(13)=bmiss

         call openmb(iunit_bufr,subset,ldate)
         call ufbint(iunit_bufr,hdrabi,13,1,iret,hdrstr)
             
         !call ufbint(iunit_bufr,timeinfo,6,1,iret,&
         !   'YEAR MNTH DAYS HOUR MINU SECO')
         !call ufbint(iunit_bufr,latlon,2,1,iret,&
         !   'CLAT CLON')
         !call ufbint(iunit_bufr,satinfo,10,1,iret,&
         !   'SAID GCLONG SCLF SSNX SSNY NPPR NPPC SAZA SOZA LSQL')
         !hcl call ufbrep(iunit_bufr,radout,13,nband,iret,&
         !hcl    'SIDP RDTP RDCM SCCF SCBW SPRD RDNE TMBRST CLDMNT NCLDMNT CLTP SDTB PCCF')
        
         call ufbrep(iunit_bufr,cminfo,1,nband,iret,'NCLDMNT')
         call ufbrep(iunit_bufr,radout(2,:),1,nband,iret,'TMBRST')
         call ufbrep(iunit_bufr,sdtbout(:),1,nband,iret,'SDTB')         

         !call ufbrep(iunit_bufr,radout,2,nband,iret,'CHNM TMBR')
         !call ufbint(iunit_bufr,cminfo,1,1,iret,'CLDMNT')
         call WRITSB(iunit_bufr)
      end do
   end do

   call closbf(iunit_bufr)
   close(iunit_bufr)
   close(iunit_dx)

end subroutine write_bufr

subroutine calc_solar_zenith_angle(nx, ny, xlat, xlon, xtime, julian, solzen)

! the calulcation is adapted from subroutines radconst and calc_coszen in
! WRF phys/module_radiation_driver.F

   implicit none

   integer(i_kind),   intent(in)    :: nx, ny, julian
   real(r_kind),      intent(in)    :: xlat(nx,ny), xlon(nx,ny)
   character(len=22), intent(in)    :: xtime
   real(r_kind),      intent(inout) :: solzen(nx,ny)

   real(r_kind) :: obliq = 23.5
   real(r_kind) :: deg_per_day = 360.0/365.0
   real(r_kind) :: slon   ! longitude of the sun
   real(r_kind) :: declin ! declination of the sun
   real(r_kind) :: hrang, da, eot, xt, tloctm, rlat
   integer(i_kind) :: gmt, minute, i, j

   ! calculate longitude of the sun from vernal equinox
   if ( julian >= 80 ) slon = (julian - 80 ) * deg_per_day
   if ( julian <  80 ) slon = (julian + 285) * deg_per_day

   declin = asin(sin(obliq*deg2rad)*sin(slon*deg2rad)) ! in radian

   read(xtime(12:13), '(i2)') gmt
   read(xtime(15:16), '(i2)') minute

   da = 6.2831853071795862*(julian-1)/365.
   eot = (0.000075+0.001868*cos(da)-0.032077*sin(da) &
          -0.014615*cos(2.0*da)-0.04089*sin(2.0*da))*(229.18)
   xt = gmt + (minute + eot)/60.0

   do j = 1, ny
      do i = 1, nx
         tloctm = xt + xlon(i,j)/15.0
         hrang = 15.0*(tloctm-12.0) * deg2rad
         rlat = xlat(i,j) * deg2rad
         solzen(i,j) = acos( sin(rlat)*sin(declin) + &
                             cos(rlat)*cos(declin)*cos(hrang) )
         solzen(i,j) = solzen(i,j) * rad2deg
      enddo
   enddo

   return
end subroutine calc_solar_zenith_angle

end program Goes_ReBroadcast_nc2bufr

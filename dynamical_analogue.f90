!***********************************************************************
!  DESCRIPTION:
!       Main Program to predict short-term climate using dynamical
!       analogue method based on NCEP CFSv2 operational products
!
!  PRECONDITIONS REQUIRED:
!       NCEP CFSv2 raw data should be processed on Unix-like platform
!       or in Windows/Cygwin environment first
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Li Zhenkun, SCC
!       Revised 10/2018 by Li Zhenkun, SCC
!
!***********************************************************************
  program dynamical_analogue

     use params_mod
     use inout_mod
     use eof_mod
     use date_mod, only :  new_date
     use utils_mod

     implicit none

     integer, parameter                 ::  stn_num = 523 ! observation station number
     real(4)                            ::  time_begin, time_end
     character(256)                     ::  progname, namelistfile
     integer                            ::  ierr
     integer                            ::  nlon, nlat
     integer                            ::  grid_num
     integer                            ::  fct_year, fct_mon, fct_day
     integer                            ::  samp_sgl_num
     integer                            ::  samp_max ! the maximum lenth of historical samples
     real(4), allocatable               ::  curt(:), hist_single(:), hist_valid(:, :), hist(:, :)
     integer, allocatable               ::  samp_dates(:, :), samp_dates_valid(:, :)
     integer                            ::  samp_num
     integer                            ::  incre
     integer                            ::  year, mon, day
     real(4), allocatable               ::  mean(:), sdev(:)
     integer                            ::  abnorm_num
     real(4)                            ::  minima = 1.e-25
     real(4)                            ::  undef = -999.0
     integer                            ::  min_dim
     integer                            ::  trans_opt = 2 !trans_opt = 0 for raw, 1 for anomaly, 2 for normalized anomoly
     real(4), allocatable               ::  eigenvalue(:)
     real(4), allocatable               ::  eigenvector(:, :)
     real(4), allocatable               ::  pc(:, :)
     real(4)                            ::  eigenvalue_sum
     integer                            ::  lead_num
     real(4), allocatable               ::  similar_val(:)
     real(4), allocatable               ::  proj_coeff_hist(:, :)
     real(4), allocatable               ::  proj_coeff_curt(:)
     real(4), allocatable               ::  prec_fct(:, :)
     integer                            ::  iyear
     integer                            ::  igrid, isamp, inum

     call cpu_time(time_begin)

!----------------------------------------------------------------------
!    Get the names of this program and namelist control file
!----------------------------------------------------------------------
     call getarg(0, progname)
     call getarg(1, namelistfile)
     call init_params(namelistfile, ierr)
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'Parameter initialization not completed'
        write( 6, '(a)' ) 'Usage : '
        write( 6, '(3a)' ) '        ', trim(progname), ' namelist'
        write( 6, '(a)' ) ' '
        write( 6, '(a)' ) 'Check argument and namelist syntax'
        stop
     end if

!----------------------------------------------------------------------
!    Check the length of leading days and average days are valid or not
!    according to the forecast period of NCEP CFSv2 (43 days), that is,
!    lead_days + ave_days should be NO more than 45 days
!----------------------------------------------------------------------
     if ( lead_days > 44 ) then
        write ( 6, '(a)' ) 'Leading days setting error, check namelist argument'
        stop
     end if
     if ( ave_days + lead_days > 45 ) then
        write ( 6, '(a)' ) 'Average days setting error, check namelist argument'
        stop
     end if

!----------------------------------------------------------------------
!    Check the domain settings are valid or not
!----------------------------------------------------------------------
     if ( domain_west < 0 )    domain_west = domain_west + 360
     if ( domain_west >= 360 ) domain_west = domain_west - 360
     if ( domain_east < 0 )    domain_east = domain_east + 360
     if ( domain_east >= 360 ) domain_east = domain_east - 360
     if ( domain_west > domain_east .or. domain_south > domain_north ) then
        write( 6, '(a)' ) 'Domain settings error, check namelist argument'
        stop
     end if
     nlon = domain_east - domain_west + 1
     nlat = domain_north - domain_south + 1
     grid_num = nlon * nlat * var_num
     write( 6, '(a, i3, a, i3)' ) 'Longitude grid num = ', nlon, '; Latitude grid num = ', nlat

!----------------------------------------------------------------------
!    Get varibles of current forecast from pre-processed binary file
!----------------------------------------------------------------------
     fct_year = fct_date/10000
     fct_mon  = (fct_date - fct_year*10000)/100
     fct_day  = mod(fct_date, 100)
     write( 6, '(a, i4, 1x, i2, 1x, i2)' ) 'The current forecast date is ', fct_year, fct_mon, fct_day

     allocate( curt(grid_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: curt'
        stop
     end if

     write( 6, '(a)' ) 'STEP 1: acquire current forecast data ...'
     call get_cfs_field(input_dir_curt, path_separator, var_num, var_names(1:var_num), fct_year, fct_mon, fct_day, &
                        lead_days, ave_days, domain_west, domain_south, nlon, nlat, curt, ierr)
     if ( ierr /= 0 ) then
        write( 6, '(a, i4, 1x, i2, 1x, i2, a)' ) 'ERROR: current forecast data of ', fct_year, fct_mon, fct_day, ' is incomplete'
        stop
     end if

!----------------------------------------------------------------------
!    Get varibles of historical forecast from pre-processed binary file
!----------------------------------------------------------------------
     allocate( hist_single(grid_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: hist_single'
        stop
     end if

     samp_sgl_num = samp_sgl_len / samp_interval
     samp_max = (samp_sgl_num * 2 + 1) * (end_year - start_year + 1)
     allocate( hist_valid(grid_num, samp_max), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: hist_valid'
        stop
     end if

     allocate( samp_dates_valid(3, samp_max), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: samp_dates_valid'
        stop
     end if

     write( 6, '(a)' ) 'STEP 2: acquire historical forecast data ...'
     samp_num = 0
     do iyear = start_year, end_year
        do isamp = -samp_sgl_num, samp_sgl_num

           incre = isamp * samp_interval
           call new_date(iyear, fct_mon, fct_day, incre, year, mon, day)
           if ( year < start_year .or. year > end_year ) cycle
           call get_cfs_field(input_dir_hist, path_separator, var_num, var_names(1:var_num), year, mon, day, &
                              lead_days, ave_days, domain_west, domain_south, nlon, nlat, hist_single, ierr)
           if ( ierr /= 0 ) then
              write( 6, '(a, i4, 1x, i2, 1x, i2, a)' ) 'WARNING: forecast data of ', year, mon, day, ' is incomplete'
              cycle
           else
              samp_num = samp_num + 1
              hist_valid(:, samp_num) = hist_single
              samp_dates_valid(1, samp_num) = year
              samp_dates_valid(2, samp_num) = mon
              samp_dates_valid(3, samp_num) = day
           end if

        end do
     end do
     write( 6, '(a, i)' ) 'Valid number of historical forecast samples is ', samp_num

     if ( samp_num < anal_fct_num ) anal_fct_num = samp_num

     deallocate( hist_single )

     allocate( hist(grid_num, samp_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: hist'
        stop
     end if

     allocate( samp_dates(3, samp_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: samp_dates'
        stop
     end if

     hist(:, :) = hist_valid(:, 1:samp_num)
     samp_dates(:, :) = samp_dates_valid(:, 1:samp_num)

     deallocate( hist_valid )
     deallocate( samp_dates_valid )

!----------------------------------------------------------------------
!    Normalize current forecast data using mean and standard deviation
!    of historical forecast samples
!----------------------------------------------------------------------
     write( 6, '(a)' ) 'STEP 3: normalize current forecast data ...'
     allocate( mean(grid_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: mean'
        stop
     end if

     allocate( sdev(grid_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: sdev'
        stop
     end if

     do igrid = 1, grid_num
        call statistic(samp_num, hist(igrid, :), mean(igrid), sdev(igrid))
     end do

     abnorm_num = count(sdev <= minima)
     if ( abnorm_num /= 0 ) then
        write( 6, '(a)' ) 'ERROR: Undefind value found in matrix standrad deviation'
        stop
     end if
     curt = (curt - mean) / sdev

     deallocate( mean, sdev )

!----------------------------------------------------------------------
!    Do Multivariate EOF decomposition to historical forecast to obtain
!    eigenvalues and eigenvectors
!----------------------------------------------------------------------
     write( 6, '(a)' ) 'STEP 4: do EOF decomposition ...'
     min_dim = min(grid_num, samp_num)

     allocate( eigenvalue(min_dim), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: eigenvalue'
        stop
     end if

     allocate( eigenvector(grid_num, min_dim), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: eigenvector'
        stop
     end if

     allocate( pc(min_dim, samp_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: pc'
        stop
     end if

     call eof(grid_num, samp_num, min_dim, trans_opt, undef, hist, eigenvalue, eigenvector, pc, ierr)

     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: EOF decomposition can NOT be rightly done'
        stop
     end if

     deallocate ( hist )

!----------------------------------------------------------------------
!    Seek the number of the first several leading EOFs whose cumulative
!    variance contribute rate exceed a certain level
!----------------------------------------------------------------------
     write( 6, '(a)' ) 'STEP 5: seek suitable number of EOFs ...'
     eigenvalue_sum = sum(eigenvalue)

     lead_num = 10
     write( 6, '(a, i2, a)' ) 'The cumulative variance contribute rates of the first ', lead_num, ' leading EOFs are: '
     do inum = 1, lead_num
        write( 6, '(a, i2, 2x, f6.3, a)' ) '        ', inum, 100.0 * sum(eigenvalue(1:inum)) / eigenvalue_sum, '%'
     end do

     deallocate( eigenvalue )

     allocate( proj_coeff_hist(lead_num, samp_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: proj_coeff_hist'
        stop
     end if

     allocate( proj_coeff_curt(lead_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: proj_coeff_curt'
        stop
     end if

     proj_coeff_hist(:, :) = pc(1:lead_num, :)

     deallocate( pc )

!----------------------------------------------------------------------
!    Calculate projection coefficients using dot product of current
!    forecast and selected eigenvectors after current forecast is
!    normalized
!----------------------------------------------------------------------
     write( 6, '(a)' ) 'STEP 6: calculate projection coefficients ...'
     do inum = 1, lead_num
        proj_coeff_curt(inum) = dot_product(eigenvector(:, inum), curt(:))
     end do

     deallocate( eigenvector )
     deallocate( curt )

!----------------------------------------------------------------------
!    Calculate similarity coefficients between current forecast and
!    historical forecasts
!----------------------------------------------------------------------
     write( 6, '(a)' ) 'STEP 7: calculate similarity coefficients ...'
     allocate( similar_val(samp_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: similar_val'
        stop
     end if

     do isamp = 1, samp_num
        call simi_disc_val(lead_num, proj_coeff_curt, proj_coeff_hist(:, isamp), similar_val(isamp))
     end do

     deallocate( proj_coeff_curt )
     deallocate( proj_coeff_hist )

!----------------------------------------------------------------------
!    Reorder similarity coefficients in ascending order to obtain specific
!    dates of the most similar historical forecasts
!----------------------------------------------------------------------
     write( 6, '(a)' ) 'STEP 8: reorder similarity coefficients ...'
     call sort(samp_num, 3, similar_val, samp_dates)
     write( 6, '(a, i3, a)' ) 'The first ', anal_fct_num, ' most similar discrete values are: '
     do isamp = 1, anal_fct_num
        write( 6, '(a, i2, 1x, f)' ) '        ', isamp, similar_val(isamp)
     end do
     write( 6, '(a, i3, a)' ) 'The first ', anal_fct_num, ' most similar forecasts dates are: '
     do isamp = 1, anal_fct_num
        write( 6, '(a, i2, 2x, i4, 1x, i2, 1x, i2)' ) '        ', isamp, samp_dates(1, isamp), samp_dates(2, isamp), samp_dates(3, isamp)
     end do
     deallocate( similar_val )

!----------------------------------------------------------------------
!    Calculate precipitation anomaly percentage of the most similar
!    historical forecasts
!----------------------------------------------------------------------
     write( 6, '(a)' ) 'STEP 9: calculate precipitation anomaly percentage ...'
     allocate( prec_fct(stn_num, anal_fct_num), stat = ierr )
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT allocate memory space to varible: prec_fct'
        stop
     end if

     do isamp = 1, anal_fct_num
        call get_prec_anom(input_dir_prec, path_separator, stn_num, samp_dates(1, isamp), samp_dates(2, isamp), &
                           samp_dates(3, isamp), lead_days, ave_days, prec_fct(:, isamp), ierr)
        if ( ierr /= 0 ) then
           write( 6, '(a)' ) 'ERROR: can NOT compute precipitation anomaly percentage'
           stop
        end if
     end do

     deallocate( samp_dates )

!----------------------------------------------------------------------
!    Output the forecast results for plotting using NCL software
!---------------------------------------------------------------------
     write( 6, '(a)' ) 'STEP 10: output forecast results ...'
     call write_out(output_dir, path_separator, fct_date, lead_days, ave_days, stn_num, anal_fct_num, prec_fct, ierr)
     if ( ierr /= 0 ) then
        write( 6, '(a)' ) 'ERROR: can NOT write results to ASCII file'
        stop
     end if

     deallocate( prec_fct )
     write( 6, '(a)' ) '=============== SUCCESSFUL TERMINATION OF DYNAMICAL ANALOGUE PREDICTION ==============='
     call cpu_time(time_end)
     write( 6, '(a, f9.3, a)' ) '=============== compute time (seconds)      =   ', time_end - time_begin, ' ==============='

  end program dynamical_analogue

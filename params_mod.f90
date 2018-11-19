!***********************************************************************
!  DESCRIPTION:
!       Module to load some necessary parameters to control the whole
!       behavior of dynamical analogue method
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Li Zhenkun, SCC
!       Revised 10/2018 by Li Zhenkun, SCC
!
!***********************************************************************
  module params_mod

     implicit none

!----------------------------------------------------------------------
!    Parameter definitions for namelist file
!
!    fct_date        = the date when CFSv2 begin to forecast
!    lead_days       = leading days for dynamical analogue prediction
!    ave_days        = period length for avarage
!    domain_west     = west edge of interested area
!    domain_east     = east edge of interested area
!    domain_south    = south edge of interested area
!    domain_north    = north edge of interested area
!    start_year      = start year to search for analogous historical forecast
!    end_year        = end year to search for analogous historical forecast
!    samp_interval   = the interval of sampling (days)
!    samp_sgl_len    = the sampling length of single side (days)
!    anal_fct_num    = the most analogous forecast number
!    var_num         = varible number
!    var_names       = varible names
!    input_dir_curt  = directory of current CFSv2 forecast
!    input_dir_hist  = directory of historical CFSv2 forecast
!    input_dir_prec  = directory of situ precipitation observations
!    output_dir      = directory for output
!    path_separator  = path separator for cross-platform, that is "\" for Windows and "/" for Unix-like system
!
!----------------------------------------------------------------------
     integer, parameter                 ::  max_num = 10
     integer                            ::  fct_date
     integer                            ::  lead_days, ave_days
     integer                            ::  domain_west, domain_east, domain_south, domain_north
     integer                            ::  start_year, end_year
     integer                            ::  samp_interval
     integer                            ::  samp_sgl_len
     integer                            ::  anal_fct_num ! selected forecast number
     integer                            ::  var_num ! for further research purpose
     character(len=256)                 ::  var_names(max_num)
     character(len=256)                 ::  input_dir_curt, input_dir_hist, input_dir_prec, output_dir
     character(len=1)                   ::  path_separator

     public          ::  init_params   !initialize the parameters listed above

     contains

!----------------------------------------------------------------------
!    Public subroutine to initialize control parameters
!----------------------------------------------------------------------
     subroutine init_params(filename, ierr)

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        character(len=*), intent(in)       ::  filename   ! namelist file name
        integer, intent(out)               ::  ierr       ! error identifier

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        character(len=*), parameter        ::  subname = '(init_params)'
        integer                            ::  iunit = 255 ! file unit

        namelist /params/ fct_date, lead_days, ave_days, domain_west, domain_east, domain_south, domain_north, start_year, end_year, &
                          samp_interval, samp_sgl_len, anal_fct_num, var_num, var_names, input_dir_curt, input_dir_hist, input_dir_prec, &
                          output_dir, path_separator

!------------------------------------------------------------------
!    Open the namelist file and obtain parameter's value
!------------------------------------------------------------------
        open( unit = iunit, file = filename, status = 'old', action = 'read', err = 100 )

        read( iunit, params, err = 200 )

        ierr = 0
        return

100     write( 6, '(3a)' ) subname, 'can NOT read namelist file ', trim(filename)
        ierr = 1
        return

200     write( 6, '(3a)' ) subname, 'can NOT read namelist stanza: params  ',  trim(filename)
        ierr = 1
        close( iunit )

     end subroutine init_params

  end module params_mod

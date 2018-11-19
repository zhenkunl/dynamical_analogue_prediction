!***********************************************************************
!  DESCRIPTION:
!       Module for some useful utilities
!
!  REVISION  HISTORY:
!       Prototype 01/2015 by Li Zhenkun, SCC
!
!***********************************************************************
  module utils_mod

     implicit none

     public      ::  simi_disc_val
     public      ::  sort

     contains

!------------------------------------------------------------------
!    Public subroutine to compute Similar Discrete Values between
!    two vectors
!------------------------------------------------------------------
     subroutine simi_disc_val(n, x, y, value)

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  n
        real(4), intent(in)                ::  x(n), y(n)
        real(4), intent(out)               ::  value

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        integer                            ::  i
        real(4)                            ::  dx(n)
        real(4)                            ::  s, d, e

        do i = 1, n
           dx(i) = x(i) - y(i)
        end do

        e = 0.
        d = 0.
        do i = 1, n
           e = e + dx(i)
           d = d + abs(dx(i))
        end do
        e = e / real(n)
        d = d / real(n)

        s = 0.
        do i = 1, n
           s = s + abs(dx(i) - e)
        end do
        s = s / real(n)

        value = (s + d) / 2.0

     end subroutine simi_disc_val

!------------------------------------------------------------------
!    Public subroutine to reorder an array in ascending order and
!    its corresponding two-dimension array
!------------------------------------------------------------------
     subroutine sort(n, m, x, y)

        implicit none

!------------------------------------------------------------------
!    dummy arguments
!------------------------------------------------------------------
        integer, intent(in)                ::  n, m
        real(4), intent(inout)             ::  x(n)
        integer, intent(inout)             ::  y(m, n)

!------------------------------------------------------------------
!    local variables
!------------------------------------------------------------------
        integer                            ::  l, ir
        real(4)                            ::  rra
        integer                            ::  rrb(3)
        integer                            ::  i, j

        l = n/2 + 1
        ir = n
        do
           if ( l > 1 ) then
              l = l - 1
              rra = x(l)
              rrb(:) = y(:, l)
           else
              rra = x(ir)
              rrb(:) = y(:, ir)
              x(ir) = x(1)
              y(:, ir) = y(:, 1)
              ir = ir - 1
              if ( ir == 1 ) then
                 x(1) = rra
                 y(:, 1) = rrb(:)
                 return
              end if
           end if
           i = l
           j = l + l
           do while ( j <= ir )
              if ( j .lt. ir ) then
                 if ( x(j) < x(j+1) ) j = j + 1
              end if
              if ( rra < x(j) ) then
                 x(i) = x(j)
                 y(:, i) = y(:, j)
                 i = j
                 j = j + j
              else
                 j = ir + 1
              endif
           end do
           x(i) = rra
           y(:, i) = rrb(:)
        end do

     end subroutine sort

  end module utils_mod

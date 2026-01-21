module extinction_module
   USE constants
   implicit none
!   integer, parameter :: dp = selected_real_kind(15, 307)
  contains

  subroutine extcurve_obs(wave, R_V, NH_EBV, model, extinction_curves)
    implicit none
    ! Inputs
    real(dp), dimension(:) :: wave ! Wavelength vector in microns
    real(dp), optional :: R_V             ! Ratio of visual extinction to reddening
    real(dp), optional :: NH_EBV ! Gas-to-dust ratio, default 5.8e21
    character(len=*), optional :: model ! Model selection, default 'ODonnell94'

    ! Outputs
    real(dp), dimension(2, size(wave)), intent(out) :: extinction_curves

    ! Local variables
    real(dp), dimension(size(wave)) :: x, a, b
    real(dp), dimension(:), allocatable :: c1, c2
    real(dp) :: NH_EBV_default, R_V_default
    real(dp) :: y
    integer :: i

    ! Assign default values
    if (present(NH_EBV)) then
        NH_EBV_default = NH_EBV
    else
        NH_EBV_default = 5.8d21
    end if

    if (.not. present(model)) then
        model = trim('ODonnell94')
    end if

    if (.not. present(R_V)) then
        R_V_default = 4.d0
    else
        R_V_default = R_V
    end if

    ! Convert wavelength to inverse microns
    x = 1.0d0 / wave

    ! Initialize a and b arrays
    a = 0.0d0
    b = 0.0d0

    ! Far-Infrared (x < 0.3)
    do i = 1, size(wave)
      if (x(i) < 0.3d0) then
        a(i) = 0.574d0 * x(i)**1.61d0
        b(i) = -0.527d0 * x(i)**1.61d0
      end if
    end do

    ! Infrared (0.3 <= x < 1.1)
    do i = 1, size(wave)
      if (x(i) >= 0.3d0 .and. x(i) < 1.1d0) then
        a(i) = 0.574d0 * x(i)**1.61d0
        b(i) = -0.527d0 * x(i)**1.61d0
      end if
    end do

    ! Optical/NIR (1.1 <= x < 3.3)
    do i = 1, size(wave)
      if (x(i) >= 1.1d0 .and. x(i) < 3.3d0) then
        y = x(i) - 1.82d0

        if (trim(model) == 'CCM89') then
          c1 = (/ 1.0d0, 0.17699d0, -0.50447d0, -0.02427d0, 0.72085d0, &
                0.01979d0, -0.77530d0, 0.32999d0 /)
          c2 = (/ 0.0d0, 1.41338d0, 2.28305d0, 1.07233d0, -5.38434d0, &
                -0.62251d0, 5.30260d0, -2.09002d0 /)
        else
          c1 = (/ 1.0d0, 0.104d0, -0.609d0, 0.701d0, 1.137d0, &
                -1.718d0, -0.827d0, 1.647d0, -0.505d0 /)
          c2 = (/ 0.0d0, 1.952d0, 2.908d0, -3.989d0, -7.985d0, &
                11.102d0, 5.491d0, -10.805d0, 3.347d0 /)
        end if

        a(i) = poly(c1, y)
        b(i) = poly(c2, y)
      end if
    end do

    ! Mid-UV (3.3 <= x < 8.0)
    do i = 1, size(wave)
      if (x(i) >= 3.3d0 .and. x(i) < 8.0d0) then
        y = x(i)

        a(i) = 1.752d0 - 0.316d0 * y - (0.104d0 / ((y - 4.67d0)**2 + 0.341d0))
        b(i) = -3.090d0 + 1.825d0 * y + (1.206d0 / ((y - 4.62d0)**2 + 0.263d0))
      end if
    end do

    ! Far-UV (8.0 <= x <= 11.0)
    do i = 1, size(wave)
      if (x(i) >= 8.0d0 .and. x(i) <= 11.0d0) then
        y = x(i) - 8.0d0

        c1 = (/ -1.073d0, -0.628d0, 0.137d0, -0.070d0 /)
        c2 = (/ 13.670d0, 4.257d0, -0.420d0, 0.374d0 /)

        a(i) = poly(c1, y)
        b(i) = poly(c2, y)
      end if
    end do

    ! Compute A_lambda/AV and A_lambda/NH
    extinction_curves(1, :) = a + b / R_V_default
    extinction_curves(2, :) = (a + b / R_V_default) * (R_V_default / NH_EBV_default)

  end subroutine extcurve_obs

  function poly(coeff, x) result(value)
    implicit none
    real(8), dimension(:), intent(in) :: coeff
    real(8), intent(in) :: x
    real(8) :: value
    integer :: i

    value = 0.0d0
    do i = 1, size(coeff)
      value = value + coeff(i) * x**(i - 1)
    end do
  end function poly

end module extinction_module

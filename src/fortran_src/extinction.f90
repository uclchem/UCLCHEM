module extinction_module
   use constants, only: dp
   implicit none

   private
   public :: extcurve_obs

  contains

  pure subroutine extcurve_obs(wave, R_V, NH_EBV, model, extinction_curves)
    real(dp), intent(in), dimension(:) :: wave       ! Wavelength vector in microns
    real(dp), intent(in), optional :: R_V            ! Ratio of visual extinction to reddening
    real(dp), intent(in), optional :: NH_EBV         ! Gas-to-dust ratio, default 5.8e21
    character(len=*), intent(inout), optional :: model  ! Model selection, default 'ODonnell94'

    real(dp), dimension(2, size(wave)), intent(out) :: extinction_curves

    real(dp), dimension(size(wave)) :: x, a, b
    real(dp), dimension(:), allocatable :: c1, c2
    real(dp) :: NH_EBV_default, R_V_default
    real(dp) :: y
    integer :: i

    ! Assign default values
    if (present(NH_EBV)) then
        NH_EBV_default = NH_EBV
    else
        NH_EBV_default = 5.8e21_dp
    end if

    if (.not. present(model)) then
        model = trim("ODonnell94")
    end if

    if (.not. present(R_V)) then
        R_V_default = 4.0_dp
    else
        R_V_default = R_V
    end if

    ! Convert wavelength to inverse microns
    x = 1.0_dp / wave

    ! Initialize a and b arrays
    a = 0.0_dp
    b = 0.0_dp

    ! Far-Infrared (x < 0.3)
    do i = 1, size(wave)
      if (x(i) < 0.3_dp) then
        a(i) = 0.574_dp * x(i)**1.61_dp
        b(i) = -0.527_dp * x(i)**1.61_dp
      end if
    end do

    ! Infrared (0.3 <= x < 1.1)
    do i = 1, size(wave)
      if (x(i) >= 0.3_dp .and. x(i) < 1.1_dp) then
        a(i) = 0.574_dp * x(i)**1.61_dp
        b(i) = -0.527_dp * x(i)**1.61_dp
      end if
    end do

    ! Optical/NIR (1.1 <= x < 3.3)
    do i = 1, size(wave)
      if (x(i) >= 1.1_dp .and. x(i) < 3.3_dp) then
        y = x(i) - 1.82_dp

        if (trim(model) == "CCM89") then
          c1 = (/ 1.0_dp, 0.17699_dp, -0.50447_dp, -0.02427_dp, 0.72085_dp, &
                0.01979_dp, -0.77530_dp, 0.32999_dp /)
          c2 = (/ 0.0_dp, 1.41338_dp, 2.28305_dp, 1.07233_dp, -5.38434_dp, &
                -0.62251_dp, 5.30260_dp, -2.09002_dp /)
        else
          c1 = (/ 1.0_dp, 0.104_dp, -0.609_dp, 0.701_dp, 1.137_dp, &
                -1.718_dp, -0.827_dp, 1.647_dp, -0.505_dp /)
          c2 = (/ 0.0_dp, 1.952_dp, 2.908_dp, -3.989_dp, -7.985_dp, &
                11.102_dp, 5.491_dp, -10.805_dp, 3.347_dp /)
        end if

        a(i) = poly(c1, y)
        b(i) = poly(c2, y)
      end if
    end do

    ! Mid-UV (3.3 <= x < 8.0)
    do i = 1, size(wave)
      if (x(i) >= 3.3_dp .and. x(i) < 8.0_dp) then
        y = x(i)

        a(i) = 1.752_dp - 0.316_dp * y - (0.104_dp / ((y - 4.67_dp)**2 + 0.341_dp))
        b(i) = -3.090_dp + 1.825_dp * y + (1.206_dp / ((y - 4.62_dp)**2 + 0.263_dp))
      end if
    end do

    ! Far-UV (8.0 <= x <= 11.0)
    do i = 1, size(wave)
      if (x(i) >= 8.0_dp .and. x(i) <= 11.0_dp) then
        y = x(i) - 8.0_dp

        c1 = (/ -1.073_dp, -0.628_dp, 0.137_dp, -0.070_dp /)
        c2 = (/ 13.670_dp, 4.257_dp, -0.420_dp, 0.374_dp /)

        a(i) = poly(c1, y)
        b(i) = poly(c2, y)
      end if
    end do

    ! Compute A_lambda/AV and A_lambda/NH
    extinction_curves(1, :) = a + b / R_V_default
    extinction_curves(2, :) = (a + b / R_V_default) * (R_V_default / NH_EBV_default)

  end subroutine extcurve_obs

  pure function poly(coeff, x) result(value)
    real(dp), dimension(:), intent(in) :: coeff
    real(dp), intent(in) :: x
    real(dp) :: value
    integer :: i

    value = 0.0_dp
    do i = 1, size(coeff)
      value = value + coeff(i) * x**(i - 1)
    end do
  end function poly

end module extinction_module

!> For accumulating and writing time means of surface and flux fields.
!
!  The fields are those of the 2-D time-mean output of SPEEDY ver41. They are
!  accumulated every time step and written at the end of each averaging period
!  (monthly or every nsteps_mean time steps) to mean_yyyymmddhhmm.nc, named after
!  the start of the period.
module mean_output
    use types, only: p, sp
    use params
    use netcdf
    use input_output, only: check

    implicit none

    private
    public update_means

    integer, parameter :: nfields = 29

    character(len=6), parameter :: names(nfields) = (/ character(len=6) :: &
        & 'ps', 'mslp', 'ts', 'tskin', 'soilw', 'alb', 'u0', 'v0', 't0', 'rh0', &
        & 'cloudc', 'clstr', 'cltop', 'prtop', 'lst', 'sst', 'ssta', &
        & 'precls', 'precnv', 'evap', 'ustr', 'vstr', 'tsr', 'olr', 'ssr', 'slr', 'shf', &
        & 'lshf', 'sshf' /)

    character(len=60), parameter :: long_names(nfields) = (/ character(len=60) :: &
        & 'surface pressure', &
        & 'mean-sea-level pressure', &
        & 'surface temperature', &
        & 'skin temperature', &
        & 'soil wetness availability', &
        & 'surface albedo', &
        & 'near-surface u-wind', &
        & 'near-surface v-wind', &
        & 'near-surface air temperature', &
        & 'relative humidity at the lowest model level', &
        & 'cloud cover (deep clouds)', &
        & 'cloud cover (stratiform clouds)', &
        & 'pressure at cloud top', &
        & 'highest precipitation level index', &
        & 'land-surface temperature', &
        & 'sea-surface temperature', &
        & 'SST anomaly w.r.t. observed climatology', &
        & 'large-scale precipitation', &
        & 'convective precipitation', &
        & 'evaporation', &
        & 'u-stress (downward)', &
        & 'v-stress (downward)', &
        & 'top-of-atmosphere shortwave radiation (downward)', &
        & 'outgoing longwave radiation (upward)', &
        & 'surface shortwave radiation (net downward)', &
        & 'surface longwave radiation (net upward)', &
        & 'sensible heat flux (upward)', &
        & 'net heat flux into land surface (downward, x land fraction)', &
        & 'net heat flux into sea surface (downward, x sea fraction)' /)

    character(len=6), parameter :: units(nfields) = (/ character(len=6) :: &
        & 'Pa', 'Pa', 'K', 'K', '1', '1', 'm/s', 'm/s', 'K', '1', &
        & '1', '1', 'Pa', '1', 'K', 'K', 'K', &
        & 'mm/day', 'mm/day', 'mm/day', 'N/m^2', 'N/m^2', 'W/m^2', 'W/m^2', 'W/m^2', 'W/m^2', &
        & 'W/m^2', 'W/m^2', 'W/m^2' /)

    ! Indices of the fields defined over land or sea only
    integer, parameter :: ilst = 15, isst = 16, issta = 17

    real(p) :: sums(ix,il,nfields) = 0.0 !! Accumulated fields
    integer :: nsum = 0                  !! Number of accumulated time steps
    integer :: period_start_step = 0     !! Number of time steps at the start of the period
    integer :: period_start(5) = -1      !! Start date of the period (year, month, day, hour, minute)

contains
    !> Adds the fields of the last time step to the time means, and writes the
    !  means at the end of an averaging period.
    subroutine update_means(steps_done)
        use date, only: model_datetime, start_datetime

        integer, intent(in) :: steps_done !! Number of time steps done since the start of the run

        logical :: end_of_period

        if (nsteps_mean == 0) return

        if (period_start(1) < 0) then
            period_start = (/ start_datetime%year, start_datetime%month, start_datetime%day, &
                & start_datetime%hour, start_datetime%minute /)
        end if

        call accumulate

        ! The model datetime is the end of the last time step
        if (nsteps_mean > 0) then
            end_of_period = mod(steps_done, nsteps_mean) == 0
        else
            end_of_period = model_datetime%day == 1 .and. &
                & model_datetime%hour*60 + model_datetime%minute < 24*60/nsteps
        end if

        if (end_of_period) then
            call write_means(steps_done)

            sums = 0.0
            nsum = 0
            period_start_step = steps_done
            period_start = (/ model_datetime%year, model_datetime%month, model_datetime%day, &
                & model_datetime%hour, model_datetime%minute /)
        end if
    end subroutine

    !> Adds the current surface and flux fields to the sums.
    subroutine accumulate
        use auxiliaries, only: psg, ts, tskin, u0, v0, t0, rh0, cloudc, clstr, cltop, prtop, &
            & precls, precnv, evap, ustr, vstr, tsr, olr, ssr, slr, shf, hfluxn
        use land_model, only: stl_am, soilw_am, fmask_l
        use sea_model, only: sst_am, sstan_am
        use mod_radcon, only: albsfc
        use boundaries, only: phis0
        use physical_constants, only: p0, grav, rgas

        real(p) :: tsg(ix,il)
        real(p) :: gam0

        ! Mean-sea-level pressure (as in SPEEDY ver41)
        gam0 = 0.006/grav
        tsg = 0.5*(t0 + max(255.0_p, min(295.0_p, t0)))

        sums(:,:,1) = sums(:,:,1) + p0*psg
        sums(:,:,2) = sums(:,:,2) + p0*psg*(1.0 + gam0*phis0/tsg)**(1.0/(rgas*gam0))
        sums(:,:,3) = sums(:,:,3) + ts
        sums(:,:,4) = sums(:,:,4) + tskin
        sums(:,:,5) = sums(:,:,5) + soilw_am
        sums(:,:,6) = sums(:,:,6) + albsfc
        sums(:,:,7) = sums(:,:,7) + u0
        sums(:,:,8) = sums(:,:,8) + v0
        sums(:,:,9) = sums(:,:,9) + t0
        sums(:,:,10) = sums(:,:,10) + rh0
        sums(:,:,11) = sums(:,:,11) + cloudc
        sums(:,:,12) = sums(:,:,12) + clstr
        sums(:,:,13) = sums(:,:,13) + p0*cltop
        sums(:,:,14) = sums(:,:,14) + prtop
        sums(:,:,15) = sums(:,:,15) + stl_am
        sums(:,:,16) = sums(:,:,16) + sst_am
        sums(:,:,17) = sums(:,:,17) + sstan_am

        ! Fluxes of water (g/(m^2 s) -> mm/day), momentum and energy
        sums(:,:,18) = sums(:,:,18) + precls*86.4
        sums(:,:,19) = sums(:,:,19) + precnv*86.4
        sums(:,:,20) = sums(:,:,20) + evap(:,:,3)*86.4
        sums(:,:,21) = sums(:,:,21) - ustr(:,:,3)
        sums(:,:,22) = sums(:,:,22) - vstr(:,:,3)
        sums(:,:,23) = sums(:,:,23) + tsr
        sums(:,:,24) = sums(:,:,24) + olr
        sums(:,:,25) = sums(:,:,25) + ssr
        sums(:,:,26) = sums(:,:,26) + slr
        sums(:,:,27) = sums(:,:,27) + shf(:,:,3)
        sums(:,:,28) = sums(:,:,28) + hfluxn(:,:,1)*fmask_l
        sums(:,:,29) = sums(:,:,29) + hfluxn(:,:,2)*(1.0 - fmask_l)

        nsum = nsum + 1
    end subroutine

    !> Writes the time means of the current averaging period to a NetCDF file.
    subroutine write_means(steps_done)
        use date, only: start_datetime
        use geometry, only: radang
        use land_model, only: fmask_l
        use sea_model, only: fmask_s

        integer, intent(in) :: steps_done !! Number of time steps done since the start of the run

        character(len=20) :: file_name = 'mean_yyyymmddhhmm.nc'
        character(len=32) :: time_template = 'hours since yyyy-mm-dd hh:mm:0.0'
        real(sp) :: means(ix,il), bounds(2)
        integer :: ncid, timedim, londim, latdim, nvdim, timevar, boundsvar, lonvar, latvar
        integer :: varids(nfields), n, k

        write (file_name(6:17),'(i4.4,4i2.2)') period_start

        write (time_template(13:16),'(i4.4)') start_datetime%year
        write (time_template(18:19),'(i2.2)') start_datetime%month
        write (time_template(21:22),'(i2.2)') start_datetime%day
        write (time_template(24:25),'(i2.2)') start_datetime%hour
        write (time_template(27:28),'(i2.2)') start_datetime%minute

        print '(A,A)', 'Write time-mean dataset: ', file_name

        call check(nf90_create(file_name, nf90_clobber, ncid))

        ! Time is the middle of the averaging period, with the period as bounds
        call check(nf90_def_dim(ncid, "time", nf90_unlimited, timedim))
        call check(nf90_def_dim(ncid, "nv", 2, nvdim))
        call check(nf90_def_var(ncid, "time", nf90_real4, timedim, timevar))
        call check(nf90_put_att(ncid, timevar, "units", time_template))
        call check(nf90_put_att(ncid, timevar, "bounds", "time_bnds"))
        call check(nf90_def_var(ncid, "time_bnds", nf90_real4, (/ nvdim, timedim /), boundsvar))

        call check(nf90_def_dim(ncid, "lon", ix, londim))
        call check(nf90_def_dim(ncid, "lat", il, latdim))
        call check(nf90_def_var(ncid, "lon", nf90_real4, londim, lonvar))
        call check(nf90_put_att(ncid, lonvar, "long_name", "longitude"))
        call check(nf90_def_var(ncid, "lat", nf90_real4, latdim, latvar))
        call check(nf90_put_att(ncid, latvar, "long_name", "latitude"))

        do n = 1, nfields
            call check(nf90_def_var(ncid, trim(names(n)), nf90_real4, &
                & (/ londim, latdim, timedim /), varids(n)))
            call check(nf90_put_att(ncid, varids(n), "long_name", trim(long_names(n))))
            call check(nf90_put_att(ncid, varids(n), "units", trim(units(n))))
            call check(nf90_put_att(ncid, varids(n), "cell_methods", "time: mean"))
            if (n == ilst .or. n == isst .or. n == issta) then
                call check(nf90_put_att(ncid, varids(n), "_FillValue", nf90_fill_real))
            end if
        end do

        call check(nf90_enddef(ncid))

        bounds = (/ period_start_step, steps_done /)*24.0/real(nsteps,sp)
        call check(nf90_put_var(ncid, timevar, (/ sum(bounds)/2.0 /), (/ 1 /)))
        call check(nf90_put_var(ncid, boundsvar, bounds, (/ 1, 1 /)))
        call check(nf90_put_var(ncid, lonvar, (/ (3.75*k, k = 0, ix-1) /)))
        call check(nf90_put_var(ncid, latvar, (/ (radang(k)*90.0/asin(1.0), k = 1, il) /)))

        do n = 1, nfields
            means = real(sums(:,:,n)/nsum, sp)

            ! Mask fields defined over land or sea only
            if (n == ilst) where (fmask_l == 0.0) means = nf90_fill_real
            if (n == isst .or. n == issta) where (fmask_s == 0.0) means = nf90_fill_real

            call check(nf90_put_var(ncid, varids(n), means, (/ 1, 1, 1 /)))
        end do

        call check(nf90_close(ncid))
    end subroutine
end module

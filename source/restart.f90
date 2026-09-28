!> For writing and reading restart files.
!
!  A restart file contains both time levels of the prognostic spectral variables
!  and the state of the land and sea/ice models, as in SPEEDY ver41. A run
!  restarted from a file written at 00:00 reproduces the continuous run exactly.
module restart
    use types, only: p
    use params
    use netcdf
    use input_output, only: check

    implicit none

    private
    public write_restart, read_restart

    integer, parameter :: spec3d(5) = (/ 2, mx, nx, kx, 2 /)      !! Shape of 3D spectral fields
    integer, parameter :: spec2d(4) = (/ 2, mx, nx, 2 /)          !! Shape of 2D spectral fields
    integer, parameter :: spectr(6) = (/ 2, mx, nx, kx, 2, ntr /) !! Shape of tracer fields

contains
    !> Writes the model state to the file restart_yyyymmddhhmm.nc.
    subroutine write_restart
        use date, only: model_datetime
        use prognostics, only: vor, div, t, ps, tr
        use land_model, only: stl_lm
        use sea_model, only: sst_om, tice_om, sice_om

        character(len=23) :: file_name = 'restart_yyyymmddhhmm.nc'
        integer :: ncid, ridim, mdim, ndim, levdim, stepdim, trdim, londim, latdim
        integer :: vorvar, divvar, tvar, psvar, trvar, stlvar, sstvar, ticevar, sicevar

        write (file_name(9:20),'(i4.4,4i2.2)') model_datetime%year, model_datetime%month, &
            & model_datetime%day, model_datetime%hour, model_datetime%minute

        print '(A,A)', 'Write restart file: ', file_name

        call check(nf90_create(file_name, nf90_clobber, ncid))

        ! Date of the model state
        call check(nf90_put_att(ncid, nf90_global, "year", model_datetime%year))
        call check(nf90_put_att(ncid, nf90_global, "month", model_datetime%month))
        call check(nf90_put_att(ncid, nf90_global, "day", model_datetime%day))
        call check(nf90_put_att(ncid, nf90_global, "hour", model_datetime%hour))
        call check(nf90_put_att(ncid, nf90_global, "minute", model_datetime%minute))

        ! Spectral fields are stored as real and imaginary parts (dimension ri)
        call check(nf90_def_dim(ncid, "ri", 2, ridim))
        call check(nf90_def_dim(ncid, "m", mx, mdim))
        call check(nf90_def_dim(ncid, "n", nx, ndim))
        call check(nf90_def_dim(ncid, "lev", kx, levdim))
        call check(nf90_def_dim(ncid, "time_level", 2, stepdim))
        call check(nf90_def_dim(ncid, "tracer", ntr, trdim))
        call check(nf90_def_dim(ncid, "lon", ix, londim))
        call check(nf90_def_dim(ncid, "lat", il, latdim))

        call check(nf90_def_var(ncid, "vor", nf90_double, &
            & (/ ridim, mdim, ndim, levdim, stepdim /), vorvar))
        call check(nf90_def_var(ncid, "div", nf90_double, &
            & (/ ridim, mdim, ndim, levdim, stepdim /), divvar))
        call check(nf90_def_var(ncid, "t", nf90_double, &
            & (/ ridim, mdim, ndim, levdim, stepdim /), tvar))
        call check(nf90_def_var(ncid, "ps", nf90_double, (/ ridim, mdim, ndim, stepdim /), psvar))
        call check(nf90_def_var(ncid, "tr", nf90_double, &
            & (/ ridim, mdim, ndim, levdim, stepdim, trdim /), trvar))
        call check(nf90_def_var(ncid, "stl_lm", nf90_double, (/ londim, latdim /), stlvar))
        call check(nf90_def_var(ncid, "sst_om", nf90_double, (/ londim, latdim /), sstvar))
        call check(nf90_def_var(ncid, "tice_om", nf90_double, (/ londim, latdim /), ticevar))
        call check(nf90_def_var(ncid, "sice_om", nf90_double, (/ londim, latdim /), sicevar))

        call check(nf90_enddef(ncid))

        call put_spectral(ncid, vorvar, spec3d, vor)
        call put_spectral(ncid, divvar, spec3d, div)
        call put_spectral(ncid, tvar, spec3d, t)
        call put_spectral(ncid, psvar, spec2d, ps)
        call put_spectral(ncid, trvar, spectr, tr)
        call check(nf90_put_var(ncid, stlvar, stl_lm))
        call check(nf90_put_var(ncid, sstvar, sst_om))
        call check(nf90_put_var(ncid, ticevar, tice_om))
        call check(nf90_put_var(ncid, sicevar, sice_om))

        call check(nf90_close(ncid))
    end subroutine

    !> Reads the model state from restart_file. The date of the restart file must
    !  be the start date of the run.
    subroutine read_restart
        use date, only: start_datetime
        use prognostics, only: vor, div, t, ps, tr
        use land_model, only: set_land_state
        use sea_model, only: set_sea_state

        integer :: ncid, year, month, day, hour, minute
        real(p), dimension(ix,il) :: stl, sst, tice, sice

        print '(A,A)', 'Read restart file: ', trim(restart_file)

        call check(nf90_open(trim(restart_file), nf90_nowrite, ncid))

        ! Check the date of the model state
        call check(nf90_get_att(ncid, nf90_global, "year", year))
        call check(nf90_get_att(ncid, nf90_global, "month", month))
        call check(nf90_get_att(ncid, nf90_global, "day", day))
        call check(nf90_get_att(ncid, nf90_global, "hour", hour))
        call check(nf90_get_att(ncid, nf90_global, "minute", minute))

        if (year /= start_datetime%year .or. month /= start_datetime%month .or. &
            & day /= start_datetime%day .or. hour /= start_datetime%hour .or. &
            & minute /= start_datetime%minute) then
            print '(A,I4.4,4I2.2)', 'Date of restart file (yyyymmddhhmm): ', &
                & year, month, day, hour, minute
            stop 'Restart file date does not match the start date in the namelist'
        end if

        ! Atmospheric state
        call get_spectral(ncid, "vor", spec3d, vor)
        call get_spectral(ncid, "div", spec3d, div)
        call get_spectral(ncid, "t", spec3d, t)
        call get_spectral(ncid, "ps", spec2d, ps)
        call get_spectral(ncid, "tr", spectr, tr)

        ! Land and sea/ice model states
        call get_grid(ncid, "stl_lm", stl)
        call get_grid(ncid, "sst_om", sst)
        call get_grid(ncid, "tice_om", tice)
        call get_grid(ncid, "sice_om", sice)

        call check(nf90_close(ncid))

        call set_land_state(stl)
        call set_sea_state(sst, tice, sice)
    end subroutine

    !> Writes a spectral field of any shape (passed by sequence association) as
    !  real and imaginary parts.
    subroutine put_spectral(ncid, varid, counts, field)
        integer, intent(in)    :: ncid, varid
        integer, intent(in)    :: counts(:)                   !! Shape of the variable in the file
        complex(p), intent(in) :: field(product(counts(2:)))

        real(p) :: parts(2,size(field))

        parts(1,:) = real(field, p)
        parts(2,:) = aimag(field)

        call check(nf90_put_var(ncid, varid, parts, count=counts))
    end subroutine

    !> Reads a spectral field of any shape (passed by sequence association) from
    !  its real and imaginary parts.
    subroutine get_spectral(ncid, name, counts, field)
        integer, intent(in)       :: ncid
        character(len=*), intent(in) :: name
        integer, intent(in)       :: counts(:)                   !! Shape of the variable in the file
        complex(p), intent(inout) :: field(product(counts(2:)))

        real(p) :: parts(2,size(field))
        integer :: varid

        call check(nf90_inq_varid(ncid, name, varid))
        call check(nf90_get_var(ncid, varid, parts, count=counts))

        field = cmplx(parts(1,:), parts(2,:), kind=p)
    end subroutine

    !> Reads a grid-point field.
    subroutine get_grid(ncid, name, field)
        integer, intent(in)          :: ncid
        character(len=*), intent(in) :: name
        real(p), intent(out)         :: field(ix,il)

        integer :: varid

        call check(nf90_inq_varid(ncid, name, varid))
        call check(nf90_get_var(ncid, varid, field))
    end subroutine
end module

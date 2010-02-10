
program simple

  use iso_c_binding

  real(kind=c_double) :: reltol_f, t_initial_f, t_final_f
  type(c_ptr) :: x_f, abstol_f
  integer(kind=c_int) :: neq
  real(kind=c_double), target :: x(1), abstol(1)

  interface
     subroutine do_run(neq, x_f, abstol_f, reltol_f, &
          t_initial_f, t_final_f) bind(c)
       use iso_c_binding
       integer(kind=c_int), value :: neq
       type(c_ptr), value :: x_f
       type(c_ptr), value :: abstol_f
       real(kind=c_double), value :: reltol_f
       real(kind=c_double), value :: t_initial_f
       real(kind=c_double), value :: t_final_f
     end subroutine do_run
  end interface

  neq = int(1, kind=c_int)
  x(1) = real(1d0, kind=c_double)
  x_f = c_loc(x)
  abstol(1) = real(1d-8, kind=c_double)
  abstol_f = c_loc(abstol)
  reltol_f = real(1d-8, kind=c_double)
  t_initial_f = real(0, kind=c_double)
  t_final_f = real(3, kind=c_double)
  call do_run(neq, x_f, abstol_f, reltol_f, t_initial_f, t_final_f)

  write(*,*) 'x = ', x(1)
  write(*,*) 'exp(-t_f) = ', exp(-t_final_f)

end program simple

  subroutine vf_f(neq, t_f, y_f, ydot_f) bind(c)
    use iso_c_binding
    integer(kind=c_int), value :: neq
    real(kind=c_double), value :: t_f
    type(c_ptr), value :: y_f
    type(c_ptr), value :: ydot_f

    real(kind=c_double), pointer :: y(:)
    real(kind=c_double), pointer :: ydot(:)

    call c_f_pointer(y_f, y, (/ 1 /))
    call c_f_pointer(ydot_f, ydot, (/ 1 /))

    ydot(1) = real(-1d0, kind=c_double) * y(1)

    write(*,*) 'vf_f: t, y, ydot = ', t_f, y(1), ydot(1)

  end subroutine vf_f

  subroutine jac_f(neq, t_f, y_f, dy_f) bind(c)
    use iso_c_binding
    integer(kind=c_int), value :: neq
    real(kind=c_double), value :: t_f
    type(c_ptr), value :: y_f
    type(c_ptr), value :: dy_f

    real(kind=c_double), pointer :: y(:)
    real(kind=c_double), pointer :: dy(:)

    call c_f_pointer(y_f, y, (/ 1 /))
    call c_f_pointer(dy_f, dy, (/ 1 /))

    dy(1) = real(-1d0, kind=c_double)

    write(*,*) 'jac_f: t, y, dy = ', t_f, y(1), dy(1)

  end subroutine jac_f

  subroutine jtimes_f(neq, t_f, y_f, v_f, Jv_f) bind(c)
    use iso_c_binding
    integer(kind=c_int), value :: neq
    real(kind=c_double), value :: t_f
    type(c_ptr), value :: y_f
    type(c_ptr), value :: v_f
    type(c_ptr), value :: Jv_f

    real(kind=c_double), pointer :: y(:)
    real(kind=c_double), pointer :: v(:)
    real(kind=c_double), pointer :: Jv(:)

    call c_f_pointer(y_f, y, (/ 1 /))
    call c_f_pointer(v_f, v, (/ 1 /))
    call c_f_pointer(Jv_f, Jv, (/ 1 /))

    Jv(1) = real(-1d0, kind=c_double) * v(1)

    write(*,*) 'jtimes_f: t, y, v, Jv = ', t_f, y(1), v(1), Jv(1)

  end subroutine jtimes_f

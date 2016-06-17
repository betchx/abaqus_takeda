
!
subroutine K_UEXTERNALDB(LOP,LRESTART, CTIME, DTIME,KSTEP,KINC)
  use debug
  implicit none
  integer, intent(in)  :: LOP
  integer, intent(in)  :: LRESTART
  real, intent(in)  :: CTIME(2)
  real, intent(in)  :: DTIME
  integer, intent(in)  :: KSTEP
  integer, intent(in)  :: KINC

  character(256) :: outdir, jobname
  integer :: lenoutdir, lenjobname

  integer, save :: last_inc = 0


  select case (LOP)
  case(0)
    !begining of analysis
    call getoutdir(outdir, lenoutdir)
    call getjobname(jobname, lenjobname)
    call init_debug(outdir(1:lenoutdir)//'/'//jobname(1:lenjobname)//'-debug.xml')
    call debug_on()
    call enter('Analysis')
  case(1)
    !begining of increment
    if(last_inc.ne.KINC) call enter('Increment', KINC)
  case(2)
    call leave('Increment', KINC)
    last_inc = 0
  case(3)
    if(last_inc.ne.0) call leave('Increment',last_inc)
    call leave('Analysis')
    call terminate_debug()
  end select
  last_inc = KINC
end subroutine


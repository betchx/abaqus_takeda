
module debug
    implicit none

    integer, parameter :: nio = (105)
    logical, private, save:: enable = (.false.)
    logical, private, save:: is_opened = (.false.)
    integer, private, save :: indent = (0)
    integer, private, save :: sub = (0) !< subsection number, 0 means omit
    integer, private, save :: sec = (1) !< section number

    interface enter
        module procedure enter1, enter_i
    end interface

    interface leave
        module procedure leave1, leave_i
    end interface

    interface trace
        module procedure trace1, trace2, trace3, trace_sis, trace_si, trace_sr
    end interface

    interface CheckPoint
        module procedure CheckPoint1, CheckPoint2, CheckPoint0
    end interface

    interface attr
        module procedure attr_i, attr_s, attr_r
    end interface

    interface to_s
      module procedure i2s,r2g,b2s
    end interface

contains

    subroutine debug_on()
        enable = is_opened
    end subroutine debug_on

    subroutine debug_off()
        enable = .false.
    end subroutine debug_off

    subroutine init_debug(basename)
        implicit none
        character(*), intent(in):: basename

        open(nio,file=basename,BUFFERED='no')
        is_opened = .true.
        call debug_on
        indent = 0
        call trace('<?xml version="1.0" encoding="Shift_JIS" ?>')
        call trace('<Debug>')
    end subroutine init_debug

    subroutine terminate_debug()
        implicit none
        if(enable) then
            if(indent .ne. 0) then
                write(nio,'(a)') 'DEBUG: XMLの対応が取れていません．要確認'
                write(nio,'(a,i0)') 'indent:',indent
                write(0,*) 'DEBUG: XMLの対応が取れていません．要確認'
                write(0,'(a,i0)') 'indent:',indent
            endif
            call trace('</Debug>')
            enable = .false.
            close(nio)
            is_opened = .false.
        endif
    end subroutine terminate_debug

    !> デバッグ時にはtxtをエラー出力に表示する
    subroutine msg(txt)
      implicit none
      character(*),intent(in):: txt
      if(enable) write(0,*) txt
      call called('msg' // attr("txt",txt))
    end subroutine msg

    !> txtをエラー出力に表示する
    subroutine echo(txt)
      implicit none
      character(*),intent(in):: txt
      write(0,*) txt
      call called('echo' // attr("txt",txt))
    end subroutine echo

    subroutine trace1(txt)
        implicit none
        character(*), intent(in):: txt
        integer i
        i = 0
        if(enable) then
            if( indent > 100 )then
                stop 'デバッグ対応メッセージにエラー'
            elseif( indent > 0 )then
            !        do i=1, indent
            !          write(nio,'(" ")',advance='no')
            !        enddo
            endif
            write(nio,'(a)') trim(txt)
            flush(nio)
        endif
    end subroutine trace1

    subroutine called(txt)
        implicit none
        character(*), intent(in)::txt
        call trace3('<',txt,'/>')
    end subroutine called
    subroutine called2(txt)
        implicit none
        character(*), intent(in)::txt
        call trace3('<',txt,'/>')
    end subroutine called2

    !  subroutine called_with(txt, n, arr)
    !    implicit none
    !    character(*),intent(in)::txt, arr(*)
    !    integer,intent(in) :: n
    !    integer:: i, L
    !    character,allocatable :: buf(:)
    !    L=len_trim(txt)+1
    !    do i=1, n
    !      L = L + len_trim(arr(2*i-1))
    !      L = L + 4 !  ="",
    !      L = L + len_trim(arr(2*i))
    !    enddo
    !    allocate( buf(L+2))
    !    write(buf, 77) txt,(trim(arr(i)),i=1,n)
    !    77 format(a,' ',(a,'="',a,'",'))
    !    call called(buf(1:(len_trim(buf)-1))
    !    deallocate(buf)
    !  end subroutine

    subroutine trace2(a1,a2)
        implicit none
        character(*), intent(in):: a1,a2
        call trace1( trim(a1) // trim(a2) )
    end subroutine trace2

    subroutine trace3(a1,a2,a3)
        implicit none
        character(*), intent(in):: a1,a2,a3
        call trace1( trim(a1) // trim(a2) // trim(a3) )
    end subroutine trace3

    subroutine trace_si(a1,i1)
        implicit none
        character(*),intent(in):: a1
        integer, intent(in)::i1
        call trace2(a1, i2s(i1))
    end subroutine trace_si

    subroutine trace_sis(a1, i1, a2)
        implicit none
        character(*),intent(in):: a1,a2
        integer, intent(in)::i1
        call trace3(a1, i2s(i1),a2)
    end subroutine trace_sis

    subroutine trace_sr(a1,r1)
        implicit none
        character(*),intent(in):: a1
        real, intent(in)::r1
        call trace2(a1, r2g(r1))
    end subroutine trace_sr

    subroutine enter1(tag)
        implicit none
        character(*),intent(in):: tag
        call trace3('<',tag,'>')
        indent = indent + 1
    end subroutine enter1

    !> @brief 整数情報つき
    subroutine enter_i(tag,i)
        implicit none
        character(*),intent(in):: tag
        integer,intent(in)::i
        call enter1(tag // i2s(i))
    end subroutine enter_i

    subroutine leave1(tag)
        implicit none
        character(*),intent(in):: tag
        indent = indent - 1
        if(indent < 0) indent = 0
        call trace3('</',tag,'>')
    end subroutine leave1

    subroutine leave_i(tag,i)
        implicit none
        character(*),intent(in)::tag
        integer, intent(in):: i
        call leave1( tag // i2s(i) )
    end subroutine leave_i

    pure function i2s(i)
        integer,intent(in) :: i
        character(20)::i2s
        write(i2s,'(i0)') i
    end function i2s

    pure function r2f(r)
        character(10) :: r2f
        real*8,intent(in):: r
        write(r2f,'(f10.4)') r
    end function r2f

    pure function r2e(r)
        character(15) :: r2e
        real(8),intent(in):: r
        write(r2e,'(e15.6)') r
    end function r2e

    pure function r2g(r)
        character(15) :: r2g
        real(8),intent(in):: r
        write(r2g,'(g15.6)') r
    end function r2g

    pure function b2s(b)
      character(15) :: b2s
      logical, intent(in):: b
      write(b2s,*) b
    end function b2s


    logical function is_debug()
        is_debug = enable
    end function is_debug

    subroutine CheckPoint1(a)
        integer,intent(in):: a
        sec = a
        sub = 0
        call CheckPoint0
    end subroutine CheckPoint1

    subroutine CheckPoint2(a,b)
        integer, intent(in):: a,b
        sec = a
        sub = b
        call CheckPoint0
    end subroutine CheckPoint2

    subroutine CheckPoint0()
        character(30) :: buf

        if(sub==0) then
            write(buf,'(i0)') sec
            sec = sec + 1
        else
            write(buf,'(i0,"-",i0)') sec, sub
            sub = sub + 1
        endif
        call echo('Check Point' // trim(buf))
        !call enter1('CheckPoint')
        !call trace1(buf)
        !call leave1('CheckPoint')
    end subroutine CheckPoint0

    pure integer function section()
        section = sec
    end function section
    pure integer function subsection()
        subsection = sub
    end function subsection

    pure integer function attr_len(a,b)
        character(*),intent(in)::a,b
        attr_len = len_trim(a) + len_trim(b) + 4
    end function attr_len

    function attr_s(key, val)
        character(*),intent(in):: key,val
        character(attr_len(key,val)) :: attr_s
        attr_s = " " // trim(adjustl(key)) // '="' // trim(adjustl(val)) // '"'
    end function attr_s

    function attr_i(key, val)
        character(*),intent(in):: key
        integer,intent(in)::val
        character(attr_len(key,i2s(val))) :: attr_i
        attr_i = attr_s(key, i2s(val) )
    end function attr_i

    function attr_r(key, val)
        character(*),intent(in):: key
        real(8),intent(in)::val
        character(attr_len(key,r2g(val))) :: attr_r
        attr_r = attr_s(key, r2g(val) )
    end function attr_r

    function tee(txt)
        character(*),intent(in)::txt
        character(len(txt)) ::tee
        call trace(txt)
        tee = txt
    end function tee
    subroutine breakpoint(txt)
      character(*),intent(in) :: txt
      call msg(txt)
      if(enable) pause
    end subroutine breakpoint
    subroutine bp(txt)
      character(*),intent(in) ::txt
      call breakpoint(txt)
    end subroutine bp
end module debug


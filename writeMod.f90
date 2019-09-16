module writeMod
  use paramMod
  implicit none

  contains

    subroutine openOutputFile(name_ad, isVertical)
      character (len=*) :: name_ad
      character (len=34),parameter :: path='/home/ecaas/decomposition/output/'
      logical, intent(in) :: isVertical
      open(unit=1,file = path//trim(name_ad)//"_pool.txt",   form="formatted", action="write", status='replace', iostat=ios)
      open(unit=15,file = path//trim(name_ad)//"_an.txt",   form="formatted", action="write", status='replace', iostat=ios)
      open(unit=2,file = path//trim(name_ad)//"_change.txt", form="formatted", action="write", status='replace', iostat=ios)


      if (isVertical) then
        open(unit=10,file = path//trim(name_ad)//"_vertical.txt", form="formatted", action="write", status='replace', iostat=ios)
        Write(unit=10, fmt=*) 'time,depth_level,pool_nr,net_vertical_transport, transport_upper, transport_lower'!header
        open(unit=3,file = path//trim(name_ad)//"_LITflux.txt", form="formatted", action="write", status='replace', iostat=ios)
        open(unit=4,file = path//trim(name_ad)//"_SAPSOMflux.txt", form="formatted", action="write", status='replace', iostat=ios)
        open(unit=7,file = path//trim(name_ad)//"_MYCSAPflux.txt", form="formatted", action="write", status='replace', iostat=ios)
        open(unit=8,file = path//trim(name_ad)//"_MYCSOMflux.txt", form="formatted", action="write", status='replace', iostat=ios)
        open(unit=9,file = path//trim(name_ad)//"_SOMflux.txt", form="formatted", action="write", status='replace', iostat=ios)
        write(unit=3, fmt=*) 'time, depth_level,LITmSAPr,LITmSAPk,LITsSAPr,LITsSAPk'
        write(unit=4, fmt=*) 'time, depth_level,SAPrSOMp,SAPrSOMa,SAPrSOMc,SAPkSOMp,SAPkSOMa,SAPkSOMc'
        write(unit=7, fmt=*) 'time, depth_level,EcM_SAPr,EcM_SAPk,ErM_SAPr,ErM_SAPk,AM_SAPr,AM_SAPk'
        write(unit=8, fmt=*) 'time,depth_level,EcM_SOMp, EcMSOMa,EcMSOMc, ErMSOMp,ErMSOMa,ErMSOMc,AMSOMp,AMSOMa,AMSOMc'
        Write(unit=9, fmt=*) 'time,depth_level,SOMaSAPr,SOMaSAPr,SOMpSOMa,SOMcSOMa'!header
      end if

      if( ios/=0) then
        write(6,*) 'Error opening file for writing'
        stop
      endif
      write(unit=1, fmt=*) 'time, depth_level, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header
      write(unit=15, fmt=*) 'time, depth_level, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header
      write(unit=2, fmt=*) 'time, depth_level, HR, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header

    end subroutine openOutputFile

    subroutine closeFiles(isVertical)
      logical, intent(in) :: isVertical
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=7)
      close(unit=8)
      close(unit=9)
      close(unit=15)
      if (isVertical) then
        close(unit=10)
      end if
    end subroutine closeFiles

end module writeMod

module writeMod
  use paramMod
  implicit none

  contains

    subroutine openOutputFile(name_ad, isVertical)
      character (len=*) :: name_ad
      logical, intent(in) :: isVertical
      open(unit=1,file = "pool_matrix_"//trim(name_ad)//".txt",   form="formatted", action="write", status='replace', iostat=ios)
      open(unit=15,file = "a_matrix_"//trim(name_ad)//".txt",   form="formatted", action="write", status='replace', iostat=ios)
      open(unit=2,file = "change_matrix_"//trim(name_ad)//".txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=3,file = "LITflux_"//trim(name_ad)//".txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=4,file = "SAPSOMflux_"//trim(name_ad)//".txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=7,file = "MYCSAPflux_"//trim(name_ad)//".txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=8,file = "MYCSOMflux_"//trim(name_ad)//".txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=9,file = "SOMflux_"//trim(name_ad)//".txt", form="formatted", action="write", status='replace', iostat=ios)

      if (isVertical) then
        open(unit=10,file = "vertical_"//trim(name_ad)//".txt", form="formatted", action="write", status='replace', iostat=ios)
        Write(unit=10, fmt=*) 'time,depth_level,pool_nr,net_vertical_transport, transport_upper, transport_lower'!header
      end if

      if( ios/=0) then
        write(6,*) 'Error opening file for writing'
        stop
      endif
      write(unit=1, fmt=*) 'time, depth_level, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header
      write(unit=15, fmt=*) 'time, depth_level, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header
      write(unit=2, fmt=*) 'time, depth_level, HR, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header
      write(unit=3, fmt=*) 'time, depth_level,LITmSAPr,LITmSAPk,LITsSAPr,LITsSAPk'
      write(unit=4, fmt=*) 'time, depth_level,SAPrSOMp,SAPrSOMa,SAPrSOMc,SAPkSOMp,SAPkSOMa,SAPkSOMc'
      write(unit=7, fmt=*) 'time, depth_level,EcM_SAPr,EcM_SAPk,ErM_SAPr,ErM_SAPk,AM_SAPr,AM_SAPk'
      write(unit=8, fmt=*) 'time,depth_level,EcM_SOMp, EcMSOMa,EcMSOMc, ErMSOMp,ErMSOMa,ErMSOMc,AMSOMp,AMSOMa,AMSOMc'
      Write(unit=9, fmt=*) 'time,depth_level,SOMaSAPr,SOMaSAPr,SOMpSOMa,SOMcSOMa'!header
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

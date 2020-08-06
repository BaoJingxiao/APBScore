!APBScore v0.1 (by Bao Jingxiao)

!读取计算参数文件
module ReadConfig
    !计算参数文件样式：(不包括行首的!)
    !#configures for calculate APBScore
    !#file path (should end with /):
    !/file/path/
    !#complex topology file name:
    !XXX.top
    !#complex restart file name (should be NetCDF format):
    !XXX_min.rst
    !#first receptor residue ID:
    !1
    !#last receptor residue ID:
    !100
    !#first ligand residue ID:
    !101
    !#last ligand residue ID:
    !101
    !#fitting coefficient data file path and name:
    !/file/path/FitCoeff.dat
    !#end
    !申明变量
    implicit none
    character*1024 :: TopFile                                                 !TOP文件路径与名称
    character*1024 :: RstFile                                                 !RST文件路径与名称
    integer*4 :: RecResStart,RecResStop,LigResStart,LigResStop                !受体起始、终止残基号，配体起始、终止残基号
    character*1024 :: FitFile                                                 !拟合系数文件路径与名称
    !子程序
    contains
        !读取计算参数文件
        subroutine ReadConfFileConfig(ConfFile)
            !申明变量
            implicit none
            character*1024 :: ConfFile                                        !计算参数文件路径与名称
            character*1024 :: FilePath                                        !文件所在文件夹的绝对路径
            integer*4 :: FileID                                               !读取文件时的ID号
            !读取参数
            FileID=10
            open(FileID,file=trim(ConfFile),status="old",action="read")
            read(FileID,*)
            read(FileID,*);read(FileID,'(a)')FilePath
            read(FileID,*);read(FileID,'(a)')TopFile
            read(FileID,*);read(FileID,'(a)')RstFile
            read(FileID,*);read(FileID,'(i10)')RecResStart
            read(FileID,*);read(FileID,'(i10)')RecResStop
            read(FileID,*);read(FileID,'(i10)')LigResStart
            read(FileID,*);read(FileID,'(i10)')LigResStop
            read(FileID,*);read(FileID,'(a)')FitFile
            close(FileID)
            TopFile=trim(FilePath)//trim(TopFile)
            RstFile=trim(FilePath)//trim(RstFile)
        end subroutine ReadConfFileConfig
end module ReadConfig

!读取TOP文件数据
module ReadTop
    !调用模块
    use ReadConfig
    !申明变量
    implicit none
    integer*4 :: NATOM                                                        !体系原子总数
    integer*4 :: NATOMTYPE                                                    !体系原子类型总数
    integer*4 :: NBONH                                                        !体系中有氢原子参与的键数量
    integer*4 :: NBONA                                                        !体系中无氢原子参与的键数量
    integer*4 :: NRES                                                         !体系残基总数
    character*4,allocatable,dimension(:) :: ATOMNAME                          !体系各原子的名称
    real*8,allocatable,dimension(:) :: CHARGE                                 !体系各原子的电荷
    integer*4,allocatable,dimension(:) :: ATOMICNUM                           !体系各原子的元素序号
    integer*4,allocatable,dimension(:) :: ATOMTYPE                            !体系各原子的原子类型序号
    integer*4,allocatable,dimension(:) :: NBPARM                              !体系各原子类型对间的非键作用参数编号
    character*4,allocatable,dimension(:) :: RESNAME                           !体系各残基的名称
    integer*4,allocatable,dimension(:) :: RESPOINT                            !体系各残基的起始原子序号
    real*8,allocatable,dimension(:) :: LJA                                    !体系各原子类型对间的范德华作用中的A值
    real*8,allocatable,dimension(:) :: LJB                                    !体系各原子类型对间的范德华作用中的B值
    integer*4,allocatable,dimension(:) :: BONDSH                              !体系各有氢原子参与的键的信息（(原子1序号-1)*3,(原子2序号-1)*3,键类型号）
    integer*4,allocatable,dimension(:) :: BONDSA                              !体系各无氢原子参与的键的信息（(原子1序号-1)*3,(原子2序号-1)*3,键类型号）
    integer*4 :: RecAtomStart,RecAtomStop,RecAtomTot                          !受体起始,终止原子序号，以及原子总数
    integer*4 :: LigAtomStart,LigAtomStop,LigAtomTot                          !配体起始,终止原子序号，以及原子总数
    !子程序
    contains
        !读取TOP文件参数指针
        subroutine ReadTopPointer()
            !申明变量
            implicit none
            integer*4 :: FileID                                               !读取文件时的ID号
            character*80 :: Line                                              !文件中的一行
            integer*4 :: FlagEOF                                              !用于判断是否读到文件尾
            integer*4 :: POINTERS(31)                                         !TOP文件参数的指针
            !读取TOP文件参数指针
            FileID=10
            open(FileID,file=trim(TopFile),status="old",action="read")
            do
                read(FileID,'(a80)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                !POINTERS
                if (trim(Line)=="%FLAG POINTERS") then
                    read(FileID,*)
                    read(FileID,'(3(10(i8),/),1(i8))')POINTERS
                    NATOM=POINTERS(1)
                    NATOMTYPE=POINTERS(2)
                    NBONH=POINTERS(3)
                    NBONA=POINTERS(4)
                    NRES=POINTERS(12)
                    exit
                end if
            end do
            close(FileID)
        end subroutine ReadTopPointer
        !读取TOP文件参数
        subroutine ReadTopData()
            !申明变量
            implicit none
            integer*4 :: FileID                                               !读取文件时的ID号
            character*80 :: Line                                              !文件中的一行
            integer*4 :: FlagEOF                                              !用于判断是否读到文件尾
            integer*4 :: n                                                    !循环变量
            !初始化数据数组
            allocate(ATOMNAME(NATOM));ATOMNAME=""
            allocate(CHARGE(NATOM));CHARGE=0
            allocate(ATOMICNUM(NATOM));ATOMICNUM=0
            allocate(ATOMTYPE(NATOM));ATOMTYPE=0
            allocate(NBPARM(NATOMTYPE*NATOMTYPE));NBPARM=0
            allocate(RESNAME(NRES));RESNAME=""
            allocate(RESPOINT(NRES));RESPOINT=0
            allocate(LJA(NATOMTYPE*(NATOMTYPE+1)/2));LJA=0
            allocate(LJB(NATOMTYPE*(NATOMTYPE+1)/2));LJB=0
            allocate(BONDSH(NBONH*3));BONDSH=0
            allocate(BONDSA(NBONA*3));BONDSA=0
            !读取TOP文件参数
            FileID=10
            open(FileID,file=trim(TopFile),status="old",action="read")
            do
                read(FileID,'(a80)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                select case (trim(Line))
                    !ATOM_NAME
                    case ("%FLAG ATOM_NAME")
                        read(FileID,*);read(FileID,'(20(a4))')(ATOMNAME(n),n=1,NATOM);cycle
                    !CHARGE
                    case ("%FLAG CHARGE")
                        read(FileID,*);read(FileID,'(5(e16.8))')(CHARGE(n),n=1,NATOM);cycle
                    !ATOMIC_NUMBER
                    case ("%FLAG ATOMIC_NUMBER")
                        read(FileID,*);read(FileID,'(10(i8))')(ATOMICNUM(n),n=1,NATOM);cycle
                    !ATOM_TYPE_INDEX
                    case ("%FLAG ATOM_TYPE_INDEX")
                        read(FileID,*);read(FileID,'(10(i8))')(ATOMTYPE(n),n=1,NATOM);cycle
                    !NONBONDED_PARM_INDEX
                    case ("%FLAG NONBONDED_PARM_INDEX")
                        read(FileID,*);read(FileID,'(10(i8))')(NBPARM(n),n=1,NATOMTYPE*NATOMTYPE);cycle
                    !RESIDUE_LABEL
                    case ("%FLAG RESIDUE_LABEL")
                        read(FileID,*);read(FileID,'(20(a4))')(RESNAME(n),n=1,NRES);cycle
                    !RESIDUE_POINTER
                    case ("%FLAG RESIDUE_POINTER")
                        read(FileID,*);read(FileID,'(10(i8))')(RESPOINT(n),n=1,NRES);cycle
                    !LENNARD_JONES_ACOEF
                    case ("%FLAG LENNARD_JONES_ACOEF")
                        read(FileID,*);read(FileID,'(5(e16.8))')(LJA(n),n=1,NATOMTYPE*(NATOMTYPE+1)/2);cycle
                    !LENNARD_JONES_BCOEF
                    case ("%FLAG LENNARD_JONES_BCOEF")
                        read(FileID,*);read(FileID,'(5(e16.8))')(LJB(n),n=1,NATOMTYPE*(NATOMTYPE+1)/2);cycle
                    !BONDS_INC_HYDROGEN
                    case ("%FLAG BONDS_INC_HYDROGEN")
                        read(FileID,*);read(FileID,'(10(i8))')(BONDSH(n),n=1,NBONH*3);cycle
                    !BONDS_WITHOUT_HYDROGEN
                    case ("%FLAG BONDS_WITHOUT_HYDROGEN")
                        read(FileID,*);read(FileID,'(10(i8))')(BONDSA(n),n=1,NBONA*3);cycle
                end select
            end do
            close(FileID)
        end subroutine ReadTopData
        !确定受体配体起始终止原子序号
        subroutine GetRecLigData()
            !申明变量
            implicit none
            !确定受体起始终止原子序号
            if (RecResStart<1) then
                write(*,'("ERROR: Receptor start residue ID should not be less than 1")')
                stop
            else
                RecAtomStart=RESPOINT(RecResStart)
            end if
            if (RecResStop<NRES) then
                RecAtomStop=RESPOINT(RecResStop+1)-1
            else if (RecResStop==NRES) then
                RecAtomStop=NATOM
            else
                write(*,'("ERROR: Receptor stop residue ID should not be more than the total residue number of the complex")')
                stop
            end if
            RecAtomTot=RecAtomStop-RecAtomStart+1
            if (RecAtomTot<=0) then
                write(*,'("ERROR: Receptor stop residue ID should not be less than its start residue ID")')
                stop
            end if
            !确定配体起始终止原子序号
            if (LigResStart<1) then
                write(*,'("ERROR: Ligand start residue ID should not be less than 1")')
                stop
            else
                LigAtomStart=RESPOINT(LigResStart)
            end if
            if (LigResStop<NRES) then
                LigAtomStop=RESPOINT(LigResStop+1)-1
            else if (LigResStop==NRES) then
                LigAtomStop=NATOM
            else
                write(*,'("ERROR: Ligand stop residue ID should not be more than the total residue number of the complex")')
                stop
            end if
            LigAtomTot=LigAtomStop-LigAtomStart+1
            if (LigAtomTot<=0) then
                write(*,'("ERROR: Ligand stop residue ID should not be less than its start residue ID")')
                stop
            end if
            if ((LigAtomStart<=RecAtomStop).and.(RecAtomStart<=LigAtomStop)) then
                write(*,'("ERROR: Residues should not be receptor and ligand at the same time")')
                stop
            end if
        end subroutine GetRecLigData
end module ReadTop

!读取RST文件数据
module ReadRst
    !调用模块
    use netcdf
    use ReadConfig
    use ReadTop
    !申明变量
    implicit none
    real*8,allocatable,dimension(:) :: Crd_All                                !体系中各原子的XYZ坐标
    !子程序
    contains
        !读取RST文件体系坐标数据
        subroutine ReadRstCrd()
            !申明变量
            implicit none
            integer*4 :: NCFlag                                               !NetCDF文件操作的返回值
            integer*4 :: NCFileID                                             !NC文件号
            integer*4 :: NCVarID_Crd                                          !NC文件中坐标数据的变量号
            integer*4 :: NCStart_Crd(3),NCCount_Crd(3)                        !读取NC中坐标数据的起始位置和数据数量
            integer*4 :: NCDimID_Atom                                         !NC文件中原子数维度的变量号
            integer*4 :: NCTotAtom                                            !RST文件中的原子总数
            character*4 :: Temp1                                              !临时变量
            character*8 :: Temp2,Temp3                                        !临时变量
            !初始化数据数组
            allocate(Crd_All(NATOM*3));Crd_All=0.0D0
            NCStart_Crd=(/1,1,1/)
            NCCount_Crd=(/3,NATOM,1/)
            !打开轨迹文件
            NCFlag=nf90_open(trim(RstFile),NF90_NOWRITE,NCFileID)
            if (NCFlag/=nf90_noerr) then
                write(*,'("ERROR: Cannot read NetCDF format restart file: ",a)')trim(RstFile)
                stop
            end if
            NCFlag=nf90_inq_varid(NCFileID,"coordinates",NCVarID_Crd)
            !确定轨迹文件中的原子总数
            NCFlag=nf90_inq_dimid(NCFileID,"atom",NCDimID_Atom)
            NCFlag=nf90_inquire_dimension(NCFileID,NCDimID_Atom,Temp1,NCTotAtom)
            !判断TOP与RST文件中的原子总数是否相等
            if (NCTotAtom/=NATOM) then
                write(Temp2,'(i8)')NATOM
                write(Temp3,'(i8)')NCTotAtom
                write(*,'("ERROR: Missmatch atom number between topology (",a,") and restart (",a,") file")')&
                     &trim(adjustl(Temp2)),trim(adjustl(Temp3))
                stop
            end if
            !读取轨迹坐标数据
            NCFlag=nf90_get_var(NCFileID,NCVarID_Crd,Crd_All,NCStart_Crd,NCCount_Crd)
            !关闭轨迹文件
            NCFlag=nf90_close(NCFileID)
        end subroutine ReadRstCrd
end module ReadRst

!读取拟合系数文件数据
module ReadFit
    !调用模块
    use ReadConfig
    !申明变量
    implicit none
    integer*4 :: RecElementTot,LigElementTot                                  !要计算范德华相互作用的受体和配体元素类型总数
    integer*4,allocatable,dimension(:) :: RecElementIndex,LigElementIndex     !要计算范德华相互作用的受体和配体元素序数
    integer*4 :: FlagElementRec(99)                                           !用于记录1到99号各元素是(>=1)否(0)会被用作范德华计算时的受体元素，以及在kVDW中对应的序号
    integer*4 :: FlagElementLig(99)                                           !用于记录1到99号各元素是(>=1)否(0)会被用作范德华计算时的配体元素，以及在kVDW中对应的序号
    real*8,allocatable,dimension(:,:) :: kVDW                                 !各元素类型对的范德华拟合系数(Rec,Lig)
    real*8 :: kELE,kHB,C0                                                     !复合物间总静电能和总氢键能的拟合系数，以及拟合常数C0
    !子程序
    contains
        !读取拟合系数文件数据
        subroutine ReadFitData()
            !申明变量
            implicit none
            integer*4 :: TotLine                                              !拟合系数文件总行数
            character*9,allocatable,dimension(:) :: CoeffName                 !拟合系数文件中各系数的名称
            real*8,allocatable,dimension(:) :: CoeffData                      !拟合系数文件中各系数的数值
            integer*4 :: FileID                                               !读取文件时的ID号
            character*22 :: Line                                              !文件中的一行
            integer*4 :: FlagEOF                                              !用于判断是否读到文件尾
            integer*4,allocatable,dimension(:,:) :: FlagkVDW                  !用于记录相应参数是否读取于拟合系数文件
            integer*4 :: FlagkELE,FlagkHB,FlagC0                              !用于记录相应参数是否读取于拟合系数文件
            integer*4 :: n,m                                                  !循环变量
            integer*4 :: Temp1,TempR,TempL                                    !临时变量
            character*9 :: TempName                                           !临时变量
            character*2 :: Temp2,Temp3                                        !临时变量
            !确定拟合系数文件总行数
            FileID=10
            TotLine=0
            open(FileID,file=trim(FitFile),status="old",action="read")
            do
                read(FileID,'(a22)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                if (trim(Line)/="") TotLine=TotLine+1
            end do
            close(FileID)
            !读取拟合系数文件各系数的名称和数值
            allocate(CoeffName(TotLine));CoeffName=""
            allocate(CoeffData(TotLine));CoeffData=0.0D0
            Temp1=0
            open(FileID,file=trim(FitFile),status="old",action="read")
            do
                read(FileID,'(a22)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                if (trim(Line)/="") then
                    Temp1=Temp1+1
                    read(Line(1:9),'(a9)')CoeffName(Temp1)
                    read(Line(11:22),'(f12.8)')CoeffData(Temp1)
                end if
            end do
            close(FileID)
            !确定要计算范德华作用的受体和配体元素类型
            FlagElementRec=0
            FlagElementLig=0
            do n=1,TotLine,1
                TempName=adjustl(CoeffName(n))
                if (TempName(1:3)=="VDW") then
                    do m=5,len(TempName),1
                        if (TempName(m:m)=="_") then
                            read(TempName(5:(m-1)),'(i2)')Temp1
                            FlagElementRec(Temp1)=1
                            read(TempName((m+1):len(TempName)),'(i2)')Temp1
                            FlagElementLig(Temp1)=1
                            exit
                        end if
                    end do
                end if
            end do
            !判断拟合系数文件中是否记录了常见受体（H,C,N,O,S,Zn）和配体原子(H,C,N,O,S,P,F,Cl,Br,I)的VDW系数
            if (FlagElementRec(1)==0) write(*,'("Warning: VDW interaction between receptor atom H to ligand will be ingored")')
            if (FlagElementRec(6)==0) write(*,'("Warning: VDW interaction between receptor atom C to ligand will be ingored")')
            if (FlagElementRec(7)==0) write(*,'("Warning: VDW interaction between receptor atom N to ligand will be ingored")')
            if (FlagElementRec(8)==0) write(*,'("Warning: VDW interaction between receptor atom O to ligand will be ingored")')
            if (FlagElementRec(16)==0) write(*,'("Warning: VDW interaction between receptor atom S to ligand will be ingored")')
            if (FlagElementRec(30)==0) write(*,'("Warning: VDW interaction between receptor atom Zn to ligand will be ingored")')
            if (FlagElementLig(1)==0) write(*,'("Warning: VDW interaction between ligand atom H to receptor will be ingored")')
            if (FlagElementLig(6)==0) write(*,'("Warning: VDW interaction between ligand atom C to receptor will be ingored")')
            if (FlagElementLig(7)==0) write(*,'("Warning: VDW interaction between ligand atom N to receptor will be ingored")')
            if (FlagElementLig(8)==0) write(*,'("Warning: VDW interaction between ligand atom O to receptor will be ingored")')
            if (FlagElementLig(9)==0) write(*,'("Warning: VDW interaction between ligand atom F to receptor will be ingored")')
            if (FlagElementLig(15)==0) write(*,'("Warning: VDW interaction between ligand atom P to receptor will be ingored")')
            if (FlagElementLig(16)==0) write(*,'("Warning: VDW interaction between ligand atom S to receptor will be ingored")')
            if (FlagElementLig(17)==0) write(*,'("Warning: VDW interaction between ligand atom Cl to receptor will be ingored")')
            if (FlagElementLig(35)==0) write(*,'("Warning: VDW interaction between ligand atom Br to receptor will be ingored")')
            if (FlagElementLig(53)==0) write(*,'("Warning: VDW interaction between ligand atom I to receptor will be ingored")')
            !确定要计算范德华作用的受体和配体元素类型总数以及各元素在kVDW中对应的编号
            RecElementTot=0
            LigElementTot=0
            do n=1,99,1
                if (FlagElementRec(n)>0) then
                    RecElementTot=RecElementTot+1
                    FlagElementRec(n)=RecElementTot
                end if
                if (FlagElementLig(n)>0) then
                    LigElementTot=LigElementTot+1
                    FlagElementLig(n)=LigElementTot
                end if
            end do
            !确定要计算范德华相互作用的受体和配体元素序数
            allocate(RecElementIndex(RecElementTot));RecElementIndex=0
            allocate(LigElementIndex(LigElementTot));LigElementIndex=0
            TempR=0
            TempL=0
            do n=1,99,1
                if (FlagElementRec(n)>0) then
                    TempR=TempR+1
                    RecElementIndex(TempR)=n
                end if
                if (FlagElementLig(n)>0) then
                    TempL=TempL+1
                    LigElementIndex(TempL)=n
                end if
            end do
            !确定各元素类型对的范德华拟合系数
            allocate(kVDW(RecElementTot,LigElementTot));kVDW=0.0D0
            allocate(FlagkVDW(RecElementTot,LigElementTot));FlagkVDW=0
            do n=1,TotLine,1
                TempName=adjustl(CoeffName(n))
                if (TempName(1:3)=="VDW") then
                    do m=5,len(TempName),1
                        if (TempName(m:m)=="_") then
                            read(TempName(5:(m-1)),'(i2)')TempR
                            read(TempName((m+1):len(TempName)),'(i2)')TempL
                            kVDW(FlagElementRec(TempR),FlagElementLig(TempL))=CoeffData(n)
                            FlagkVDW(FlagElementRec(TempR),FlagElementLig(TempL))=1
                            exit
                        end if
                    end do
                end if
            end do
            !判断各元素类型对的范德华拟合系数是否读取于拟合系数文件
            do n=1,RecElementTot,1
                do m=1,LigElementTot,1
                    if (FlagkVDW(n,m)==0) then
                        write(Temp2,"(i2)")RecElementIndex(n)
                        write(Temp3,"(i2)")LigElementIndex(m)
                        write(*,'("Warning: No coefficient for VDW_",a,"_",a," in data file")')&
                             &trim(adjustl(Temp2)),trim(adjustl(Temp3))
                    end if
                end do
            end do
            !确定复合物间总静电能和总氢键能的拟合系数，以及拟合常数C0
            kELE=0.0D0;FlagkELE=0
            kHB=0.0D0;FlagkHB=0
            C0=0.0D0;FlagC0=0
            do n=1,TotLine,1
                TempName=adjustl(CoeffName(n))
                select case (trim(TempName))
                    case ("ELE_com")
                        kELE=CoeffData(n)
                        FlagkELE=1
                    case ("HB_com")
                        kHB=CoeffData(n)
                        FlagkHB=1
                    case ("C0")
                        C0=CoeffData(n)
                        FlagC0=1
                end select
            end do
            !判断总静电能和总氢键能的拟合系数，以及拟合常数C0是否读取于拟合系数文件
            if (FlagkELE==0) write(*,'("Warning: No coefficient for ELE_com in data file")')
            if (FlagkHB==0) write(*,'("Warning: No coefficient for HB_com in data file")')
            if (FlagC0==0) write(*,'("Warning: No coefficient for C0 in data file")')
        end subroutine ReadFitData
end module ReadFit

!准备氢键能量计算中的各种参数
module PrepHB
    !调用模块
    use ReadTop
    !申明变量
    implicit none
    integer*4,allocatable,dimension(:,:) :: ATOMBONDS                         !体系中各原子的成键情况，(n,0)记录第n个原子成键的数量，(n,1-XXX)记录与第n个原子成键的原子序号
    real*8,allocatable,dimension(:) :: ATOMRADII                              !体系中各原子在XScore中的原子半径
    integer*4,allocatable,dimension(:) :: ATOMHBD                             !体系中各原子能提供的最大氢键数量
    integer*4,allocatable,dimension(:) :: ATOMHBA                             !体系中各原子能接受的最大氢键数量
    integer*4 :: MaxHBPerAtom                                                 !计算过程中每个原子最多形成的临时氢键数量
    !子程序
    contains
        !确定体系中各原子间的成键情况
        subroutine GetAtomBonds(MaxAtomBond)
            !申明变量
            implicit none
            integer*4 :: MaxAtomBond                                          !计算过程中每个原子最多能相连的原子的数量
            integer*4 :: n                                                    !循环变量
            integer*4 :: Temp1,Temp2                                          !临时变量
            !初始化数据数组
            allocate(ATOMBONDS(NATOM,0:MaxAtomBond));ATOMBONDS=0
            !确定体系中各原子间的成键情况
            do n=1,NBONH,1
                Temp1=BONDSH(n*3-2)/3+1
                Temp2=BONDSH(n*3-1)/3+1
                ATOMBONDS(Temp1,0)=ATOMBONDS(Temp1,0)+1
                if (ATOMBONDS(Temp1,0)>MaxAtomBond) ATOMBONDS(Temp1,0)=MaxAtomBond
                ATOMBONDS(Temp2,0)=ATOMBONDS(Temp2,0)+1
                if (ATOMBONDS(Temp2,0)>MaxAtomBond) ATOMBONDS(Temp2,0)=MaxAtomBond
                ATOMBONDS(Temp1,ATOMBONDS(Temp1,0))=Temp2
                ATOMBONDS(Temp2,ATOMBONDS(Temp2,0))=Temp1
            end do
            do n=1,NBONA,1
                Temp1=BONDSA(n*3-2)/3+1
                Temp2=BONDSA(n*3-1)/3+1
                ATOMBONDS(Temp1,0)=ATOMBONDS(Temp1,0)+1
                if (ATOMBONDS(Temp1,0)>MaxAtomBond) ATOMBONDS(Temp1,0)=MaxAtomBond
                ATOMBONDS(Temp2,0)=ATOMBONDS(Temp2,0)+1
                if (ATOMBONDS(Temp2,0)>MaxAtomBond) ATOMBONDS(Temp2,0)=MaxAtomBond
                ATOMBONDS(Temp1,ATOMBONDS(Temp1,0))=Temp2
                ATOMBONDS(Temp2,ATOMBONDS(Temp2,0))=Temp1
            end do
        end subroutine GetAtomBonds
        !确定体系中各原子在XScore中的原子半径
        subroutine GetAtomHBRadii()
            !申明变量
            implicit none
            integer*4 :: n                                                    !循环变量
            !初始化数据数组
            allocate(ATOMRADII(NATOM));ATOMRADII=0.0D0
            !确定体系中各原子在XScore中的原子半径
            do n=1,NATOM,1
                select case (ATOMICNUM(n))
                    case (6)                                                  !C
                        ATOMRADII(n)=1.9D0
                    case (7)                                                  !N
                        ATOMRADII(n)=1.8D0
                    case (8)                                                  !O
                        ATOMRADII(n)=1.7D0
                    case (9)                                                  !F
                        ATOMRADII(n)=1.5D0
                    case (15)                                                 !P
                        ATOMRADII(n)=2.1D0
                    case (16)                                                 !S
                        ATOMRADII(n)=2.0D0
                    case (17)                                                 !Cl
                        ATOMRADII(n)=1.8D0
                    case (35)                                                 !Br
                        ATOMRADII(n)=2.0D0
                    case (53)                                                 !I
                        ATOMRADII(n)=2.2D0
                    case (12)                                                 !Mg
                        ATOMRADII(n)=1.2D0
                    case (25)                                                 !Mn
                        ATOMRADII(n)=1.2D0
                    case (26)                                                 !Fe
                        ATOMRADII(n)=1.2D0
                    case (30)                                                 !Zn
                        ATOMRADII(n)=1.2D0
                end select
            end do
        end subroutine GetAtomHBRadii
        !设置计算过程中每个原子最多形成的临时氢键数量,并确定体系中各原子能提供和接受的最大氢键数量
        subroutine GetAtomHBMaxNum(MaxHBNum)
            !申明变量
            implicit none
            integer*4 :: MaxHBNum                                             !计算过程中每个原子最多形成的临时氢键数量
            integer*4 :: n,m                                                  !循环变量
            !设置计算过程中每个原子最多形成的临时氢键数量
            MaxHBPerAtom=MaxHBNum
            !初始化数据数组
            allocate(ATOMHBD(NATOM));ATOMHBD=0
            allocate(ATOMHBA(NATOM));ATOMHBA=0
            !确定体系中各原子能提供的最大氢键数量
            do n=1,NATOM,1
                select case (ATOMICNUM(n))
                    case (7,8)                                                !N,O
                        do m=1,ATOMBONDS(n,0),1
                            if (ATOMICNUM(ATOMBONDS(n,m))==1) then
                                ATOMHBD(n)=ATOMHBD(n)+1
                            end if
                        end do
                    case (12)                                                 !Mg
                        ATOMHBD(n)=6
                    case (25)                                                 !Mn
                        ATOMHBD(n)=6
                    case (26)                                                 !Fe
                        ATOMHBD(n)=6
                    case (30)                                                 !Zn
                        ATOMHBD(n)=4
                end select
            end do
            !确定体系中各原子能接受的最大氢键数量
            do n=1,NATOM,1
                select case (ATOMICNUM(n))
                    case (7)                                                  !N
                        if (ATOMBONDS(n,0)<4) then
                            ATOMHBA(n)=1
                        end if
                    case (8)                                                  !O
                        ATOMHBA(n)=2
                    case (16)                                                 !S
                        if (ATOMBONDS(n,0)<2) then
                            ATOMHBA(n)=2
                        end if
                end select
            end do
        end subroutine GetAtomHBMaxNum
end module PrepHB

!主程序
program main
    !调用模块
    use ReadConfig
    use ReadTop
    use ReadRst
    use ReadFit
    use PrepHB
    !申明变量
    implicit none
    character*1024 :: ConfFile                                                !计算参数文件路径与名称
    real*8 :: ScoreVDW,ScoreELE,ScoreHB                                       !打分函数中各分项的分值
    real*8 :: APBScore                                                        !打分函数的最终打分值
    !输入计算参数文件
    call getarg(1,ConfFile)
    !ConfFile="/media/Data/Dynamics/Test_APBScore/APBScore.config"
    !读取计算参数文件
    call ReadConfFileConfig(ConfFile)
    !读取TOP文件数据
    call ReadTopPointer()
    call ReadTopData()
    !确定受体配体起始终止原子序号
    call GetRecLigData()
    !读取RST文件体系坐标数据
    call ReadRstCrd()
    !读取拟合系数文件数据
    call ReadFitData()
    !计算打分函数中的范德华项分值和静电项分值
    ScoreVDW=0.0D0
    ScoreELE=0.0D0
    call CalScoreVDWELE(ScoreVDW,ScoreELE)
    !计算打分函数中的氢键项分值
    ScoreHB=0.0D0
    if (kHB/=0.0D0) then
        !确定体系中各原子间的成键情况
        call GetAtomBonds(16)
        !确定体系中各原子在XScore中的原子半径
        call GetAtomHBRadii()
        !设置计算过程中每个原子最多形成的临时氢键数量,并确定体系中各原子能提供和接受的最大氢键数量
        call GetAtomHBMaxNum(32)
        !计算受体配体间的氢键项分值
        call CalScoreHB(ScoreHB)
    end if
    !计算最终的打分值
    APBScore=ScoreVDW+ScoreELE+ScoreHB+C0
    if (APBScore>=10.0D0) then
        APBScore=10.0D0
    end if
    write(*,'("APBScore = ",f8.4)')APBScore
end

!计算受体配体间的范德华项分值和静电项分值
subroutine CalScoreVDWELE(ScoreVDW,ScoreELE)
    !调用模块
    use ReadTop
    use ReadRst
    use ReadFit
    !申明变量
    implicit none
    real*8 :: ScoreVDW                                                        !打分函数中范德华项的分值
    real*8 :: ScoreELE                                                        !打分函数中静电项的分值
    real*8,allocatable,dimension(:,:) :: RECLIG_VDWA                          !受体配体原子对间范德华的A值
    real*8,allocatable,dimension(:,:) :: RECLIG_VDWB                          !受体配体原子对间范德华的B值
    real*8 :: Dist2,Dist6,Dist12                                              !两原子间距离的平方、6次方、12次方
    real*8 :: VDWEnergy,ELEEnergy                                             !两原子间的范德华能量和静电能量
    integer*4 :: n,m                                                          !循环变量
    integer*4 :: Temp1,Temp2                                                  !临时变量
    !确定受体配体原子对间的范德华能中的A和B值
    allocate(RECLIG_VDWA(RecAtomTot,LigAtomTot));RECLIG_VDWA=0.0D0
    allocate(RECLIG_VDWB(RecAtomTot,LigAtomTot));RECLIG_VDWB=0.0D0
    do n=LigAtomStart,LigAtomStop,1
        Temp1=n-LigAtomStart+1
        do m=RecAtomStart,RecAtomStop,1
            Temp2=m-RecAtomStart+1
            if (n<m) then
                RECLIG_VDWA(Temp2,Temp1)=LJA(NBPARM(NATOMTYPE*(ATOMTYPE(n)-1)+ATOMTYPE(m)))
                RECLIG_VDWB(Temp2,Temp1)=LJB(NBPARM(NATOMTYPE*(ATOMTYPE(n)-1)+ATOMTYPE(m)))
            else
                RECLIG_VDWA(Temp2,Temp1)=LJA(NBPARM(NATOMTYPE*(ATOMTYPE(m)-1)+ATOMTYPE(n)))
                RECLIG_VDWB(Temp2,Temp1)=LJB(NBPARM(NATOMTYPE*(ATOMTYPE(m)-1)+ATOMTYPE(n)))
            end if
        end do
    end do
    !判断不参与打分的受体元素序号
    do n=RecAtomStart,RecAtomStop,1
        if ((ATOMICNUM(n)<=0).or.(ATOMICNUM(n)>=100)) then
            write(*,'("Error: Atomic number of receptor atoms should not be less than 1 nor larger than 99")')
            stop
        else
            if (FlagElementRec(ATOMICNUM(n))==0) then
                FlagElementRec(ATOMICNUM(n))=-1
            end if
        end if
    end do
    do n=1,99,1
        if (FlagElementRec(n)==-1) then
            write(*,'("Notice: VDW interaction between ligand and receptor atom with atomic number ",i2," will be ignored")')n
        end if
    end do
    !判断不参与打分的配体元素序号
    do n=LigAtomStart,LigAtomStop,1
        if ((ATOMICNUM(n)<=0).or.(ATOMICNUM(n)>=100)) then
            write(*,'("Error: Atomic number of ligand atoms should not be less than 1 nor larger than 99")')
            stop
        else
            if (FlagElementLig(ATOMICNUM(n))==0) then
                FlagElementLig(ATOMICNUM(n))=-1
            end if
        end if
    end do
    do n=1,99,1
        if (FlagElementLig(n)==-1) then
            write(*,'("Notice: VDW interaction between receptor and ligand atom with atomic number ",i2," will be ignored")')n
        end if
    end do
    !计算受体配体原子对间的范德华项分值和静电项的分值
    ScoreVDW=0.0D0
    ScoreELE=0.0D0
    do n=LigAtomStart,LigAtomStop,1
        Temp1=n-LigAtomStart+1
        do m=RecAtomStart,RecAtomStop,1
            Temp2=m-RecAtomStart+1
            Dist2=sum((Crd_All((m*3-2):(m*3))-Crd_All((n*3-2):(n*3)))**2)
            if ((FlagElementRec(ATOMICNUM(m))>0).and.(FlagElementLig(ATOMICNUM(n))>0)) then
                Dist6=Dist2**3
                Dist12=Dist6**2
                VDWEnergy=RECLIG_VDWA(Temp2,Temp1)/Dist12-RECLIG_VDWB(Temp2,Temp1)/Dist6
                ScoreVDW=ScoreVDW+VDWEnergy*kVDW(FlagElementRec(ATOMICNUM(m)),FlagElementLig(ATOMICNUM(n)))
            end if
            if (kELE/=0.0D0) then
                ELEEnergy=CHARGE(m)*CHARGE(n)/sqrt(Dist2)
                ScoreELE=ScoreELE+ELEEnergy*kELE
            end if
        end do
    end do
end subroutine CalScoreVDWELE

!计算受体配体间的氢键项分值
subroutine CalScoreHB(ScoreHB)
    !调用模块
    use ReadTop
    use ReadRst
    use ReadFit
    use PrepHB
    !申明变量
    implicit none
    real*8 :: ScoreHB                                                         !受体配体间的氢键项分值
    real*8 :: CrdDB(NATOM*3)                                                  !与各氢键供体原子相连的重原子的平均XYZ坐标
    real*8 :: CrdAB(NATOM*3)                                                  !与各氢键受体原子相连的重原子的平均XYZ坐标
    integer*4 :: D,A                                                          !氢键供体重原子、氢键受体重原子
    integer*4 HBIndex                                                         !氢键序号
    integer*4 :: HBAtomDNum(NATOM)                                            !记录各氢键供体原子形成的氢键数量
    integer*4 :: HBAtomDA(MaxHBPerAtom,NATOM)                                 !记录与各氢键供体原子形成氢键的受体原子号
    integer*4 :: HBIndexD(MaxHBPerAtom,NATOM)                                 !记录与各氢键供体原子形成氢键的氢键序号
    real*8 :: HBEnergyD(MaxHBPerAtom,NATOM)                                   !记录与各氢键供体原子形成氢键的能量
    integer*4 :: HBAtomANum(NATOM)                                            !记录各氢键受体原子形成的氢键数量
    integer*4 :: HBAtomAD(MaxHBPerAtom,NATOM)                                 !记录与各氢键受体原子形成氢键的供体原子号
    integer*4 :: HBIndexA(MaxHBPerAtom,NATOM)                                 !记录与各氢键受体原子形成氢键的氢键序号
    real*8 :: HBEnergyA(MaxHBPerAtom,NATOM)                                   !记录与各氢键受体原子形成氢键的能量
    real*8 :: TempCrd(3)                                                      !临时原子坐标
    integer*4 :: Temp1                                                        !临时变量
    integer*4 :: n,m                                                          !循环变量
    !计算与各氢键供体或受体原子相连的重原子的平均XYZ坐标
    CrdDB=0.0D0
    CrdAB=0.0D0
    do n=1,NATOM,1
        if (ATOMHBD(n)>0) then
            TempCrd=0.0D0
            Temp1=0
            do m=1,ATOMBONDS(n,0),1
                if (ATOMICNUM(ATOMBONDS(n,m))>1) then
                    TempCrd=TempCrd+Crd_All((ATOMBONDS(n,m)*3-2):(ATOMBONDS(n,m)*3))
                    Temp1=Temp1+1
                end if
            end do
            CrdDB((3*n-2):(3*n))=TempCrd/Temp1
        end if
        if (ATOMHBA(n)>0) then
            TempCrd=0.0D0
            Temp1=0
            do m=1,ATOMBONDS(n,0),1
                if (ATOMICNUM(ATOMBONDS(n,m))>1) then
                    TempCrd=TempCrd+Crd_All((ATOMBONDS(n,m)*3-2):(ATOMBONDS(n,m)*3))
                    Temp1=Temp1+1
                end if
            end do
            CrdAB((3*n-2):(3*n))=TempCrd/Temp1
        end if
    end do
    !初始化数组
    HBIndex=0
    HBAtomDNum=0;HBAtomANum=0
    HBAtomDA=0;HBAtomAD=0
    HBIndexD=0;HBIndexA=0
    HBEnergyD=0.0D0;HBEnergyA=0.0D0
    !计算受体中的各氢键供体与配体中的各氢键受体间形成的氢键的能量
    do D=RecAtomStart,RecAtomStop,1
        if (ATOMHBD(D)>0) then
            do A=LigAtomStart,LigAtomStop,1
                if (ATOMHBA(A)>0) then
                    !计算D与A间的氢键能量
                    call CalHB(D,A,CrdDB,CrdAB,HBIndex,HBAtomDNum,HBAtomANum,HBAtomDA,HBAtomAD,&
                               HBIndexD,HBIndexA,HBEnergyD,HBEnergyA)
                end if
            end do
        end if
    end do
    !计算配体中的各氢键供体与受体中的各氢键受体间形成的氢键的能量
    do D=LigAtomStart,LigAtomStop,1
        if (ATOMHBD(D)>0) then
            do A=RecAtomStart,RecAtomStop,1
                if (ATOMHBA(A)>0) then
                    !计算D与A间的氢键能量
                    call CalHB(D,A,CrdDB,CrdAB,HBIndex,HBAtomDNum,HBAtomANum,HBAtomDA,HBAtomAD,&
                               HBIndexD,HBIndexA,HBEnergyD,HBEnergyA)
                end if
            end do
        end if
    end do
    !对于每个氢键供体原子，仅保留氢键能量最高的前最大允许氢键数量个氢键
    do n=1,NATOM,1
        if (ATOMHBD(n)>0) then
            call KeepBestHB(HBEnergyD(:,n),HBAtomDA(:,n),HBIndexD(:,n),HBEnergyA,HBAtomAD,HBIndexA,ATOMHBD(n))
        end if
    end do
    !对于每个氢键受体原子，仅保留氢键能量最高的前最大允许氢键数量个氢键
    do n=1,NATOM,1
        if (ATOMHBA(n)>0) then
            call KeepBestHB(HBEnergyA(:,n),HBAtomAD(:,n),HBIndexA(:,n),HBEnergyD,HBAtomDA,HBIndexD,ATOMHBA(n))
        end if
    end do
    !计算最终保留的氢键贡献的打分值
    ScoreHB=0.0D0
    do n=1,NATOM,1
        if (ATOMHBD(n)>0) then
            do m=1,MaxHBPerAtom,1
                if (HBAtomDA(m,n)>0) then
                    ScoreHB=ScoreHB+HBEnergyD(m,n)
                end if
            end do
        end if
        if (ATOMHBA(n)>0) then
            do m=1,MaxHBPerAtom,1
                if (HBAtomAD(m,n)>0) then
                    ScoreHB=ScoreHB+HBEnergyA(m,n)
                end if
            end do
        end if
    end do
    ScoreHB=ScoreHB*0.5*kHB
    !子程序
    contains
        !计算两原子间的氢键能量
        subroutine CalHB(D,A,CrdDB,CrdAB,HBIndex,HBAtomDNum,HBAtomANum,HBAtomDA,HBAtomAD,HBIndexD,HBIndexA,HBEnergyD,HBEnergyA)
            !申明变量
            implicit none
            integer*4 :: D,A                                                  !氢键供体重原子、氢键受体重原子
            real*8 :: CrdDB(NATOM*3)                                          !与各氢键供体原子相连的重原子的平均XYZ坐标
            real*8 :: CrdAB(NATOM*3)                                          !与各氢键受体原子相连的重原子的平均XYZ坐标
            integer*4 HBIndex                                                 !氢键序号
            integer*4 :: HBAtomDNum(NATOM)                                    !记录各氢键供体原子形成的氢键数量
            integer*4 :: HBAtomDA(MaxHBPerAtom,NATOM)                         !记录与各氢键供体原子形成氢键的受体原子号
            integer*4 :: HBIndexD(MaxHBPerAtom,NATOM)                         !记录与各氢键供体原子形成氢键的氢键序号
            real*8 :: HBEnergyD(MaxHBPerAtom,NATOM)                           !记录与各氢键供体原子形成氢键的能量
            integer*4 :: HBAtomANum(NATOM)                                    !记录各氢键受体原子形成的氢键数量
            integer*4 :: HBAtomAD(MaxHBPerAtom,NATOM)                         !记录与各氢键受体原子形成氢键的供体原子号
            integer*4 :: HBIndexA(MaxHBPerAtom,NATOM)                         !记录与各氢键受体原子形成氢键的氢键序号
            real*8 :: HBEnergyA(MaxHBPerAtom,NATOM)                           !记录与各氢键受体原子形成氢键的能量
            real*8 :: Dist,Dist2                                              !两原子间的距离和距离的平方
            real*8 :: Angle1,Angle2                                           !三个原子间的角度
            real*8 :: FDist,FAngle1,FAngle2,FHB                               !计算氢键能量时的临时变量
            real*8 :: PI=3.1415926536D0                                       !圆周率
            real*8 :: Rad60,Rad120,Rad180                                     !角度60、120、180对应的弧度
            !计算角度60、120、180对应的弧度
            Rad60=PI/3.0D0
            Rad120=PI*2.0D0/3.0D0
            Rad180=PI
            !计算D与A间的氢键情况
            Dist2=sum((Crd_All((D*3-2):(D*3))-Crd_All((A*3-2):(A*3)))**2)
            if (Dist2<(ATOMRADII(D)+ATOMRADII(A))**2) then
                !计算FDist
                Dist=sqrt(Dist2)
                if (Dist<(ATOMRADII(D)+ATOMRADII(A)-0.7D0)) then
                    FDist=1.0D0
                else
                    FDist=(ATOMRADII(D)+ATOMRADII(A)-Dist)/0.7D0
                end if
                !计算DB,D,A间的夹角和FAngle1
                if (ATOMBONDS(D,0)==0) then
                    Angle1=Rad180
                else
                    call CalAngle(CrdDB((D*3-2):(D*3)),Crd_All((D*3-2):(D*3)),Crd_All((A*3-2):(A*3)),Angle1)
                end if
                if (Angle1>=Rad120) then
                    FAngle1=1.0D0
                else if (Angle1>=Rad60) then
                    FAngle1=(Angle1-Rad60)/Rad60
                else
                    FAngle1=0.0D0
                end if
                !计算D,A,AB间的夹角和FAngle2
                call CalAngle(Crd_All((D*3-2):(D*3)),Crd_All((A*3-2):(A*3)),CrdAB((A*3-2):(A*3)),Angle2)
                if (Angle2>=Rad120) then
                    FAngle2=1.0D0
                else if (Angle2>=Rad60) then
                    FAngle2=(Angle2-Rad60)/Rad60
                else
                    FAngle2=0.0D0
                end if
                !计算氢键能量
                FHB=FDist*FAngle1*FAngle2
                if (FHB>1.0D-3) then
                    HBIndex=HBIndex+1
                    HBAtomDNum(D)=HBAtomDNum(D)+1
                    if (HBAtomDNum(D)>MaxHBPerAtom) HBAtomDNum(D)=MaxHBPerAtom
                    HBAtomANum(A)=HBAtomANum(A)+1
                    if (HBAtomANum(A)>MaxHBPerAtom) HBAtomANum(A)=MaxHBPerAtom
                    HBAtomDA(HBAtomDNum(D),D)=A
                    HBAtomAD(HBAtomANum(A),A)=D
                    HBEnergyD(HBAtomDNum(D),D)=FHB
                    HBEnergyA(HBAtomANum(A),A)=FHB
                    HBIndexD(HBAtomDNum(D),D)=HBIndex
                    HBIndexA(HBAtomANum(A),A)=HBIndex
                end if
            end if
        end subroutine CalHB
        !计算三原子间(向量B->A和B->C间)的夹角
        subroutine CalAngle(CrdA,CrdB,CrdC,Angle)
            !申明变量
            implicit none
            real*8 :: CrdA(3),CrdB(3),CrdC(3)                                 !ABC原子的XYZ坐标
            real*8 :: Angle                                                   !向量B->A和B->C间的夹角
            real*8 :: Temp1,Temp2,Temp3                                       !临时变量
            !进行计算
            Temp1=(CrdA(1)-CrdB(1))*(CrdC(1)-CrdB(1))+(CrdA(2)-CrdB(2))*(CrdC(2)-CrdB(2))+(CrdA(3)-CrdB(3))*(CrdC(3)-CrdB(3))
            Temp2=sum((CrdA-CrdB)**2)
            Temp3=sum((CrdC-CrdB)**2)
            Angle=acos(Temp1/sqrt(Temp2*Temp3))
        end subroutine CalAngle
        !对氢键供体或受体原子，依据其最多能参与的氢键数量，删去能量较低的多余的氢键
        subroutine KeepBestHB(HBEnergy1,HBAtom1,HBIndex1,HBEnergy2,HBAtom2,HBIndex2,KeepHBNum)
            !申明变量
            implicit none
            real*8 :: HBEnergy1(MaxHBPerAtom)                                 !原子A作为供体或受体参与的各个氢键的能量
            integer*4 :: HBAtom1(MaxHBPerAtom)                                !原子A作为供体或受体参与的各个氢键中相应的受体或供体原子序号
            integer*4 :: HBIndex1(MaxHBPerAtom)                               !原子A作为供体或受体参与的各个氢键的氢键编号
            real*8 :: HBEnergy2(MaxHBPerAtom,NATOM)                           !原子A作为供体或受体时，体系中各个原子作为受体或供参与的各个氢键的能量
            integer*4 :: HBAtom2(MaxHBPerAtom,NATOM)                          !原子A作为供体或受体时，体系中各个原子作为受体或供参与的各个氢键的供体或受体原子序号
            integer*4 :: HBIndex2(MaxHBPerAtom,NATOM)                         !原子A作为供体或受体时，体系中各个原子作为受体或供参与的各个氢键的氢键编号
            integer*4 :: KeepHBNum                                            !对于原子A要保留的能量最高的氢键数量
            integer*4 :: KeepFlag(MaxHBPerAtom)                               !用于存储原子A的各个氢键是否要保留
            real*8 :: TempEnergy(MaxHBPerAtom)                                !临时变量
            integer*4 :: Temp1,Temp2                                          !临时变量
            real*8 :: Temp3                                                   !临时变量
            integer*4 :: n,m                                                  !循环变量
            !确定HBEnergy1中要保留的能量最高的前KeepHBNum个氢键
            KeepFlag=0
            TempEnergy=HBEnergy1
            Temp1=0
            do while (Temp1<KeepHBNum)
                Temp2=1
                Temp3=TempEnergy(1)
                do n=2,MaxHBPerAtom,1
                    if (TempEnergy(n)>Temp3) then
                        Temp2=n
                        Temp3=TempEnergy(n)
                    end if
                end do
                if (Temp3>0.0D0) then
                    KeepFlag(Temp2)=1
                    TempEnergy(Temp2)=-1.0D0
                    Temp1=Temp1+1
                else
                    exit
                end if
            end do
            !判断是否要删去多余的氢键
            do n=1,MaxHBPerAtom,1
                if ((KeepFlag(n)==0).and.(HBIndex1(n)>0)) then
                    !删去多余的氢键
                    do m=1,MaxHBPerAtom,1
                        if (HBIndex2(m,HBAtom1(n))==HBIndex1(n)) then
                            HBEnergy2(m,HBAtom1(n))=0.0D0
                            HBAtom2(m,HBAtom1(n))=0
                            HBIndex2(m,HBAtom1(n))=0
                            exit
                        end if
                    end do
                    HBEnergy1(n)=0.0D0
                    HBAtom1(n)=0
                    HBIndex1(n)=0
                end if
            end do
        end subroutine KeepBestHB
end subroutine CalScoreHB

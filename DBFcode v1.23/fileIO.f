
!     File input and outputs
      module fileIO
      implicit none      
      
      contains
      
      
      
!     find the number of lines in a file
      subroutine  fileread
      use errors, only: error
      implicit none
      character(len=255) :: wd1
      character(len=255) :: wd2
!
      integer ierr
      character(len=:), allocatable :: ext
      character(len=:), allocatable :: filename
      character(len=:), allocatable :: fname
!
!     Find working directory
      call finddirectory(wd1,wd2)
!
!
!
!     Search for .INP file in both directories
      call searchinpfile(filename,ierr)
!
      if (ierr==0) then
          write(*,*) 'Located the *.INP file!'
      end if
!
!      
!
!     Read element number and element type from .INP file
      call read_abqinpfile(filename,wd1,wd2,ierr)
!     
!     
!
!
!
      end subroutine  fileread
      
      
      
      
      
!     find the working directory
      subroutine  finddirectory(wd1,wd2)
      implicit none
      character(len=255), intent(out) :: wd1
      character(len=255), intent(out) :: wd2
!
      character(len=255) :: cwd
      character(len=:), allocatable :: dir
      integer :: pos1, pos2


      call getcwd(cwd)
      dir = trim(cwd)
      pos1 = scan(dir, '/\', .true.)
      wd1 = dir(1:pos1)
      write(*,*) 'working directory: '
      write(*,*) trim(wd1)
      
!     The optional directory - Linux system may have two subfolders (ARC)
      pos2 = scan(wd1(1:pos1-1), '/\', .true.)
      wd2 = wd1(1:pos2)
      write(*,*) 'sub-working directory: '
      write(*,*) trim(wd2)
      

      end subroutine  finddirectory
    
          
      
      
      
!     find the working directory
      subroutine  searchinpfile(filename,ierr)
      implicit none

      character(len=:), allocatable, intent(out) :: filename
      integer, intent(out) :: ierr
!     
      integer :: num_files, i, pos1, pos2, pos3, pos4
      character(len=256), allocatable :: filenames(:)
      character(len=256) :: str1
      character(len=:), allocatable :: str2

      num_files = command_argument_count()
      
      
      allocate(filenames(num_files))
      
      do i = 1, num_files
          call get_command_argument(i, filenames(i))
      end do

     
      
      
!     Set error flag
      ierr=1
      
      do i = 1, num_files
          

          
          pos1 = scan(filenames(i), '.', .true.)
          pos2 = scan(filenames(i), 's', .true.)
          pos3 = scan(filenames(i), 'i', .true.)
          pos4 = scan(filenames(i), 'm', .true.)
          
          
          
          if ((pos1>0).and.(pos2>0).and.(pos3>0).and.(pos4>0)) then
              
              if ((pos1<pos2).and.(pos2<pos3).and.(pos3<pos4)) then
              
              
                  str1 = filenames(i)
              
                  str2 = str1(1:pos4-4)
          
                  filename = str2 // '.inp'
          
                  ierr=0
              
                  exit
              
              endif
              
                  
          endif
      

      end do
      

      write(*,*) '*.INP filename: '
      write(*,*) filename
      
      end subroutine  searchinpfile
          
      
      
      
      
      
!     Read ABAQUS .INP file to find:
!     1. ABAQUS element type
!     2. Total number of elements in the mesh
      subroutine read_abqinpfile(finame,wd1,wd2,ierr)
      use globalvariables, only: eltyp, numel,
     + filename, foldername 
      use userinputs, only: defaultfilename, defaultfoldername
      use errors, only: error
      implicit none
      character(len=:), allocatable, intent(in) :: finame
      character(len=255), intent(in) :: wd1
      character(len=255), intent(in) :: wd2
      integer, intent(in) :: ierr
!     Variables
      character(len=256) :: line1, eltyp1
      character(len=:), allocatable :: line2
      character(len=:), allocatable :: fname
!     logical there
      integer  n, i, sline, eline, sum, pos, pos1
      integer nnpel, my_iostat
      
   
      
!     If the *.INP file is located earlier
!     ierr = 0 means --> located the file!
      if (ierr==0) then
          fname = trim(finame)

          write(*,*) 'folder+file: ', fname

          
      
!         Read from the getcwd command
          open(unit=11, file=fname, status='old', iostat=my_iostat)
          
!     Use the working directory with the defaultfilename
      else
          
          write(*,*) 'Could not locate the *.INP file!'
          write(*,*) 'Using the default filename with
     + the working directory!'
      
          

          
!         If not working - try working directory
          fname = trim(wd1) // defaultfilename
          open(unit=11, file=fname, status='old', iostat=my_iostat)

      
!         If not working - try sub-working directory
          if  (my_iostat/=0) then
              close(11)
              write(*,*) 'Could not locate the *.INP file again!'
              write(*,*) 'Using the default filename with
     + the sub-working directory this time!'
              fname = trim(wd2) // defaultfilename
              open(unit=11, file=fname, status='old', iostat=my_iostat)
           endif    

!         If still not working - try the default definition
          if  (my_iostat/=0) then
              close(11)
              write(*,*) 'Could not locate the *.INP file
     + again and again!'
              write(*,*) 'Using the default working directory 
     + with the default filename as the final trial!'
          
              fname = defaultfoldername // defaultfilename
              write (*,*) fname
              open(unit=11, file=fname, status='old', iostat=my_iostat)
          endif
      
 
      

      
      
!         If nothing works - call an error message
          if (my_iostat/=0) then
              close(11)
              write(*,*) 'Could not locate the *.INP file!'
              call error(9)
          endif
      
          write(*,*) '*.INP file location:'
          write(*,*) fname      
      
      endif
!     Reading *.INP file completed!
      
      
!     Assign the global variables
      pos1 = scan(fname, '/\', .true.)
      foldername = fname(1:pos1)
      filename = fname(pos1+1:len(fname))      
      
!     Output Global variables
      write(*,*) 'Global variables for file I/O'
      write(*,*) 'foldername: ', foldername
      write(*,*) 'filename: ', filename
      
      
      
!     Get the total number of lines in the file
      call linesfile(fname,len(fname),n)
      write(*,*) 'number of lines in the .inp file: ', n      
            
      
!     Starting line in the code
      sline = 0
      do i=1, n
    
          read(11,'(A)') line1
     
          if (line1(1:15) == '*Element, type=') then
              eltyp1 = line1(16:256)
              sline = i + 1
              exit
          end if

      end do
    

    
      eltyp = trim(eltyp1)
    
      write(*,*) 'eltyp (*.INP file): ', eltyp
    
   
      nnpel = 0
      call findnumberofnodesperelement(eltyp,len(eltyp),nnpel)

      write(*,*) 'nnpel for displacement (*.INP file)', nnpel
    
!     
      if (eltyp=='') then
          call error(9)
      end if
      
    
    
      sum = 0
      do i=sline, n
    


          read(11,'(A)') line1
        
          !write(*,*) 'line: ', line1

          pos=1
!         Read also the next line
          if (nnpel/=20) then
            
                    
        
              sum =0
            
          else
            
              if (mod(i-sline+1,2)==1) then
                
                  sum =0
                
              end if
            
          end if
        
                
        
        
        
        
          line2 = line1
          do while (pos>0)
            if (sline>0) then

                pos= scan(line2,',')
                if (pos>0) then
                    line2 = line2(pos+1:len(line2))
                    sum = sum +1

                else 
                    pos = 0
                end if
                    
            else
                pos = 0
            end if
        end do
        
            

        
            

        
        !write(*,*) 'sum: ', sum
        
!       Read also the next line
        if (nnpel/=20) then
            
            if (sum/=nnpel) then
                eline = i-1
                exit
            end if
            
        else
        
            
            if (mod(i-sline+1,2)==0) then
                
                if (sum/=nnpel) then
                    eline = i-1
                    exit
                end if
                
            end if
            
            
            
        endif
        
    
      end do
    
    
    !write(*,*) 'end line: ', eline
      if (nnpel/=20) then
          numel = eline - sline + 1
      else
          numel = (eline - sline)/2
      endif
    
      write(*,*) 'numel (*.INP file): ', numel 
    
    
      close(11)
      
      return
      end subroutine read_abqinpfile
      
      
      
      
      
      !     find the number of lines in a file
      subroutine  linesfile(filein,filesize,n)
      implicit none
      integer, intent(in) :: filesize
      character(len =filesize), intent(in) :: filein
      integer, intent(out) :: n
      ! locals
      integer :: iostatus, unit_read
      real    :: dummy


      unit_read = 9
      open(unit=unit_read, file=filein)

      n =0
      do
        read(unit_read, *, iostat=iostatus) dummy
        if (iostatus < 0) exit
        n = n + 1
      end do


      close(unit_read)

      end subroutine  linesfile
    
    
    

    
    
!     Find the number of nodes for displacement analysis (not for GND analysis)
      subroutine findnumberofnodesperelement(eltyp,length,nnpel)  
      implicit none
      character(len=length), intent(in) :: eltyp
      integer, intent(in) :: length
      integer, intent(out) :: nnpel
    
!     Determine the element type from the string
!     nnpel: number of nodes per element
!     numpt: number of integration points per element

      select case(eltyp)
!
!     2D - 3 node linear triangular
      case('CPS3')
          
        nnpel = 3
        
!     2D - 3 node linear triangular
      case('CPE3')
          
        nnpel = 3        
        
      
!     2D - 4 node linear quadrilateral
      case('CPS4')
          
        nnpel = 4
        
!     2D - 4 node linear quadrilateral
      case('CPE4')
          
        nnpel = 4
          
!     2D - 6 node quadratic triangular
      case('CPS6')
          
        nnpel = 6
        
!     2D - 6 node quadratic triangular
      case('CPE6')
          
        nnpel = 6        
          
!     2D - 6 node quadratic triangular reduced integration
      case('CPS6R')
          
        nnpel = 6
          
      case('CPE6R')
          
        nnpel = 6
          
!     2D - 8 node quadratic quadrilateral
      case('CPS8')
          
        nnpel = 8
        
!     2D - 8 node quadratic quadrilateral
      case('CPE8')
          
        nnpel = 8        
      
!     2D - 8 node quadratic quadrilateral with reduced integration
      case('CPS8R')
          
          nnpel = 8

!     2D - 8 node quadratic quadrilateral with reduced integration
      case('CPE8R')
          
          nnpel = 8
      
!     3D - 4 node linear tetrahedron
      case('C3D4')
          
          nnpel = 4
          
!     3D - 6 node linear prismatic
      case('C3D6')
          
         nnpel = 6 
          
        
!     3D - 8 node linear hexahedron
      case('C3D8')
          
          nnpel = 8
                  
  
!     3D - 10 node quadratic tetrahedron
      case('C3D10')
          
          nnpel = 10
          
!     3D - 10 node quadratic tetrahedron reduced integration
      case('C3D10R')
          
         nnpel = 10
          
!     3D - 15 node quadratic prismatic
      case('C3D15')
          
         nnpel = 15
          
!     3D - 15 node quadratic prismatic reduced integration
      case('C3D15R')
          
          nnpel = 15
          
!     3D - 20 node hexahedron
      case('C3D20')
          
          nnpel = 20
          
!     3D - 20 node hexahedron reduced integration
      case('C3D20R')
          
          nnpel = 20
          
      end select
    
    
      end subroutine findnumberofnodesperelement
      
      
      end module fileIO
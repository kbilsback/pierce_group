; BPCH_RENUMBER copies a binary punch file to another punch file,
; but changes one of many fields (see below).

pro Bpch_Sep_Sal, InFile, OutFile, Quiet=Quiet, Tau0=ThisTau0

   ; Close input & output file if I/O error happens
   On_IoError, Quit
            
   ; Arguments
   if ( N_Elements( InFile  ) ne 1 ) then Message, 'INFILE not passed!'
   if ( N_Elements( OutFile ) ne 1 ) then OutFile = InFile + '.new'
 
   ; Keywords
   Verbose = 1L - Keyword_Set( Quiet )

   ; External functions
   FORWARD_FUNCTION Str2Byte, Little_Endian

   ; Open the file
   Open_File, InFile, Ilun, /F77, /GET_LUN, SWAP_ENDIAN=little_endian()

   ; Open the file
   Open_File, OutFile, Ilun_OUT, /F77, /GET_LUN, /Write, $
              SWAP_ENDIAN=little_endian()

   ; Define some variables
   fti        = bytarr(40)
   toptitle   = bytarr(80)
   modelname  = bytarr(20)
   modelres   = fltarr(2)
   mhalfPolar = -1L
   mcenter180 = -1L
   unit       = bytarr(40)
   reserved   = bytarr(40)
   dim        = lonarr(6)
   skip       = -1L

   ; Get the file-top string
   readu, ilun, fti
   writeu, ilun_out, fti
   if ( Verbose ) then begin
      print, '---------------------------------------------'
      print, 'FTI       : ', string(fti)
   endif

   ; read the top title
   readu,ilun,toptitle
   writeu, ilun_Out, toptitle
   if ( Verbose ) then print, 'Title:', string(toptitle)

   ; Count of data blocks
   count = 0L

   ; Loop thru file
   while ( not EOF( ilun ) ) do begin 

      category = bytarr(40)
      tracer   = 0L
      tau0     = 0D
      tau1     = 0D
      skip     = 0L

      ; This line is for the new file format
      readu,ilun,modelname,modelres,mhalfpolar,mcenter180
      readu,ilun,category,tracer,unit,tau0,tau1,reserved,dim,skip

      ; Create data array 
      print, Dim   
      Data = FltArr( Dim[0], Dim[1], Dim[2] )
      readu, Ilun, Data

      if ( Tau0 eq ThisTau0 ) then begin
       ; print modified things

         if ( Verbose ) then begin
            print, '---------------------------------------------'
            print, 'Ilun      : ', Ilun
            print, 'ModelName : ', String( ModelName )
            print, 'ModelRes  : ', String( ModelRes  )
            print, 'MHalfPolar: ', MHalfPolar
            print, 'MCenter180: ', MCenter180
            print, 'Category  : ', String( Category )
            print, 'Tracer    : ', Fix( tracer )
            print, 'Unit      : ', String( Unit )
            print, 'Reserved  : ', String( Reserved )
            print, 'TAU0, TAU1: ', Tau0, Tau1
            print, 'Dim       : ', Dim
            print, 'Skip      : ', Skip
        endif

         ; This line is for the new file format
         WriteU, ilun_OUT, modelname, modelres, mhalfpolar, mcenter180
         WriteU, ilun_OUT, category, tracer,    unit,       tau0,       $
                           tau1,     reserved,  dim,        skip
         WriteU, Ilun_OUT, data

        

      endif

      ;pause

      ; Increment count of data blocks
      count = count + 1L

      ; Undefine data array for safety's sake
      UnDefine, Data
    endwhile
 
Quit:

    On_IoError, Null

    ; Close input file
    Close,    Ilun
    Free_lun, Ilun

    ; Close output file
    Close,    Ilun_OUT
    Free_LUN, Ilun_OUT
end
 
 
 

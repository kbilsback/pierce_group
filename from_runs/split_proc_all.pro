pro split_proc_all

   ; Splits into monthly bpch files (change as necesary)

   ; Output directory
   Dir = '../'

   locs=['ts']
   orunname=['china'] ; output name
   years=[2016]
   ;months=[1,2,3,4,5,6,7,8,9,10,11,12]
   months=[2]
   
   lenl=size(locs,/DIMENSIONS)-1
   for l = 0,lenl[0],1 do begin
     loc=locs[l]
     leny=size(years,/DIMENSIONS)-1
     for y = 0,leny[0],1 do begin
       year=years[y]
       lenm=size(months,/DIMENSIONS)-1
       for m = 0,lenm[0],1 do begin
         month=months[m]
           date=year*long(10000)+month*long(100)+long(1)
           if month lt 12 then begin
            date2=year*long(10000)+(month+1)*long(100)+long(1)
           endif else begin
            date2=(year+1)*long(10000)+(1)*long(100)+long(1)
           endelse
           print,date,date2
           datestr = STRTRIM(date, 2)
           datestr2 = STRTRIM(date2, 2)
;           daystr = STRTRIM(day, 2)
           ; datestr = STRING(date, FORMAT='(I8)')
;           print,datestr
           ;Bpch_Sep_Sal,'ctm.bpch','ctm.'+datestr+'.bpch',Tau0=nymd2tau(date)
           Bpch_Sep_Sal, Dir+'trac_avg.geosfp_4x5_TOMAS15.201602010000','ctm.'+datestr+'.bpch',Tau0=nymd2tau(date)
           InFile = 'ctm.'+datestr+'.bpch'
           print,InFile
;           OutFile = loc+datestr+'_'+datestr+'.nc'
           ;OutFile = runname+'_'+'%DATE%.nc'
           OutFile = orunname+'_'+datestr+'.nc'
           print,OutFile
           ctm_cleanup
           bpch2coards, InFile, OutFile, DIAGN='IJ-AVG-$'
           ;OutFile2 = orunname+'_'+datestr+'_nuc.nc'
           ;bpch2coards, InFile, OutFile2, DIAGN='TOMAS-3D'
           ;OutFile3 = orunname+'_'+datestr+'_PORL-L.nc'
           ;bpch2coards, InFile, OutFile3, DIAGN='PORL-L=$'
;           OutFile4 = orunname+'_'+datestr+'_PL-SUL.nc'
;           bpch2coards, InFile, OutFile4, DIAGN='PL-SUL=$'
;           OutFile5 = orunname+'_'+datestr+'_WETDCV.nc'
;           bpch2coards, InFile, OutFile5, DIAGN='WETDCV-$'
;           OutFile6 = orunname+'_'+datestr+'_WETDLS.nc'
;           bpch2coards, InFile, OutFile6, DIAGN='WETDLS-$'
;           OutFile7 = orunname+'_'+datestr+'_WETDLS.nc'
;           bpch2coards, InFile, OutFile7, DIAGN='WETDLS-$'
;           OutFile8 = orunname+'_'+datestr+'_DRYD.nc'
;           bpch2coards, InFile, OutFile8, DIAGN='DRYD-FLX'
           ;OutFile9 = orunname+'_'+datestr+'_CHEM-L.nc'
           ;bpch2coards, InFile, OutFile9, DIAGN='CHEM-L=$'
;           bpch2nc, InFile, OutFile
           ctm_cleanup
;         endfor
       endfor
     endfor
   endfor

end

;
; NAME:
;    TAP
;
; PURPOSE:
;    GUI driven software package for the analysis of extrasolar
;    transit light curves using wavelet likelihood and the analytic
;    model of Mandel and Agol 2002 (EXOFAST implementation)
;
; AUTHORs:
;    J. Zachary Gazak  (zgazak@ifa.hawaii.edu)
;    John A. Johnson
;    John Tonry
;
;
; CALLING SEQUENCE:
;    IDL> tap

function TAP_MPFit_function,param,time=time,dfit=dfit,finefit=finefit,flux=flux $
                            ,redfit=redfit,tdat=tdat,fdat=fdat,efit=efit,rebin=rebin $
                            ,smooth=smooth,deconLD=deconLD,adv_opt=adv_opt
  
;; Elements of Param: 
;;  param[0] = Period
;;  param[1] = Inclination 
;;  param[2] = a/R*
;;  param[3] = Rp/R*
;;  param[4] = Mid Transit
;;  param[5] = Linear LD
;;  param[6] = Quad LD
;;  param[7] = Eccentricity
;;  param[8] = Omega
;;  param[9] = OOT t^0
;;  param[10] = OOT t^1
;;  param[11] = OOT t^3
;;  param[12] = Sigma Red
;;  param[13] = Sigma White
;;  param[14] = delta light parameter (advanced)

  if n_elements(flux) eq 0 then flux = time*0+1d0
  inc = param[1]
  if inc gt 90d0 then inc = 90d0 - (inc mod 90d0)
   
  if 1-keyword_set(deconLD) then deconLD = 0
  case deconLD of
     0: begin
        u1 = param[5]
        u2 = param[6]
     end
     1: begin
        u1 = (param[6] + 2d0*param[5])/5d0
        u2 = (param[5] - 2d0*param[6])/5d0
     end
     2: begin
        u1 = 2*sqrt(param[5])*param[6]
        u2 = sqrt(param[5])*(1-2*param[6])
     end
  endcase
  
  tap_transite,time,[param[0],inc,1d0/param[2],param[3:4],u1,u2,param[7:n_elements(param)-1]],finefit ;,/trend
  

  if adv_opt[1] then begin
     ;; adjust for 3rd light correction
     ;; see johnson+ 2011ApJ...730...79J equation 9
     finefit += 10d0^(-0.4d0*param[14])
     finefit /=(1d0+10d0^(-0.4d0*param[14]))
  endif

  if rebin eq 0 then dfit = finefit
  if rebin eq 1 then dfit = ((dblarr(smooth[1])+smooth[1]^(-1d0))#reform(finefit,smooth[1],n_elements(tdat)))[0,*]
  dfit *= poly(tdat-min(tdat),param[9:11])
  
  if keyword_set(redfit) then begin
     redfit = (filterredwv((fdat-dfit),param[12],param[13],/zeropad))[0:n_elements(fdat)-1]
  endif
  inc = 0L
  if keyword_set(redfit) then return,(fdat-dfit-redfit) else return,(fdat-dfit)
end


;;; THe following functions control the delta light parameter



;; DEL Kp CODE

pro TAP::adjustdelkp_button,event
  widget_control, event.id, GET_UVALUE= uvalue
  widget_control, /Hourglass
  transit = (*self.transits)[uvalue.value]->get()
  case uvalue.param of
     0: transit.adv_opt[1] = event.select
     4: transit.params[14].fixed = event.select
  endcase
  (*self.transits)[uvalue.value]->set,transit
  transit=0L
  self->delupdate,uvalue.value
end

pro TAP::delupdate,changed
  widget_control, /Hourglass
  
 ; print,changed
  for i=0,n_elements(self.fld2[0,*])-1,1 do $
     for j=0,n_elements(self.fld2[*,0])-1,1 do $
        if j ne changed then $
           if self.fld2[j,i] ne 0 then begin
;     stop
    ; print,j,i
     case i of 
        0: widget_control,self.fld2[j,i],set_button=((*self.transits)[j]->get()).adv_opt[1] 
        1: widget_control,self.fld2[j,i],set_value=((*self.transits)[j]->get()).params[14].prior[i]
        2: widget_control,self.fld2[j,i],set_value=((*self.transits)[j]->get()).params[14].prior[i]
        3: widget_control,self.fld2[j,i],set_value=round(((*self.transits)[j]->get()).params[14].set) 
        4: widget_control,self.fld2[j,i],set_button=((*self.transits)[j]->get()).params[14].fixed
     endcase
   ;;  print,j,((*self.transits)[j]->get()).params[14].prior,((*self.transits)[j]->get()).params[14].value
  endif


  self->updatemodel
  ;self->lcplot
  
end

pro TAP::adv_oneset_all,event
  widget_control, /Hourglass

  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.params[14].set = 1
     transit.params[14].prior[0] = 2
     transit.params[14].prior[1] = 0.75
     transit.params[14].prior[2] = 0.1
     transit.params[14].value = transit.params[14].prior[1] 
     transit.params[14].fixed = 0
     (*self.transits)[i]->set,transit
     transit=0L
     widget_control,self.fld2[i,3],set_value=round(((*self.transits)[i]->get()).params[14].set) 
     widget_control,self.fld2[i,1],set_value=((*self.transits)[i]->get()).params[14].prior[1]
     widget_control,self.fld2[i,2],set_value=((*self.transits)[i]->get()).params[14].prior[2]
     widget_control,self.fld2[i,4],set_button=((*self.transits)[i]->get()).params[14].fixed
  endfor

  self->updatemodel
  self->lcplot

  self->setup,'gaussianpriors'
  
  
end


pro TAP::adv_activate_all,event
  widget_control, /Hourglass
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.adv_opt[1] = 1
     (*self.transits)[i]->set,transit
     transit=0L
     widget_control,self.fld2[i,0],set_button=((*self.transits)[i]->get()).adv_opt[1] 
  endfor  
end

pro TAP::adjustdelkp_event,event
  if n_elements(*event.value) gt 0 then begin
     widget_control, event.id, GET_UVALUE= uvalue
     widget_control, /Hourglass
     transit = (*self.transits)[uvalue.value]->get()
     case uvalue.param of
        0: transit.adv_opt[1] = *event.value
        1: begin
           transit.params[14].prior = [transit.adv_opt[1]*2, *event.value, transit.params[14].prior[2]]
           transit.params[14].value = *event.value
           for i=0,self.num_transits-1,1 do begin
              if i ne uvalue.value then begin
                 trans = (*self.transits)[i]->get()
                 if trans.params[14].set eq transit.params[14].set then begin
                    trans.params[14].prior = transit.params[14].prior
                    trans.params[14].value = transit.params[14].value
                    (*self.transits)[i]->set,trans
                 endif
                 trans = 0L
              endif
           endfor
        end
        2: begin
           transit.params[14].prior =  [transit.adv_opt[1]*2, transit.params[14].prior[1], *event.value]
           for i=0,self.num_transits-1,1 do begin
              if i ne uvalue.value then begin
                 trans = (*self.transits)[i]->get()
                 if trans.params[14].set eq transit.params[14].set then begin
                    trans.params[14].prior = transit.params[14].prior
                    (*self.transits)[i]->set,trans
                 endif
                 trans = 0L
              endif
           endfor
        end
        3: begin
           transit.params[14].set = *event.value
        end
        else: stop
     endcase
   ;;  transit.params[14].value = transit.params[14].prior[1]
   ;  if transit.adv_opt[1] and self.adv_opt[1] then transit.params[14].fixed = 0 else transit.params[14].fixed = 1
    (*self.transits)[uvalue.value]->set,transit

    
 ;    print,transit.params[14].set
 ;    print,transit.params[14].prior
 ;    print,transit.params[14].value
 ;    print,transit.params[14].fixed
     transit = 0L
     
                                ;widget_control,self.fld[uvalue.value,0],set_value = transit.adv_opt[1]
                                ;widget_control,self.fld[uvalue.value,1],set_value = transit.params[14].prior[1]
                                ;widget_control,self.fld[uvalue.value,2],set_value
                                ;= transit.params[14].prior[2]
                                ;          stop
     self->delupdate,uvalue.value
                                ;  stop
  endif
end






;;;; helper functions for repeated math
function tap_diff,ptr
  return,max(*ptr)-min(*ptr)
end

function tap_mm,ptr,adjust=adjust,shift=shift
  if 1-keyword_set(shift) then shift = 0d0
  if 1-keyword_set(adjust) then adjust = 0d0
  return,[min(*ptr)-shift-tap_diff(ptr)*adjust,max(*ptr)-shift+tap_diff(ptr)*adjust]
end

function tap_colors
  colors= {white:  -1 ,$
           red:     0 ,$
           blue:    1 ,$  
           black:   2 ,$
           violet:  3 ,$
           green:   4 ,$
           gray:    5 ,$
           orange:  6 ,$
           yellow:  7 $
          }
  tvlct,0,0,0,colors.black       ; black
  tvlct,255,255,255,colors.white ; white
  tvlct,255,0,0,colors.red       ; red
  tvlct,0,180,100,colors.green   ; green
  tvlct,0,110,220,colors.blue    ; blue
  tvlct,250,110,0,colors.orange  ; orange
  tvlct,120,120,120,colors.gray  ; gray
  tvlct,130,60,170,colors.violet ; indigo_violet_1
  tvlct,255,255,0,colors.yellow  ; indigo_violet_1
  return,colors
end


pro tap::write_setup
  openw,lun,'transit_setup.ascii',width=2500,bufsize=0,/get_lun
  printf,lun,'# TAP combination of setup parameters and light curves.  Adjust parameter setup matrix to'
  printf,lun,'#  change the setup before loading this file as a "Transit File" in a new instance of TAP'
  printf,lun,'#'
  printf,lun,string("#",'transit','long_int','int_length_min','cadence_multiplier',format='(a1,a9,a10,a18,a20)')
  for i=0,self.num_transits-1,1 do begin
     transit = (*self.transits)[i]->get()
     if transit.rebin eq 0 then vals = [1,1] else vals = transit.t_int
     string= string('#',i+1,transit.rebin,vals[0],vals[1],format='(a1,i9,i10,d15.4,i15)')
     printf,lun,string
     transit = 0L
     vals=0L
  endfor
  printf,lun,'#'
  printf,lun,string("#",'transit','set','param','value','lock',$
                    'low_lim','hi_lim','limited','prior=1_penalty=2',$
                    'val','sig',format='(a1,a9,a5,a20,a29,a5,a29,a29,a9,a20,a15,a15)')
  for i=0,self.num_transits-1,1 do begin
     transit = (*self.transits)[i]->get()
     for j=0,n_elements(transit.params)-1,1 do begin
        string = string('#',i+1,transit.params[j].set,strjoin(strsplit(transit.params[j].param,' ',/extract),'_'),transit.params[j].value,transit.params[j].fixed,$
                        transit.params[j].limits[0],transit.params[j].limits[1],transit.params[j].limited[0],transit.params[j].prior[0],transit.params[j].prior[1],$
                        transit.params[j].prior[2],$
                        format='(a1,i9,i5,a20,d29.10,i5,d29.10,d29.10,i9,i20,d15.5,d15.5)')
        printf,lun,string
   ;     print,string
   ;     stop
        string=0L
     endfor
     transit=0L
  endfor
  for i=0,self.num_transits-1,1 do begin
     transit = (*self.transits)[i]->get()
     printf,lun,'# transit '+string(i+1,format='(i5)')
     for j=0,n_elements(transit.transit_tday)-1,1 do printf,lun,transit.transit_tday[j],transit.transit_flux[j],format='(d25.10,d15.10)'
     if i ne self.num_transits-1 then printf,lun,-1,-1,format='(i10,i10)'
  endfor
  close,lun
end



pro tap::execute_mcmc,restart=restart
  self.message = 'Writing ascii setup and transit file'
  self->message
  self->write_setup

  self.message = 'Running Markov Chain Monte Carlo'
  if keyword_set(restart) then self.message = 'Resuming Markov Chain Monte Carlo'
  self->message
  
  transit = (*self.transits)[0]->get()
  transit.diff = self.diff
  (*self.transits)[0]->set,transit
  transit = 0L
  
  if keyword_set(restart) then $
     x = tapmcmc(base_parent=self,input_transits=self.transits,restart=restart,input_parameterization=self.parameterize_id) else $
        x = tapmcmc(base_parent=self,input_transits=self.transits,input_parameterization=self.parameterize_id)
  
                                ;x = obj_new('TAPmcmc',base_parent=self, input_transits=self.transits)
  x->start
  
  self.message = 'MCMC Complete.'
  self->message
  
  ptr_free,self.mcmc_stuff
  self.mcmc_stuff = ptr_new(x->info())
  cd,current=thisdir
  (*self.mcmc_stuff).savefile =  thisdir+'/'+ (*self.mcmc_stuff).savefile
  widget_control,self.mcmc_fld[1],set_value = (*self.mcmc_stuff).savefile
  widget_control,self.buttons[0],sensitive=1
  x->destroy
  self->loadmcmc
end

pro TAP_Event,event
  Widget_Control, event.id, Get_UValue=info
  Call_Method, info.method, info.object, event
end


function tap::userquery,message
 D = dialog_message(message,/question,/center)
 if strcmp(d,'Yes') then return,1
 return,0
end


function TAP::stripspecial,string
  
  if strcmp(string,'tau_o') then return,'$\tau_o$'

  string = strjoin(strsplit(string,"_",/extract),"\_")
  
  if strcmp(string,strjoin(strsplit(string,"*",/extract),'*')) eq 0 then string += ' '
  string = strjoin(strsplit(string,"*",/extract),"$_*$")
  
;  print,string
  if n_elements(strsplit(string,"^")) eq 2 then string = strmid(string,0,(strsplit(string,"^"))[1]-1)+'$^'+strmid(string,(strsplit(string,"^"))[1])+'$'
  
  return,string
end

pro tap::remake_t,noreset=noreset
  transit = (*self.transits)[self.active_transit-1]->get() 
  case transit.rebin of
     0: t = ptr_new(transit.transit_tday)
     1: begin 
        if transit.t_int[0] eq 0 then begin
           transit.t_int[0] = median(transit.transit_tday[1:n_elements(transit.transit_tday)-1]-transit.transit_tday[0:n_elements(transit.transit_tday)-2]) 
           transit.t_int[1] = 1
        endif
        
        t1 = dblarr(transit.t_int[1])+1d0
        t2 =(findgen(transit.t_int[1])+1-5d-1*(transit.t_int[1]+1d0))*(transit.t_int[0]/transit.t_int[1])/1440d0
        t3 = (dblarr(n_elements(transit.transit_tday))+1)
   
        t = ptr_new((reform((t1)#transit.transit_tday +t2#t3,transit.t_int[1]*n_elements(transit.transit_tday),1))[*,0]) 
        t1 = 0
        t2 = 0
        t3 = 0
     end
  endcase
  ptr_free,transit.model_tfine,transit.model_ffine
  transit.model_tfine = ptr_new(*t)
  transit.model_ffine = ptr_new((*t)*0d0)
  ptr_free,t
  (*self.transits)[self.active_transit-1]->set,transit
  transit = 0L
  if 1-keyword_set(noreset) then begin
     widget_control,self.settings[2],set_value=((*self.transits)[self.active_transit-1]->get()).t_int[0]
     widget_control,self.settings[3],set_value=round(((*self.transits)[self.active_transit-1]->get()).t_int[1])
  endif
end

pro TAP::ButtonEvent,event
  widget_control, event.id, GET_UVALUE= uvalue
  widget_control, /Hourglass
  
  case uvalue.value of
     'advopts': begin
        if self.adv_opt[uvalue.which] then self.adv_opt[uvalue.which] = 0 else self.adv_opt[uvalue.which] = 1
        widget_control,self.settings[4],sensitive=self.adv_opt[1]
        
       ; stop
                                ;     widget_control,self.settings[5],sensitive=self.adv_opt[0]
        
        case uvalue.which of
           0: print,'bad code!!'
           1: begin
              for i=0,self.num_transits-1,1 do begin
                 transit = (*self.transits)[i]->get()
                 transit.adv_opt[1] = self.adv_opt[1]
                 if self.adv_opt[1] then transit.params[14].fixed = 0 else begin
                    transit.params[14].fixed =1
                    transit.params[14].value = 0d0
                 endelse
                 (*self.transits)[i]->set,transit
                 transit = 0L
              endfor
           end
           else: stop
        endcase

        self->prepcol
        
       ;; if uvalue.which eq 0 then $
       ;;    for i=0,self.num_transits-1,1 do begin
       ;;    transit = (*self.transits)[i]->get()
       ;;    transit.adv_opt[0] = self.adv_opt[0]
       ;;    if self.adv_opt[0] then transit.params[14].fixed = 0 else  transit.params[14].fixed =1
       ;;    (*self.transits)[i]->set,transit
       ;;    transit = 0L
       ;; endfor        
                                ;print,self.adv_opt[1]
     end
     'Adjust Advanced Parameters':begin
        self->setup,'adjustadvanced'
        centertlb,(*self.extra_windows)[5]
        widget_control,(*self.extra_windows)[5],/realize
     end
     'create_plots': if self.create_plots then self.create_plots = 0 else self.create_plots = 1 
     'create_mcmcascii': begin
                                ;  print,self.create_mcmcascii
        if self.create_mcmcascii then self.create_mcmcascii = 0 else self.create_mcmcascii = 1
                                ;  print,self.create_mcmcascii
     end
     'calc_gr': if self.calc_gr then self.calc_gr = 0 else self.calc_gr = 1
     'calc_effl': if self.calc_effl then self.calc_effl = 0 else self.calc_effl = 1
     'create_2d': if self.create_2d then self.create_2d=0 else self.create_2d=1
     'diffset': begin
        self.diff = event.value
        self->lcplot
     end
     'Save Current Setup Button': begin
        path = dialog_pickfile(dialog_parent=(*self.bases)[0],title='Select / Type save file')
        if path ne '' then begin
           widget_control,self.transitfile2_fld[1],set_value = path
        endif 
     end
     'Delete Active':begin        
        transit = (*self.transits)[self.active_transit-1]->get()
        for k=0,n_elements(transit.params)-1,1 do begin
           ptr_free,transit.params[k].mcmc_chain
           ptr_free,transit.params[k].refined_mcmc_chain
        endfor
        ptr_free,transit.model_t,transit.model_f,transit.model_tfine,transit.model_ffine,transit.mcmc_files
        transit = 0L
        (*self.transits)[self.active_transit-1]->destroy

        if self.num_transits eq 1 then (*self.transits)=0L else $
           if self.active_transit eq self.num_transits then (*self.transits)=(*self.transits)[0:self.active_transit-2] else $
              if self.active_transit eq 1 then (*self.transits)=(*self.transits)[self.active_transit:self.num_transits-1] else $
                 (*self.transits)=[(*self.transits)[0:self.active_transit-2],(*self.transits)[self.active_transit:self.num_transits-1]]
        
        self.num_transits-=1
        self.active_transit = max([1,self.active_transit-1])  
        
        ;widget_control,(*self.bases)[2],MAP=0
        self->setup,'multi'
        widget_control,(*self.bases)[2],MAP=1
        self->setup,'cross lc links'
        self->setup,'fit'
        
        self->lcplot
        self.message = 'light curve deleted.'
        self->message
     end
     'Saved Setup File Button': begin
        path = dialog_pickfile(dialog_parent=(*self.bases)[0],title='Select Save File',filter='TAP_setup.idlsav')
        if path ne '' then begin
           widget_control,self.transitfile1_fld[1],set_value = path
           widget_control,self.buttons[1],sensitive=1
        endif 
     end
     'Load Setup Button': begin
        self.message = 'Loading existing setup...'
        self->message
        widget_control,self.transitfile1_fld[1],get_value=path
        path = path[0]
        restore,path
        ptr_free,self.transits
        self.transits = ptr_new(TAP_state)
        TAP_state = 0L
        transit = (*self.transits)[0]->get() 
        
        self.parameterize_id = transit.parameterize_id
       ;; self.parameterize_num = transit.parameterize_num
;;           self.setup_ld = transit.ldtype
        self->setup_parameterization
        
        self.setup_LD = transit.LDtype 
        transit = 0L
        
        self.diff = ((*self.transits)[0]->get()).diff
  
        self.adv_opt = ((*self.transits)[0]->get()).adv_opt
        self.active_transit =  n_elements(*self.transits) 

        self.num_transits = n_elements(*self.transits)
        self->setup,'fit'      

        self->setup,'multi'      
        self->setup,'adjustparams'
        self->setup,'cross lc links'
        self->setup,'adjustadvanced'
        self->setup,'adjustlimlocks'

        self->lcplot
        self.numcol=2
        self->prepcol
        resume = 0
        if strcmp('-1',(*((*self.transits)[0]->get()).mcmc_files)[0]) eq 0 then $
           if n_elements(*((*self.transits)[0]->get()).mcmc_files) lt (((*self.transits)[0]->get()).mcmc_params)[0] then $
              resume = self->userquery('TAP has detected an imcomplete MCMC execution.  Continue it?')
        self.message = 'Load complete.'
        self->message
        if resume then self->execute_mcmc,restart='/'+strjoin((strsplit(path,'/',/extract))[0:n_elements(strsplit(path,'/',/extract))-2],'/')+'/'  else begin
           transit = (*self.transits)[0]->get()
           *transit.mcmc_files = '-1'
           (*self.transits)[0]->set,transit
           transit = 0L
        endelse
      end
     'Transit File Button': begin
        path = dialog_pickfile(dialog_parent=(*self.bases)[0],title='Select Transit File',/must_exist)
        if path ne '' then begin
           widget_control,self.transitfile_fld[1],set_value = path
           widget_control,self.buttons[2],sensitive=1
        endif
     end
     'FileType': self.filetype = event.value
     'ParamType': begin
        case event.value of
           'Basic [i, a/R*]': self.parameterize_id = 'basic'
           'T, b':self.parameterize_id = 'adv2'
           'tau_o, b': self.parameterize_id = 'adv1'    
           else: self.parameterize_id = 'basic' 
        endcase
        self.message = 'Switching Parameterization'
        self->message
        self->setup_parameterization
        self.numcol = 1
        self->PrepCol
     end
     'RebinType': begin
        widget_control,self.setup_smooth,sensitive=0
        case event.value of
           'None': self.setup_rebin = 0
          ; 'Rebin to data cadence': self.setup_rebin = 1
           'Rebin to "Input Integration"': begin
              self.setup_rebin = 1
              widget_control,self.setup_smooth,sensitive=1
           end
        endcase
     end
     'LDType':begin
        widget_control,self.setup_smooth,sensitive=0
        case event.value of
           '2u1+u2, u1-2u2': self.setup_LD = 1
           'u1 and u2': self.setup_LD = 0
        endcase
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get() 
           transit.LDtype = self.setup_LD
           (*self.transits)[i]->set,transit
           transit = 0L
        endfor
     end   
     'RebinTypeA': begin
        transit = (*self.transits)[self.active_transit-1]->get() 
        case event.value of
           'None': begin
              transit.rebin = 0
              widget_control,self.settings[2],sensitive=0
              widget_control,self.settings[3],sensitive=0
           end
           'Rebin to "Input Integration"': begin
              transit.rebin = 1
              widget_control,self.settings[2],sensitive=1
              widget_control,self.settings[3],sensitive=1
;              widget_control,self.setup_smooth,sensitive=1
           end
        endcase
        (*self.transits)[self.active_transit-1]->set,transit
        transit = 0L
        self->remake_t
        self->updatemodel
        self->lcplot
        
     end
     'settint': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        self.smooth_val[uvalue.set] = *event.value
;event.value
     end
     'settint2': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 1
        transit = (*self.transits)[self.active_transit-1]->get() 
        transit.t_int[uvalue.set] = *event.value
        (*self.transits)[self.active_transit-1]->set,transit
        transit = 0L
        self->remake_t,/noreset
        self->updatemodel
        self->lcplot
        
;event.value
     end
     'Load Transit Button': begin
        self->loadtransit,event
        ;for i=1,n_elements(*self.bases)-1,1 do widget_control,(*self.bases)[i],MAP=0
        ;self->setup,'load'
        ;widget_control,(*self.bases)[1],MAP=1
     end
     'Clear Data Path': begin
        widget_control,self.transitfile_fld[1],set_value = ''
        widget_control,self.buttons[2],sensitive=0
     end
     'Clear Setup Path': begin
        widget_control,self.transitfile1_fld[1],set_value = ''
        widget_control,self.buttons[1],sensitive=0
     end
     'ActiveTransit': begin
        self.active_transit = event.value
        widget_control,self.settings[0],set_value=string(self.active_transit,format='(i2.2)')+": "+((*self.transits)[self.active_transit-1]->get()).fname
        widget_control,self.settings[1],set_value=((*self.transits)[self.active_transit-1]->get()).rebin   
        widget_control,self.settings[2],set_value=((*self.transits)[self.active_transit-1]->get()).t_int[0]
        widget_control,self.settings[3],set_value=round(((*self.transits)[self.active_transit-1]->get()).t_int[1])
        if ((*self.transits)[self.active_transit-1]->get()).rebin eq 1 then begin
           widget_control,self.settings[2],sensitive=1
           widget_control,self.settings[3],sensitive=1
        endif else begin
           widget_control,self.settings[2],sensitive=0
           widget_control,self.settings[3],sensitive=0
        endelse
        self->lcplot
        self->setup,'adjustparams'
        self->setup,'cross lc links'
        self->prepcol
     end
     'Setup Cross LC Locks': begin
        centertlb,(*self.extra_windows)[0]
        widget_control,(*self.extra_windows)[0],/realize
     end
     'Manual Parameter Adjustment': begin
        self->setup,'adjustparams'
        centertlb,(*self.extra_windows)[1]
        widget_control,(*self.extra_windows)[1],/realize
     end
     'Adjust Limits and Locks': begin
        self->setup,'adjustlimlocks'
        centertlb,(*self.extra_windows)[2]
        widget_control,(*self.extra_windows)[2],/realize
     end 
     'MCMC Parameters': begin
        self->setup,'mcmcparams'
        centertlb,(*self.extra_windows)[3]
        widget_control,(*self.extra_windows)[3],/realize
     end
     'Gaussian Priors': begin
        self->setup,'gaussianpriors'
        centertlb,(*self.extra_windows)[4]
        widget_control,(*self.extra_windows)[4],/realize
     end
     'Execute MCMC': begin
        self->execute_mcmc
     end
     'MCMC File Button': begin
        path = dialog_pickfile(dialog_parent=self.tap_base,title='Select MCMC SETUP file',/must_exist,filter='TAP_setup.idlsav')
        if path ne '' then begin
           widget_control,self.mcmc_fld[1],set_value = path
           widget_control,self.buttons[0],sensitive=1
        endif
     end
     'Load MCMC Button': begin
        self->loadmcmc
    
     end
     'Clear MCMC Path': begin
        widget_control,self.mcmc_fld[1],set_value=''
        widget_control,self.buttons[0],sensitive=0
                                ; self->setup,'inference'
     end 
     'xlock_freeall': begin
        widget_control, event.id, GET_UVALUE= uvalue
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get()
           transit.params[uvalue.param].set = i+1
           widget_control,self.fld[i,uvalue.param],set_value = round(transit.params[uvalue.param].set)
           (*self.transits)[i]->set,transit
           transit = 0L
        endfor
    
        self->updatemodel
     end
     'xlock_lockall': begin
        widget_control, event.id, GET_UVALUE= uvalue
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get()
           transit.params[uvalue.param].set = 1
           widget_control,self.fld[i,uvalue.param],set_value = round(transit.params[uvalue.param].set)
           (*self.transits)[i]->set,transit
           transit = 0L
        endfor
    
        self->updatemodel
     end
     else: print,'Unknown Button Event "'+uvalue.value+'"'
  endcase
end

pro TAP::QuitWidget,event
  widget_control, event.id, GET_UVALUE= uvalue

  ;help,/heaps
;  stop
  widget_control,(*self.extra_windows)[uvalue.wid],map=0
  if uvalue.wid eq 0 then self.fld[*,*] = 0
;     for i=0,n_elements(self.fld[*,0])-1,1 do $
;        for j=0,n_elements(self.fld[0,*])-1 do if self.fld[i,j] ne 0 then begin
;     widget_control,self.fld[i,j],/destroy
;     self.fld[i,j] = 0
;  endif
  widget_control,(*self.extra_windows)[uvalue.wid],/destroy
  
  case uvalue.wid of
     0: self->setup,'cross lc links'
     1: self->setup,'adjustparams'
     2: self->setup,'adjustlimlocks'
     3: self->setup,'mcmcparams'
     4: self->setup,'gaussianpriors'
     5: self->setup,'adjustadvanced'
     else: print,string(uvalue.wid,format='(i)')+' is an unknown setup'
  endcase 
  self->prepcol
  uvalue = 0L
end

;; pro TAP::QuitLinks,event
;;   widget_control,(*self.extra_windows)[0],/destroy
;;   self->setup,'cross lc links'
;; end

;; pro TAP::QuitManualAdjust,event
;;   widget_control,(*self.extra_windows)[1],/destroy
;;   self->setup,'adjustparams'
;; end

;; pro TAP::QuitAdjustLL,event
;;   widget_control,(*self.extra_windows)[2],/destroy
;;  ; self->setup,'adjustlimlocks'
;; end

;; pro tap::quitmcmcparams,event
;;   widget_control,(*self.extra_windows)[3],/destroy
;; end

pro TAP::MCMC_inference
  self.message = 'Conducting Bayesian inference.'
  self->message

  for i=0,n_elements(*self.transits)-1,1 do begin
     self.message = 'Transit '+string(i+1,format='(i2.2)')+' of '+string(self.num_transits,format='(i2.2)')
     self->message
     self.active_transit = i+1
     transit = (*self.transits)[i]->get()
     for k=0,n_elements(transit.params)-1,1 do begin
        bigarr = *transit.params[k].refined_mcmc_chain
        if strcmp(transit.params[k].param,'Inclination') then if (where(bigarr ge 90d0))[0] ne -1 then bigarr[where(bigarr ge 90)] = 90d0 - (bigarr[where(bigarr ge 90d0)] mod 90d0)
        
        sorted = bigarr[sort(bigarr)]
        range = n_elements(bigarr)/100d0
        transit.params[k].mcmc_val = [sorted[50d0*range],$
                                      sorted[84.135d0*range]-sorted[50d0*range],$
                                      sorted[50d0*range]-sorted[15.865d0*range]]
        if transit.params[k].fixed then transit.params[k].mcmc_val[1:2] = -1d0
        transit.params[k].value =  transit.params[k].mcmc_val[0]
        
        sorted = 0L
        range = 0L
        bigarr = 0L
     endfor
     (*self.transits)[i]->set,transit
     self.message = 'Combined '+string(self.mcmc_complete,format='(i2.2)')+' MCMC chains.'
     if transit.mcmc_params[0] gt 1 then self->message
     transit=0L
     self->updatemodel
     self.numcol=4
     self->prepcol
     self->lcplot
  endfor
  self.message = 'MCMC inference complete.'
  self->message
end

pro TAP::MCMC_inference_basic
  self.message = 'Conducting Bayesian inference.. basic param set.'
  self->message

  for i=0,n_elements(*self.transits)-1,1 do begin
     self.message = 'Transit '+string(i+1,format='(i2.2)')+' of '+string(self.num_transits,format='(i2.2)')
     self->message
     self.active_transit = i+1
     transit = (*self.transits)[i]->get()
     for k=0,n_elements(transit.basic_params)-1,1 do begin
        bigarr = *transit.basic_params[k].refined_mcmc_chain
        if strcmp(transit.basic_params[k].param,'Inclination') then if (where(bigarr ge 90d0))[0] ne -1 then bigarr[where(bigarr ge 90)] = 90d0 - (bigarr[where(bigarr ge 90d0)] mod 90d0)
        
        sorted = bigarr[sort(bigarr)]
        range = n_elements(bigarr)/100d0
        transit.basic_params[k].mcmc_val = [sorted[50d0*range],$
                                      sorted[84.135d0*range]-sorted[50d0*range],$
                                      sorted[50d0*range]-sorted[15.865d0*range]]
        if transit.basic_params[k].fixed then transit.basic_params[k].mcmc_val[1:2] = -1d0
        transit.basic_params[k].value =  transit.basic_params[k].mcmc_val[0]
        
        sorted = 0L
        range = 0L
        bigarr = 0L
     endfor
     (*self.transits)[i]->set,transit
     self.message = 'Combined '+string(self.mcmc_complete,format='(i2.2)')+' MCMC chains.'
     if transit.mcmc_params[0] gt 1 then self->message
     transit=0L
     self->updatemodel
     self.numcol=4
     self->prepcol
     self->lcplot
  endfor
  
  if self.create_2d then self->create_2d
  if self.create_ascii then self->create_ascii
  if self.create_plots then self->create_plots
    
  widget_control,self.mcmc_fld[1],get_value=path
  path=path[0]
  tap_state = (*self.transits)
  for i=0,n_elements(tap_state)-1,1 do begin
     transit = tap_state[i]->get()
     for j=0,n_elements(transit.basic_params)-1,1 do begin
        *transit.basic_params[j].mcmc_chain = -1
        *transit.basic_params[j].refined_mcmc_chain = -1  
        *transit.params[j].mcmc_chain = -1
        *transit.params[j].refined_mcmc_chain = -1              
     endfor
     tap_state[i]->set,transit     
  endfor
  save,TAP_state,filename=path
;;   for i=0,n_elements(tap_state)-1,1 do begin
;;      st = TAP_state[i]->get()
;;      for j=0,n_elements(st.params)-1,1 do begin
;;         ptr_free,st.params[j].mcmc_chain
;;         ptr_free,st.params[j].refined_mcmc_chain
;;      endfor
;;      ptr_free,st.model_t,st.model_f,st.model_tfine,st.model_ffine,st.mcmc_files
;;      TAP_state[i]->destroy
;;      st = 0L
;;   endfor 
  TAP_state = 0L
  

  self.message = 'MCMC inference complete.'
  self->message
end







pro TAP::plot_open, filename, $
                      LANDSCAPE=landscape, $
                      XSIZE=xsize, $
                      YSIZE=ysize, $
                      INCHES=inches, $
                      COLOR=color, $
                      ENCAPSULATED=encapsulated, $
                      BITS_PER_PIXEL=bits_per_pixel, $
                      _REF_EXTRA=_extra

  self.base_plot = !d.name
  set_plot, 'PS', COPY=keyword_set(COLOR), INTERPOLATE=keyword_set(COLOR)
  
  device, FILENAME  = filename, $
          LANDSCAPE = keyword_set(LANDSCAPE), $
          XSIZE     = xsize, $
          YSIZE     = ysize, $
          XOFFSET   = xoffset, $
          YOFFSET   = yoffset, $
          INCHES    = keyword_set(INCHES), $
          COLOR     = keyword_set(COLOR), $
          BITS_PER_PIXEL = 8, $
          ENCAPSULATED   = keyword_set(ENCAPSULATED), $
          _EXTRA    = _extra
end




pro TAP::create_plots
  maxthick = 12

  self.message = 'Creating MCMC plots...'
  self->message
  
  tempP = !p
  tempX = !x
  tempY = !y

  !p.charsize=1.1
  !p.charthick=4
  !x.thick = 5
  !y.thick = 5
  !p.thick= 5
  !p.font = 0
  ;; first plot--ALL LCs together
  tap_readcol,self.plot_dir+"ascii_phased_data.ascii",p,f,rf,t,m,r,format='(d,d,d,d,d,d)'
  tap_readcol,self.plot_dir+"ascii_phased_model.ascii",t2,p2,mf2,format='(d,d,d)'

  all_xr = mm(p*24d0)
  
  ysize=4
  self->plot_open,self.plot_dir+'plot_alltransit_phased_lc.eps',xsize=6,ysize=ysize,/inches,/encapsulated,/color
  yr = [min(f)-(15d0*stdev(r)),max(f)+(3d0*stdev(r))]
 
  if 0 then begin
     plot,p2[sort(p2)]*24d0,(mf2+r2)[sort(p2)],psym=8,color=(*self.colors).black,background=(*self.colors).white,yrange=yr, $
          /ys,/xs,xrange=mm(p2*24d0),symsize=.5,xthick=4,ythick=4,xtitle='Hours from Mid Transit',ytitle='Relative Flux'
                                ; oplot,p[sort(p)],m[sort(p)],color=(*self.colors).blue,thick=4
     oplot,p2[sort(p2)]*24d0,r2[sort(p2)]+(min(mf2+r2)-(8d0*stdev(r2))),color=(*self.colors).black,psym=8,symsize=.5
                                ; hline,(min(f)-(8d0*stdev(r))),color=(*self.colors).blue,thick=4
  endif else begin
     plot,p[sort(p)]*24d0,f[sort(p)],psym=8,color=(*self.colors).black,background=(*self.colors).white,yrange=yr, $
          /ys,/xs,xrange=mm(p*24d0),symsize=.5,xthick=4,ythick=4,xtitle='Hours from Mid Transit',ytitle='Relative Flux'
;;                                 ; oplot,p[sort(p)],m[sort(p)],color=(*self.colors).blue,thick=4
     oplot,p[sort(p)]*24d0,r[sort(p)]+(min(f)-(8d0*stdev(r))),color=(*self.colors).black,psym=8,symsize=.5
;;                                 ; hline,(min(f)-(8d0*stdev(r))),color=(*self.colors).blue,thick=4
  endelse
     
  sharpcorners,thick=4,color=(*self.colors).black
  self->plot_close

    model = interpol(mf2[sort(p2)],p2[sort(p2)],p[sort(p)])
  
  ysize=4
  self->plot_open,self.plot_dir+'plot_alltransit_phase_modandresid.eps',xsize=6,ysize=ysize,/inches,/encapsulated,/color
  yr = [min(f)-(15d0*stdev(r)),max(f)+(3d0*stdev(r))]
  plot,p[sort(p)]*24d0,model+r[sort(p)],psym=8,color=(*self.colors).black,background=(*self.colors).white,yrange=yr, $
       /ys,/xs,xrange=mm(p*24d0),symsize=.5,xthick=4,ythick=4,xtitle='Hours from Mid Transit',ytitle='Relative Flux'
  oplot,p[sort(p)]*24d0,model,color=(*self.colors).blue,thick=5
  oplot,p[sort(p)]*24d0,r[sort(p)]+(min(f)-(8d0*stdev(r))),color=(*self.colors).black,psym=8,symsize=.5
  hline,(min(f)-(8d0*stdev(r))),color=(*self.colors).blue,thick=5
  
  sharpcorners,thick=4,color=(*self.colors).black
  self->plot_close
 
  
  !p.multi=[0,1,1]
  ysize = min([8,3d0 + 1d0*n_elements(*self.transits)])
  self->plot_open,self.plot_dir+'plot_alltransit_lightcurve.eps',xsize=6,ysize=ysize,/inches,/encapsulated,/color
  
  
  for i=0,n_elements(*self.transits)-1,1 do $
     (*self.plot_windows)[0].xrange = [$
     min([(*self.plot_windows)[0].xrange[0],min(((*self.transits)[i]->get()).transit_tday- $
                                                (((*self.transits)[i]->get()).basic_params[where(strcmp(((*self.transits)[i]->get()).basic_params.param,'Mid Transit')  $
                                                                                           eq 1)].value-(((*self.transits)[i]->get()).basic_params[0].value*((*self.transits)[i]->get()).epoch)))]),$
                                                                 max([(*self.plot_windows)[0].xrange[1],$
                                                                      max(((*self.transits)[i]->get()).transit_tday-$
                                                                          (((*self.transits)[i]->get()).basic_params[where(strcmp(((*self.transits)[i]->get()).basic_params.param,'Mid Transit') eq 1)].value $
                                                                           - (((*self.transits)[i]->get()).basic_params[0].value*((*self.transits)[i]->get()).epoch)))])]
  



  tot = self.num_transits-1
  lci = (*self.transits)[0]->get()
  lcf = (*self.transits)[self.num_transits-1]->get()
 ; yrange = [.995*min(lci.transit_flux)-tot*.01,1.00*max(lcf.transit_flux)+((tot)*.019)]
  diff = self.diff
  if diff eq 0 then diff = 10d0*max([lci.residuals,lcf.residuals])
  yrange = [min(lci.transit_flux)-tot*(diff/2d0)-diff,max(lcf.transit_flux)+((tot+.5)*diff)]

  for i=0,self.num_transits-1,1 do begin
   lc = (*self.transits)[i]->get() 
  ;   diff = 5d0*max(lc.residuals)
     midt = lc.basic_params[where(strcmp(((*self.transits)[i]->get()).basic_params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.basic_params[0].value)
     if i eq 0 then $
        plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,yrange=yrange,$
color=(*self.colors).black,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,$
             /xs,/ys,xtitle='Hours from Mid Transit',ytitle='Relative Flux + Constant',title=title,/nodata
     
     oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,color=(*self.colors).black
     
     oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(diff*i),color=lc.modcol,thick=5
     oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)*(i+1)),psym=8,symsize=.6,color=(*self.colors).black
     hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=5
     oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=5
     oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=5
     xyouts,!x.crange[0]+.05d0*(!x.crange[1]-!x.crange[0]),1d0+diff*i+diff/6,$
            string(i+1,format='(i2.2)')+': '+lc.fname,color=(*self.colors).black,charsize=1,charthick=2,align=0
     lc = 0L   
     midt = 0L
  endfor
  lci = 0L
  lcf = 0L
  tot = 0L
  yrange = 0L
  self->plot_close
  tap_readcol,self.plot_dir+"ascii_phased_data.ascii",p,f,rf,t,m,r,format='(d,d,d,d,d,d)'
  tap_readcol,self.plot_dir+"ascii_phased_model.ascii",t2,p2,mf2,format='(d,d,d)'
  
  yr = [min(f)-(10d0*stdev(r)),max(f)+2.5d0*stdev(r)]
  for i=0,n_elements(*self.transits)-1,1 do begin
     ysize=4
     lc = (*self.transits)[i]->get() 
     
     compare = 'OOT t^0'
     compare2 = 'OOT t^1'
     compare3 = 'OOT t^2'
     self->plot_open,self.plot_dir+'plot_transit_'+string(i+1,format='(i2.2)')+'_'+lc.fname+'.eps',xsize=6,ysize=ysize,/inches,/encapsulated,/color
     
                                ;   diff = 5d0*max(lc.residuals)
     midt = lc.basic_params[where(strcmp(((*self.transits)[i]->get()).basic_params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.basic_params[0].value)
     plot,(lc.transit_tday-midt)*24d0,lc.transit_flux,psym=8,symsize=.6,yrange=yr,$
          color=(*self.colors).black,background=(*self.colors).white,xrange=all_xr,$
          /xs,/ys,xtitle='Hours from Mid Transit',ytitle='Relative Flux + Constant',title=title,/nodata
     
     oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux,psym=8,symsize=.6,color=(*self.colors).black
     

     oplot,p2[sort(p2)]*24d0,mf2[sort(p2)]*poly(p2[sort(p2)]-(min(lc.transit_tday)-midt),[lc.basic_params[where(strcmp(lc.basic_params.param,compare) eq 1)].value,lc.basic_params[where(strcmp(lc.basic_params.param,compare2) eq 1)].value,lc.basic_params[where(strcmp(lc.basic_params.param,compare3) eq 1)].value]),color=lc.modcol,thick=5
;     oplot,(lc.transit_tday-midt)*24d0,(*lc.model_ffine)*poly((lc.transit_tday-min(lc.transit_tday)),[lc.basic_params[where(strcmp(lc.basic_params.param,compare) eq 1)].value,lc.basic_params[where(strcmp(lc.basic_params.param,compare2) eq 1)].value,lc.basic_params[where(strcmp(lc.basic_params.param,compare3) eq 1)].value]),color=lc.modcol,thick=5
     oplot,p2[sort(p2)]*24d0,mf2[sort(p2)],color=(*self.colors).gray,thick=5
     
     oplot,(lc.transit_tday-midt)*24d0,lc.residuals+min(f)-(5d0*stdev(r)),psym=8,symsize=.6,color=(*self.colors).black
     hline,min(f)-(5d0*stdev(r)),color=lc.modcol,thick=5
     oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+lc.rednoise,color=lc.redcol,thick=5
     oplot,(lc.transit_tday-midt)*24d0,min(f)-(5d0*stdev(r))+lc.rednoise,color=lc.redcol,thick=5
     
     sharpcorners,color=(*self.colors).black
     
    
     self->plot_close
    ; stop
  endfor

  ;spawn,'open plot_alltransit_lightcurve.eps'


  self->plot_open,self.plot_dir+'plot_alltransit_mcmc_params3.eps',xsize=7,ysize=8,/inches,/encapsulated,/color
  !y.tickname=' '   
  !p.multi = [0,2,3]
!p.charsize = 1.4
 ; !p.charsize *= (2d0/3)
  
  bin =30
 xtitle= ['Parameter Value','Parameter Value','Parameter Value','Parameter Value','Parameter Value','Parameter Value']
     order = [12,13,9,10,11]
     for j=0,n_elements(order)-1,1 do begin
     all_lock = 1
     xr = [1d50,-1d50]
     yr = [0,0]
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
       if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) ne 0 then begin ;   if transit.basic_params[order[j]].mcmc_val[1] ge 0 then begin
           r = [.1*n_elements(*transit.basic_params[order[j]].mcmc_chain),n_elements(*transit.basic_params[order[j]].mcmc_chain)-1]
           txr = tap_mm(transit.basic_params[order[j]].refined_mcmc_chain,shift=modify,adjust=.05)
           if all_lock eq 1 then xr = txr else  xr = [min([xr[0],txr[0]]),max([xr[1],txr[1]])]
           transit = 0L
           r = 0L
           all_lock = 0
        endif
     endfor
     if all_lock eq 0 then $
        for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        r = [.1*n_elements(*transit.basic_params[order[j]].mcmc_chain),n_elements(*transit.basic_params[order[j]].mcmc_chain)-1]
        if stddev((*transit.basic_params[order[j]].mcmc_chain)[r[0]:r[1]]) ne 0 then begin
           plothist,(*transit.basic_params[order[j]].mcmc_chain)[r[0]:r[1]],tx,ty,bin=(max(xr)-min(xr))/bin,/noplot
                                ;  stop
           yr = [0,max([yr[1],max(ty)])]
        endif
                                ;  print,yr
        tx = 0L
        ty = 0L
        transit = 0L
        r = 0L
     endfor
     
     ;print,yr
    
                 ; print,yr
     yr[1]*=1.1

     ;print,yr
     ;stop
     first = 1
  for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        
        if all_lock then begin
           
           if first then begin
              !x.tickname=' '
              plot,[0,2],[0,2],/nodata,title=transit.basic_params[order[j]].param,background=(*self.colors).white,thick=0,$
                   color=(*self.colors).black,charthick=4
              xyouts,1,.9,'Fixed',color=(*self.colors).black,charsize=3d,charthick=3d,align=0.5
              !x.tickname=''
              sharpcorners,color=(*self.colors).black
              
              first = 0L
           endif
           if j eq 0 then xyouts,!x.crange[1]-.3d0*(!x.crange[1]-!x.crange[0]),yr[1]-(.05*(i+1)*yr[1]),transit.fname,color=transit.color,charsize=0.8d0
           
        endif else begin
           modify = 0d0
           if strcmp(transit.basic_params[order[j]].param,'Mid Transit') then modify = transit.basic_params[order[j]].mcmc_val[0]
           r = [.1*n_elements(*transit.basic_params[order[j]].mcmc_chain),n_elements(*transit.basic_params[order[j]].mcmc_chain)-1]
        
           if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) ne 0  then begin ;transit.basic_params[order[j]].mcmc_val[1] ge 0 then begin
              if first eq 1 then begin
                 plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify,$
                          bin=(max(xr)-min(xr))/bin,$
                          color=transit.color, $
                          background=(*self.colors).white,axiscolor=(*self.colors).black, $
                          xtitle=xtitle[j],title=transit.basic_params[order[j]].param, $
                          xrange=xr,/xs, $
                          thick=maxthick,xticks=3,yrange=yr,/ys 
                 first = 0
              endif else $
                 plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify, $
                          bin=(max(xr)-min(xr))/bin,$
                                ;  bin=tap_diff(transit.basic_params[order[j]].refined_mcmc_chain)/bin,$
                          color=transit.color,/overplot, $
                          background=(*self.colors).white,axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                          title=transit.basic_params[order[j]].param,xrange=xr,/xs, $
                          thick=maxthick*(n_elements(*self.transits)-i)/n_elements(*self.transits),xticks=3,xthick=4,ythick=4,yrange=yr,/ys
              if j eq 0 then xyouts,!x.crange[1]-.3d0*(!x.crange[1]-!x.crange[0]),yr[1]-(.05*(i+1)*yr[1]),transit.fname,color=transit.color,charsize=0.8d0
           endif 
        endelse
              transit = 0L
           
           endfor
  if all_lock eq 0 then begin
     !p.multi=[6-j,2,3]
     plothist,[1,0],$
              color=(*self.colors).black, $
              background=(*self.colors).white,axiscolor=(*self.colors).black, $
              xrange=xr,/xs, $
              thick=0,xticks=3,yrange=yr,/ys,/nodata
  endif
     sharpcorners,color=(*self.colors).black
  endfor 
  !y.tickname=''   



  self->plot_close

  self->plot_open,self.plot_dir+'plot_alltransit_mcmc_params2.eps',xsize=7,ysize=2.67*2,/inches,/encapsulated,/color
  !y.tickname=' '   
  !p.multi = [0,2,2]
  !p.charsize *= (2d0/3)
  
  bin =30
  
  xtitle= ['Parameter Value','Parameter Value','Parameter Value','Parameter Value']
  order = [5,6,7,8]
  for j=0,n_elements(order)-1,1 do begin
     all_lock = 1
     xr = [1d50,-1d50]
     yr = [0,0]
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) ne 0 then begin ;transit.basic_params[order[j]].mcmc_val[1] ge 0 then begin
           r = [.1*n_elements(*transit.basic_params[order[j]].mcmc_chain),n_elements(*transit.basic_params[order[j]].mcmc_chain)-1]
           txr = tap_mm(transit.basic_params[order[j]].refined_mcmc_chain,shift=modify,adjust=.05)
           if all_lock eq 1 then xr = txr else  xr = [min([xr[0],txr[0]]),max([xr[1],txr[1]])]
           transit = 0L
           r = 0L
           all_lock = 0
        endif
     endfor
     if all_lock eq 0 then $
        for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        r = [.1*n_elements(*transit.basic_params[order[j]].mcmc_chain),n_elements(*transit.basic_params[order[j]].mcmc_chain)-1]
        if stddev((*transit.basic_params[order[j]].mcmc_chain)[r[0]:r[1]]) ne 0 then begin
           plothist,(*transit.basic_params[order[j]].mcmc_chain)[r[0]:r[1]],tx,ty,bin=(max(xr)-min(xr))/bin,/noplot
                                ;  stop
           yr = [0,max([yr[1],max(ty)])]
        endif
                                ;  print,yr
        tx = 0L
        ty = 0L
        transit = 0L
        r = 0L
     endfor
                                ; print,yr
     yr[1]*=1.1
     first = 1
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        r = [.1*n_elements(*transit.basic_params[order[j]].mcmc_chain),n_elements(*transit.basic_params[order[j]].mcmc_chain)-1]

        if all_lock then begin
                                ;   yr = [0,2]
           if first then begin
              !x.tickname=' '
              plot,[0,2],[0,2],/nodata,title=transit.basic_params[order[j]].param,background=(*self.colors).white,thick=0,color=(*self.colors).black,charthick=4
              xyouts,1,.9,'Fixed',color=(*self.colors).black,charsize=3d,charthick=3d,align=0.5
                                ;     sharpcorners,color=(*self.colors).black
              
              !x.tickname=''
              first = 0L
           endif
           if j eq 0 then xyouts,!x.crange[1]-.3d0*(!x.crange[1]-!x.crange[0]),yr[1]-(.05*(i+1)*yr[1]),transit.fname,color=transit.color,charsize=0.8d0
           
        endif else begin
           modify = 0d0
           if strcmp(transit.basic_params[order[j]].param,'Mid Transit') then modify = transit.basic_params[order[j]].mcmc_val[0]
          ; stop
            if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) ne 0 then begin
              if first eq 1 then begin
                 plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify,$
                          bin=tap_diff(( (*self.transits)[0]->get()).basic_params[order[j]].refined_mcmc_chain)/bin,$
                          color=transit.color, $
                          background=(*self.colors).white,axiscolor=(*self.colors).black, $
                          xtitle=xtitle[j],title=transit.basic_params[order[j]].param, $
                          xrange=xr,/xs, $
                          thick=maxthick,xticks=3,yrange=yr,/ys 
                 first = 0
              endif else $
                 plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify, $
                          bin=tap_diff(( (*self.transits)[0]->get()).basic_params[order[j]].refined_mcmc_chain)/bin,$
                                ;  bin=tap_diff(transit.basic_params[order[j]].refined_mcmc_chain)/bin,$
                          color=transit.color,/overplot, $
                          background=(*self.colors).white,axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                          title=transit.basic_params[order[j]].param,xrange=xr,/xs, $
                          thick=maxthick*(n_elements(*self.transits)-i)/n_elements(*self.transits),xticks=3,xthick=4,ythick=4,yrange=yr,/ys
              if j eq 0 then xyouts,!x.crange[1]-.3d0*(!x.crange[1]-!x.crange[0]),yr[1]-(.05*(i+1)*yr[1]),transit.fname,color=transit.color,charsize=0.8d0
           endif 
     


        endelse
        transit = 0L
        
     endfor
     if all_lock eq 0 then begin
        !p.multi=[4-j,2,2]
        plothist,[1,0],$
                 color=(*self.colors).black, $
                 background=(*self.colors).white,axiscolor=(*self.colors).black, $
                 xrange=xr,/xs, $
                 thick=0,xticks=3,yrange=yr,/ys,/nodata
     endif
     sharpcorners,color=(*self.colors).black
  endfor 
  !y.tickname=''   



  self->plot_close
  !p.charsize*=(3d0/2)

  
  self->plot_open,self.plot_dir+'plot_alltransit_mcmc_params1.eps',xsize=7,ysize=8,/inches,/encapsulated,/color
  !y.tickname=' '   
  !p.multi = [0,2,3]
  !p.charsize = 1.4d0
  
  bin =30
  
  xtitle= ['Days','Degrees','Parameter Value','Parameter Value','Days since Mid Transit']
  order = [0,1,2,3,4]
  !x.thick = 6
  !y.thick = 6

  for j=0,n_elements(order)-1,1 do begin
     yr = [0,0]
     all_lock = 1
  
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) ne 0 then begin ;;if transit.basic_params[order[j]].mcmc_val[1] ge 0 then begin
                                ;  xmin = floor(min(*transit.basic_params[where(strcmp(((*self.transits)[0]->get()).basic_params.param,'Mid Transit') eq 1)].refined_mcmc_chain))    
           modify = 0d0
           if strcmp(transit.basic_params[order[j]].param,'Mid Transit') then modify = transit.basic_params[order[j]].mcmc_val[0]
           plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify,tx,ty,$
                    bin=tap_diff(( (*self.transits)[0]->get()).basic_params[order[j]].refined_mcmc_chain)/bin,/noplot
           txr = tap_mm(transit.basic_params[order[j]].refined_mcmc_chain,shift=modify,adjust=.05)
           if all_lock eq 1 then xr = txr else xr = [min([xr[0],txr[0]]),max([xr[1],txr[1]])]
           if strcmp(transit.basic_params[order[j]].param,'Inclination') then xr[1] = min([90d0,xr[1]]) 
           
                                ;  stop
           yr = [0,max([yr[1],max(ty)])]
                                ;  print,yr
           tx = 0L
           ty = 0L
           transit = 0L
           r = 0L
           all_lock = 0
        endif
     endfor
     ;; plot the hist
     first = 1
     yr[1]*= 1.1d0
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        
        if all_lock then begin
           yr = [0,2]
           if first then begin
              !x.tickname=' '
              plot,[0,2],[0,2],/nodata,title=transit.basic_params[order[j]].param,background=(*self.colors).white,$
                   color=(*self.colors).black,thick=4,charthick=4
              xyouts,1,.9,'Fixed',color=(*self.colors).black,charsize=3d,charthick=3d,align=0.5
              !x.tickname=''
              first = 0L
              sharpcorners,color=(*self.colors).black
           endif
           if j eq 0 then xyouts,!x.crange[1]-.3d0*(!x.crange[1]-!x.crange[0]),yr[1]-(.05*(i+1)*yr[1]),transit.fname,color=transit.color,charsize=0.8d0
           
        endif else begin
           modify = 0d0
           if strcmp(transit.basic_params[order[j]].param,'Mid Transit') then modify = transit.basic_params[order[j]].mcmc_val[0]
           if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) ne 0  then begin ; if transit.basic_params[order[j]].mcmc_val[1] ge 0 then begin
              if first eq 1 then begin
                 plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify,$
                          bin=tap_diff(( (*self.transits)[0]->get()).basic_params[order[j]].refined_mcmc_chain)/bin,$
                          color=transit.color, $
                          background=(*self.colors).white,axiscolor=(*self.colors).black, $
                          xtitle=xtitle[j],title=transit.basic_params[order[j]].param, $
                          xrange=xr,/xs, $
                          thick=maxthick,xticks=3,yrange=yr,/ys 
                 first = 0
              endif else $
                 plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify, $
                          bin=tap_diff(( (*self.transits)[0]->get()).basic_params[order[j]].refined_mcmc_chain)/bin,$
                                ;  bin=tap_diff(transit.basic_params[order[j]].refined_mcmc_chain)/bin,$
                          color=transit.color,/overplot, $
                          background=(*self.colors).white,axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                          title=transit.basic_params[order[j]].param,xrange=xr,/xs, $
                          thick=maxthick*(n_elements(*self.transits)-i)/n_elements(*self.transits),xticks=3,xthick=4,ythick=4,yrange=yr,/ys
              if j eq 0 then xyouts,!x.crange[1]-.3d0*(!x.crange[1]-!x.crange[0]),yr[1]-(.05*(i+1)*yr[1]),transit.fname,color=transit.color,charsize=0.8d0
           endif 
        endelse
              transit = 0L
              
           endfor
     if all_lock eq 0 then begin
        !p.multi=[6-j,2,3]
        plothist,[1,0],$
                 color=(*self.colors).black, $
                 background=(*self.colors).white,axiscolor=(*self.colors).black, $
                 xrange=xr,/xs, $
                 thick=0,xticks=3,yrange=yr,/ys,/nodata
     endif
     sharpcorners,color=(*self.colors).black
  endfor 
  !y.tickname=''   
  self->plot_close
  




  
  ;;order = [0,1]
  
  ;;for j=0,n_elements(order)-1,1 do $
  ;;   for k=0,n_elements(*self.transits)-1,1 do begin
  ;;   transit  = (*self.transits)[k]->get() 
  ;;   if 
  ;;  endfor
    

  for i=0,n_elements(*self.transits)-1,1 do begin
     !y.tickname=' '
  !p.charsize = 1.4d0
  
     
     ;; plot system variables
     transit  = (*self.transits)[i]->get() 
     self->plot_open,self.plot_dir+'plot_'+string(i+1,format='(i2.2)')+"_"+transit.fname+'_mcmc_params1.eps',xsize=7,ysize=8,/inches,/encapsulated,/color
     !p.multi = [0,2,3]
     bin =30


  xmin = floor(min(*transit.basic_params[where(strcmp(((*self.transits)[0]->get()).basic_params.param,'Mid Transit') eq 1)].refined_mcmc_chain))    
  
     xtitle= ['Days','Degrees','Parameter Value','Parameter Value','Days since '+sigfig(xmin,8)]
     order = [0,1,2,3,4]
     for j=0,n_elements(order)-1,1 do begin
        if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) eq 0  then begin ; if transit.basic_params[order[j]].fixed eq 1 then begin
           !x.tickname=' '
           plot,[0,2],[0,2],/nodata,title=transit.basic_params[order[j]].param,background=(*self.colors).white,$
                color=(*self.colors).black,thick=4,charthick=4
           xyouts,1,1,'Fixed',color=(*self.colors).black,charsize=4d,charthick=3d,align=0.5
           sharpcorners,xthick=!x.thick,ythick=!y.thick,color=(*self.colors).black

           !x.tickname=''

        endif else begin
           modify = 0d0
           if strcmp(transit.basic_params[order[j]].param,'Mid Transit') then modify = xmin
           xrange = tap_mm(transit.basic_params[order[j]].refined_mcmc_chain,adjust=.1,shift=modify)
           if order[j] eq 1 then xrange[1] = min([xrange[1],90d0])
           plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify, $
                    bin=tap_diff(transit.basic_params[order[j]].refined_mcmc_chain)/bin, $
                    color=(*self.colors).black,background=(*self.colors).white, $
                    axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                    title=transit.basic_params[order[j]].param, $
                    xrange=xrange,/xs, $
                    thick=6,xticks=3,charthick=4
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify,color=(*self.colors).blue,thick=6  
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify+transit.basic_params[order[j]].mcmc_val[1], $
                 color=(*self.colors).blue,thick=6,linestyle=1
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify-transit.basic_params[order[j]].mcmc_val[2], $
                 color=(*self.colors).blue,thick=6,linestyle=1
           sharpcorners,xthick=!x.thick,ythick=!y.thick,color=(*self.colors).black
        
endelse
     endfor
     self->plot_close
;; was ysize=2.7
     self->plot_open,self.plot_dir+'plot_'+string(i+1,format='(i2.2)')+"_"+transit.fname+'_mcmc_params2.eps',xsize=7,ysize=2.67*2,/inches,/encapsulated,/color
     !p.multi = [0,2,2]
     !p.charsize *= (1d0/2)
     
     xtitle= ['Parameter Value','Parameter Value','Parameter Value','Parameter Value']
     order = [5,6,7,8]
     for j=0,n_elements(order)-1,1 do begin
        if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) eq 0  then begin ; if transit.basic_params[order[j]].fixed eq 1 then begin
           !x.tickname=' '
           plot,[0,2],[0,2],/nodata,title=transit.basic_params[order[j]].param,background=(*self.colors).white,$
                color=(*self.colors).black,thick=4,charthick=4
           xyouts,1,1,'Fixed',color=(*self.colors).black,charsize=4d,charthick=3d,align=0.5
           sharpcorners,xthick=!x.thick,ythick=!y.thick,color=(*self.colors).black
          
           !x.tickname=''
        endif else begin
           modify = 0d0
           if strcmp(transit.basic_params[order[j]].param,'Mid Transit') then modify = xmin
           plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify, $
                    bin=tap_diff(transit.basic_params[order[j]].refined_mcmc_chain)/bin, $
                    color=(*self.colors).black,background=(*self.colors).white, $
                    axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                    title=transit.basic_params[order[j]].param, $
                    xrange=tap_mm(transit.basic_params[order[j]].refined_mcmc_chain,adjust=.1,shift=modify),/xs, $
                    thick=6,xticks=3,charthick=4
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify,color=(*self.colors).blue,thick=6  
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify+transit.basic_params[order[j]].mcmc_val[1],$
                 color=(*self.colors).blue,thick=6,linestyle=1
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify-transit.basic_params[order[j]].mcmc_val[2], $
                 color=(*self.colors).blue,thick=6,linestyle=1
           sharpcorners,xthick=!x.thick,ythick=!y.thick,color=(*self.colors).black

        endelse
     endfor
     self->plot_close
     !p.charsize *= 2d0

     self->plot_open,self.plot_dir+'plot_'+string(i+1,format='(i2.2)')+"_"+transit.fname+'_mcmc_params3.eps',xsize=7,ysize=8,/inches,/encapsulated,/color
     !p.multi = [0,2,3]
     xtitle= ['Parameter Value','Parameter Value','Parameter Value','Parameter Value','Parameter Value','Parameter Value']
     order = [12,13,9,10,11]
     for j=0,n_elements(order)-1,1 do begin
        if stddev(*transit.basic_params[order[j]].refined_mcmc_chain) eq 0  then begin ; if transit.basic_params[order[j]].fixed eq 1 then begin
           !x.tickname=' '
           plot,[0,2],[0,2],/nodata,title=transit.basic_params[order[j]].param,background=(*self.colors).white,$
                color=(*self.colors).black,thick=4,charthick=4
           xyouts,1,1,'Fixed',color=(*self.colors).black,charsize=4d,charthick=3d,align=0.5
           sharpcorners,xthick=!x.thick,ythick=!y.thick,color=(*self.colors).black
           !x.tickname=''
        endif else begin
           modify = 0d0
           if strcmp(transit.basic_params[order[j]].param,'Mid Transit') then modify = xmin
           plothist,(*transit.basic_params[order[j]].refined_mcmc_chain)-modify, $
                    bin=tap_diff(transit.basic_params[order[j]].refined_mcmc_chain)/bin, $
                    color=(*self.colors).black,background=(*self.colors).white, $
                    axiscolor=(*self.colors).black,xtitle=xtitle[j], $
                    title=transit.basic_params[order[j]].param, $
                    xrange=tap_mm(transit.basic_params[order[j]].refined_mcmc_chain,adjust=.1,shift=modify),/xs, $
                    thick=6,xticks=3,charthick=4
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify,color=(*self.colors).blue,thick=6       
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify+transit.basic_params[order[j]].mcmc_val[1], $
                 color=(*self.colors).blue,thick=6,linestyle=1
           vline,transit.basic_params[order[j]].mcmc_val[0]-modify-transit.basic_params[order[j]].mcmc_val[2], $
                 color=(*self.colors).blue,thick=6,linestyle=1
           sharpcorners,xthick=!x.thick,ythick=!y.thick,color=(*self.colors).black

        endelse
     endfor
     self->plot_close
     !y.tickname=''      


     ;spawn,'open '+self.plot_dir+'plot_'+transit.fname+'_mcmc_LD.eps'
     
     
     
  endfor
  !y.tickname=''    
                                ; stop
  
  
  
  
  
  !p = tempp
  !x = tempx
  !y = tempy
  tempp = 0L
  tempx = 0L
  tempy = 0L

 ; stop
end


;; more memory efficient
function TAP::gelmanrubin
 ;;;; stop

  self.message = '...Calculating Gelman-Rubin statistic'
  self->message
  spawn,'ls '+self.plot_dir+'/MCMC_chains/',chainlist
  totlink = 0d
  
  for i=0,n_elements(chainlist)-1,1 do begin
     restore,self.plot_dir+'/MCMC_chains/'+chainlist[i]
     n_transit = n_elements(tap_state)
                                ; help,tap_state
     
     for j=0,n_elements(tap_state)-1,1 do begin
        x = tap_state[j]->get()
        range = fillarr(1,.1*n_elements(*x.params[0].mcmc_chain),n_elements(*x.params[0].mcmc_chain)-1)
        
        if j eq 0 then mcmcs = dblarr(n_elements(tap_state),$
                                      n_elements(x.params),$
                                      n_elements(range))+sqrt(-1)
        
        for k=0,n_elements(x.params)-1,1 do mcmcs[j,k,0:n_elements(range)-1] = (*x.params[k].mcmc_chain)[range]
     endfor
     
     if i eq 0 then begin r = dblarr((size(mcmcs))[2],(size(mcmcs))[1])
        z_c = replicate({params: dblarr((size(mcmcs))[2])},n_elements(chainlist))
        w_c = z_c
     endif
     
     for k=0,(size(mcmcs))[1]-1,1 do begin ;; loop over transits
        
                                ;  nlink = (size(mcmcs))[4]
                                ;  totlink += (size(mcmcs))[4]
        ;; for i=0,n_elements(z_c)-1,1 do begin  ;; loop over nchains
        good = where(finite(mcmcs[0,0,*]))
        nlink = n_elements(good)
        totlink += nlink
        for j=0,n_elements(z_c[0].params)-1,1 do $
           z_c[i].params[j] = $
           (1d0/nlink)*total((mcmcs[0,j,*])[good])
        
        
        ;;    for i=0,n_elements(w_c)-1,1 do begin
        good = where(finite(mcmcs[0,0,*]))
        nlink = n_elements(good)
        for j=0,n_elements(w_c[0].params)-1,1 do $
           w_c[i].params[j] = $
           (1d0/(nlink-1))*total(((mcmcs[0,j,*])[good]-z_c[i].params[j])^2d0)
        
     endfor
     mcmcs = 0L
     
     for z=0,n_elements(tap_state)-1,1 do begin
        transit = tap_state[z]->get()   
        for k=0,n_elements(transit.params)-1,1 do begin
           ptr_free,transit.params[k].mcmc_chain
           ptr_free,transit.params[k].refined_mcmc_chain
        endfor
        ptr_free,transit.model_t,transit.model_f,transit.model_tfine,transit.model_ffine,transit.mcmc_files
        transit = 0L
        tap_state[z]->destroy
     endfor
     tap_state=0L
  endfor
  
  for k=0,n_elements(*self.transits)-1,1 do begin ;; loop over transits
     w_z = dblarr(n_elements(z_c[0].params))
     z_dd= dblarr(n_elements(z_c[0].params))
     b_z = dblarr(n_elements(z_c[0].params))
     vpz = dblarr(n_elements(z_c[0].params))
     
     for i=0,n_elements(w_z)-1,1 do begin
        w_z[i] = (1d0/n_elements(chainlist))*total(w_c.params[i])
        z_dd[i]= (1d0/n_elements(chainlist))*total(z_c.params[i])
        b_z[i] = (totlink/(n_elements(chainlist)-1))*total((z_c.params[i]-z_dd[i])^2d0)
        vpz[i] = (((totlink-1)/totlink)*w_z[i]) + (1d0/totlink)*b_z[i]
        r[i,k] = sqrt(vpz[i]/w_z[i])
     endfor
     
     if (where(x.params.fixed))[0] ne -1 then r[where(x.params.fixed),k] = 0d0
     
     w_z = 0L
     z_dd = 0L
     b_z = 0L
     vpz = 0L
  endfor
  z_c = 0L
  w_c = 0L
  ;; print,x.params.param,format=format
  ;; print,'Median GR:'
  ;; print,median(r,dimension=2)
  ;; print,'Average GR:'
  ;; print,average(r,2)
  if (size(r))[0] gt 1 then maxr = max(average(r,2)) else maxr = max(r) 
  if maxr gt 1.1d0 then self.message = '...G-R indicates lack of convergence, check outputs.' else $
     self.message = '...G-R shows no evidence of non-convergence.'
  self->message
  return,r
end


pro TAP::plot_close
  device, /CLOSE_FILE
  set_plot, self.base_plot
end

function TAP::makeline2,array
 ; print,array
 ; stop
  if array[1] eq -1 then line = string(array[0],format='(d20.10)')+' [Fixed]' else $
     line = string(array[0],format='(d20.10)')+' +'+strtrim(string(array[1],format='(d15.10)'),2)+' -'+strtrim(string(array[2],format='(d15.10)'),2)
  return,line
end

function TAP::makeline3,val

 ; stop

  
  line = ' & '
 ; i=1
;  while (sigfig(val,i) ne val) and i lt 15 do i++
  ;while (val*(10d0^i) lt 100) do i++
  
 ; f = "(d"+string(i+5,format='(i2.2)')+"."+string(i,format='(i3.3)')+")"  
  f = '(e12.3)'
  line += string(val,format=f)
  
  return,line
end

function TAP::makeline,array
  line = ' & '
  if 0 then begin
     print,array
     stop
  endif
  if array[1] eq -1 then begin
     i=1
     while (sigfig(array[0],i) ne array[0]) and i lt 15 do i++
     line += sigfig(array[0],i)
     line += '\tablenotemark{a}' 
  endif else begin
     if abs(array[0]) ge 1 then begin
        line += sigfig(array[0],strlen(sigfig(array[1],2))-2+strlen(strtrim(floor(array[0]),2)))
        line +=  ' $^{+'+sigfig(array[1],2)+'}_{-'+sigfig(array[2],2)+'}$'
     endif else begin
     ;   print,'needs fix'
        i = 0d0
        while (abs(array[0])*(10^i) lt 1) do i++
        if array[0] lt 0 then i++
        line += sigfig(array[0],strlen(sigfig(array[1],2))-2-i+strlen(strtrim(floor(array[0]),2)))
        line +=  ' $^{+'+sigfig(array[1],2)+'}_{-'+sigfig(array[2],2)+'}$'
     endelse   
  endelse
  if 0 then print,line
  array = 0L
  return,line
end  

function TAP::effective_length,final_corr_len=final_corr_len
  self.message = '...Calculating Effective Lengths'
  self->message
  spawn,'ls '+self.plot_dir+'/MCMC_chains/',chainlist
  
  final_eff_len = dblarr(self.num_transits, n_elements(((*self.transits)[0]->get()).params),n_elements(chainlist))
  final_corr_len = dblarr(self.num_transits, n_elements(((*self.transits)[0]->get()).params),n_elements(chainlist))
    
  for i=0,n_elements(chainlist)-1,1 do begin
     restore,self.plot_dir+'/MCMC_chains/'+chainlist[i]
     n_transit = n_elements(tap_state)
     
     corr_len= dblarr(self.num_transits, n_elements(((*self.transits)[0]->get()).params))
     eff_len = corr_len+sqrt(-1)
     
     fillsize = 1d10
     for j=0,n_elements(tap_state)-1,1 do begin
        x = tap_state[j]->get()
        range = self->burnrange(tnum=j)
        range = fillarr(1,range[0],range[1])
        fillsize = min([fillsize,n_elements(range)])
        x = 0L
     endfor
     
     for j=0,n_elements(tap_state)-1,1 do begin
        x = tap_state[j]->get()
        range = self->burnrange(tnum=j)
        range = fillarr(1,range[0],range[1])
        if j eq 0 then mcmcs = dblarr(n_elements(tap_state),$
                                      n_elements(x.params),$
                                      fillsize)
        for k=0,n_elements(x.params)-1,1 do mcmcs[j,k,*] = (*x.params[k].mcmc_chain)[range[n_elements(range)-fillsize:*]]
        x = 0L
     endfor

    ; stop

     for l=0,n_elements(tap_state)-1,1 do $
        for k=0,n_elements((tap_state[0]->get()).params)-1,1 do begin
        if (tap_state[l]->get()).params[k].fixed eq 0 then begin
           c_j  = [1d0]
           jact = [0d0]
           tot_links = n_elements(range)*10d0
           while c_j[n_elements(c_j)-1] gt 0.48d0 and jact[n_elements(jact)-1] le tot_links/2d0 do begin              
              if n_elements(c_j) ge 2 then begin
                 del =  c_j[n_elements(c_j)-2]-c_j[n_elements(c_j)-1] 
                 step = max([jact[n_elements(c_j)-1]-jact[n_elements(c_j)-2],10])^(.1d0/del)
                 step -= (step mod 10)
                 if step lt 10 then step = 10
                 if step gt 100 then step = 20
                 if step gt 1000 then step = 100
                 if finite(step) ne 1 then step = 100
                 jact = [jact,jact[n_elements(jact)-1]+step]
              endif else jact = [jact,jact[n_elements(jact)-1]+10d0]
              juse = jact[n_elements(jact)-1]/10d0
              avg0 = mean(mcmcs[l,k,*],/double)
              avg1 = mcmcs[l,k,0:n_elements(mcmcs[0,0,*])-(juse+1)]-avg0
              avg2 = mcmcs[l,k,(juse):n_elements(mcmcs[0,0,*])-1]-avg0   
              avg3 = mean((mcmcs[l,k,*]-avg0)^2d0,/double)
              c_j = [c_j,mean(avg1*avg2,/double)/avg3]
           endwhile
           corr_len[l,k]= interpol(jact,c_j,0.5d0)
           eff_len[l,k] = tot_links/corr_len[l,k]
        endif

      ;  stop
     endfor
     final_eff_len[*,*,i] = eff_len
     final_corr_len[*,*,i] = corr_len

     for z=0,n_elements(tap_state)-1,1 do begin
        transit = tap_state[z]->get()   
        for k=0,n_elements(transit.params)-1,1 do begin
           ptr_free,transit.params[k].mcmc_chain
           ptr_free,transit.params[k].refined_mcmc_chain
        endfor
        ptr_free,transit.model_t,transit.model_f,transit.model_tfine,transit.model_ffine,transit.mcmc_files
        transit = 0L
        tap_state[z]->destroy
     endfor
     tap_state=0L
  endfor
  
  mcmcs = 0L
  eff_len = 0L
  corr_len = 0L
  juse=0L
  avg0 = 0L
  avg2 = 0L
  avg3 = 0L
  avg1 = 0L
  c_j = 0L
;  print,max(final_eff_len)
;  stop
  
  return,final_eff_len
end


pro TAP::plotter2d,filename,x,y
  plot_all = 1
  
  tempP = !p
  tempX = !x
  tempY = !y
  
  !p.charsize=1.1
  !p.charthick=4
  !x.thick = 5
  !y.thick = 5
  !p.thick= 5
 ; !p.font = 0

  ymod = 0
 
  xmod=0
 
  for b=0,self.num_transits-1,1 do begin
     if plot_all or b eq 0 then begin 
        
        transit = (*self.transits)[b]->get()
     if transit.params[where(strcmp(((*self.transits)[b]->get()).params.param,x))].fixed eq 0 and $
        transit.params[where(strcmp(((*self.transits)[b]->get()).params.param,y))].fixed eq 0 then begin
        
        ps_start,filename=filename+'.ps'

  case y of  
     'Inclination': yt = 'Inclination [degrees]'
     'a/R*': yt = 'a/R!ds!n'
     'Rp/Rs': yt = 'R!dp!n/R!ds!n'
     'tau_o': yt = textoidl('\tau')+'!do!n'
     'b': yt = 'b'
     else: yt = y
  endcase
  case x of
     'Inclination': xt = 'Inclination [degrees]'
     'a/R*': xt = 'a/R!ds!n'
     'Rp/Rs': xt = 'R!dp!n/R!ds!n'
     'tau_o': xt = textoidl('\tau')+'!do!n'
     'b': xt = 'b'
     else: xt = x
  endcase


        ncont = 40d0
        colors=tap_colors()
        
        !p.charthick = 1.5
        !p.charsize = 1.2
        
        z = tap_hist2d(*transit.params[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain,$
                   *transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain,$
                   binsize1=(max(*transit.params[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain)-$
                             min(*transit.params[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain))/ncont,$
                   binsize2=(max(*transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain)-$
                             min(*transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain))/ncont,$
                   obin1=xv,obin2=yv,binedge1=0,binedge2=0)
        z/=max(z)

        ;stop

;;  xv -= (max(*transit.paramsnnn[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain)-min(*transit.params[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain))/(2*ncont)
;;  yv -= (max(*transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain)-min(*transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain))/(2*ncont)
        contour,z,xv,yv,levels=[.0001],color=colors.black,/nodata,/fill,/path_double,xtitle=xt,ytitle=yt,xrange=mm(*transit.params[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain),yrange=mm(*transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain),xticks=4,/xs,/ys
        oplot,*transit.params[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain,*transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain,psym=8,color=colors.gray,symsize=.2
        contour,z,xv,yv,levels=[.2,.4,.6,.8],color=colors.black,thick=3,/overplot,/path_double ;,/closed
        oploterror,[transit.params[where(strcmp(transit.params.param,x) eq 1)].mcmc_val[0]],[transit.params[where(strcmp(transit.params.param,y) eq 1)].mcmc_val[0]],[transit.params[where(strcmp(transit.params.param,x) eq 1)].mcmc_val[1]],[transit.params[where(strcmp(transit.params.param,y) eq 1)].mcmc_val[1]],errcolor=colors.black,/hibar,psym=8,symsize=1.2,errthick=5,color=colors.black
        oploterror,[transit.params[where(strcmp(transit.params.param,x) eq 1)].mcmc_val[0]],[transit.params[where(strcmp(transit.params.param,y) eq 1)].mcmc_val[0]],[transit.params[where(strcmp(transit.params.param,x) eq 1)].mcmc_val[2]],[transit.params[where(strcmp(transit.params.param,y) eq 1)].mcmc_val[2]],errcolor=colors.black,/lobar,psym=8,symsize=1.2,errthick=5,color=colors.black
        
        sharpcorners,color=colors.black,thick=!x.thick
        ps_end,/png,resize=50
        
        spawn,'mv '+filename+'.png '+self.plot_dir+filename+'_transit'+string(b+1,format='(i2.2)')+'_'+transit.fname+'.png'
        spawn,'rm '+filename+'.ps'
        
       if b eq 0 then   if (*self.plots2d)[0] eq -1 then *self.plots2d = [filename+'_transit'+string(b+1,format='(i2.2)')+'_'+transit.fname] else *self.plots2d = [*self.plots2d,filename+'_transit'+string(b+1,format='(i2.2)')+'_'+transit.fname]
        ;;if strcmp(y,'Eccentricity') or strcmp(x,'Eccentricity') then begin
        ;;   spawn,'open '+self.plot_dir+filename+'.png'
        ;;   stop
        ;;endif
     endif
     
     transit = 0L
     z = 0L
  endif
  endfor
  xv=0L
  yv=0L 
  x = 0L
  y = 0L
  
end
  

pro TAP::plotter2d_basic,filename,x,y
  plot_all = 1
  
  tempP = !p
  tempX = !x
  tempY = !y
  
  !p.charsize=1.1
  !p.charthick=4
  !x.thick = 5
  !y.thick = 5
  !p.thick= 5
 ; !p.font = 0

  ymod = 0
  case y of  
     'Inclination': yt = 'Inclination [degrees]'
     'a/R*': yt = 'a/R!ds!n'
     'Rp/Rs': yt = 'R!dp!n/R!ds!n'
     else: yt = y
  endcase
  
  xmod=0
  case x of
     'Inclination': xt = 'Inclination [degrees]'
     'a/R*': xt = 'a/R!ds!n'
     'Rp/Rs': xt = 'R!dp!n/R!ds!n'
     else: xt = x
  endcase

  for b=0,self.num_transits-1,1 do begin
     if plot_all or b eq 0 then begin 
        
        transit = (*self.transits)[b]->get()

       ; stop

        if stddev(*transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].refined_mcmc_chain) ne 0 and $
           stddev(*transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].refined_mcmc_chain) ne 0 then begin
           
           ps_start,filename=filename+'.ps'
           ncont = 40d0
        colors=tap_colors()
        
        !p.charthick = 1.5
        !p.charsize = 1.2
        
        z = tap_hist2d(*transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].refined_mcmc_chain,$
                   *transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].refined_mcmc_chain,$
                   binsize1=(max(*(transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].refined_mcmc_chain))-$
                             min(*(transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].refined_mcmc_chain)))/ncont,$
                   binsize2=(max(*(transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].refined_mcmc_chain))-$
                             min(*(transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].refined_mcmc_chain)))/ncont,$
                   obin1=xv,obin2=yv,binedge1=0,binedge2=0)
        z/=max(z)
;;  xv -= (max(*transit.paramsnnn[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain)-min(*transit.params[where(strcmp(transit.params.param,x) eq 1)].refined_mcmc_chain))/(2*ncont)
;;  yv -= (max(*transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain)-min(*transit.params[where(strcmp(transit.params.param,y) eq 1)].refined_mcmc_chain))/(2*ncont)
        contour,z,xv,yv,levels=[.0001],color=colors.black,/nodata,/fill,/path_double,xtitle=xt,ytitle=yt,$
                xrange=mm(*(transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].refined_mcmc_chain)),$
                yrange=mm(*(transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].refined_mcmc_chain)),$
                xticks=4,/xs,/ys
        oplot,*(transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].refined_mcmc_chain),$
              *(transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].refined_mcmc_chain),psym=8,color=colors.gray,symsize=.2
        contour,z,xv,yv,levels=[.2,.4,.6,.8],color=colors.black,thick=3,/overplot,/path_double ;,/closed
        oploterror,[transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].mcmc_val[0]],$
                   [transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].mcmc_val[0]],$
                   [transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].mcmc_val[1]],$
                   [transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].mcmc_val[1]],errcolor=colors.black,$
                   /hibar,psym=8,symsize=1.2,errthick=5,color=colors.black
        oploterror,[transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].mcmc_val[0]],$
                   [transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].mcmc_val[0]],$
                   [transit.basic_params[where(strcmp(transit.basic_params.param,x) eq 1)].mcmc_val[2]],$
                   [transit.basic_params[where(strcmp(transit.basic_params.param,y) eq 1)].mcmc_val[2]],errcolor=colors.black,$
                   /lobar,psym=8,symsize=1.2,errthick=5,color=colors.black
        
        sharpcorners,color=colors.black,thick=!x.thick
        ps_end,/png,resize=50
        
        spawn,'mv '+filename+'.png '+self.plot_dir+filename+'_transit'+string(b+1,format='(i2.2)')+'_'+transit.fname+'.png'
        spawn,'rm '+filename+'.ps'
        
       if b eq 0 then if (*self.plots2d)[0] eq -1 then *self.plots2d = [filename+'_transit'+string(b+1,format='(i2.2)')+'_'+transit.fname] else *self.plots2d = [*self.plots2d,filename+'_transit'+string(b+1,format='(i2.2)')+'_'+transit.fname]
        ;;if strcmp(y,'Eccentricity') or strcmp(x,'Eccentricity') then begin
        ;;   spawn,'open '+self.plot_dir+filename+'.png'
        ;;   stop
        ;;endif
     endif
     
     transit = 0L
     z = 0L
  endif
  endfor
  xv=0L
  yv=0L 
  x = 0L
  y = 0L
  
  end
  
pro TAP::create_2d
  
  self.message = 'Creating 2-D MCMC distribution plots...'
  self->message
  ptr_free,self.plots2d
  self.plots2d = ptr_new(-1)
  ;; numplot=0

  x = 'Inclination'
  y = 'a/R*'
  filename = 'p2D_inc_aRs'
  self.message = '  1/5: '+x+' vs '+y+'... '+filename
  self->message
  self->plotter2d_basic,filename,x,y
  ;if numplot eq 0 then *self.plots2d = [filename] else *self.plots2d = [*self.plots2d,filename]
  ;numplot++
  
  x = 'Inclination'
  y = 'Rp/R*'
  filename = 'p2D_inc_RpRs'
  self.message = '  2/5: '+x+' vs '+y+'... '+filename
  self->message
  self->plotter2d_basic,filename,x,y
 ; if numplot eq 0 then *self.plots2d = [filename] else *self.plots2d = [*self.plots2d,filename]
 ; numplot++
  
  x = 'Inclination'
  y = 'Eccentricity'
  filename = 'p2D_inc_ecc'
  self.message = '  3/5: '+x+' vs '+y+'... '+filename
  self->message
  self->plotter2d_basic,filename,x,y
 ; if numplot eq 0 then *self.plots2d = [filename] else *self.plots2d = [*self.plots2d,filename]
 ; numplot++
  
  x = 'Eccentricity'
  y = 'Omega'
  filename = 'p2D_ecc_omega'
  self.message = '  4/5: '+x+' vs '+y+'... '+filename
  self->message
  self->plotter2d_basic,filename,x,y
 ;  if numplot eq 0 then *self.plots2d = [filename] else *self.plots2d = [*self.plots2d,filename]
 ;  numplot++
   
   x = 'Rp/R*'
   y = 'a/R*'
   filename = 'p2D_RpRs_aRs'
   self.message = '  5/5: '+x+' vs '+y+'... '+filename
   self->message
  self->plotter2d_basic,filename,x,y


   x = 'Linear LD'
   y = 'Quadratic LD'
   filename = 'p2D_u1_u2'
   self.message = ' 6/5: '+x+' vs '+y+'... '+filename
   self->message
  self->plotter2d_basic,filename,x,y
;   if numplot eq 0 then *self.plots2d = [filename] else *self.plots2d = [*self.plots2d,filename]
;  numplot++
    
  if strcmp(self.parameterize_id,'adv1') then begin
     x = 'b'
     y = 'tau_o'
     filename = 'p2D_b_tauo'
     self.message = ' 7/5: '+x+' vs '+y+'... '+filename
     self->message
     self->plotter2d,filename,x,y
  endif

  if strcmp(self.parameterize_id,'adv2') then begin
     x = 'b'
     y = 'T'
     filename = 'p2D_b_T'
     self.message = ' 7/5: '+x+' vs '+y+'... '+filename
     self->message
     self->plotter2d,filename,x,y
  endif

end


pro TAP::create_ascii
  self.message = 'Writing MCMC results to ascii files...'
  self->message

  ;help,r
  ;stop
 ; help,r
 ; help,eff_len

;  stop
  openw,strlun,self.plot_dir+'MCMC_tables.txt',/get_lun,width=2500,bufsize=0
  printf,strlun,"::Overview of TAP MCMC Parameters:"
  printf,strlun,"TAP version "+self.version  
  printf,strlun,"TAPmcmc version "+((*self.transits)[0]->get()).mcmc_version
  printf,strlun,"MCMC Chains "+string(self.mcmc_complete,format='(i)')
  printf,strlun,"Initial Chain Length "+string((((*self.transits)[self.active_transit-1])->get()).mcmc_params[1],format='(i)')
  printf,strlun,"Total Inference Links "+string(10d0*(n_elements(*((*self.transits)[0]->get()).params[0].refined_mcmc_chain)),format='(i)')
  
  transit = (*self.transits)[0]->get()
  case transit.ldtype of 
     2: begin
        ldt1 = 'Drawn from ($\mu_{1}$+$\mu_{2}$)$^2$ and 0.5$\mu_{1}(\mu_{1} + \mu_{2})^{-1}$' 
        ldt = 'Drawn from (u1 + u2)^2 and 0.5u1*(u1+u1)^-1'
     end
     1: begin
        ldt1 = 'Drawn from 2$\mu_{1}$+$\mu_{2}$ and $\mu_{1}$-2$\mu_{2}$' 
        ldt = 'Drawn from 2u1+u2 and u1-2u2' 
     end
     0: begin
        ldt1 = 'Drawn from $\mu_{1}$ and $\mu_{2}$' 
        ldt = 'Drawn from u1 and u2'
     end
  endcase
  
  transit =  0L
  printf,strlun,"LD MCMC Distribution"+ldt
 
  rebin=[-1]
  texp = [-1]
  texp1 = [-1]
  for i=0,n_elements(*self.transits)-1,1 do begin
     rebin = [rebin,((*self.transits)[i]->get()).rebin]
     texp  = [texp, ((*self.transits)[i]->get()).t_int[0]]
     texp1  = [texp1, ((*self.transits)[i]->get()).t_int[1]]
  endfor

  rebin_txt1 = ''
  rebin_txt2 = ''
  rebin_tbl = 0
  ;stop

  if n_elements(uniq((rebin[1:self.num_transits])[sort(rebin[1:self.num_transits])])) eq 1 then begin
     case rebin[1] of 
        0: rebin_txt = 'No Re-sampling'
        1: begin
           rebin_txt = 'Resample and Rebin Mandel \& Agol'
           case robust_sigma(texp[1:self.num_transits]) of
              0: case robust_sigma(texp1[1:self.num_transits]) of
                 0: rebin_txt1 = string(texp[1],format='(d9.6)')+' minutes, '+string(texp1[1],format='(i2.2)')+' Samples'
              endcase
              else: begin
                 rebin_txt1 = 'Table \ref{tbl:rebin} values.'
                 rebin_tbl = 1
              end
           endcase
        end
     endcase
  endif else begin
     rebin_txt = 'Varied, see Table \ref{tbl:rebin}'
     rebin_tbl=1
  endelse
  
  printf,strlun,"Long Int Mode:  "+rebin_txt
  if strcmp(rebin_txt1,'') eq 0 then printf,strlun,"Resample technique  "+rebin_txt1

  for i=0,n_elements(*self.transits)-1,1 do printf,strlun,"Transit "+ $
     string(i+1,format='(i2.2)')+': '+((*self.transits)[i]->get()).fname

  printf,strlun,""
  printf,strlun,""

  ;print,rebin_tbl
  ;stop
  if rebin_tbl then begin

     printf,strlun,"::Overview of TAP MCMC Rebin Parameters"
     printf,strlun,"                 Transit                                    Mode                               Technique"
     
     for i=0,n_elements(*self.transits)-1,1 do begin
        transit = string(i+1,format='(i2.2)')+': '+((*self.transits)[i]->get()).fname
        case rebin[i+1] of
           0: begin
              mode = 'No Re-sampling'
              technique = '...'
           end
           1: begin
              mode = 'Resample and Rebin Mandel & Agol'
              technique = string(texp[i+1],format='(d9.6)')+' minutes, '+string(texp1[i+1],format='(i2.2)')+' Samples'
           end
        endcase
        printf,strlun,transit,mode,technique,format='(a40,a30,a30)'
     endfor
 
  printf,strlun,""
  printf,strlun,""
endif

  for i=0,n_elements(*self.transits)-1,1 do begin
     if i eq 0 then begin
        header = [string(i+1,format='(i2.2)')+": "+strmid(((*self.transits)[i]->get()).fname,0,15)]
        format = '(a15,a50'
     endif else begin
        header = [header,string(i+1,format='(i2.2)')+": "+strmid(((*self.transits)[i]->get()).fname,0,15)]
        format += ',a50'
     endelse
  endfor
  format += ')'
  printf,strlun,"::Wavelet basis red noise MCMC analysis"
  printf,strlun,"  each column is [ Parameter value +1sig -1sig) ]"
  printf,strlun,'Parameter',header,format=format
  for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
     line = string(((*self.transits)[0]->get()).params[i].param,format='(A15)')
     for j=0,n_elements(*self.transits)-1,1 do begin
        line += string(self->makeline2(((*self.transits)[j]->get()).params[i].mcmc_val),format='(a50)')
     endfor
     printf,strlun,line  
  endfor

  printf,strlun,""
  printf,strlun,""  
  printf,strlun,"::Multi Curve MCMC Parameter Jump Rates"
  printf,strlun,"  each column is [ Jump rate (requested rate) ]"
  for i=0,n_elements(*self.transits)-1,1 do begin
     if i eq 0 then begin
        header = [string(i+1,format='(i2.2)')+": "+strmid(((*self.transits)[i]->get()).fname,0,15)]
        format = '(a15,a20'
     endif else begin
        header = [header,string(i+1,format='(i2.2)')+": "+strmid(((*self.transits)[i]->get()).fname,0,15)]
        format += ',a20'
     endelse
  endfor
  format += ')'
  printf,strlun,'Parameter',header,format=format
  for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
     line = string(((*self.transits)[0]->get()).params[i].param,format='(A15)')
     for j=0,n_elements(*self.transits)-1,1 do begin
        if ((*self.transits)[j]->get()).params[i].mcmc_val[1] eq -1 then $
           line += string('[fixed]',format='(a20)') else $
              line +=  string(string((((*self.transits)[j]->get()).params[i].jumpct/((*self.transits)[j]->get()).params[i].jumptot),format='(d4.2)') +$
                      ' ('+string(((*self.transits)[j]->get()).params[i].accept,format='(d4.2)')+')',format='(a20)')
     endfor
     printf,strlun,line  
  endfor
  
  printf,strlun,""
  printf,strlun,""  
  printf,strlun,"::Multi Curve MCMC Parameter Lock Matrix"
  printf,strlun,"  each column is parameter set"
  for i=0,n_elements(*self.transits)-1,1 do begin
     if i eq 0 then begin
        header = ['Transit '+string(i+1,format='(i2.2)')]
        format = '(a15,a12'
     endif else begin
        header = [header,string('Transit '+string(i+1,format='(i2.2)'),format='(a12)')]
        format += ',a12'
     endelse
  endfor
  format += ')'
  printf,strlun,'Parameter',header,format=format
  for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
     line = string(((*self.transits)[0]->get()).params[i].param,format='(A15)')
     for j=0,n_elements(*self.transits)-1,1 do begin
        if ((*self.transits)[j]->get()).params[i].mcmc_val[1] eq -1 then $
           line += string(string(((*self.transits)[j]->get()).params[i].set,format='(i3)')+' [fixed]',format='(a12)') else $
              line +=  string(((*self.transits)[j]->get()).params[i].set,format='(i12)')
     endfor
     printf,strlun,line  
  endfor

  close,strlun
  
  

  openw,strlun,self.plot_dir+'MCMC_tables.tex',/get_lun,width=2500,bufsize=0
  printf,strlun,'\documentclass[preprint,10pt]{aastex}'
  printf,strlun,'\usepackage{graphicx}'
  printf,strlun,'\usepackage{epstopdf}'
  printf,strlun,'\usepackage{morefloats}'
  printf,strlun,'\begin{document}'
  printf,strlun,""
  printf,strlun,"\begin{deluxetable}{lc}"  
  printf,strlun,"\tablewidth{0pt}"
  printf,strlun,"\tablecaption{Overview of TAP MCMC Parameters}"
  printf,strlun,"\tablehead{"
  printf,strlun,"\colhead{Parameter} & \colhead{Value}}"
  printf,strlun,"\startdata"
  printf,strlun,"TAP version & "+self.version+'\\'  
  printf,strlun,"TAPmcmc version & "+((*self.transits)[0]->get()).mcmc_version+'\\'
  printf,strlun,'\hline'
  printf,strlun,"MCMC Chains & "+string(self.mcmc_complete,format='(i)')+'\\'
  printf,strlun,"Chain Length & "+string((((*self.transits)[self.active_transit-1])->get()).mcmc_params[1],format='(i)')+'\\'
  printf,strlun,"Total Inference Links & "+string(10d0*(n_elements(*((*self.transits)[0]->get()).params[0].refined_mcmc_chain)), $
                                                  format='(i)')+'\\ \hline'
  printf,strlun,"LD MCMC Distribution &"+ldt1+'\\'

  rebin=[-1]
  texp = [-1]
  texp1 = [-1]
  for i=0,n_elements(*self.transits)-1,1 do begin
     rebin = [rebin,((*self.transits)[i]->get()).rebin]
     texp  = [texp, ((*self.transits)[i]->get()).t_int[0]]
     texp1  = [texp1, ((*self.transits)[i]->get()).t_int[1]]
  endfor
  
  rebin_txt1 = ''
  rebin_txt2 = ''
  rebin_tbl = 0
  if n_elements(uniq((rebin[1:self.num_transits])[sort(rebin[1:self.num_transits])])) eq 1 then begin
     case rebin[1] of 
        0: rebin_txt = 'No Re-sampling'
        1: begin
           rebin_txt = 'Resample and Rebin Mandel \& Agol'
           case robust_sigma(texp[1:self.num_transits]) of
              0: case robust_sigma(texp1[1:self.num_transits]) of
                 0: rebin_txt1 = string(texp[1],format='(d9.6)')+' minutes, '+string(texp1[1],format='(i2.2)')+' Samples'
              endcase
              else: begin
                 rebin_txt1 = 'Table \ref{tbl:rebin} values.'
                 rebin_tbl = 1
              end
           endcase
        end
     endcase
  endif else begin
     rebin_txt = 'Varied, see Table \ref{tbl:rebin}'
     rebin_tbl=1
  endelse
  
  printf,strlun,"Long Int Mode: & "+rebin_txt+'\\'
  if strcmp(rebin_txt1,'') eq 0 then printf,strlun,"Resample technique & "+rebin_txt1+'\\' 

  printf,strlun,'\hline'

  for i=0,n_elements(*self.transits)-1,1 do printf,strlun,"Transit "+ $
     string(i+1,format='(i2.2)')+' & '+self->stripspecial(((*self.transits)[i]->get()).fname)+'\\'
  printf,strlun,"\enddata"  
  printf,strlun,"\tablecomments{TAP MCMC "+self.version+"}"
  printf,strlun,"\label{tbl:mcmcpar}"
  printf,strlun,"\end{deluxetable}"
  printf,strlun,""

  if rebin_tbl then begin
     printf,strlun,""
     printf,strlun,"\begin{deluxetable}{lcc}"  
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Overview of TAP MCMC Rebin Parameters}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Transit} & \colhead{Mode} & \colhead{Technique}}"
     printf,strlun,"\startdata"
     for i=0,n_elements(*self.transits)-1,1 do begin
        transit = string(i+1,format='(i2.2)')+': '+self->stripspecial(((*self.transits)[i]->get()).fname)
        case rebin[i+1] of
           0: begin
              mode = 'No Re-sampling'
              technique = '\nodata'
           end
           1: begin
              mode = 'Resample and Rebin Mandel \& Agol'
              technique = string(texp[i+1],format='(d9.6)')+' minutes, '+string(texp1[i+1],format='(i2.2)')+' Samples'
           end
        endcase
        printf,strlun,transit+' & '+mode+' & '+technique+' \\'
     endfor
   ;  printf,strlun,"\hline"
     

   ;  stop
     printf,strlun,"\enddata"  
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+"}"
     printf,strlun,"\label{tbl:rebin}"
     printf,strlun,"\end{deluxetable}"
     printf,strlun,""
  endif
  mode = 0L
  technique = 0L
  transit = 0L
  rebin_txt = 0L
  rebin_txt1 = 0L
  rebin = 0L
  texp =  0L
  
  
  transit = (*self.transits)[0]->get()

  
  printf,strlun,""
  printf,strlun,"\begin{deluxetable}{cl}"  
  printf,strlun,"\tablewidth{0pt}"
  printf,strlun,"\tablecaption{TAP MCMC Parameterization}"
  printf,strlun,"\tablehead{"
  printf,strlun,"\colhead{Parameter} & \colhead{Notes}}"
  printf,strlun,"\startdata"
  printf,strlun,"Period  &   [days]  \\"
  case transit.parameterize_id of
     'basic': begin
        printf,strlun,"Inclination & [degrees] \\" 
        printf,strlun," /n $\frac{a}{R_*}$  &  \\"
     end
     'adv2': begin
        printf,strlun,"b &  $= \left[ \frac{a}{R_*}  cos(i) \frac{1-e^2}{1+esin(\omega)}  \right] $ \\"
        printf,strlun," T  &  $ = 2\frac{R_*}{a} \frac{P}{\pi} \frac{1-e^2}{1+esin(\omega)} \sqrt{1-b^2} $ \\"
     end
  endcase
  printf,strlun,"$\frac{R_p}{R_*}$  &   \\"
  printf,strlun,"Mid-Transit (T$_{\rm{mid}}$)  &  [days] \\"
  printf,strlun,"    &      \\"
  printf,strlun,"    &      \\"
  printf,strlun,"    &      \\"
  printf,strlun,"    &      \\"
  printf,strlun,"    &      \\"
  printf,strlun,"    &      \\"
  printf,strlun,"    &      \\"
  printf,strlun,"    &      \\"
  printf,strlun,"    &      \\"
  printf,strlun,"\enddata"  
  printf,strlun,"\tablecomments{TAP MCMC "+self.version+"}"
  printf,strlun,"\label{tbl:rebin}"
  printf,strlun,"\end{deluxetable}"
  printf,strlun,""
  
transit = 0L


  if self.calc_gr then begin
     r = self->gelmanrubin()
     
   ;;;  stop

     write_table = [-1,-1]
     while write_table[1] lt n_elements(*self.transits)-1 do begin
        
        format = ''
        fixed = 1
        
        write_table[0] = write_table[1]+1
        write_table[1] = min([write_table[0]+5,n_elements(*self.transits)-1])
        for i=write_table[0],write_table[1],1 do format = format+'c'
        printf,strlun,""
        printf,strlun,"\begin{deluxetable}{l"+format+"}"  
                                ; if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
        printf,strlun,"\tablewidth{0pt}"
        printf,strlun,"\tablecaption{Gelman-Rubin statistic for non-convergence.}"
        printf,strlun,"\tablehead{"
        printf,strlun,"\colhead{Parameter}"
        printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{Value} \\"
        for i=write_table[0],write_table[1],1 do printf,strlun,$
           "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,5)+" }"
        printf,strlun,"}"
        printf,strlun,"\startdata"
        for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(self->stripspecial(((*self.transits)[0]->get()).params[i].param),format='(A15)')
           for j=write_table[0],write_table[1],1 do begin
              if r[i,j] eq 0 then begin
                 line += ' & \nodata\tablenotemark{a}' 
              endif else begin
                 if r[i,j] gt 1.1d0 then extra = ['\textbf{','}'] else extra = ['','']
                 line += ' & '+extra[0]+string(r[i,j],format='(d6.3)')+extra[1]
              endelse
           endfor
           line += '\\'
           printf,strlun,line  
        endfor
        
        printf,strlun,"\enddata"
        if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
        if max(r) gt 1.10d0 then converge_string = 'Gelman-Rubin statistic shows evidence for non-convergence.' $
        else converge_string = 'Gelman-Rubin statistic shows no evidence for non-convergence.'
        
        
        printf,strlun,"\tablecomments{TAP MCMC "+self.version+"  "+converge_string+" }"
        printf,strlun,"\label{tbl:tapmcmc1}"
        printf,strlun,"\end{deluxetable}"
        
     endwhile
     r=0L
  endif
  
  if self.calc_effl then begin
     eff_len = self->effective_length(final_corr_len=corr_len)
     selector = [1,10,100,1000,10000,100000,1000000]
                                ; stop
     
     write_table = [-1,-1]
     while write_table[1] lt n_elements(*self.transits)-1 do begin    
        format = ''
        fixed = 1
        
        write_table[0] = write_table[1]+1
        write_table[1] = min([write_table[0]+5,n_elements(*self.transits)-1])
        for i=write_table[0],write_table[1],1 do format = format+'c'
        printf,strlun,""
        printf,strlun,"\begin{deluxetable}{l"+format+"}"  
                                ; if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
        printf,strlun,"\tablewidth{0pt}"
        printf,strlun,"\tablecaption{MCMC Chain Effective Lengths (Correlation Length)}"
        printf,strlun,"\tablehead{"
        printf,strlun,"\colhead{Parameter}"
        printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{Value} \\"
        for i=write_table[0],write_table[1],1 do printf,strlun,$
           "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,5)+" }"
        printf,strlun,"}"
        printf,strlun,"\startdata"
        for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(self->stripspecial(((*self.transits)[0]->get()).params[i].param),format='(A15)')
           for j=write_table[0],write_table[1],1 do begin
              if finite(total(eff_len[j,i,*])) eq 0 then begin
                 line += ' & \nodata\tablenotemark{a}' 
              endif else begin
                 ;;   if r[i,j] gt 1.1d0 then extra = ['\textbf{','}'] else extra = ['','']
                 extra = ['','']
            ;;     doit = (where((corr_len[j,i,*] mod selector) eq corr_len[j,i,*]))[0]
                 format = '(i8)'
                 line += ' & '+extra[0]+string(total(eff_len[j,i,*]),format='(i10)')+" "+$
                         strjoin(strsplit("("+string(total(corr_len[j,i,*]),format=format)+")",' ',/extract),'')
              endelse
           endfor
           line += '\\'
           printf,strlun,line  
        endfor
        printf,strlun,"\enddata"
        if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
        string = 'Effective MCMC lengths.  Lengths $\gg$ 1 assure that the MCMC is well mixed.'
        
        printf,strlun,"\tablecomments{TAP MCMC "+self.version+" Effective MCMC lengths.  Lengths $\gg$ 1 assure that the MCMC is well mixed.}"
        printf,strlun,"\label{tbl:tapmcmc1}"
        printf,strlun,"\end{deluxetable}"
        
     endwhile
     eff_len = 0L
  endif
                                ;stop
  
  if self.create_plots then begin
     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_phase_modandresid.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+'.  Global model (blue solid line) and all data overplotted.  For multiple light curves, data points represent the residuals from each light curves model.}'
     printf,strlun,'\end{figure*}'

     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_phased_lc.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+'.  All data plotted as (global model) + residuals}'
     printf,strlun,'\end{figure*}'

     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_lightcurve.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+'}'
     printf,strlun,'\end{figure*}'

     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_mcmc_params1.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+' parameter distributions.}'
     printf,strlun,'\end{figure*}'

     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_mcmc_params2.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+' parameter distributions.}'
     printf,strlun,'\end{figure*}'

     printf,strlun,'\begin{figure*}[h]'
     printf,strlun,'\centering'
     printf,strlun,'\includegraphics[]{'+self.plot_dir+'plot_alltransit_mcmc_params3.eps}'
     printf,strlun,'\caption{TAP MCMC '+self.version+' parameter distributions.}'
     printf,strlun,'\end{figure*}'
  endif

  if self.create_2d then begin
    ; stop
     for i=0,n_elements(*self.plots2d)-1,1 do begin
        printf,strlun,'\begin{figure*}[h]'
        printf,strlun,'\centering'
        printf,strlun,'\includegraphics[width=6in]{'+self.plot_dir+(*self.plots2d)[i]+'.png}'
        printf,strlun,'\caption{TAP MCMC '+self.version+': plot file is '+self->stripspecial((*self.plots2d)[i])+'.png}'
        printf,strlun,'\end{figure*}'
        printf,strlun,''
     endfor
  endif

  write_table = [-1,-1]
  while write_table[1] lt n_elements(*self.transits)-1 do begin

     format = ''
     fixed = 1
     
     write_table[0] = write_table[1]+1
     write_table[1] = min([write_table[0]+2,n_elements(*self.transits)-1])
     for i=write_table[0],write_table[1],1 do format = format+'c'
     printf,strlun,""
     printf,strlun,"\begin{deluxetable}{l"+format+"}"  
                                ; if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Wavelet basis red noise MCMC analysis}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Parameter}"
     printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{Value} \\"
     for i=write_table[0],write_table[1],1 do printf,strlun,$
        "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,15)+" }"
     printf,strlun,"}"
     printf,strlun,"\startdata"
     for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(self->stripspecial(((*self.transits)[0]->get()).params[i].param),format='(A15)')
        for j=write_table[0],write_table[1],1 do begin
           line += self->makeline(((*self.transits)[j]->get()).params[i].mcmc_val)
        endfor
        line += '\\'
        printf,strlun,line  
     endfor
     basic = 0
     if strcmp(self.parameterize_id,'basic') eq 0 then begin
        basic = 1
        printf,strlun,'\\ \hline \hline \\'
        i = where(strcmp('Inclination',((*self.transits)[0]->get()).basic_params.param))
        line = string(self->stripspecial(((*self.transits)[0]->get()).basic_params[i].param),format='(A15)')+ '\tablenotemark{b} ' 
        for j=write_table[0],write_table[1],1 do begin
           line += self->makeline(((*self.transits)[j]->get()).basic_params[i].mcmc_val) 
        endfor
        line += '\\'
        printf,strlun,line  
        i = where(strcmp('a/R*',((*self.transits)[0]->get()).basic_params.param))
        line = string(self->stripspecial(((*self.transits)[0]->get()).basic_params[i].param),format='(A15)') + ' \tablenotemark{b} ' 
        for j=write_table[0],write_table[1],1 do begin
           line += self->makeline(((*self.transits)[j]->get()).basic_params[i].mcmc_val) 
        endfor
        line += '\\'
        printf,strlun,line  
      
     endif

     printf,strlun,"\enddata"
     if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
     if basic then printf,strlun,"\tablenotetext{b}{Value Inferred from Parameterization.}"
     
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+".  Parameters calcualted the MCMC run, or inferred from MCMC parameters.}"
     printf,strlun,"\label{tbl:tapmcmc1}"
     printf,strlun,"\end{deluxetable}"
     
  endwhile
  
  if 0 then $
     if strcmp(self.parameterize_id,'basic') eq 0 then begin

     write_table = [-1,-1]
     while write_table[1] lt n_elements(*self.transits)-1 do begin
        format = ''
        fixed = 1
        
        write_table[0] = write_table[1]+1
        write_table[1] = min([write_table[0]+2,n_elements(*self.transits)-1])
        for i=write_table[0],write_table[1],1 do format = format+'c'
        printf,strlun,""
        printf,strlun,"\begin{deluxetable}{l"+format+"}"  
                                ; if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
        printf,strlun,"\tablewidth{0pt}"
        printf,strlun,"\tablecaption{Wavelet basis red noise MCMC analysis}"
        printf,strlun,"\tablehead{"
        printf,strlun,"\colhead{Parameter}"
        printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{Value} \\"
        for i=write_table[0],write_table[1],1 do printf,strlun,$
           "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,15)+" }"
        printf,strlun,"}"
        printf,strlun,"\startdata"
        for i=0,n_elements(((*self.transits)[0]->get()).basic_params)-1,1 do begin 
           line = string(self->stripspecial(((*self.transits)[0]->get()).basic_params[i].param),format='(A15)')
           for j=write_table[0],write_table[1],1 do begin
              line += self->makeline(((*self.transits)[j]->get()).basic_params[i].mcmc_val) 
           endfor
           line += '\\'
           printf,strlun,line  
        endfor
        
        printf,strlun,"\enddata"
        if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
        
        printf,strlun,"\tablecomments{TAP MCMC "+self.version+".  System parameters inferred from parameterization.}"
        printf,strlun,"\label{tbl:tapmcmc1}"
        printf,strlun,"\end{deluxetable}"
        
     endwhile
  endif


;;  stop

  spawn,'ls '+self.plot_dir+'/MCMC_chains/',chainlist
  restore,self.plot_dir+'/MCMC_chains/'+chainlist[0]
 
  write_table = [-1,-1]
  while write_table[1] lt n_elements(*self.transits)-1 do begin
     
     format = ''
     fixed = 1
     
     write_table[0] = write_table[1]+1
     write_table[1] = min([write_table[0]+5,n_elements(*self.transits)-1])
     for i=write_table[0],write_table[1],1 do format = format+'c'
     printf,strlun,""
     ;; PARAMETER LOCKS
     printf,strlun,"\begin{deluxetable}{l"+format+"}"    ;;if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Multi Curve MCMC Parameter Lock Matrix}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Parameter}"
     printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{$\beta$ (Characteristic Jump Size)} \\"
     for i=write_table[0],write_table[1],1 do printf,strlun,$
  "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,5)+"... }"
     printf,strlun,"}"
    
     printf,strlun,"\startdata"
     for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(self->stripspecial(((*self.transits)[0]->get()).params[i].param),format='(A15)')
        for j=write_table[0],write_table[1],1 do begin
           if ((*self.transits)[j]->get()).params[i].mcmc_val[1] eq -1 then $
              line += '& \nodata' else $
                 line += self->makeline3((tap_state[j]->get()).params[i].beta)
        endfor
        line += '\\'
        printf,strlun,line  
     endfor
     
     printf,strlun,"\enddata"
  ;   if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
     
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+".  Characteristic step sizes for the first MCMC chain.}"
     printf,strlun,"\label{tbl:tapmcmc_betas}"
     printf,strlun,"\end{deluxetable}"
     
  endwhile

  for z=0,n_elements(tap_state)-1,1 do begin
     transit = tap_state[z]->get()   
     for k=0,n_elements(transit.params)-1,1 do begin
        ptr_free,transit.params[k].mcmc_chain
        ptr_free,transit.params[k].refined_mcmc_chain
     endfor
     ptr_free,transit.model_t,transit.model_f,transit.model_tfine,transit.model_ffine,transit.mcmc_files
     transit = 0L
     tap_state[z]->destroy
  endfor
  tap_state=0L
  



  write_table = [-1,-1]
  while write_table[1] lt n_elements(*self.transits)-1 do begin
     
     format = ''
     fixed = 1
     
     write_table[0] = write_table[1]+1
     write_table[1] = min([write_table[0]+5,n_elements(*self.transits)-1])
     for i=write_table[0],write_table[1],1 do format = format+'c'
     printf,strlun,""
     ;; PARAMETER LOCKS
     printf,strlun,"\begin{deluxetable}{l"+format+"}"    ;;if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Multi Curve MCMC Parameter Lock Matrix}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Parameter}"
     printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{Jump Rate (Requested)} \\"
     for i=write_table[0],write_table[1],1 do printf,strlun,$
  "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,5)+"... }"
     printf,strlun,"}"
     printf,strlun,"\startdata"
     for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(self->stripspecial(((*self.transits)[0]->get()).params[i].param),format='(A15)')
        for j=write_table[0],write_table[1],1 do begin
           if ((*self.transits)[j]->get()).params[i].mcmc_val[1] eq -1 then $
              line += '& \nodata' else $
                 line += ' & ' + string((((*self.transits)[j]->get()).params[i].jumpct/((*self.transits)[j]->get()).params[i].jumptot),format='(d4.2)') +$
                         ' ('+string(((*self.transits)[j]->get()).params[i].accept,format='(d4.2)')+')'
        endfor
        line += '\\'
        printf,strlun,line  
     endfor
     
     printf,strlun,"\enddata"
  ;   if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
     
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+".  Jump rates (with request).}"
     printf,strlun,"\label{tbl:tapmcmc2}"
     printf,strlun,"\end{deluxetable}"
     
  endwhile


  
  write_table = [-1,-1]
  while write_table[1] lt n_elements(*self.transits)-1 do begin
     
     format = ''
     fixed = 1
     
     write_table[0] = write_table[1]+1
     write_table[1] = min([write_table[0]+5,n_elements(*self.transits)-1])
     for i=write_table[0],write_table[1],1 do format = format+'c'
     printf,strlun,""
     ;; PARAMETER LOCKS
     printf,strlun,"\begin{deluxetable}{l"+format+"}"    ;;if n_elements(*self.transits) gt 4 then printf,strlun,"\rotate"
     printf,strlun,"\tablewidth{0pt}"
     printf,strlun,"\tablecaption{Multi Curve MCMC Parameter Lock Matrix}"
     printf,strlun,"\tablehead{"
     printf,strlun,"\colhead{Parameter}"
     printf,strlun,"& \multicolumn{"+string(write_table[1]-write_table[0]+1,format='(i2.2)')+"}{c}{MCMC Parameter Set} \\"
     for i=write_table[0],write_table[1],1 do printf,strlun,$
  "   & \colhead{"+string(i+1,format='(i2.2)')+": "+strmid(self->stripspecial(((*self.transits)[i]->get()).fname),0,5)+"... }"
     printf,strlun,"}"
     printf,strlun,"\startdata"
     for i=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin 
        line = string(self->stripspecial(((*self.transits)[0]->get()).params[i].param),format='(A15)')
        for j=write_table[0],write_table[1],1 do begin
           line += ' & ' + string(((*self.transits)[j]->get()).params[i].set,format='(i3)')
           if ((*self.transits)[j]->get()).params[i].mcmc_val[1] eq -1 then line += '\tablenotemark{a}' 
        endfor
        line += '\\'
        printf,strlun,line  
     endfor
     
     printf,strlun,"\enddata"
     if fixed then printf,strlun,"\tablenotetext{a}{Value Fixed in MCMC Analysis.}"
     
     printf,strlun,"\tablecomments{TAP MCMC "+self.version+".  Transits with the same values for any parameter row are locked together in the MCMC analysis.}"
     printf,strlun,"\label{tbl:tapmcmc2}"
     printf,strlun,"\end{deluxetable}"
     
     endwhile




  printf,strlun,"\end{document}"
  close,strlun
  
  openw,strlun,self.plot_dir+'ascii_phased_data.ascii',/get_lun,width=2500,bufsize=0
  printf,strlun,'#  Phase [Days from Tmid]    Flux_corr      Raw_Flux  Orig_T_days     model_flux  Residual'
  printf,strlun,"#   tap_readcol,'"+self.plot_dir+"ascii_phased_data.ascii',p,f,rf,t,m,r,format='(d,d,d,d,d,d)'"
  
  compare = 'OOT t^0'
  compare2 = 'OOT t^1'
  compare3 = 'OOT t^2'

  for i=0,self.num_transits-1,1 do begin
     lc = (*self.transits)[i]->get() 
     midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
     fit = ptr_new(poly((lc.transit_tday-min(lc.transit_tday)),[lc.params[where(strcmp(lc.params.param,compare) eq 1)].value,lc.params[where(strcmp(lc.params.param,compare2) eq 1)].value,lc.params[where(strcmp(lc.params.param,compare3) eq 1)].value]))
     for j=0,n_elements(lc.transit_tday)-1,1 do begin
        
        printf,strlun,(lc.transit_tday[j]-midt),lc.transit_flux[j]/(*fit)[j],lc.transit_flux[j],lc.transit_tday[j],(*lc.model_f)[j],lc.transit_flux[j]-(*lc.model_f)[j],format='(d20.8,d20.8,d20.8,d25.8,d20.8,d20.8)'
     endfor
     ptr_free,fit
  endfor
  close,strlun
  compare = 0L
  compare2 = 0L
  


  
  transit = (*self.transits)[0]->get()
  tpar = transit.params.value
  tpar[9] = 1d0
  tpar[10] = 0d0
  tpar[11] = 0d0

  dum = TAP_MPFit_function(tpar,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,redfit=redfit,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int,adv_opt=transit.adv_opt)
  
  openw,strlun,self.plot_dir+'ascii_phased_model.ascii',/get_lun,width=2500,bufsize=0
  printf,strlun,"#   tap_readcol,'"+self.plot_dir+"ascii_phased_model.ascii',t,p,mf,format='(d,d,d)'"
  printf,strlun,'# Phased from '+transit.fname+' Parameters'
  printf,strlun,'#      T_days         Phase [hours from Tmid]     Model Flux    '
  
  tpar = ((*self.transits)[0]->get()).params.value

  midt = transit.params[where(strcmp(transit.params.param,'Mid Transit') eq 1)].value - (transit.epoch*transit.params[0].value)
  for j=0,n_elements(ffit)-1,1 do printf,strlun,(*transit.model_tfine)[j],((*transit.model_tfine)[j]-midt),ffit[j],format='(d20.8,d17.8,d22.8)'
  transit=0L
  for k=1,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[k]->get()
     tpar[9] = 1d0
     tpar[10] = 0d0
     tpar[11] = 0d0
     
     dum = TAP_MPFit_function(tpar,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,redfit=redfit,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int,adv_opt=transit.adv_opt)
     printf,strlun,'# Phased from '+transit.fname+' Parameters'
     printf,strlun,'#      T_days         Phase [hours from Tmid]     Model Flux'
     
     midt = transit.params[where(strcmp(transit.params.param,'Mid Transit') eq 1)].value - (transit.epoch*transit.params[0].value)

     for j=0,n_elements(ffit)-1,1 do printf,strlun,(*transit.model_tfine)[j],((*transit.model_tfine)[j]-midt),ffit[j],format='(d20.8,d17.8,d22.8)'
     
     
  endfor
  
  
  close,strlun
  midt=0L
  dum = 0L
  ffit=0L
  dfit=0L
  transit = 0L
  tpar = 0
  

 ; stop

  if self.create_mcmcascii then begin
     self.message = '   ...Writing MCMC chains to ascii files'
     self->message
     


     for i=0,n_elements(*self.transits)-1,1 do begin
        transit = (*self.transits)[i]->get()
        openw,strlun,self.plot_dir+'ascii_'+transit.fname+'_'+string(i+1,format='(i2.2)')+'_MCMC.txt',/get_lun,width=2500,bufsize=0
        
        if strcmp(self.parameterize_id,'basic') then begin
           head = '## TAP '+self.version+' Analysis of '+string(transit.mcmc_params[0],format='(i2.2)')+' MCMC Chains of '+transit.fname+' '+string(i+1,format='(i2.2)')
           printf,strlun,head
           ;; head = '##   '+string(self.num_combined,format='(i)')+' MCMC chains combined.'
           ;; printf,strlun,head
           head = '##   Burn-in stripped, every 10th link saved, leaving '+ $
                  string(n_elements(*transit.params[0].refined_mcmc_chain),format='(i)')+ $
                  ' links in this file, representative of '$
                  +string(n_elements(*transit.params[0].refined_mcmc_chain)*10d0,format='(i)')+ ' total MCMC links.'
           printf,strlun,head
           
           head = '##'
           format='('
           for ii=0,n_elements(transit.params.param)-1 do begin
              if ii ne 0 then format+=','
              if strcmp(transit.params[ii].param,'Mid Transit') then head += string(transit.params[ii].param,format='(a30)') else $
                 if ii eq 0 then head+= string(transit.params[ii].param,format='(a18)') else $
                    head+= string(transit.params[ii].param,format='(a20)') 
              if strcmp(transit.params[ii].param,'Mid Transit') then format += 'd30.15' else format+='d20.15'
           endfor
           
           head+=string('Likelihood',format='(a30)')
           
           format+=',d30.10)'
           
           printf,strlun,head      

           for j=0d0,n_elements(*transit.params[0].refined_mcmc_chain)-1d0,1d0 do begin
              printf,strlun,(*transit.params[0].refined_mcmc_chain)[j], $
                     (*transit.params[1].refined_mcmc_chain)[j], $          
                     (*transit.params[2].refined_mcmc_chain)[j], $
                     (*transit.params[3].refined_mcmc_chain)[j], $          
                     (*transit.params[4].refined_mcmc_chain)[j], $
                     (*transit.params[5].refined_mcmc_chain)[j], $          
                     (*transit.params[6].refined_mcmc_chain)[j], $
                     (*transit.params[7].refined_mcmc_chain)[j], $          
                     (*transit.params[8].refined_mcmc_chain)[j], $
                     (*transit.params[9].refined_mcmc_chain)[j], $          
                     (*transit.params[10].refined_mcmc_chain)[j], $
                     (*transit.params[11].refined_mcmc_chain)[j], $          
                     (*transit.params[12].refined_mcmc_chain)[j], $         
                     (*transit.params[13].refined_mcmc_chain)[j], $      
                     (*transit.params[14].refined_mcmc_chain)[j], $
                     (*transit.refined_likelihood_chain)[j],$
                     format=format 
           endfor
           
           close,strlun
           transit = 0L
           
     endif else begin
        head = '## TAP '+self.version+' Analysis of '+string(transit.mcmc_params[0],format='(i2.2)')+' MCMC Chains of '+transit.fname+' '+string(i+1,format='(i2.2)')
        printf,strlun,head
        head = '##   Burn-in stripped, every 10th link saved, leaving '+ $
               string(n_elements(*transit.params[0].refined_mcmc_chain),format='(i)')+ $
               ' links in this file, representative of '$
               +string(n_elements(*transit.params[0].refined_mcmc_chain)*10d0,format='(i)')+ ' total MCMC links.'
        printf,strlun,head
        
        head = '##'
        format='('
        for ii=0,n_elements(transit.params.param)-1 do begin
           if ii ne 0 then format+=','
           if strcmp(transit.params[ii].param,'Mid Transit') then head += string(transit.params[ii].param,format='(a30)') else $
              if ii eq 0 then head+= string(transit.params[ii].param,format='(a18)') else $
                 head+= string(transit.params[ii].param,format='(a20)') 
           if strcmp(transit.params[ii].param,'Mid Transit') then format += 'd30.15' else format+='d20.15'
        endfor
   
        head+=string('Inclination[inf]','a/R*[inf]',format='(a20,a20)')
        head+=string('Likelihood',format='(a30)')
        format+=',d20.15,d20.15,d30.10)'
        printf,strlun,head
        
        for j=0d0,n_elements(*transit.params[0].refined_mcmc_chain)-1d0,1d0 do begin
           printf,strlun,(*transit.params[0].refined_mcmc_chain)[j], $
                  (*transit.params[1].refined_mcmc_chain)[j], $          
                  (*transit.params[2].refined_mcmc_chain)[j], $
                  (*transit.params[3].refined_mcmc_chain)[j], $          
                  (*transit.params[4].refined_mcmc_chain)[j], $
                  (*transit.params[5].refined_mcmc_chain)[j], $          
                  (*transit.params[6].refined_mcmc_chain)[j], $
                  (*transit.params[7].refined_mcmc_chain)[j], $          
                  (*transit.params[8].refined_mcmc_chain)[j], $
                  (*transit.params[9].refined_mcmc_chain)[j], $          
                  (*transit.params[10].refined_mcmc_chain)[j], $
                  (*transit.params[11].refined_mcmc_chain)[j], $          
                  (*transit.params[12].refined_mcmc_chain)[j], $         
                  (*transit.params[13].refined_mcmc_chain)[j], $
                  (*transit.params[14].refined_mcmc_chain)[j], $
                  (*transit.basic_params[where(strcmp(transit.basic_params.param,'Inclination'))].refined_mcmc_chain)[j],$
                  (*transit.basic_params[where(strcmp(transit.basic_params.param,'a/R*'))].refined_mcmc_chain)[j],$
                  (*transit.refined_likelihood_chain)[j],$
                  format=format 
        endfor
        
        close,strlun
        transit = 0L
        
        
        
     endelse
  endfor
     
  endif
  
  
end


pro TAP::LoadTransit,event
  widget_control,self.transitfile_fld[1],get_value=path
  path = path[0]

 ; stop
  case self.filetype of
     'IDL Save File': BEGIN
        self.num_transits+=1
        self.active_transit+=1
        
        restore,path
        self->preptransit,input=lc,type='struct'
     end
     'ASCII File': Begin
        tap_readcol,path,temp1,temp2,format='d,d'
        tap_readcol,path,flag,epoch,sets,param,value,lock,llim,ulim,limited,pripen,ppval,ppsig,$
                    format='(a,i,i,a,d,i,d,d,d,d,d,d)'
        tap_readcol,path,flag,tnum,lint1,lint2,lint3,format='(a,i,i,d,i)'
        
        bounds = where(temp1 eq -1 and temp2 eq -1)
        if bounds[0] eq -1 then begin
           ntransit = 1
           bounds=[-1,n_elements(temp1)]
        endif else begin
           ntransit = n_elements(where(temp1 eq -1 and temp2 eq -1))+1
           bounds=[-1,bounds,n_elements(temp1)]
        endelse
        self.message = 'Detected '+string(ntransit,format='(i2.2)')+' light curves.'
        self->message
        
        for i=0,ntransit-1,1 do begin   
           self.num_transits+=1
           self.active_transit+=1
           
           t1 = temp1[bounds[i]+1:bounds[i+1]-1]
           t2 = temp2[bounds[i]+1:bounds[i+1]-1]
     
           struct= replicate({hjd:0d0, f: 0d0, e:0d0},n_elements(t1))
           struct.hjd=t1[sort(t1)]
           struct.f = t2[sort(t1)]

           setpars = where(tnum eq self.active_transit)
           if setpars[0] ne -1 then begin
              int_set = [lint1[setpars],lint2[setpars],lint3[setpars]]           
           endif else int_set = [self.setup_rebin,self.smooth_val]

         ;  stop
           setpars = where(epoch eq self.active_transit)
           if setpars[0] ne -1 then begin
              set = replicate({epoch: 0d0,$
                               set: 0,$
                               param: '',$
                               value: 0d0,$
                               lock: 0L,$
                               llim: 0d0,$
                               ulim:0d0,$
                               limited: 0d0,$
                               pripen: 0,$
                               ppval: 0d0,$
                               ppsig: 0d0},n_elements(setpars))
              set.epoch = epoch[setpars]
              set.set = sets[setpars]
              set.param = param[setpars]
              set.value = value[setpars]
              set.lock = lock[setpars]
              set.llim = llim[setpars]
              set.ulim = ulim[setpars]  
              set.limited = limited[setpars]   
              set.pripen = pripen[setpars]   
              set.ppval = ppval[setpars]   
              set.ppsig = ppsig[setpars]  
           endif else set = {param: ''}
           
           
           self->preptransit,input=struct,type='text',set=set,int_set=int_set
           int_set=0L
           set = 0L
           line=0L
           struct = 0L
           wc = 0L
        endfor
        temp1=0L
        temp2=0L
     end
  endcase
  self->setup,'multi'
  self->setup,'cross lc links'
  self->setup,'fit'

  self.message = 'light curve(s) loaded.'
  self->message

  ;; self->setup,'inference'
end

pro TAP::OrganizeTransits,input,_extra=_extra
  if self.num_transits eq 1 then *self.transits = obj_new('transit') else *self.transits = [*self.transits,obj_new('transit')]
  
  transit = ptr_new(self->prepTransitStruct(input,_extra=_extra))
  (*transit).rebin = self.setup_rebin
  (*transit).t_int = self.smooth_val
  (*transit).epoch = 0

  (*self.transits)[self.active_transit-1]->set,*transit

  input = 0L
  
  ptr_free,transit
  self->setupmodel
end

pro TAP::setupmodel
  transit = (*self.transits)[self.active_transit-1]->get()
 
  self->setparinfo
  temp = *self.parinfo

  (*self.parinfo)[0].fixed = 1
  (*self.parinfo)[10].fixed = 1
  (*self.parinfo)[11].fixed = 1
  (*self.parinfo)[5].fixed = 1
  (*self.parinfo)[6].fixed = 1
  (*self.parinfo)[8].fixed = 1
  (*self.parinfo)[9].fixed = 1
  
 ;; if (where(transit.params.fixed eq 2))[0] ne -1 then (*self.parinfo)[where(transit.params.fixed eq 2)].fixed = 1
  
 ;; x = mpfit('TAP_MPFit_function',$
 ;;           transit.params.value,$
 ;;           functargs={time:*transit.model_tfine ,$
 ;;                      flux:*transit.model_ffine ,$
 ;;                      tdat:transit.transit_tday, $   
 ;;                      fdat:transit.transit_flux, rebin:transit.rebin, smooth:transit.t_int $
 ;;                     },parinfo=(*self.parinfo),quiet=1)
  
  *self.parinfo = temp
 ;; temp = 0L
                                ;(*self.parinfo)[10].fixed = 0
                                ;(*self.parinfo)[11].fixed = 0
                                ;(*self.parinfo)[5].fixed = 0
                                ;(*self.parinfo)[6].fixed = 0
  
  if (where(transit.params.fixed eq 2))[0] ne -1 then transit.params[where(transit.params.fixed eq 2)].fixed = 0
 ;; transit.params.value = x
 ;; x = 0L
  (*self.transits)[self.active_transit-1]->set,transit
  
  transit = 0L
  self->updateModel
end


pro TAP::adjustslider,event
  widget_control, event.id, GET_UVALUE= uvalue
  transit = (*self.transits)[self.active_transit-1]->get()
  transit.params[where(strcmp(transit.params.param,uvalue.value) eq 1)].value = event.value
  (*self.transits)[self.active_transit-1]->set,transit
  transit = 0L
  self->updatemodel
end


function tap::disentangle_parameterization,transit
  case self.parameterize_id of
     'basic': begin
        transit.basic_params.value = transit.params.value
     end
     'adv1': begin
        print,"Not using this Parameterization!!!!"
        stop
     end
  'adv2': begin
     ;;   stop
        gamma1 = 1d0+transit.params[where(strcmp(transit.params.param,'Eccentricity'))].value*sin(transit.params[where(strcmp(transit.params.param,'Omega'))].value)
        gamma2 = sqrt(1-transit.params[where(strcmp(transit.params.param,'Eccentricity'))].value^2d0)
        ;; first just do the unchanged ones
        for i=0,n_elements(transit.basic_params)-1 do begin
           param = transit.basic_params[i].param
           if max(strcmp(param,['Period','Rp/R*','Mid Transit','Linear LD','Quadratic LD','Eccentricity',$
                                'Omega','OOT t^0','OOT t^1','OOT t^2','Sigma Red','Sigma White','Delta Light'])) then $
                                   transit.basic_params[i].value = transit.params[where(strcmp(transit.params.param,param))].value 

           ;;                 b         pi*T
           ;;  i = acos[ ----------- *  ---- ]
           ;;            sqrt(1-b^2)      P         
           if strcmp(param,'Inclination') then $
              transit.basic_params[i].value = (180/!dpi)*acos((transit.params[where(strcmp(transit.params.param,'b'))].value/$
                                                               sqrt(1-(transit.params[where(strcmp(transit.params.param,'b'))].value)^2)) * $
                                                              (!dpi*transit.params[where(strcmp(transit.params.param,'T'))].value)/$
                                                              transit.params[where(strcmp(transit.params.param,'Period'))].value)
           ;;   a        P       gamma2^2 
           ;;  --- =  ------- * ---------- * sqrt(1-b^2)
           ;;   R*     pi* T      gamma1
           if strcmp(param,'a/R*') then $
              transit.basic_params[i].value = ( (gamma2^2/gamma1) * $
                                                (transit.params[where(strcmp(transit.params.param,'Period'))].value / $
                                                 (!dpi*transit.params[where(strcmp(transit.params.param,'T'))].value)) * $
                                                sqrt(1-(transit.params[where(strcmp(transit.params.param,'b'))].value)^2))
        endfor
     end
endcase
  return,transit
end



pro TAP::updateModel
  transit = (*self.transits)[self.active_transit-1]->get()
  
  ;;  
  links = transit.params.set
  params= transit.params.value
  transit=0L
  
  for i=0,self.num_transits-1,1 do begin
     transit = (*self.transits)[i]->get()
     if (where(transit.params.set eq links))[0] ne -1 then begin
        transit.params[where(transit.params.set eq links)].value = params[where(transit.params.set eq links)]
        redfit = dblarr(n_elements(transit.rednoise))
        
     ;   stop
        transit = self->disentangle_parameterization(transit)
     ;   stop
        dum = TAP_MPFit_function(transit.basic_params.value,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,redfit=redfit,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int,adv_opt=transit.adv_opt)
        transit.rednoise = redfit
        *transit.model_f = dfit*1d0
        transit.residuals = transit.transit_flux - dfit
        *transit.model_ffine = ffit
        redfit=0L
        ffit = 0L
        dfit = 0L
        dum = 0L
        
        if abs((transit.params[where(strcmp(transit.params.param,'Mid Transit'))].value - median(transit.transit_tday))-transit.epoch*transit.params[where(strcmp(transit.params.param,'Period'))].value) gt .5*transit.params[where(strcmp(transit.params.param,'Period'))].value then begin
           epoch = fillarr(1,-5000,5000)
           transit.epoch = epoch[where(abs((transit.params[where(strcmp(transit.params.param,'Mid Transit'))].value - median(transit.transit_tday))-epoch*transit.params[where(strcmp(transit.params.param,'Period'))].value) eq min(abs((transit.params[where(strcmp(transit.params.param,'Mid Transit'))].value - median(transit.transit_tday))-epoch*transit.params[where(strcmp(transit.params.param,'Period'))].value)))]
           epoch = 0L
        endif
        
        (*self.transits)[i]->set,transit
     endif
     transit = 0L
  endfor
  links = 0L
  params = 0L
  self->lcplot
  self.numcol=2
  self->PrepCol
end

pro TAP::adjustxlock_event,event
  widget_control, event.id, GET_UVALUE= uvalue
  if n_elements(*event.value) eq 0 then (*event.value) = 0
  
  case uvalue.type of
     'val':begin
        transit = (*self.transits)[uvalue.value]->get()
        transit.params[uvalue.param].set = (*event.value)
        (*self.transits)[uvalue.value]->set,transit
        transit=0L
        self.active_transit = uvalue.value+1
        self->updatemodel
     end
     else: print,"unknown adjustxlock_event"
  endcase
end

pro TAP::SetParInfo
  ptr_free,self.parinfo

  self.parinfo = ptr_new(replicate({fixed: 0, $
                                    limited: [0,0], $
                                    limits: [0.,0.], $
                                    relstep: 0.d0},$
                                   n_elements(((*self.transits)[self.active_transit-1]->get()).params.fixed)))
  
  (*self.parinfo).fixed = ((*self.transits)[self.active_transit-1]->get()).params.fixed
  if (where((*self.parinfo).fixed eq 2))[0] ne -1 then (*self.parinfo)[where((*self.parinfo).fixed eq 2)].fixed = 1
  
  (*self.parinfo).limited = ((*self.transits)[self.active_transit-1]->get()).params.limited
  (*self.parinfo).limits = ((*self.transits)[self.active_transit-1]->get()).params.limits
  
  transit = (*self.transits)[self.active_transit-1]->get()
  if transit.rebin gt 0 then begin
     (*self.parinfo)[5].fixed = 1 
     (*self.parinfo)[6].fixed = 1 
     (*self.parinfo)[9].fixed = 1 
     (*self.parinfo)[10].fixed = 1 
  endif
  transit = 0L
  
end

function tap::window
  return,(*self.plot_windows)[0]
end

function TAP::PrepTransitStruct,input,set=set,int_set=int_set
  
  ;help,set
  ;help,set,/struct
 ; stop
  
  scale = (input.hjd)[1]-(input.hjd)[0]
                                ; t_final = ptr_new(makearr(((mm(double(input.hjd)))[1]-(mm(double(input.hjd)))[0])*(24d0*60d0), mm(double(input.hjd))+[-10d0*scale,10d0*scale]))
  t_final = ptr_new(fillarr((1d0/(24d0*60d0)),mm(double(input.hjd))+[-5d0*scale,5d0*scale]))
  
 ; stop
  self.setup_rebin = int_set[0]
  self.smooth_val = int_set[1:2]

  case self.setup_rebin of
     0: t = ptr_new(input.hjd)
     1: t = ptr_new((reform((dblarr(self.smooth_val[1])+1)#input.hjd $
                            +(findgen(self.smooth_val[1])+1-5d-1*(self.smooth_val[1]+1d0))*$
                            (self.smooth_val[0]/self.smooth_val[1])/1440d0 $
                            #(dblarr(n_elements(input.hjd))+1),self.smooth_val[1]*n_elements(input.hjd),1))[*,0])
  endcase  

       ;  parameterize_params: replicate({param: '', limits: [0d0,0d0], fixed: 0, limited: [0,0], $
           ;                                  value: 0d0, accept: 0.44d0, beta: 44d-2, set:self.active_transit, $
           ;                                  b_last: 0d0, sval: 2d0, prior: [0d0,0d0,0d0], $
           ;                                  mcmc_chain:ptr_new([-1]), jumpct: 0d0, jumptot: 0d0, $
           ;                                  curr_link: 0d0, new_link: 0d0, mcmc_val: dblarr(3), $
           ;                                  refined_mcmc_chain: ptr_new(), runval: dblarr(3),likelihoods:ptr_new([-1]),$
           ;                                  active: 0L},20) ,$
           ;;  parameterize_type: self.parameterization,$
  transit = {transit_tday: input.hjd  , $
             transit_flux: input.f    , $ 
             params: replicate({param: '', limits: [0d0,0d0], fixed: 0, limited: [0,0], $
                                value: 0d0, accept: 0.44d0, beta: 44d-2, set:self.active_transit, $
                                b_last: 0d0, sval: 2d0, prior: [0d0,0d0,0d0], $
                                mcmc_chain:ptr_new([-1]), jumpct: 0d0, jumptot: 0d0, $
                                curr_link: 0d0, new_link: 0d0, mcmc_val: dblarr(3), $
                                refined_mcmc_chain: ptr_new(), runval: dblarr(3),likelihoods:ptr_new([-1])},n_elements(*self.label_indices)-1) ,$
             basic_params:replicate({param: '', limits: [0d0,0d0], fixed: 0, limited: [0,0], $
                                     value: 0d0, accept: 0.44d0, beta: 44d-2, set:self.active_transit, $
                                     b_last: 0d0, sval: 2d0, prior: [0d0,0d0,0d0], $
                                    mcmc_chain:ptr_new([-1]), jumpct: 0d0, jumptot: 0d0, $
                                     curr_link: 0d0, new_link: 0d0, mcmc_val: dblarr(3), $
                                     refined_mcmc_chain: ptr_new(), runval: dblarr(3),likelihoods:ptr_new([-1])},n_elements(*self.basic_indices)-1) ,$
             likelihood_chain: ptr_new(),$
             refined_likelihood_chain: ptr_new(),$
             parameterize_id: '',$
             epoch: 0,$
             new_redl: 0d0,$
             curr_redl: 0d0,$
             mcmc_links: 0d0,$
             model_t: ptr_new(*t_final) ,$
             model_f: ptr_new(*t_final*0d0),$
             model_tfine: ptr_new(*t),$
             model_ffine: ptr_new((*t)*0d0),$
             rebin: 0L, $     
             ldtype: 0L,$
             t_int: dblarr(2),$
             residuals: dblarr(n_elements(input.hjd)) ,$
             rednoise:  dblarr(n_elements(input.hjd)) ,$
             color: 0L ,$
             modcol: 0L ,$
             redcol: 0L ,$
             fname: '',$
             mcmc_params: [10,1d5,1d5],$
             mcmc_complete: 0L ,$
             mcmc_files: ptr_new('-1') ,$
             mcmc_version: '' ,$
             adv_opt: lonarr(2), $
             diff: 0d0, $
             min_effective_length: dblarr(2), $
             expansion_ptr1: ptr_new(),$
             expansion_ptr2: ptr_new(),$
             expansion_ptr3: ptr_new(),$
             expansion_ptr4: ptr_new(),$
             expansion_ptr5: ptr_new()$
            } 
  
  
  widget_control,self.transitfile_fld[1],get_value=path
  path = path[0]
                                ; if self.active_transit+2 gt n_tags(*self.colors)
  transit.color = (*self.colors).((self.active_transit+2)-(n_tags(*self.colors)-2)*fix(self.active_transit/(n_tags(*self.colors)-2)))
                                ;print,transit.color
  transit.modcol = (*self.colors).blue
  transit.redcol = (*self.colors).red
  transit.ldtype = self.setup_ld
  if n_elements(strsplit(path,'/')) gt 1 then $
     transit.fname = (strsplit(((strsplit(path,'/',/extract))[n_elements(strsplit(path,'/',/extract))-1]),'.',/extract))[0] else $
        transit.fname = (strsplit(((strsplit(path,'\',/extract))[n_elements(strsplit(path,'\',/extract))-1]),'.',/extract))[0]
  transit.min_effective_length = self.eff_len
;;print,transit.fname
  ptr_free,transit.expansion_ptr1,transit.expansion_ptr2,transit.expansion_ptr3,transit.expansion_ptr4,transit.expansion_ptr5
  
  for i=1,n_elements(*self.label_indices)-1,1 do transit.params[i-1].param = strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)
  for i=1,n_elements(*self.basic_indices)-1,1 do transit.basic_params[i-1].param = strmid((*self.basic_indices)[i],0,strlen((*self.basic_indices)[i])-1)

  transit.parameterize_id = self.parameterize_id
 ;; transit.parameterize_num = self.parameterize_num

  input = 0L

;;  compare = ''
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].value = 
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].min =
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].max = 
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].locked = 
;;  transit.params[where(strcmp(transit.params.param,compare) eq 1)].limited = 

  

  comparr = transit.params.param
 ; print,comparr
 ; stop
  
  for i=0,n_elements(transit.params.param)-1,1 do comparr[i]=strjoin(strsplit(strlowcase(comparr[i]),' ',/extract),'_')
  for i=1,n_elements(*self.label_indices)-1,1 do begin
     compare = strjoin(strsplit(strlowcase(strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)),' ',/extract),'_')
     comp0 = where(strcmp(comparr,compare))
     comp = where(strcmp(STRLOWCASE(compare),STRLOWCASE(set.param)))
     comp2 = where(strcmp((*self.defaults).param,strlowcase(compare)))
     
 ;    print,compare
  ;   print,comp0,comp,comp2
     
     transit.params[comp0].value   =  (*self.defaults)[comp2].val
     transit.params[comp0].limits  = ((*self.defaults)[comp2].limit)[1:2]
     transit.params[comp0].fixed   =  (*self.defaults)[comp2].lock
     transit.params[comp0].limited = intarr(2)+((*self.defaults)[comp2].limit)[0]
     transit.params[comp0].prior   =  (*self.defaults)[comp2].prior

;     print,transit.params[comp0].value
  ;  stop

     if transit.params[comp0].value eq -99 and strcmp(comparr[i-1],'mid_transit') then begin
        transit.params[comp0].value  = transit.transit_tday[where(min(transit.transit_flux) eq transit.transit_flux)]
        transit.params[comp0].limits = transit.params[comp0].value+((*self.defaults)[comp2].limit)[1:2]
       ; stop
    ;   print, transit.params[comp0].value, transit.params[comp0].limits
    ;;  ;  stop
     endif
   
     if comp[0] ne -1 then begin
        transit.params[comp0].value = set[comp].value
        if set[comp].lock then begin
           transit.params[comp0].fixed = 1 
           transit.params[comp0].limits = [set[comp].value*.1,(set[comp].value+.01)*10]
        endif else begin
           transit.params[comp0].fixed = 2
           transit.params[comp0].limits = [set[comp].llim,set[comp].ulim]
        endelse
        if set[comp].limited then transit.params[comp0].limited = 1
        if set[comp].set then transit.params[comp0].set = set[comp].set
        transit.params[comp0].prior = [set[comp].pripen,set[comp].ppval,set[comp].ppsig]
     endif

  ;   stop
     
  endfor

  
                    ; print_struct,transit.params
  scale = 0L
  ptr_free,t_final,t
  return,transit
end


function tap::burnrangebasic,tnum=tnum
  range = findgen(n_elements(*((*self.transits)[0]->get()).params[0].mcmc_chain),2)
  range[*,1] = 1d0
  i0 = 0
  i1 = self.num_transits-1
  if keyword_set(tnum) then begin
     i0 = tnum
     i1 = tnum
  endif
  for i=i0,i1,1 do begin
     transit = (*self.transits)[i]->get()
     ;;; replaced basic_params with params.... check hsould be same
     for j=0,n_elements(transit.params)-1,1 do $
        if transit.params[j].fixed eq 0 then begin
        med = median(*transit.params[j].mcmc_chain)
        first = (*transit.params[j].mcmc_chain)[0]
        if first-med le 0 then cross = (where(*transit.params[j].mcmc_chain-med gt 0))[0] else $
           cross = (where(*transit.params[j].mcmc_chain-med le 0))[0]
        if cross[0] ne -1 then range[where(range lt cross),*] = 0d0
     endif 
  endfor
  
  ;;print,mm(range[*,0])
  ;;print,mm(range[where(range[*,1])])
  range = mm(range[where(range[*,1])])
  return,range
end

function tap::burnrange,tnum=tnum
  range = findgen(n_elements(*((*self.transits)[0]->get()).params[0].mcmc_chain),2)
  range[*,1] = 1d0
  i0 = 0
  i1 = self.num_transits-1
  if keyword_set(tnum) then begin
     i0 = tnum
     i1 = tnum
  endif
  for i=i0,i1,1 do begin
     transit = (*self.transits)[i]->get()
     for j=0,n_elements(transit.params)-1,1 do $
        if transit.params[j].fixed eq 0 then begin
        med = median(*transit.params[j].mcmc_chain)
        first = (*transit.params[j].mcmc_chain)[0]
        if first-med le 0 then cross = (where(*transit.params[j].mcmc_chain-med gt 0))[0] else $
           cross = (where(*transit.params[j].mcmc_chain-med le 0))[0]
        if cross[0] ne -1 then range[where(range lt cross),*] = 0d0
     endif 
  endfor
  
  ;;print,mm(range[*,0])
  ;;print,mm(range[where(range[*,1])])
  range = mm(range[where(range[*,1])])
  return,range
end



pro TAP::Message,event,message=message
  if keyword_set(message) then self.message = message
  widget_control,self.message_window,/append,$
                 set_value='('+curr_date(format='hh:mm:ss')+') '+self.message
                                ;print,self.message
  self.message=''
end

pro TAP::PrepTransit,event,input=input,type=type,_extra=_extra
  case type of
     'struct' :BEGIN
        input.hjd = input[sort(input.hjd)].hjd
        input.f = input[sort(input.hjd)].f
        input.e = input[sort(input.hjd)].e
        self->organizetransits,input,_extra=_extra
        self.message='Imported existing Transit Structure'
        self->message
     END
     'text' : BEGIN
        self->organizetransits,input,_extra=_extra
 ;       self.message='Imported ASCII Transit File'
     END
  endcase
  


  self->lcplot
end


;function TAP::PrepModel,event
;  setp=ptr_new([-1, 0.d0])
;  
;  return,model
;end

pro TAP::MenuMap,event
  Widget_Control, event.id, Get_UValue=info
  
  case info.type of
     'main': begin
        for i=1,n_elements(*self.bases)-1,1 do widget_control,(*self.bases)[i],MAP=0
        case info.value of
          ; 'General': begin
          ;    widget_control,self.general_base,/map
          ; end
           'Load Transit': begin
             ; self->setup,'load'
              widget_control,   (*self.bases)[1],/map
           end
           'Manage Transits': begin
             ; self->setup,'multi'
              widget_control,   (*self.bases)[2],/map
           end
           'Fit': begin
             ; self->setup,'fit'
              widget_control,   (*self.bases)[3],/map
           end
           'MCMC Inference': begin
             ; self->setup,'inference'
              widget_control,   (*self.bases)[4],/map
           end
           'Parameterization': begin
              widget_control,   (*self.bases)[5],/map
           end
           else:  print,'Unknown MenuMap Event "'+info.value+'"'
        endcase
     end
  endcase
end

pro TAP::start
  centertlb,self.tap_base
  widget_control,self.tap_base,/realize
  widget_control,self.tap_base, XOffset=0, YOffset=0

  widget_control, (*self.plot_windows)[0].window, Get_Value=wid
  (*self.plot_windows)[0].w_id=wid
  
  window, xsize=(*self.plot_windows)[0].x $
          , ysize=(*self.plot_windows)[0].y $
          , /pixmap,/free
  (*self.plot_windows)[0].pix_window = !d.window
  self->lcplot
end

pro TAP::destroy
  self.message = 'Destroying Widget'
  self->message

  for i=0,n_elements(self.buttons)-1,1 do if (self.buttons)[i] then widget_control,(self.buttons)[i],/destroy
  for i=0,n_elements(*self.slider)-1,1 do  if (*self.slider)[i].id then widget_control,(*self.slider)[i].id,/destroy
  for i=0,n_elements(self.settings)-1,1 do if self.settings[i] then widget_control,self.settings[i],/destroy
  for i=n_elements(*self.bases)-1,0,-1 do if (*self.bases)[i] then widget_control,(*self.bases)[i],/destroy
  if self.tap_base then widget_control,self.tap_base,/destroy
  ptr_free,self.colors
  ptr_free,self.slider
  ptr_free,self.plot_windows
  ptr_free,self.parinfo
  
  ptr_free,self.mcmc_stuff

  for i=0,n_elements(*self.extra_windows)-1,1 do if (*self.extra_windows)[i] then widget_control,(*self.extra_windows)[i],map=0
  for i=0,n_elements(self.fld[*,0])-1,1 do $
     for j=0,n_elements(self.fld[0,*])-1 do if self.fld[i,j] ne 0 then begin
     widget_control,self.fld[i,j],/destroy
     self.fld[i,j] = 0
  endif
  for i=0,n_elements(*self.extra_windows)-1,1 do if (*self.extra_windows)[i] then widget_control,(*self.extra_windows)[i],/destroy


  for i=0,self.num_transits-1,1 do begin
     transit = (*self.transits)[i]->get()   
     for k=0,n_elements(transit.params)-1,1 do begin
        ptr_free,transit.params[k].mcmc_chain
        ptr_free,transit.params[k].refined_mcmc_chain
     endfor
     ptr_free,transit.model_t,transit.model_f,transit.model_tfine,transit.model_ffine,transit.mcmc_files
     transit = 0L
     (*self.transits)[i]->destroy
  endfor
  
  ptr_free,self.transits
 
  ptr_free,self.extra_windows
  ptr_free,self.label_indices
  ptr_free,self.label
  ptr_free,self.menus
  ptr_free,self.bases

  obj_destroy,self
end


pro TAP_Cleanup,event

end


pro tap::PrepCol,runval=runval,es=es
  cols = self.numcol

  if 1-keyword_set(es) then et = strarr(n_elements(*self.label)) else begin
     et = strarr(n_elements(*self.label))
     et[0] = es
  endelse
  
  if keyword_set(runval) then begin
     cols = 4
     col2 = ['Value',string(((*self.transits)[runval-1]->get()).params.runval[0],format='(d17.8)')]
     col3 = ['85.1%',string(((*self.transits)[runval-1]->get()).params.runval[1],format='(d9.5)')]
     col4 = ['15.9%',string(((*self.transits)[runval-1]->get()).params.runval[2],format='(d9.5)')] 
     
     badformat = where(strcmp(strmid(col3,0,1),'*'))
     if badformat[0] ne -1 then col3[badformat] = string(((*self.transits)[runval-1]->get()).params[badformat].runval[1],format='(e9.2)')
     badformat = 0L
     badformat = where(strcmp(strmid(col4,0,1),'*'))
     if badformat[0] ne -1 then col4[badformat] = string(((*self.transits)[runval-1]->get()).params[badformat].runval[2],format='(e9.2)')
     badformat = 0L
     
  endif else begin
     if cols gt 1 then begin
        col2 = ['Value',string(((*self.transits)[self.active_transit-1]->get()).params.value,format='(d17.8)')]
        if cols gt 2 then begin
           col3 = ['85.1%',string(((*self.transits)[self.active_transit-1]->get()).params.mcmc_val[1],format='(d9.5)')]
           col4 = ['15.9%',string(((*self.transits)[self.active_transit-1]->get()).params.mcmc_val[2],format='(d9.5)')]

           badformat = where(strcmp(strmid(col3,0,1),'*'))
           if badformat[0] ne -1 then col3[badformat] = string(((*self.transits)[self.active_transit-1]->get()).params[badformat].mcmc_val[1],format='(e9.2)')
           badformat = 0L
           badformat = where(strcmp(strmid(col4,0,1),'*'))
           if badformat[0] ne -1 then col4[badformat] = string(((*self.transits)[self.active_transit-1]->get()).params[badformat].mcmc_val[2],format='(e9.2)')
           badformat = 0L
        endif
     endif
  endelse
  
  skip = -1
  if self.adv_opt[1] eq 0 then begin
     skip = [skip,where(strcmp('Delta Light:',*self.label_indices))]
     widget_control,(*self.label)[where(strcmp('Delta Light:',*self.label_indices))],set_value=''
  endif
  
  line = ''
  case cols of 
     1: begin
        for i=0,n_elements(*self.label_indices)-1,1 do begin
           if (where(i eq skip))[0] eq -1 then begin
              string = string((*self.label_indices)[i],format='(A14)') 
              widget_control,(*self.label)[i],set_value=string
              string = 0L
           endif
        endfor
        ;widget_control,(*self.label)[i],set_value="note: [locked] {penalized} (gaussian prior)"  
     end
     2: begin
        for i=0,n_elements(*self.label_indices)-1,1 do begin
           if (where(i eq skip))[0] eq -1 then begin
              string = string((*self.label_indices)[i],col2[i],format='(A15,A17)')
              if i gt 0 then begin
                 if ((*self.transits)[self.active_transit-1]->get()).params[where(strcmp(((*self.transits)[self.active_transit-1]->get()).params.param,strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)))].prior[0] eq 1 then string = string('{'+(*self.label_indices)[i]+'}',col2[i],format='(A15,A17)')
                 if ((*self.transits)[self.active_transit-1]->get()).params[where(strcmp(((*self.transits)[self.active_transit-1]->get()).params.param,strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)))].prior[0] eq 2 then string = string('('+(*self.label_indices)[i]+')',col2[i],format='(A15,A17)')
                 if ((*self.transits)[self.active_transit-1]->get()).params[where(strcmp(((*self.transits)[self.active_transit-1]->get()).params.param,strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)))].fixed then string = string('['+(*self.label_indices)[i]+']',col2[i],format='(A15,A17)')
              endif
              widget_control,(*self.label)[i],set_value=string
              string = 0L
           endif
        endfor
        widget_control,(*self.label)[i],set_value="note: [locked] {penalized} (gaussian prior)" 
     end
     4: begin
        if (where(1d0*col3[1:n_elements(col3)-1] eq -1))[0] ne -1 then col3[where(1d0*col3[1:n_elements(col3)-1] eq -1)+1] = string('Fixed',format='(A8)')
        if (where(1d0*col4[1:n_elements(col4)-1] eq -1))[0] ne -1 then col4[where(1d0*col4[1:n_elements(col4)-1] eq -1)+1] = string('Fixed',format='(A8)')
        for i=0,n_elements(*self.label_indices)-1,1 do begin
           if (where(i eq skip))[0] eq -1 then begin
              string = string(et[i]+(*self.label_indices)[i],col2[i],col3[i],col4[i],format='(A15,A17,A10,A10)')
              if i gt 0 then begin
                 if ((*self.transits)[self.active_transit-1]->get()).params[where(strcmp(((*self.transits)[self.active_transit-1]->get()).params.param,strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)))].prior[0] eq 1 then string =string(et[i]+'{'+(*self.label_indices)[i]+'}',col2[i],col3[i],col4[i],format='(A15,A17,A10,A10)')
                 if ((*self.transits)[self.active_transit-1]->get()).params[where(strcmp(((*self.transits)[self.active_transit-1]->get()).params.param,strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)))].prior[0] eq 2 then string = string(et[i]+'('+(*self.label_indices)[i]+')',col2[i],col3[i],col4[i],format='(A15,A17,A10,A10)')
                 if ((*self.transits)[self.active_transit-1]->get()).params[where(strcmp(((*self.transits)[self.active_transit-1]->get()).params.param,strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)))].fixed then string = string(et[i]+'['+(*self.label_indices)[i]+']',col2[i],col3[i],col4[i],format='(A15,A17,A10,A10)')
              endif
              widget_control,(*self.label)[i],set_value=string
              string=0L
           endif
        endfor
        widget_control,(*self.label)[i],set_value="   *NOTE: [locked] {penalized} (gaussian prior)" 
     end
     else: print,'DEBUG: PrepCol unknown # of cols:',cols
     
  endcase
  
  cols=0L
  col1=0L
  col2=0L
  col3=0L
  col4=0L
  
end

pro tap::setup,set_this
  case set_this of 
 'adjustadvanced':begin
        (*self.extra_windows)[5] = widget_base(title='Advanced Parameters',/column,frame=3,uname='adjustadvanced')

        XManager, 'TAP' $
                  , (*self.extra_windows)[5] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[5] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:5, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[5],frame=0,/column,/base_align_center)
       
        row = widget_base(work_base,frame=1,/col,/base_align_center)
        ;;  label = widget_label(row_one,value='Delta Kp Correction Term')
        radio = widget_base(row,column=1,/nonexclusive)
        button = widget_button(radio,value='Activate 3rd Light Correction Term',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'advopts',which:1})
        widget_control,button,set_button=self.adv_opt[1]
        
        self.settings[4] = widget_base(row,frame=1,/column,sensitive=self.adv_opt[1])
        row2 = widget_base(self.settings[4],frame=0,/row,/base_align_center)
        
        button = widget_button(row2,$
                               font = self.buttonfont,$
                               value = 'Activate All',$
                               uvalue={object:self, wid:5, method:'adv_activate_all'})
        
        button = widget_button(row2,$
                               font = self.buttonfont,$
                               value = 'One Set',$
                               uvalue={object:self, wid:5, method:'adv_oneset_all'})
        
      
        for i=0,n_elements(*self.transits)-1,1 do begin
           transit = (*self.transits)[i]->get()
           row2 = widget_base(self.settings[4],frame=0,/row,/base_align_center)
           label=widget_label(row2,FONT=self.buttonfont,Value=transit.fname+' Delta:',xsize=150,/align_center)
           radio = widget_base(row2,column=1,/nonexclusive)
           self.fld2[i,0] = widget_button(radio,value='Active?',$ 
                                  uvalue={object:self,method:'adjustdelkp_button',value:i,$
                                          type: 'val', param:0 })
           widget_control,self.fld2[i,0],set_button=transit.adv_opt[1]
           
           
                                ;   self.fld2[i,0] =  coyote_field2(row2, LABELFONT=self.buttonfont,FIELDFONT=self.textfont,$
                                ;                                   TITLE='active [0 or 1]:', /integervalue, $
                                ;                                   UVALUE={object:self, method:'adjustdelkp_event',value:i,$
                                ;                                           type: 'val', param:0 },$
                                ;                                   VALUE=transit.adv_opt[1], /positive, XSIZE=2,scr_ysize=30,$
                                ;                                   event_pro='TAP_event',textid=textid)
           ;;; VALUE of delta Kp
           self.fld2[i,1] =  coyote_field2(row2, LABELFONT=self.buttonfont,FIELDFONT=self.textfont,$
                                           TITLE='Value:', /doublevalue, $
                                           UVALUE={object:self, method:'adjustdelkp_event',value:i,$
                                                   type: 'val', param:1 },$
                                           decimal=4,/cr_only,$
                                           VALUE=transit.params[14].prior[1], XSIZE=7,scr_ysize=30,$
                                           event_pro='TAP_event',textid=textid)
           ;; SIGMA of delta Kp
           self.fld2[i,2] =  coyote_field2(row2, LABELFONT=self.buttonfont,FIELDFONT=self.textfont,$
                                           TITLE='sigma:', /doublevalue, $
                                           UVALUE={object:self, method:'adjustdelkp_event',value:i,$
                                                   type: 'val', param:2 },$
                                           decimal=4,/cr_only,/positive,$
                                           VALUE=transit.params[14].prior[2], XSIZE=7,scr_ysize=30,$
                                           event_pro='TAP_event',textid=textid)

           ;;; SET value for delta kp
           self.fld2[i,3] =  coyote_field2(row2, LABELFONT=self.buttonfont,FIELDFONT=self.textfont,$
                                           TITLE='set:', /integervalue,/positive, $
                                           UVALUE={object:self, method:'adjustdelkp_event',value:i,$
                                                   type: 'val', param:3 },$
                                           VALUE=round(transit.params[14].set), XSIZE=7,scr_ysize=30,$
                                           event_pro='TAP_event',textid=textid)
           radio = widget_base(row2,column=1,/nonexclusive)
           self.fld2[i,4] = widget_button(radio,value='Lock',$ 
                                          uvalue={object:self,method:'adjustdelkp_button',value:i,$
                                                  type: 'val', param:4 })
           widget_control,self.fld2[i,4],set_button=transit.params[14].fixed           
           transit=0L
        endfor
        
     end
     'parameterize': begin
        (*self.bases)[5] = widget_base((*self.bases)[0] ,$
                                       /row,$
                                       MAP=0)
        
        col = widget_base((*self.bases)[5],/column,frame=1)
     ;;   lbl = widget_label(col,value='Select Parameterization')
        
        
        bg = cw_bgroup(col,$
                       *self.parameterization_options,$
                       /column,$
                       LABEL_top='System Parameterization:',$
                       /RETURN_NAME,$
                       /NO_RELEASE,$
                       set_value =self.parameterization_num,$
                       UVALUE={object:self, method:'ButtonEvent', value:'ParamType'},$      
                       FONT=self.buttonfont,$
                       /EXCLUSIVE)


        col = widget_base((*self.bases)[5],/column,frame=1)
     ;;   lbl = widget_label(col,value='Select Parameterization')
        
        bg = cw_bgroup(col,$
                       *self.ld_parameterization_options,$
                       /column,$
                       LABEL_top='Limb-Darkening Param:',$
                       /RETURN_NAME,$
                       /NO_RELEASE,$
                       set_value =self.setup_LD,$
                       UVALUE={object:self, method:'ButtonEvent', value:'LDType'},$      
                       FONT=self.buttonfont,$
                       /EXCLUSIVE)
                                ;self.setup_rebin = 0
                                ;   options = 0L
        
        col = widget_base((*self.bases)[5],/column,frame=0)
        lbl = widget_label(col,value='Notes')
        
        message= 'Parametrization Options:'+self.cr+$
                 '[Note that TAP will calculate inferred paramters i and a/rs if a parameterization other than "Basic" is selected.'+self.cr+self.cr+$
                 '"Basic": General system parameters.'+self.cr+self.cr+$
                 '"T, b": Transit duration and impact parameter replace a/R* and inclination in this parameterization which decreases correlation between b and a/R*'+self.cr+self.cr+$
                 '"Limb Darkening Selection:"'+self.cr+$
                 'u1 and u2 draws directly from these (correlated) parameters, while "2u1+u2, u1-2u2" draws from those orthogonalized distributions (Holman 2006).  If both parameters are locked, this is disregarded.  If one parameter is locked, TAP defaults to draw the other from either u1 or u2.  Regardless of the case, TAP outputs u1 and u2 values so no transformations are needed.'
        
        
        text = widget_text(col,font=self.textfont,value=message,/scroll,/wrap,xsize=50,ysize=11)

     end
     'load': begin
        ys = 80

        (*self.bases)[1] = widget_base((*self.bases)[0] ,$
                                        /column,$
                                        MAP=0)

        bigrow =  widget_base( (*self.bases)[1] ,$
                           /column,$
                           /Base_align_center,frame=1)
        row = widget_base( bigrow,$
                           /ROW,$
                           /Base_align_center,frame=0)
        
        button = widget_button(row,$
                               font=self.buttonfont,$
                               xsize=100,$
                               value='Transit File',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Transit File Button'})
    
        path = coyote_field2(row,$
                             labelfont=self.buttonfont,$
                             FIELDFONT=self.textfont,$
                             TITLE=':',$
                             VALUE = self.transitpath,$
                             UVALUE='Transit File Field',$
                             XSIZE=65,$
                             TEXTID=textid)
        self.transitfile_fld = [path,textid]
        
        clear = widget_button(row,$
                              FONT=self.buttonfont,$
                              VALUE='Clear',$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear Data Path'})
        
        mrow = widget_base(bigrow,$
                           /ROW,$
                           /Base_align_center)
        

        if 0 then begin
           row = widget_base(mrow,frame=1,$
                             /column,ysize=ys,$
                             /Base_align_left)
           bg = cw_bgroup(row,$
                          ['ASCII File','IDL Save File'],$
                          /column,$
                          LABEL_top='File Type:',$
                          /RETURN_NAME,$
                          /NO_RELEASE,$
                          UVALUE={object:self, method:'ButtonEvent', value:'FileType'},$      
                          FONT=self.buttonfont,$
                          /EXCLUSIVE,$
                          SET_VALUE=0)   
        endif
        self.filetype = 'ASCII File'
        

        row = widget_base(mrow,frame=1,$
                          /column,ysize=ys,$
                          /Base_align_left)
        
                                ; row = widget_base(mrow,$
                                ;                   /column,frame=1,$
                                ;                   /Base_align_center)
        
                                ; radio = widget_base(row,column=1,/nonexclusive)
                                ;  label = widget_label(row,value='Long integrations?')
      ;  options = ['None','Rebin to data cadence','Rebin to "Input Integration"']
        options = ['None','Rebin to "Input Integration"']
        bg = cw_bgroup(row,$
                       options,$
                       /column,$
                       LABEL_top='Long Integrations?',$
                       /RETURN_NAME,$
                       /NO_RELEASE,$
                       set_value =self.setup_rebin,$
                       UVALUE={object:self, method:'ButtonEvent', value:'RebinType'},$      
                       FONT=self.buttonfont,$
                       /EXCLUSIVE)
                                ;self.setup_rebin = 0
        options = 0L
        
        if self.setup_rebin eq 1 then x = 1 else x = 0
        self.setup_smooth = widget_base(mrow,frame=1,$
                                        /column,$
                                        /Base_align_left, sensitive=x,ysize=ys)
        
        lbl = widget_label(self.setup_smooth,value='Input Integration')
        
        row = widget_base(self.setup_smooth,/row)
        fld = coyote_field2(row,$
                            TITLE='Minutes:',$
                            UVALUE={object:self, method:'ButtonEvent', $
                                    value:'settint',set:0},$
                            decimal=4, digits=8, $
                            VALUE=self.smooth_val[0],$
                            XSIZE=10,/doubleValue,event_pro='TAP_event',$
                            textid=textid)    
        
     ;   lbl = widget_label(self.setup_smooth,value='Samples per Integration')
        
        fld = coyote_field2(row,$
                            TITLE='N_samp:',$
                            UVALUE={object:self, method:'ButtonEvent', $
                                    value:'settint',set:1},$
                            decimal=4, digits=8, $
                            VALUE=round(self.smooth_val[1]),$
                            XSIZE=5,/integerValue,event_pro='TAP_event',$
                            textid=textid)    
        
        
        row = widget_base(mrow,frame=0,$
                          /column,$
                          /Base_align_left)
        
                                ;button = widget_button(radio,value='Long integrations',$
                                ;uname = ,$ 
         ;                      uvalue={object:self,$
          ;                             method:'setup_rebin',$
           ;                            value:0 ,$
            ;                           type: ''$
             ;                         })
        ;widget_control,button,set_button=self.setup_rebin
        

       ; row = widget_base(mrow,$
       ;                   /ROW,$
       ;                   /Base_align_center,ysize=40)
        
        self.buttons[2] =  widget_button(row,$
                                font=self.buttonfont,$
                                xsize=100,$
                                value='Load Transit',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Load Transit Button'},sensitive=0)
       ; widget_control,   (*self.bases)[1],/map
              

        row = widget_base( (*self.bases)[1] ,$
                           /ROW,$
                           /Base_align_center,frame=1)
        
        button = widget_button(row,$
                               font=self.buttonfont,$
                               xsize=100,$
                               value='Setup File',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Saved Setup File Button'})
    
        path = coyote_field2(row,$
                             TITLE=':',$
                             VALUE = self.transitpath1,$
                             UVALUE='Transit Setup Field',$
                             XSIZE=65,$
                             TEXTID=textid)
        self.transitfile1_fld = [path,textid]
        
        clear = widget_button(row,$
                              FONT=self.buttonfont,$
                              VALUE='Clear',$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear Setup Path'})
        
        self.buttons[1] =  widget_button(row,$
                                font=self.buttonfont,$
                                xsize=110,$
                                value='Load Setup',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Load Setup Button'},sensitive=0)
      ;  widget_control,   (*self.bases)[1],/map

        if 0 then begin
        row = widget_base( (*self.bases)[1] ,$
                           /ROW,$
                           /Base_align_center,frame=1)
        button = widget_button(row,$
                               font=self.buttonfont,$
                               xsize=100,$
                               value='Save Setup',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'Save Current Setup Button'})
        
        path = coyote_field2(row,$
                             TITLE=':',$
                             VALUE = self.transitpath2,$
                             UVALUE='Transit Save Setup Field',$
                             XSIZE=65,$
                             TEXTID=textid)
        self.transitfile2_fld = [path,textid]
        
        clear = widget_button(row,$
                              FONT=self.buttonfont,$
                              VALUE='Clear',$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear Data Path'})
        
        button =  widget_button(row,$
                                font=self.buttonfont,$
                                xsize=110,$
                                value='Save Setup',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Save Setup'})
     endif              
     end
     'multi': begin
        tframe = 0
        if self.num_transits ne 0 then sens=1 else sens = 0
        (*self.bases)[2] = widget_base( (*self.bases)[0],$
                                        /column,$
                                        MAP=0,sensitive=sens,frame=tframe)
        
        base = widget_base((*self.bases)[2],/row,frame=tframe)
        col = widget_base(base,/column,frame=1)
        lbl = widget_label(col,value='General')

        if self.num_transits lt 2 then sens1 = 0 else sens1 = 1
        slider = widget_slider(col,minimum=1,maximum=max([2,self.num_transits]), $
                               value=max([1,self.active_transit]), $
                               uvalue={object:self, method:'ButtonEvent', $
                                       value:'ActiveTransit'},/drag,title='Active Transit',sensitive=sens1)

        slider = cw_fslider(col,title='Plot scaling',min=0,max=.05,value=self.diff,$
                            uvalue={object:self, method:'buttonevent',value:'diffset'},drag=1,$
                            /double, /edit,xsize=150)
        
        col=widget_base(base,/column,frame=1,/base_align_center)
        
        self.settings[0] = widget_label(col,value='settings:',/dynamic_resize)
        
        
        if self.num_transits gt 0 then widget_control,self.settings[0],set_value=string(self.active_transit,format='(i2.2)')+": "+((*self.transits)[self.active_transit-1]->get()).fname+' settings:'
        
        row = widget_base(col,/row)
        col2 = widget_base(row,frame=1,/col)
      options = ['None','Rebin to "Input Integration"']
        self.settings[1] = cw_bgroup(col2,$
                                     options,$
                                     /column,$
                                     LABEL_top='Long Integrations?',$
                                     /RETURN_NAME,$
                                     /NO_RELEASE,$
                                     set_value =0,$
                                     UVALUE={object:self, method:'ButtonEvent', value:'RebinTypeA'},$      
                                     FONT=self.buttonfont,$
                                     /EXCLUSIVE)
                                ;self.setup_rebin = 0
        options = 0L
        if self.num_transits gt 0 then widget_control,self.settings[1],set_value=((*self.transits)[self.active_transit-1]->get()).rebin        
        col2 = widget_base(row,frame=1,/col)
        
        lbl = widget_label(col2,value='Input Integration')
        
      ;  row = widget_base(,/row)
        self.settings[2] = coyote_field2(col2,$
                            LABELFONT=self.buttonfont,$
                            FIELDFONT=self.textfont,$
                            TITLE='Minutes:',$
                            UVALUE={object:self, method:'ButtonEvent', $
                                    value:'settint2',set:0},$
                                         decimal=4, digits=8, $
                                         VALUE=self.smooth_val[0],$
                                         XSIZE=10,/doubleValue,event_pro='TAP_event',$
                            textid=textid)    
        
        
     ;   lbl = widget_label(self.setup_smooth,value='Samples per Integration')
        
        self.settings[3] = coyote_field2(col2,$
                            LABELFONT=self.buttonfont,$
                            FIELDFONT=self.textfont,$
                            TITLE='N_samp:',$
                            UVALUE={object:self, method:'ButtonEvent', $
                                    value:'settint2',set:1},$
                            decimal=4, digits=8, $
                            VALUE=round(self.smooth_val[1]),$
                            XSIZE=5,/integerValue,event_pro='TAP_event',$
                            textid=textid)    
        
        
        button = widget_button(col,xsize=150,value='Delete Active Transit',$
                               uvalue={object:self,$
                                       method:'ButtonEvent',$
                                       value:'Delete Active'})

        
        col=widget_base(base,/column,frame=1)    
        

        

lbl = widget_label(col,value='Inter-Transit Settings:')

        button =  widget_button(col,$
                                font=self.buttonfont,$
                                xsize=60,$
                                value='Set Links',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Setup Cross LC Locks'})    
             
        
     end
     'oldmulti': begin
        if self.num_transits ne 0 then sens=1 else sens = 0
        
   ;     if (*self.bases)[2] ne 0 then widget_control,(*self.bases)[2],/destroy
                                ;    if (*self.bases)[2] eq 0 then begin
        (*self.bases)[2] = widget_base( (*self.bases)[0],$
                                        /column,$
                                        MAP=0,sensitive=sens)
        
        row = widget_base( (*self.bases)[2] ,$
                           /column)
        
        button =  widget_button(row,$
                                font=self.buttonfont,$
                                xsize=200,$
                                value='Setup Cross LC Locks',$
                                uvalue={object:self, $
                                        method:'ButtonEvent', $
                                        value:'Setup Cross LC Locks'})
        if self.num_transits gt 0 then begin
           transits='01: '+((*self.transits)[0]->get()).fname
           for i=1,n_elements(*self.transits)-1,1 do transits = [transits,string(i+1,format='(i2.2)')+': '+((*self.transits)[i]->get()).fname]
           bg = cw_bgroup(row,$
                          transits,$
                          ROW=n_elements(transits)/4d0,$
                          LABEL_LEF='Set Active Transit:',$
                          /RETURN_NAME,$
                          /NO_RELEASE,$
                          UVALUE={object:self, method:'ButtonEvent', value:'ActiveTransit'},$      
                          FONT=self.buttonfont,$
                          /EXCLUSIVE,/scroll,$
                          x_scroll_size=500, y_scroll_size=50,$ ;,ysize=20,$
                          
                          SET_VALUE=self.active_transit-1)
           transits = 0L
        endif
                                ; endif
     end
     'fit': begin
        if self.num_transits ne 0 then sens=1 else sens = 0
     ;   if (*self.bases)[3] ne 0 then widget_control,(*self.bases)[3],/destroy
        (*self.bases)[3] = widget_base( (*self.bases)[0],$
                                        /column,$
                                        MAP=0,sensitive=sens)
        
        
        row1_base = widget_base(  (*self.bases)[3] ,$
                                  /ROW)
        
        col1_base = widget_base(row1_base,$
                                /COLUMN,$
                                /BASE_ALIGN_LEFT,$
                                FRAME=1)
        
        label = widget_label(col1_base,$
                             value='Manual Parameter Adjustments and Setup',$
                             font=self.buttonfont,$
                             /align_left)
        
;;; COL 1: manual adjust!        
        
        
        button = widget_button(col1_base,$
                               font  = self.buttonfont,$
                               value = 'Manual Parameter Adjustment',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value: 'Manual Parameter Adjustment'},$
                               /align_center)
        
        
        
        button = widget_button(col1_base,$
                               font=self.buttonfont,$
                               value='Adjust Limits and Locks',$
                               uvalue= {object:self,$
                                        method: 'ButtonEvent',$
                                        value: 'Adjust Limits and Locks'},$
                               /align_center)
        
        
        
        
        col2 = widget_base(row1_base,$
                          /column,$
                          /BASE_ALIGN_left,$
                          FRAME=2)
        

        label = widget_label(col2,$
                             value='Advanced Parameters',$
                             font=self.buttonfont,$
                             /align_center)
        button = widget_button(col2,$
                               font=self.buttonfont,$
                               value='Advanced Parameters',$
                               uvalue= {object:self,$
                                        method: 'ButtonEvent',$
                                        value: 'Adjust Advanced Parameters'},$
                               /align_center)
        
        

        
        
;;; COL1 END
        
                                ; col2_base = widget_base(row1_base,$
                                ;                         /COLUMN,$
                                ;                         /BASE_ALIGN_center,$
                                ;                         FRAME=2)
                                ; 
                                ; label = widget_label(col2_base,$
                                ;                      value='Levenberg-Markwart',$
                                ;                      font=self.buttonfont,$
                                ;                      /align_left)
  ;;; COL 2: automatic LM fit...   
        
        
  col3_base = widget_base(row1_base,$
                          /column,$
                          /BASE_ALIGN_left,$
                          FRAME=2)
  
  label = widget_label(col3_base,$
                       value='Markov Chain Monte Carlo',$
                       font=self.buttonfont,$
                       /align_center)
  ;;; col 3 MCMC
        

  button = widget_button(col3_base,$
                         value = 'MCMC Parameters',$
                         uvalue={object:self, $
                                 method:'ButtonEvent', $
                                 value: 'MCMC Parameters'},$
                         /align_center,$
                         xsize=110.,$
                                ; ysize=30.,$
                         sensitive=1)
  
  button = widget_button(col3_base,$
                         value = 'Gaussian Priors',$
                         uvalue={object:self, $
                                 method:'ButtonEvent', $
                                 value: 'Gaussian Priors'},$
                         /align_center,$
                         xsize=110.,$
                                ; ysize=30.,$
                         sensitive=1) 

  button = widget_button(col3_base,$
                         value = 'Execute Chain',$
                         uvalue={object:self, $
                                 method:'ButtonEvent', $
                                 value: 'Execute MCMC'},$
                         /align_center,$
                         xsize=110.,$
                                ; ysize=30.,$
                         sensitive=1)
  
        
        
     end
     'inference': begin
        (*self.bases)[4] = widget_base( (*self.bases)[0],$
                                        /column,$
                                        MAP=0)
      ;  print,4,(*self.bases)[4]
        row1_base = widget_base(  (*self.bases)[4] ,$
                                  /column)
        
        
        col1_base = widget_base(row1_base,$
                                /column,$
                                /BASE_ALIGN_LEFT,$
                                FRAME=1)
        label = widget_label(col1_base,value='1) Select options for output:'$
                             ,font=self.buttonfont, /align_left)

        row1 = widget_base(col1_base,/row,frame=0)
        col2 = widget_base(row1,/column,frame=0,/base_align_left,xsize=200)
        row = widget_base(col2,/row,frame=0,/base_align_center,ysize=30)
        radio = widget_base(row,column=1,/nonexclusive)
        button = widget_button(radio,value='Create .eps Plots',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'create_plots'})
        widget_control,button,set_button=self.create_plots
        

        row = widget_base(col2,/row,frame=0,/base_align_center,ysize=30)
        radio = widget_base(row,column=1,/nonexclusive)
        button = widget_button(radio,value='Plot 2D distributions (slow)',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'create_2d'})
        widget_control,button,set_button=self.create_2d
        
        

        col2 = widget_base(row1,/column,frame=0,/base_align_left,xsize=200)
        row = widget_base(col2,/row,frame=0,/base_align_center,ysize=30)
        radio = widget_base(row,column=1,/nonexclusive)
        button = widget_button(radio,value='Calculate G-R Stat',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'calc_gr'})
        widget_control,button,set_button=self.calc_gr
        
        row = widget_base(col2,/row,frame=0,/base_align_center,ysize=30)
        radio = widget_base(row,column=1,/nonexclusive)
        button = widget_button(radio,value='Calculate Eff Lengths',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'calc_effl'})
        widget_control,button,set_button=self.calc_effl
        

        col2 = widget_base(row1,/column,frame=0,/base_align_left,xsize=200)
        row = widget_base(col2,/row,frame=0,/base_align_center,ysize=30)
        radio = widget_base(row,column=1,/nonexclusive)
        button = widget_button(radio,value='MCMC chains -> ascii (big)',$ 
                               uvalue={object:self, method:'buttonevent',$
                                       value: 'create_mcmcascii'})
        widget_control,button,set_button=self.create_mcmcascii


        col1_base = widget_base(row1_base,$
                                /COLUMN,$
                                /BASE_ALIGN_LEFT,$
                                FRAME=1)
        label = widget_label(col1_base,value='2) Load a compatible TAP_setup.idlsav file'$
                             ,font=self.buttonfont, /align_left)
        row = widget_base(col1_base,/row,frame=0,/base_align_center)
        
        button = widget_button(row,$
                               font=self.buttonfont,$
                               xsize=150,$
                               value='MCMC Save File',$
                               uvalue={object:self, $
                                       method:'ButtonEvent', $
                                       value:'MCMC File Button'})
        
        if self.mcmc_fld[1] ne 0 then widget_control,self.mcmc_fld[1],get_value=pval else pval= ''
        path = coyote_field2(row,$
                             labelfont=self.buttonfont,$
                             FIELDFONT=self.textfont,$
                             TITLE=':',$
                             VALUE = pval,$
                             UVALUE='MCMC File Field',$
                             XSIZE=70,$
                             TEXTID=textid)
        self.mcmc_fld = [path,textid]
       
        pval = 0L
        
        self.buttons[0] =  widget_button(row, font=self.buttonfont, xsize=40, value='Load',$
                                         sensitive=0, uvalue={object:self, method:'ButtonEvent', $
                                                              index: 0, value:'Load MCMC Button'})
        
        clear = widget_button(row,FONT=self.buttonfont,VALUE='Clear', xsize=40,$
                              uvalue={object:self, $
                                      method:'ButtonEvent', $
                                      value:'Clear MCMC Path'})
        
        
        
        
        
        
     end
     'cross lc links': begin
        (*self.extra_windows)[0] = widget_base(title='MCMC Multi Chain Parameter Sets',/column,frame=3,uname='links')
        XManager, 'TAP' $
                  , (*self.extra_windows)[0] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[0] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:0, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[0],frame=1,/row)
        for i=0,self.num_transits-1,1 do begin
           column = widget_base(work_base,frame=1,/column)
           values = round(((*self.transits)[i]->get()).params.set)
         ;  help,values
           label=widget_label(column,font=self.buttonfont $
                              , Value=((*self.transits)[i]->get()).fname,/align_center)   
           for j=0,n_elements(((*self.transits)[i]->get()).params)-1,1 do begin
              row = widget_base(column,frame=0,/row,/align_center)
              if i eq 0 then begin
                 button = widget_button(row,value='Lock All',$
                                       uvalue={object:self,  method:'ButtonEvent', value:'xlock_lockall',param:j})
                 button = widget_button(row,value='Free',$
                                       uvalue={object:self,  method:'ButtonEvent', value:'xlock_freeall',param:j})
                 label=widget_label(row,FONT=self.buttonfont,Value=(*self.label_indices)[j+1],xsize=80,/align_center)
                 
              endif
              self.fld[i,j] = coyote_field2(row,$
                                            LABELFONT=self.buttonfont,$
                                            FIELDFONT=self.textfont,$
                                            TITLE='',$
                                            /integervalue, $
                                            UVALUE={object:self, method:'adjustxlock_event',value:i,$
                                                    type: 'val', param:j $
                                                   },$
                                            VALUE=values[j],$
                                            /positive,$
                                            XSIZE=2,$
                                            scr_ysize=30,$
                                            event_pro='TAP_event',$
                                            textid=textid)
           endfor
           values = 0L
        endfor
                     
     end
     'gaussianpriors': begin
        (*self.extra_windows)[4] = widget_base(title='MCMC Gaussian Priors',/column,frame=3,uname='mcmcgaussianpriors')
        XManager, 'TAP' $    
                  , (*self.extra_windows)[4] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[4] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:4, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[4],frame=0,/column,/base_align_center)

        transit = (*self.transits)[self.active_transit-1]->get()
        
        for i=0,n_elements(transit.params)-1,1 do begin
           row = widget_base(work_base,frame=0,/row)
           label = widget_label(row,font=self.buttonfont,value=transit.params[i].param,xsize=80,/align_right)
           fld = coyote_field2(row,$
                               LABELFONT=self.buttonfont,$
                               FIELDFONT=self.textfont,$
                               TITLE='Value:',$
                               UVALUE={object:self, method:'adjustLL_event',value:i,$
                                       type: 'prior_val'$
                                      },$
                               VALUE=((transit.params[i]).prior)[1],$
                               XSIZE=15,$
                               /doubleValue,event_pro='TAP_event',$
                               textid=textid)
           fld = coyote_field2(row,$
                               LABELFONT=self.buttonfont,$
                               FIELDFONT=self.textfont,$
                               TITLE='Value:',$
                               UVALUE={object:self, method:'adjustLL_event',value:i,$
                                       type: 'prior_sig'$
                                      },$
                               VALUE=((transit.params[i]).prior)[2],$
                               XSIZE=15,$
                               /doubleValue,event_pro='TAP_event',$
                               textid=textid)

          ;; options = ['Disabled','Penalty','Prior']
           bg = cw_bgroup(row,['Disabled','Penalty','Prior'],$
                         /return_name,/no_release,label_top='',/exclusive,$
                          set_value=(transit.params[i].prior)[0],/row,$
                          uvalue={object:self,method:'adjustLL_event',value:i,type:'prior_penalize'})
           
;           radio = widget_base(row,column=1,/nonexclusive)
                                ;          button = widget_button(radio,value='Enable',$
                                ;                              ;uname = ,$ 
   ;                               uvalue={object:self,$
    ;                                      method:'adjustLL_event',$
     ;                                     value: i,$
      ;                                    type: 'prior_penalize'$
       ;                                  })
        ;   widget_control,button,set_button= (transit.params[i].prior)[0]
          
        endfor
       ; transit = 0L
     end
     'mcmcparams': begin
        (*self.extra_windows)[3] = widget_base(title='MCMC Chain Parameters',/column,frame=3,uname='mcmcparams')
        XManager, 'TAP' $    
                  , (*self.extra_windows)[3] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[3] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:3, method:'QuitWidget'})

        base_work = widget_base((*self.extra_windows)[3],frame=1,/row,/base_align_left)
        work_base = widget_base(base_work,frame=0,/column,/base_align_left)
        work_base2 = widget_base(base_work,frame=0,/column,/base_align_right)

        base      = widget_base(work_base,frame=1,/column)
        label     = widget_label(base,FONT=self.buttonfont,Value='Number of Chains:',/align_left)
        
        fld = coyote_field2(base,$
                            LABELFONT=self.buttonfont,$
                            FIELDFONT=self.textfont,$
                            TITLE=' ',$
                            UVALUE={object:self, method:'adjustLL_event',value:'numchain',$
                                    type: 'mcmcsetup'$
                                   },$
                            VALUE=(((*self.transits)[self.active_transit-1])->get()).mcmc_params[0],$
                            XSIZE=10,$
                            
                                ;  format='(G10.2)',$
                            /doubleValue, $
                            ;/integervalue,$
                            decimal = 0,$
                            event_pro='TAP_event',$
                            textid=textid)
        
        
        ;; chain length
         ;;;  'chainmin':  self.mcmc_chainlength[0] = *event.value
         ;;;  'chainmax':  self.mcmc_chainlength[1] = *event.value
        
       ; base      = widget_base(work_base,frame=1,/column)
  label     = widget_label(base,FONT=self.buttonfont,Value='Chain Length ("links"):',/align_left)
  
  fld = coyote_field2(base,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Min. Total:',$
                      UVALUE={object:self, method:'adjustLL_event',value:'chainmin',$
                              type: 'mcmcsetup'$
                             },$
                      VALUE=(((*self.transits)[self.active_transit-1])->get()).mcmc_params[1],$
                      XSIZE=10,$
                      decimal=0,$
                                ;  format='(G10.2)',$
                      /doubleValue,event_pro='TAP_event',$
                      textid=textid)
  
  
  fld = coyote_field2(base,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Min. Effective:',$
                      UVALUE={object:self, method:'adjustLL_event',value:'effmin',$
                              type: 'mcmcsetup'$
                             },$
                      VALUE=(((*self.transits)[self.active_transit-1])->get()).min_effective_length[1],$
                      XSIZE=6,$
                      decimal=0,$
                                ;  format='(G10.2)',$
                      /doubleValue,event_pro='TAP_event',$
                      textid=textid)

  fld = coyote_field2(base,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Max. Total:',$
                      UVALUE={object:self, method:'adjustLL_event',value:'chainmax',$
                              type: 'mcmcsetup'$
                             },$
                      VALUE=(((*self.transits)[self.active_transit-1])->get()).mcmc_params[2],$
                      XSIZE=10,$
                      decimal=0,$
                                ;  format='(G10.2)',$
                      /doubleValue,event_pro='TAP_event',$
                      textid=textid)

  if 0 then begin
  fld = coyote_field2(base,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Maximum Links:',$
                      UVALUE={object:self, method:'adjustLL_event',value:'chainmax',$
                              type: 'mcmcsetup'$
                             },$
                      VALUE=(((*self.transits)[self.active_transit-1])->get()).mcmc_params[2],$
                      XSIZE=10,$
                      decimal=0,$
                                ;  format='(G10.2)',$
                      /doubleValue,event_pro='TAP_event',$
                      textid=textid)

endif

 ; base      = widget_base(work_base,frame=1,/column)
;  label     = widget_label(base,FONT=self.buttonfont,Value='Limb Darkening Selection:',/align_left)
;  options = ['u1 and u2','2u1+u2, u1-2u2']
;  bg = cw_bgroup(base,$
 ;                options,$
 ;                /column,$
 ;                LABEL_top='Draw from:',$
 ;                /RETURN_NAME,$
 ;                /NO_RELEASE,$
 ;                set_value =self.setup_LD,$
 ;                UVALUE={object:self, method:'ButtonEvent', value:'LDType'},$      
 ;                FONT=self.buttonfont,$
 ;                /EXCLUSIVE)
  
  
  message= '"Number of Chains:"'+self.cr+$
           'TAP will calculate this number of independent chains.'+self.cr+self.cr+$
           '"Minimum Total Length:"'+self.cr+$
           'The minimum number of links per chain which TAP will execute.'+self.cr+self.cr+$
           '"Minimum Effective Lengths:"'+self.cr+$
           'After reaching "Min Total Links", TAP calculates correlation length and effective chain length for each free parameter (Tegmark 2004), and extends the chains by 100,000 links at a time until all free parameters reach the minimum effective length.  For adequate mixing, that length should be >> 1, and TAP defaults to 100 per chain, or 1000 for a 10 chain run.';+self.cr+self.cr+$
  ;         '"Limb Darkening Selection:"'+self.cr+$
  ;         'u1 and u2 draws directly from these (correlated) parameters, while "2u1+u2, u1-2u2" draws from those orthogonalized distributions (Holman 2006).  If both parameters are locked, this is disregarded.  If one parameter is locked, TAP defaults to draw the other from either u1 or u2.  Regardless of the case, TAP outputs u1 and u2 values so no transformations are needed.'
  
  text = widget_text(work_base2,font=self.textfont,value=message,/scroll,/wrap,xsize=40,ysize=20)


     end
     'adjustparams': begin
        (*self.extra_windows)[1] = widget_base(title='Manual Parameter Adjustment',/column,frame=3,uname='links')
        XManager, 'TAP' $
                  , (*self.extra_windows)[1] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[1] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:1, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[1],frame=0,/column,/base_align_center)
        
        message = 'Manual Parameter Adjustment:'+self.cr+$
                  'Adjust parameters for the active LC (and any LCs linked to it) by the sliders or text entry boxes.'+self.cr+$
                  '**NOTE*** If using the text entry boxes, press enter to update the values.  Make sure that the value you have entered shows up in the main TAP widget Parameter section.'+self.cr+self.cr+$
                  'To adjust beyond the limits of the sliders, you must change the "Min" and "Max" limits in the "Adjust Limits and Locks" widget, accessable from the main TAP widget.'

        text = widget_text(work_base,font=self.textfont,value=message,/scroll,/wrap,ysize=4,xsize=70)
        
        label = widget_label(work_base,value='LC: '+((*self.transits)[self.active_transit-1]->get()).fname)



        row = widget_base(work_base,/column,/align_center,frame=1)
        
        label = widget_label(row,$
                             value='System Parameters',$
                             font=self.buttonfont,$
                             /align_left)
        
        col_params = widget_base(row,$
                                 frame=0,$
                                 column=3,$
                                 /align_left)
        for i=0,4,1 do begin
        ;   help,/heap
           (*(self.slider))[i].id =$
              CW_Fslider(col_params, $
                         title=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         min=((*self.transits)[self.active_transit-1]->get()).params[i].limits[0], $
                         max=((*self.transits)[self.active_transit-1]->get()).params[i].limits[1],$
                         value = ((*self.transits)[self.active_transit-1]->get()).params[i].value,$
                         format='(G20.15)',/double, /edit,$
                         uname=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         uvalue={object:self, method:'AdjustSlider', value:((*self.transits)[self.active_transit-1]->get()).params[i].param},$
                         drag=1)
         ;  help,/heap
         ;  stop
        endfor
        for i=7,8,1 do begin
           (*(self.slider))[i].id =$
              CW_Fslider(col_params, $
                         title=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         min=((*self.transits)[self.active_transit-1]->get()).params[i].limits[0], $
                         max=((*self.transits)[self.active_transit-1]->get()).params[i].limits[1],$
                         value = ((*self.transits)[self.active_transit-1]->get()).params[i].value,$
                         format='(G20.15)',/double, /edit,$
                         uname=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         uvalue={object:self, method:'AdjustSlider', value:((*self.transits)[self.active_transit-1]->get()).params[i].param},$
                         drag=1)
           endfor
        row = widget_base(work_base,/column,/align_left,frame=1)
        
        label = widget_label(row,$
                             value='Quadratic Limb Darkening',$
                             font=self.buttonfont,$
                             /align_left)
        
        col_params = widget_base(row,$
                                 frame=0,$
                                 column=3,$
                                 /align_left)
        
           for i=5,6,1 do begin
           (*(self.slider))[i].id =$
              CW_Fslider(col_params, $
                         title=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         min=((*self.transits)[self.active_transit-1]->get()).params[i].limits[0], $
                         max=((*self.transits)[self.active_transit-1]->get()).params[i].limits[1],$
                         value = ((*self.transits)[self.active_transit-1]->get()).params[i].value,$
                         format='(G20.15)',/double, /edit,$
                         uname=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         uvalue={object:self, method:'AdjustSlider', value:((*self.transits)[self.active_transit-1]->get()).params[i].param},$
                         drag=1)
           endfor
           row = widget_base(work_base,/column,/align_left,frame=1)
        
        label = widget_label(row,$
                             value='Data Specific Parameters',$
                             font=self.buttonfont,$
                             /align_left)
        
        col_params = widget_base(row,$
                                 frame=0,$
                                 row=2,$
                                 /align_left)
        
           for i=9,13,1 do begin
           (*(self.slider))[i].id =$
              CW_Fslider(col_params, $
                         title=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         min=((*self.transits)[self.active_transit-1]->get()).params[i].limits[0], $
                         max=((*self.transits)[self.active_transit-1]->get()).params[i].limits[1],$
                         value = ((*self.transits)[self.active_transit-1]->get()).params[i].value,$
                         format='(G20.15)',/double, /edit,$
                         uname=((*self.transits)[self.active_transit-1]->get()).params[i].param,$
                         uvalue={object:self, method:'AdjustSlider', value:((*self.transits)[self.active_transit-1]->get()).params[i].param},$
                         drag=1)
        endfor
        end
     'adjustlimlocks' : begin
        ;;   'Adjust Limits and Locks': begin
        ;;      self->setup,'adjustlimlocks'
        ;;      centertlb,(*self.extra_windows)[2]
        ;;      widget_control,(*self.extra_windows)[2],/realize
        ;;   end
        (*self.extra_windows)[2] = widget_base(title='MCMC Limits and Locks' , $
                                               /column,frame=3,uname='links')
        XManager, 'TAP' $
                  , (*self.extra_windows)[2] $
                  , /no_block       
        quit_button = widget_button((*self.extra_windows)[2] ,$
                                    font = self.buttonfont,$
                                    value = 'Quit',$
                                    uvalue={object:self, wid:2, method:'QuitWidget'})
        work_base = widget_base((*self.extra_windows)[2],frame=0,/column)
        label = widget_label(work_base,value='LC: ' + $
                             ((*self.transits)[self.active_transit-1]->get()).fname,/align_center)
        col = widget_base(work_base,/column,/align_center,frame=0)

        transit = (*self.transits)[self.active_transit-1]->get()
        
        label = widget_label(col,value='System Parameters',font=self.buttonfont,/align_left)
        col_params = widget_base(col,frame=1,/column,/align_left)
        for i=0,4,1 do self->llrow,transit,i,col_params
        for i=7,8,1 do self->llrow,transit,i,col_params
        label = widget_label(col,value='Quadratic Limb Darkening',font=self.buttonfont,/align_left)
        col_params = widget_base(col,frame=1,/column,/align_left)
        for i=5,6,1 do   self->llrow,transit,i,col_params,/disablelims
        label = widget_label(col,value='Data Specific Parameters',font=self.buttonfont,/align_left)
        col_params = widget_base(col,frame=1,/column,/align_left)
        for i=9,13,1 do self->llrow,transit,i,col_params
                           
        transit = 0L
        
        
     end
     else: print,'unknown SETUP: '+set_this
  endcase
end


pro TAP::AdjustLL_event,event
  widget_control, event.id, GET_UVALUE= uvalue
 ; print,*event.value
 ; print,uvalue.type
  transit = (*self.transits)[self.active_transit-1]->get()
  case uvalue.type of
     'min': begin
        if n_elements(*event.value) gt 0 then begin
           transit.params[uvalue.value].limits[0] = *event.value
           for i=0,self.num_transits-1,1 do begin
              if i ne self.active_transit-1 then begin
                 trans = (*self.transits)[i]->get()
                 if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                    trans.params[uvalue.value].limits[0] = *event.value
                    (*self.transits)[i]->set,trans
                 endif
                 trans = 0L
              endif
           endfor
        endif
     end
     'max': begin
        if n_elements(*event.value) gt 0 then begin
           transit.params[uvalue.value].limits[1] = *event.value
           for i=0,self.num_transits-1,1 do begin
              if i ne self.active_transit-1 then begin
                 trans = (*self.transits)[i]->get()
                 if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                    trans.params[uvalue.value].limits[1] = *event.value
                    (*self.transits)[i]->set,trans
                 endif
                 trans = 0L
              endif
           endfor
        endif
     end
     'lock': begin
        if transit.params[uvalue.value].fixed eq 0 then $
           transit.params[uvalue.value].fixed = 1 else $
              transit.params[uvalue.value].fixed = 0
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].fixed = transit.params[uvalue.value].fixed
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'limit': begin
        if transit.params[uvalue.value].limited[0] eq 0 then $
           transit.params[uvalue.value].limited = [1,1] else $
              transit.params[uvalue.value].limited = [0,0]
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].limited = transit.params[uvalue.value].limited
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'mcmcacceptrate': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        transit.params[uvalue.value].accept = (*event.value)
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].accept = transit.params[uvalue.value].accept
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'mcmcsetup': begin
        transit = 0L
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        case uvalue.value of
           'numchain':  begin
              for i=0,n_elements(*self.transits)-1,1 do begin
                 tran = (*self.transits)[i]->get()
                 tran.mcmc_params[0] = *event.value
                 (*self.transits)[i]->set,tran
                 if self.active_transit-1 eq i then transit = tran
                 tran = 0L
              endfor
           end
           'effmin': begin
              self.eff_len[1] = *event.value
              for i=0,n_elements(*self.transits)-1,1 do begin
                 tran = (*self.transits)[i]->get()
                 tran.min_effective_length[1] = *event.value
                 (*self.transits)[i]->set,tran
                 if self.active_transit-1 eq i then transit = tran
                 tran = 0L
              endfor
           end
           'chainmin':   begin
              for i=0,n_elements(*self.transits)-1,1 do begin
                 tran = (*self.transits)[i]->get()
                 tran.mcmc_params[1] = *event.value
                 (*self.transits)[i]->set,tran
                 if self.active_transit-1 eq i then transit = tran
                 tran = 0L
              endfor
           end
           'chainmax':    begin
              for i=0,n_elements(*self.transits)-1,1 do begin
                 tran = (*self.transits)[i]->get()
                 tran.mcmc_params[2] = *event.value
                 (*self.transits)[i]->set,tran
                 if self.active_transit-1 eq i then transit = tran
                 tran = 0L
              endfor
           end
        endcase
     end
     'prior_val': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        transit.params[uvalue.value].prior[1] = (*event.value)
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].prior = transit.params[uvalue.value].prior
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'prior_sig': begin
        if n_elements(*event.value) eq 0 then (*event.value) = 0
        transit.params[uvalue.value].prior[2] = (*event.value)
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].prior = transit.params[uvalue.value].prior
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     'prior_penalize': begin
        transit.params[uvalue.value].prior[0] = where(strcmp(['Disabled','Penalty','Prior'],event.value) eq 1)
        for i=0,self.num_transits-1,1 do begin
           if i ne self.active_transit-1 then begin
              trans = (*self.transits)[i]->get()
              if trans.params[uvalue.value].set eq transit.params[uvalue.value].set then begin
                 trans.params[uvalue.value].prior = transit.params[uvalue.value].prior
                 (*self.transits)[i]->set,trans
              endif
              trans = 0L
           endif
        endfor
     end
     else: begin
        print,'UNKNOWN adjust LL event'+uvalue.type
        stop
     end
  endcase
  (*self.transits)[self.active_transit-1]->set,transit
  transit=0L
                                ; endif
end


pro TAP::LLrow,transit,i,col_params,disablelims=disablelims
  row = widget_base(col_params,frame=0,/row,/base_align_left)
  label=widget_label(row,FONT=self.buttonfont, $
                     Value=transit.params[i].param+':', $
                     xsize=90,/align_right)
  fld = coyote_field2(row,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Min',$
                      UVALUE={object:self, method:'AdjustLL_event', $
                              value:i, type: 'min'},$
                      decimal=10, digits=20, $
                      VALUE=transit.params[i].limits[0],$
                      XSIZE=15,/doubleValue,event_pro='TAP_event',$
                      textid=textid)    
  fld = coyote_field2(row,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='Max',$
                      UVALUE={object:self, method:'AdjustLL_event', $
                              value:i, type: 'max'},$
                      decimal=10, digits=20, $
                      VALUE=transit.params[i].limits[1],$
                      XSIZE=15,/doubleValue,event_pro='TAP_event',$
                      textid=textid)    
  radio = widget_base(row,column=1,/nonexclusive)
  button = widget_button(radio,value='Lock',$ 
                         uvalue={object:self, method:'AdjustLL_event',$
                                 value: i, type: 'lock'})
  widget_control,button,set_button=transit.params[i].fixed
  fld = coyote_field2(row,$
                      LABELFONT=self.buttonfont,$
                      FIELDFONT=self.textfont,$
                      TITLE='MCMC Accept Rate:',$
                      UVALUE={object:self, method:'AdjustLL_event', $
                              value:i, type: 'mcmcacceptrate' },$
                      VALUE=transit.params[i].accept,$
                      XSIZE=5,$
                      decimal=2,$
                      digits=3,$
                                ;  format='(G10.2)',$
                      /doubleValue,event_pro='TAP_event',$
                      textid=textid)
  radio = widget_base(row,column=1,/nonexclusive)
  
  if keyword_set(disablelims) then sens = 0 else sens = 1
  button = widget_button(radio,value='Apply Limits to Fitting.',$
                         uvalue={object:self, method:'AdjustLL_event',$
                                 value: i, type: 'limit'},sensitive=sens)
  widget_control,button,set_button=transit.params[i].limited[0]
end

pro tap::widget_setup
  self.tap_base = widget_base(title='Transit Analysis Package '+self.version,/column)
  
  XManager, 'TAP' $
            , self.tap_base $
            , /no_block $
            , cleanup = 'TAP_cleanup'
  

  quit_button = widget_button(self.tap_base,$
                              font = self.buttonfont,$
                              value = 'Quit',$
                              uvalue={object:self, method:'Quit'})
 
  message = '('+ curr_date(format='hh:mm:ss yyyymmdd') +') TAP Tools '+self.version
  self.message_window = widget_text(self.tap_base, $
                                    font = self.textfont, $
                                    value = message, $
                                    /scroll, $
                                    ysize=4)
  
  main_base = widget_base(self.tap_base,$
                          /row,$
                          frame=5, event_pro='TAP_AllEvents')
  
  col1_base = widget_base(main_base,$
                          /COLUMN,$
                          /BASE_ALIGN_RIGHT,$
                          /FRAME)

  ysize=18.
  xsize=320.
  for i=0,n_elements(*self.label_indices),1 do begin
     (*self.label)[i] = widget_label(col1_base,$
                                     font=self.textfont, $
                                     value = ' ', $
                                     /align_left,$
                                     ysize = ysize,$
                                     xsize = xsize)
  endfor
  
 ; print,i

  self.numcol = 1
  self->PrepCol

  col2_base = widget_base(main_base,$
                          /column,$
                          /base_align_left,$
                          /frame)

  (*self.plot_windows)[0].x = 400d0
  (*self.plot_windows)[0].y = 380d0
  
  (*self.plot_windows)[0].window = widget_draw(col2_base $
                                               , xsize=(*self.plot_windows)[0].x $
                                               , ysize=(*self.plot_windows)[0].y $
                                               , uvalue='Plot Window 1')
  
 

  ;; menu
  menubar = widget_base(self.tap_base,$
                        /row)
  
  row = widget_base(menubar,$
                    /row,$
                    /toolbar,$
                    /exclusive,$
                    /base_align_center)
  
  for i=0,n_elements(*self.menus)-1 do begin
     button = widget_button(row,$
                            value=' '+(*self.menus)[i]+' ',$
                            uvalue={object:self, method:'menumap', value:(*self.menus)[i], type:'main'},$
                            /no_release,$
                            font=self.buttonfont)
     if i eq 1 then widget_control, button, /SET_BUTTON
  endfor

  ;; workspace:
  (*self.bases)[0]  = widget_base(self.tap_base,$
                                  frame=5)
  
  self->setup,'parameterize'
  self->setup,'load'
  self->setup,'multi'
  self->setup,'fit'
  self->setup,'inference'
end

pro TAP::load_defaults

 ; stop
  ptr_free,self.defaults
  
  !p.font = 0                  
  plotsym,0,1.2,/fill
  
  self.setup_ld = 2
  self.diff = 0
  self.buttonfont = ''
  self.textfont=''
  ;self.create_plots = 1
  ;self.create_ascii = 1
  ;self.create_mcmcascii=0
  self.cr = string(10B)
  self.eff_len = [1,20]
  self.slash = '/'
    
  ptr_free,self.colors
  self.colors = ptr_new(tap_colors())
  
  self.init_period = 3d0
  
  self.setup_rebin = 0  
  self.smooth_val = [29.4244d0,10d0]
  
  self.defaults = ptr_new(replicate({param: '',$
                                     val: 0d0,$
                                     lock: 0d0,$
                                     limit: dblarr(3),$
                                     prior: dblarr(3)},n_elements(*self.label_indices)-1))
  
  for i=1,n_elements(*self.label_indices)-1,1 do (*self.defaults)[i-1].param = strjoin(strsplit(strlowcase(strmid((*self.label_indices)[i],0,strlen((*self.label_indices)[i])-1)),' ',/extract),'_')
  

 ; stop
  defaults = 0
  if file_test('TAP_defaults.txt') then begin
     self.message = 'Found local TAP_defaults.txt file... implementing those defaults.'
     read = 'TAP_defaults.txt'
     defaults = 1
  endif else begin
     path = strsplit((routine_info(/source))[where(strmatch((routine_info(/SOURCE)).NAME,'tap',/fold_case))].path,$
                     self.slash,/extract)
     path = self.slash+strjoin(path[0:n_elements(path)-2],self.slash)+self.slash
     if file_test(path+'TAP_MASTER_defaults.txt') then begin
        self.message = 'No local TAP_defaults.txt file... using '+path+'TAP_MASTER_defaults.txt'        
        read = path+'TAP_MASTER_defaults.txt'
        defaults = 1
     endif else begin
        self.message = 'No TAP_defaults.txt or TAP_MASTER_defaults.txt found.  Using hard-coded default settings.'
        case self.parameterize_id of
           'basic': begin
              (*self.defaults)[1].val = 89d0
              (*self.defaults)[1].limit = [1,0,90]
              (*self.defaults)[2].val = 0.1
              (*self.defaults)[2].limit = [0,0,10]
           end
           'adv1': begin
              (*self.defaults)[1].val = 0.0   
              (*self.defaults)[1].limit = [1,0,10]
              (*self.defaults)[2].val = 0.02
              (*self.defaults)[2].limit = [0,0,1]
           end
           'adv2': begin
              (*self.defaults)[1].val = 0.01            ;; b
              (*self.defaults)[1].limit = [0,-1,1]      ;;
              (*self.defaults)[2].val = 0.02            ;; T
              (*self.defaults)[2].limit = [0,0,1]       ;; 
           end
        endcase
        (*self.defaults)[0].val = 3d0
        (*self.defaults)[0].limit = [0,1,5]
        (*self.defaults)[3].val = 0.1
        (*self.defaults)[3].limit = [0,0.001,1]
        (*self.defaults)[4].val = -99
        (*self.defaults)[4].limit = [0,-.005,+.005]
        (*self.defaults)[5].val = 0.2
        (*self.defaults)[5].limit =  [0,0,1]
        (*self.defaults)[6].val = 0.2
        (*self.defaults)[6].limit =  [0,-1,1]
        (*self.defaults)[7].val = 0
        (*self.defaults)[7].limit =  [1,0,1]
        (*self.defaults)[8].val = 0
        (*self.defaults)[8].limit =  [1,0,!dpi*2d0]
        (*self.defaults)[9].val = 1
        (*self.defaults)[9].limit = [0,0.99,1.01]
        (*self.defaults)[10].val = 0
        (*self.defaults)[10].limit = [0,-0.01, 0.01]
        (*self.defaults)[11].val = 0
        (*self.defaults)[11].limit = [0,-0.001,0.001]
        (*self.defaults)[12].val = 0
        (*self.defaults)[12].limit = [0,0,1]
        (*self.defaults)[13].val = 0.001
        (*self.defaults)[13].limit = [0,0,1]
        (*self.defaults)[14].val = 0
        (*self.defaults)[14].limit = [1,0,1]
                                ;     self->message
     endelse
  endelse 
  
  if defaults then begin
    ; readcol,read,strings,format='(a500)'
    ; print,strings
     openr,lun,read,/get_lun
     str = ''
     while ~ EOF(lun) DO BEGIN  
        readf,lun,str
        print,str
        if strmatch(str,'*;*') eq 0 and strlen(strcompress(str,/remove_all)) ne 0 then begin
           arr = strsplit(str,/extract)
           if n_elements(arr) eq 9 then begin
              ;; this is a parameter default...
              location =where(strcmp(arr[0],(*self.defaults).param))
              (*self.defaults)[location].val = arr[1]*1d0
              (*self.defaults)[location].lock = arr[2]*1d0
              (*self.defaults)[location].limit = arr[[5,3,4]]*1d0
              (*self.defaults)[location].prior = arr[6:8]*1d0
           endif
        endif
        print_struct,*self.defaults
     endwhile
     free_lun,lun

  endif

end

pro TAP::lcplot,event
  !p.multi = [0,1,1]
  device,decomposed=0
  
  wset,(*self.plot_windows)[0].pix_window
  if self.num_transits ne 0 then begin
     (*self.plot_windows)[0].xrange = [1d4,-1d4]
     for i=0,n_elements(*self.transits)-1,1 do $
        (*self.plot_windows)[0].xrange = [$
        min([(*self.plot_windows)[0].xrange[0],min(((*self.transits)[i]->get()).transit_tday- $
                                                   (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit')  $
                                                                                              eq 1)].value-(((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))]),$
        max([(*self.plot_windows)[0].xrange[1],$
             max(((*self.transits)[i]->get()).transit_tday-$
                 (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value $
                  - (((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))])]
     
     
     tot = self.num_transits-1
     lci = (*self.transits)[0]->get()
     lcf = (*self.transits)[self.num_transits-1]->get()
     diff = self.diff
     if diff eq 0 then diff = 10d0*max([lci.residuals,lcf.residuals])
     
     
     yrange = [min(lci.transit_flux)-tot*(diff/2d0)-diff,max(lcf.transit_flux)+((tot+.5)*diff)]
     
     title = string(self.active_transit,format='(i2.2)')+": "+((*self.transits)[self.active_transit-1]->get()).fname
     
     for i=0,self.num_transits-1,1 do begin
        lc = (*self.transits)[i]->get() 
                                ;   diff = 5d0*max(lc.residuals)
        midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
        if i eq 0 then $
           plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,$
                yrange=yrange,color=(*self.colors).black,background=(*self.colors).white,$
                xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',$
                ytitle='Relative Flux + Constant',title=title,/nodata
        
        oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,color=lc.color
        
        oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(diff*i),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)*(i+1)),psym=8,symsize=.6,color=lc.color
        hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=2
        lc = 0L   
        midt = 0L
        
        
        if i+1 eq self.active_transit then begin
           plotsym,7,3.0,/fill,thick=5
           oplot,[(*self.plot_windows)[0].xrange[0]*24d0],[1d0+(diff*i)],psym=8,color=(*self.colors).blue,symsize=1
           plotsym,6,3.0,/fill,thick=5
           oplot,[(*self.plot_windows)[0].xrange[1]*24d0],[1d0+(diff*i)],psym=8,color=(*self.colors).blue,symsize=1
           plotsym,0,1.2,/fill
        endif
     endfor
     lci = 0L
     lcf = 0L
     tot = 0L
     yrange = 0L
  endif else plot,[0],[0],/nodata,color=(*self.colors).gray,background=(*self.colors).gray
  wset,(*self.plot_windows)[0].w_id
  device,copy=[0,0,(*self.plot_windows)[0].x,(*self.plot_windows)[0].y,0,0,(*self.plot_windows)[0].pix_window]
  
end

pro TAP::quit,event
  self->destroy
end

pro tap::loadmcmc
  self.message = 'Restoring Savefile...'
  self->message  
  tap_state = 0L
  widget_control,self.mcmc_fld[1],get_value=path
  path = path[0]
 ; obj_destroy,self.restoreSaved
  self.restoreSaved = obj_new('IDL_Savefile',path)
;  restore,path
 ; stop
  
  if strcmp(strmid(path,strlen(path)-16,16),'TAP_setup.idlsav') then self.plot_dir = strmid(path,0,strlen(path)-16) else $
     if n_elements(strsplit(path,'/')) gt 1 then self.plot_dir = '/'+strjoin((strsplit(path,'/',/extract))[0:n_elements(strsplit(path,'/'))-2],'/')+'/' $
     else self.plot_dir = '\'+strjoin((strsplit(path,'\',/extract))[0:n_elements(strsplit(path,'\'))-2],'\')+'\'
  
  if n_elements(strsplit(path,'/')) le 1 then self.backslash = 1
  
  if where(self.restoresaved->names() eq 'TAP_STATE') eq -1 then begin
     self.message = 'Incompatible savefile.'
     self->message
  endif  else begin
     restore,path
     if (tap_state[0]->get()).mcmc_complete eq 0 then begin
        self.message = 'Savefile contains no MCMC chains.'
        self->message
        for i=0,n_elements(tap_state)-1,1 do begin
           st = TAP_state[i]->get()
           for j=0,n_elements(st.params)-1,1 do begin
              ptr_free,st.params[j].mcmc_chain
              ptr_free,st.params[j].refined_mcmc_chain
           endfor
           for j=0,n_elements(st.basic_params)-1,1 do begin
              ptr_free,st.basic_params[j].mcmc_chain
              ptr_free,st.basic_params[j].refined_mcmc_chain
           endfor
           ptr_free,st.likelihood_chain,st.refined_likelihood_chain
           ptr_free,st.model_t,st.model_f,st.model_tfine,st.model_ffine,st.mcmc_files
           TAP_state[i]->destroy
           st = 0L
        endfor
     endif else begin 
        ptr_free,self.transits
        self.transits = ptr_new(TAP_state)
        TAP_state = 0L
        self.active_transit = 1
        self.num_transits = n_elements(*self.transits)
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get()   
           transit.params.jumpct = 0
           transit.params.jumptot= 0   
           (*self.transits)[i]->set,transit
           if i eq 0 then begin
              self.parameterize_id = transit.parameterize_id
           ;;   self.parameterize_num = transit.parameterize_num
              self.setup_ld = transit.ldtype
              self->setup_parameterization
           endif
           transit=0L
        endfor
        
        self.mcmc_complete =0
        transit = (*self.transits)[0]->get()   
        new = 0
        for i=0,n_elements((*transit.mcmc_files))-1,1 do begin
           if strcmp((*transit.mcmc_files)[i],'-1') ne 1 then begin
              file = (*transit.mcmc_files)[i]
              if self.backslash then strreplace,file,'/','\'
              self.mcmc_complete++
              go = self->loadintocurr(self.plot_dir+file)
              self->addtofull,new,go
              self->addtofullbasic,new,go; self->loadintocurrbasic(self.plot_dir+file)
              file = 0
              new++
           endif
        endfor
        transit = 0L
        self->setup,'adjustparams'
        self->setup,'cross lc links'
        
        self->mcmc_inference
        self->mcmc_inference_basic
        
        self->setup,'multi'
        self->setup,'fit'
        
        
        self.diff = ((*self.transits)[0]->get()).diff
     endelse
  endelse
  obj_destroy,self.restoreSaved
  
end


pro tap::addtofullbasic,cclear,go
  if cclear eq 0 then clearit=1 else clearit = 0
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()  
     range = self->burnrangebasic(tnum=i)
     
     for j=0,n_elements(transit.basic_params)-1,1 do begin
                                ;   range = [.1*n_elements(*transit.params[j].mcmc_chain),n_elements(*transit.params[j].mcmc_chain)-1]
        if clearit then begin
           ptr_free,transit.basic_params[j].refined_mcmc_chain
           if go then transit.basic_params[j].refined_mcmc_chain = ptr_new((*transit.basic_params[j].mcmc_chain)[range[0]:range[1]])
        endif else if go then *transit.basic_params[j].refined_mcmc_chain = [*transit.basic_params[j].refined_mcmc_chain,(*transit.basic_params[j].mcmc_chain)[range[0]:range[1]]]  
     endfor
     
     if n_elements(strsplit(transit.fname,'/')) gt 1 then transit.fname = (strsplit(((strsplit(transit.fname,'/',/extract))[n_elements(strsplit(transit.fname,'/',/extract))-1]),'.',/extract))[0]  
     if n_elements(strsplit(transit.fname,'\')) gt 1 then transit.fname = (strsplit(((strsplit(transit.fname,'\',/extract))[n_elements(strsplit(transit.fname,'\',/extract))-1]),'.',/extract))[0]  

     
     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
end


pro tap::addtofull,cclear,go
  if cclear eq 0 then clearit=1 else clearit = 0
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()  
     range = self->burnrange(tnum=i)
     
     for j=0,n_elements(transit.params)-1,1 do begin
                                ;   range = [.1*n_elements(*transit.params[j].mcmc_chain),n_elements(*transit.params[j].mcmc_chain)-1]
        if clearit then begin
           ptr_free,transit.params[j].refined_mcmc_chain
           if go then transit.params[j].refined_mcmc_chain = ptr_new((*transit.params[j].mcmc_chain)[range[0]:range[1]])
        endif else if go then *transit.params[j].refined_mcmc_chain = [*transit.params[j].refined_mcmc_chain,(*transit.params[j].mcmc_chain)[range[0]:range[1]]]  
        
       ; print,transit.params[j].param
       ; help,*transit.params[j].refined_mcmc_chain
        
     endfor
     if clearit then begin
        ptr_free,transit.refined_likelihood_chain
        if go then transit.refined_likelihood_chain = ptr_new((*transit.likelihood_chain)[range[0]:range[1]])
     endif else if go then *transit.refined_likelihood_chain = [*transit.refined_likelihood_chain,(*transit.likelihood_chain)[range[0]:range[1]]]  
     
     if n_elements(strsplit(transit.fname,'/')) gt 1 then transit.fname = (strsplit(((strsplit(transit.fname,'/',/extract))[n_elements(strsplit(transit.fname,'/',/extract))-1]),'.',/extract))[0]  
     if n_elements(strsplit(transit.fname,'\')) gt 1 then transit.fname = (strsplit(((strsplit(transit.fname,'\',/extract))[n_elements(strsplit(transit.fname,'\',/extract))-1]),'.',/extract))[0]  

     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
end

function tap::loadintocurr,filename
  if file_test(filename) then begin
     restore,filename
     for i=0,n_elements(tap_state)-1,1 do begin
        tt = (*self.transits)[i]->get()  
        st = TAP_state[i]->get()
        for j=0,n_elements(tt.params)-1,1 do begin
           ptr_free,tt.params[j].mcmc_chain
           tt.params[j].mcmc_chain = ptr_new(*st.params[j].mcmc_chain)
           tt.params[j].jumpct += st.params[j].jumpct
           tt.params[j].jumptot += st.params[j].jumptot
           ptr_free,st.params[j].mcmc_chain
           ptr_free,st.params[j].refined_mcmc_chain
           ptr_free,tt.basic_params[j].mcmc_chain
           tt.basic_params[j].mcmc_chain = ptr_new(*st.basic_params[j].mcmc_chain)
           tt.basic_params[j].jumpct += st.basic_params[j].jumpct
           tt.basic_params[j].jumptot += st.basic_params[j].jumptot
           ptr_free,st.basic_params[j].mcmc_chain
           ptr_free,st.basic_params[j].refined_mcmc_chain
        endfor
        ptr_free,tt.likelihood_chain
        tt.likelihood_chain = ptr_new(*st.likelihood_chain)
        (*self.transits)[i]->set,tt
        ptr_free,st.model_t,st.model_f,st.model_tfine,st.model_ffine,st.mcmc_files
        TAP_state[i]->destroy
        tt = 0L
        st = 0L
     endfor
     return,1
  endif else return,0
end

pro tap::setup_parameterization

  ptr_free,self.parameterization_options
  self.parameterization_options = ptr_new(['Basic [i, a/R*]','T, b']);,'tau_o, b^2'])
  self.parameterization_num = -1 ;

  self.basic_indices = ptr_new(['Parameters  ',$
                                'Period:',$
                                'Inclination:',$
                                'a/R*:',$
                                'Rp/R*:',$
                                'Mid Transit:',$
                                'Linear LD:',$
                                'Quadratic LD:',$
                                'Eccentricity:',$
                                'Omega:',$
                                'OOT t^0:',$
                                'OOT t^1:',$
                                'OOT t^2:',$
                                'Sigma Red:',$
                                'Sigma White:',$
                                'Delta Light:'])
  case self.parameterize_id of
     'basic': begin
        self.parameterization_num = 0
        self.label_indices = ptr_new(['Parameters  ',$
                                      'Period:',$
                                      'Inclination:',$
                                      'a/R*:',$
                                      'Rp/R*:',$
                                      'Mid Transit:',$
                                      'Linear LD:',$
                                      'Quadratic LD:',$
                                      'Eccentricity:',$
                                      'Omega:',$
                                      'OOT t^0:',$
                                      'OOT t^1:',$
                                      'OOT t^2:',$
                                      'Sigma Red:',$
                                      'Sigma White:',$
                                      'Delta Light:'])
     end
     'adv1': begin
        self.parameterization_num = 2
        self.label_indices = ptr_new(['Parameters  ',$
                                      'Period:',$
                                      'b^2:',$
                                      'tau_o:',$
                                      'Rp/R*:',$
                                      'Mid Transit:',$
                                      'Linear LD:',$
                                      'Quadratic LD:',$
                                      'Eccentricity:',$
                                      'Omega:',$
                                      'OOT t^0:',$
                                      'OOT t^1:',$
                                      'OOT t^2:',$
                                      'Sigma Red:',$
                                      'Sigma White:',$
                                      'Delta Light:'])
     end
    'adv2': begin
        self.parameterization_num = 1
        self.label_indices = ptr_new(['Parameters  ',$
                                      'Period:',$
                                      'b:',$
                                      'T:',$
                                      'Rp/R*:',$
                                      'Mid Transit:',$
                                      'Linear LD:',$
                                      'Quadratic LD:',$
                                      'Eccentricity:',$
                                      'Omega:',$
                                      'OOT t^0:',$
                                      'OOT t^1:',$
                                      'OOT t^2:',$
                                      'Sigma Red:',$
                                      'Sigma White:',$
                                      'Delta Light:'])
     end
     else: begin
        print,'unknown parameterization type: '+self.parameterizate_id
     end
  endcase
  self->load_defaults
  ;self.label = ptr_new(lonarr(n_elements(*self.label_indices)+1))
end

     
function TAP::INIT,input_ascii=input_ascii,smooth=smooth,period=period,parameterization=parameterization,version=version
  
  self.version = version+' '

  ;self.version = 'v2.5 '
  
  self.calc_gr = 1
  self.calc_effl = 1

  ptr_free,self.label_indices,self.label
  self.label = ptr_new(lonarr(20)) 
  
  if 1-keyword_set(parameterization) then $
     self.parameterize_id = 'adv2' else $
        self.parameterize_id = parameterization
  ptr_free,self.ld_parameterization_options
  self.ld_parameterization_options = ptr_new(['u1 and u2','2u1+u2, u1-2u2','(u1+u2)^2, 0.5u1(u1+u2)^-1'])
  
  ;; self.parameterize_id = 'agol'
  

  self.create_plots = 1
  self.create_ascii = 1
  self.create_mcmcascii=0

  self->setup_parameterization

  obj_destroy,self.restoreSaved
  
  ptr_free,self.slider
  self.slider = ptr_new(replicate({id: 0},22))
    
  ptr_free,self.plot_windows
  self.plot_windows = ptr_new(replicate({x: 0d, y: 0d, w_id: 0L, pw_id: 0L, window:0L, pix_window: 0L, xrange: [0d0,0d0], yrange: [0d0,0d0]},1))
  
  ptr_free,self.menus
  ptr_free,self.bases
  self.menus = ptr_new(['Parameterization','Load Transit','Manage Transits','Fit','MCMC Inference'])
  self.bases = ptr_new(lonarr(6))
  
  ptr_free,self.extra_windows
  self.extra_windows = ptr_new(lonarr(6))
  
  ptr_free,self.transits
  self.transits = ptr_new(-1)
  temp = obj_new('transit')
  temp->destroy
  
  if keyword_set(period) then self.init_period = period 
  
  if keyword_set(smooth) then begin
     self.setup_rebin=1
     self.smooth_val = smooth
  endif
  
  self->widget_setup
  
  ptr_free,self.mcmc_stuff
  
  self.mcmc_stuff = ptr_new()
  
  if keyword_set(input_ascii) then begin
     self.filetype = 'ASCII File'
     widget_control,self.transitfile_fld[1],set_value =input_ascii
  endif
  
  for i=1,n_elements(*self.bases)-1,1 do widget_control,(*self.bases)[i],MAP=0
  
  widget_control,(*self.bases)[1],MAP=1        
  
  Widget_Control,self.tap_base, set_UVALUE=self 

  self->message
  
  self->start
  if keyword_set(input_ascii) then self->loadtransit,0
  return,1
end


pro tap__define
  struct = {tap,$
            $ ;; GENERAL STUFF
            version: '',$ 
            colors: ptr_new(), $
            buttonfont: '',$
            textfont: '',$
            cr: '',$
            restoreSaved: obj_new(),$
            defaults: ptr_new(),$
            $ ;;
            parameterization_options: ptr_new(),$
            parameterization_num: 0,$
            parameterize_id: '',$
            ld_parameterization_options: ptr_new(),$
            $ ;; MAIN WIDGET stuff
            tap_base: 0L, $
            message_window: 0L, $
            message: '',$
            plot_windows: ptr_new(),$
            diff: 0d0,$
            create_plots: 0L,$
            create_ascii: 0L,$
            calc_gr: 0L,$
            calc_effl: 0L,$
            create_2d: 0L,$
            plots2d: ptr_new(-1),$
            create_mcmcascii: 0L,$
            $ ;; EXTRA WIDGET WINDOWS
            extra_windows: ptr_new(),$
            $ ;; WIDGET pieces
            slider: ptr_new(),$
            ;;buttons: ptr_new(),$
            label_indices: ptr_new(),$
            basic_indices: ptr_new(),$
            label: ptr_new(),$
            menus: ptr_new(),$
            bases: ptr_new(),$
            settings: lonarr(6),$
            fld: lonarr(50,15),$
            fld2: lonarr(100,25),$
            numcol: 0L,$
            $ ;; TRANSIT + MODEL stuff
            setup_rebin: 0L,$
            setup_LD: 0L,$
            setup_smooth: 0L,$
            smooth_val: dblarr(2),$
            num_transits: 0d,$
            active_transit: 0d,$
          ;;  transit_colors: ptr_new(),$
            transits: ptr_new(),$
            parinfo: ptr_new(),$
            backslash: 0,$
            slash: '/',$
            base_plot: '',$
            eff_len: dblarr(2),$
            $ ;; PATHS:    
            init_period: 0d0,$
            base_path: '',$
            filetype: '',$
            transitpath:'',$
            transitpath1:'',$
            transitpath2:'',$
            transitfile_fld:[0,0], $
            transitfile1_fld:[0,0], $
            transitfile2_fld:[0,0], $
            mcmc_fld:[0,0],$
            plot_dir: '' ,$
            buttons: lonarr(5),$
            mcmc_stuff: ptr_new(),$
            mcmc_complete: 0,$
            adv_opt: lonarr(2),$
            $ ;; MCMC parameters
            $ ;; END
            exists: 0L $
}
end

pro tap,_REF_EXTRA=_extra
  version = 'v2.51'

  ;; v2.40: -added advanced parameterization b^2, T
  ;;        -added additional LD parameterization
  ;; v2.50: added advanced delta_light parameter (from johnson+ 2011ApJ...730...79J equation 9)
  ;; v2.51: -fixed advanced parameterization.... b^2->b
  ;;        -added likelihood tracking and outputting to mcmc chain files

  tap = obj_new('tap',version=version,_EXTRA=_extra)
end






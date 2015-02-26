
pro tapmcmc::destroy,event
  ptr_free,self.need_adjust
  ptr_free,self.adjust
  ptr_free,self.dec
  ptr_free,self.inc
  ptr_free,self.phi
  ptr_free,self.currlikes
  ptr_free,self.newlikes
  ptr_free,self.jumps
  ptr_free,self.sets
  
  obj_destroy,self.progbarobj1
  obj_destroy,self.progbarobj2
  
  obj_destroy,self.rand_obj
  obj_destroy,self.rand_obj2

  ptr_free,self.plot_windows

  for i=0,n_elements(*self.widget_bases)-1,1 do if (*self.widget_bases)[i] then widget_control,(*self.widget_bases)[i],/destroy
  ptr_free,self.widget_bases
  
  ptr_free,self.colors
  widget_control,self.mcmc_base,/destroy
  obj_destroy,self
 end

pro tapmcmc::storelink,event,pick=pick 
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.params[pick].jumptot++
     if (*self.jumps)[i] then begin
        transit.params.curr_link = transit.params.new_link
        transit.params[pick].jumpct++
        (*self.currlikes)[i] = (*self.newlikes)[i]
                                ;transit.curr_redl = transit.new_redl
     endif
     ;; if pick eq 4 then print,i,transit.params[pick].curr_link
     if (self.iter mod 10) eq 0 then begin
        transit = self->disentangle_parameterization(transit,/curr)
        *transit.likelihood_chain = [*transit.likelihood_chain,(*self.currlikes)[i]]
        for j=0,n_elements(transit.params)-1,1 do begin
           if strcmp(transit.params[j].param,'Linear LD') eq 1 then begin
              c1 = transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].curr_link
              c2 = transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].curr_link
              u1 = c1
              u2 = c2
              if self.deconLD eq 1 then begin
                 u1 = (c2+2d0*c1)/5d0
                 u2 = (c1-2d0*c2)/5d0
              endif
              if self.deconLD eq 2 then begin
                 u1 = 2*sqrt(c1)*c2
                 u2 = sqrt(c1)*(1-2*c2)
              endif
              *transit.params[j].mcmc_chain = [*transit.params[j].mcmc_chain,u1]
              *transit.basic_params[j].mcmc_chain = [*transit.basic_params[j].mcmc_chain,u1]
              
           endif else begin
              if strcmp(transit.params[j].param,'Quadratic LD') eq 1 then begin
                 *transit.params[j].mcmc_chain = [*transit.params[j].mcmc_chain,u2]
                 *transit.basic_params[j].mcmc_chain = [*transit.basic_params[j].mcmc_chain,u2]

                 u1 = 0L
                 u2 = 0L
                 c1 = 0L
                 c2 = 0L
              endif else begin 
                 *transit.params[j].mcmc_chain = [*transit.params[j].mcmc_chain,transit.params[j].curr_link]
                 *transit.basic_params[j].mcmc_chain = [*transit.basic_params[j].mcmc_chain,transit.basic_params[j].curr_link]
              endelse
           endelse
                                ;  *transit.params[0].likelihoods = [*transit.params[0].likelihoods,(*self.currlikes)[i]]
        endfor 
        
        self->updatemod    
        if i eq 0 then begin
           if (self.iter mod 1000) eq 0 then begin
              self->curr_analyze
              self.base_widget->prepcol,runval=self.runval,es='t'+string(self.runval,format='(i2.2)')+' '
           endif
           if (self.iter mod 100) eq 0 then self->plot
        endif
     endif
     
     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
  (*self.jumps)*=0d0
end

pro TAPmcmc::lcplot,event
  
  !p.multi = [0,1,1]
                                ;device,decomposed=0
  wset,(*self.plot_windows)[0].pix_window
  (*self.plot_windows)[0].xrange = [1d4,-1d4]
  
  if self.which_plot ge self.num_transits then begin
     for i=0,n_elements(*self.transits)-1,1 do $
        (*self.plot_windows)[0].xrange = [$
        min([(*self.plot_windows)[0].xrange[0],1d0*min(((*self.transits)[i]->get()).transit_tday- $
                                                       (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit')  $
                                                                                                  eq 1)].value-(((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))]),$
        max([(*self.plot_windows)[0].xrange[1],$
             max(((*self.transits)[i]->get()).transit_tday-$
                 (((*self.transits)[i]->get()).params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value $
                  - (((*self.transits)[i]->get()).params[0].value*((*self.transits)[i]->get()).epoch)))])]
     
     tot = self.num_transits-1
     lci = (*self.transits)[0]->get()
     lcf = (*self.transits)[self.num_transits-1]->get()
     diff = ((*self.transits)[0]->get()).diff
     if diff eq 0 then diff = 5d0*max([lci.residuals,lcf.residuals])
     
     if self.phased then begin
        yrange = [min(lci.transit_flux)-(diff/2d0)-diff,max(lcf.transit_flux)+diff]
        
                                ;  diff = 0L
        
        for i=0,self.num_transits-1,1 do begin
        lc = (*self.transits)[i]->get() 
                                ;   diff = 5d0*max(lc.residuals)
        midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
        if i eq 0 then $
           plot,(lc.transit_tday-midt)*24d0,lc.transit_flux,psym=8,symsize=.6,yrange=yrange,color=lc.color,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',ytitle='Rel. Flux + Const' else $
              oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux,psym=8,symsize=.6,color=lc.color
        
   ;     oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)),psym=8,symsize=.6,color=lc.color
    ;    hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=2
    ;    oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=2
    ;    oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=2
        lc = 0L   
        midt = 0L
     endfor
     lci = 0L
     lcf = 0L
     tot = 0L
     yrange = 0L
     
  endif else begin 
     yrange = [min(lci.transit_flux)-tot*(diff/2d0)-diff,max(lcf.transit_flux)+((tot+.5)*diff)]
     
     for i=0,self.num_transits-1,1 do begin
        
        lc = (*self.transits)[i]->get() 
                                ;   diff = 5d0*max(lc.residuals)
        midt = lc.params[where(strcmp(((*self.transits)[i]->get()).params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
        if i eq 0 then $
           plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,yrange=yrange,color=lc.color,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',ytitle='Rel. Flux + Const' else $
              oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,color=lc.color
        
     oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(diff*i),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)*(i+1)),psym=8,symsize=.6,color=lc.color
        hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=2
        oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=2
        lc = 0L   
        midt = 0L
     endfor
     lci = 0L
     lcf = 0L
     tot = 0L
     yrange = 0L
     ;; plot,((*self.transits)[0]->get()).transit_tday-((*self.transits)[0]->get()).params[where(strcmp(((*self.transits)[0]->get()).params.param,'Mid Transit') eq 1)].value,((*self.transits)[0]->get()).transit_flux,psym=8,symsize=.4,background=(*self.colors).white,color= ((*self.transits)[0]->get()).color,xrange=(*self.plot_windows)[0].xrange
     
  endelse
endif  else begin
   (*self.plot_windows)[0].xrange = [min([(*self.plot_windows)[0].xrange[0],min(((*self.transits)[self.which_plot]->get()).transit_tday-(((*self.transits)[self.which_plot]->get()).params[where(strcmp(((*self.transits)[self.which_plot]->get()).params.param,'Mid Transit') eq 1)].value-(((*self.transits)[self.which_plot]->get()).params[0].value*((*self.transits)[self.which_plot]->get()).epoch)))]),$
                                     max([(*self.plot_windows)[0].xrange[1],$
                                          max(((*self.transits)[self.which_plot]->get()).transit_tday-$
                                              (((*self.transits)[self.which_plot]->get()).params[where(strcmp(((*self.transits)[self.which_plot]->get()).params.param,'Mid Transit') eq 1)].value - (((*self.transits)[self.which_plot]->get()).params[0].value*((*self.transits)[self.which_plot]->get()).epoch)))])]
   i=0
   lc = (*self.transits)[self.which_plot]->get()
   diff = ((*self.transits)[0]->get()).diff
   if diff eq 0 then diff = 5d0*max(lc.residuals)
   midt = lc.params[where(strcmp(lc.params.param,'Mid Transit') eq 1)].value - (lc.epoch*lc.params[0].value)
   yrange = [min(lc.transit_flux)-(diff/2d0)-diff,max(lc.transit_flux)+(.5*diff)]

   plot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),psym=8,symsize=.6,yrange=yrange,color=(*self.colors).black,background=(*self.colors).white,xrange=(*self.plot_windows)[0].xrange*24d0,/xs,/ys,xtitle='Hours from Mid Transit',ytitle='Rel. Flux + Const',title=string(self.which_plot+1,format='(i2.2)')+': '+lc.fname
   oplot,(lc.transit_tday-midt)*24d0,lc.transit_flux+(diff*i),color=lc.color,psym=8,symsize=.6
   oplot,(*lc.model_tfine-midt)*24d0,*lc.model_ffine+(diff*i),color=lc.modcol,thick=2
   oplot,(lc.transit_tday-midt)*24d0,lc.residuals+yrange[0]+((diff/2d0)*(i+1)),psym=8,symsize=.6,color=lc.color
   hline,yrange[0]+(diff/2d0*(i+1)),color=lc.modcol,thick=2
   oplot,(lc.transit_tday-midt)*24d0,*lc.model_f+(diff*i)+lc.rednoise,color=lc.redcol,thick=2
   oplot,(lc.transit_tday-midt)*24d0,yrange[0]+((diff/2d0)*(i+1))+lc.rednoise,color=lc.redcol,thick=2
   lc = 0L
endelse


  wset,(*self.plot_windows)[0].w_id
  device,copy=[0,0,(*self.plot_windows)[0].x,(*self.plot_windows)[0].y,0,0,(*self.plot_windows)[0].pix_window]
  
end

pro tapmcmc::Message,event
  self.message = '('+curr_date(format='hh:mm:ss')+') '+self.message
  widget_control, (*self.widget_bases)[0],/append,$
                  set_value=self.message

  self.message=''
end

pro tapmcmc::setuprun,throw=throw
  if 1-keyword_set(throw) then throw=0d
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     tags = n_elements(transit.params.param)
     for j=0,tags-1,1 do begin
        transit.params[j].jumpct = 0d0
        transit.params[j].jumptot = 0d0
        ptr_free,transit.params[j].mcmc_chain
        transit.params[j].mcmc_chain = ptr_new(transit.params[j].value)
        if transit.params[j].limited[0] then $
           if (*transit.params[j].mcmc_chain)[0] lt transit.params[j].limits[0] then $
              (*transit.params[j].mcmc_chain)[0] = transit.params[j].limits[0]      
        if transit.params[j].limited[1] then $
           if (*transit.params[j].mcmc_chain)[0] gt transit.params[j].limits[1] then $
              (*transit.params[j].mcmc_chain)[0] = transit.params[j].limits[1]      
        transit.params[j].curr_link = (*transit.params[j].mcmc_chain)[0]
     endfor
     tags = n_elements(transit.basic_params.param)
     for j=0,tags-1,1 do begin
        transit.basic_params[j].jumpct = 0d0
        transit.basic_params[j].jumptot = 0d0
        ptr_free,transit.basic_params[j].mcmc_chain
        transit = self->disentangle_parameterization(transit,/val)
        transit.basic_params[j].mcmc_chain = ptr_new(transit.basic_params[j].value)
        transit.basic_params[j].curr_link = (*transit.basic_params[j].mcmc_chain)[0]
     endfor
     
     transit.new_redl = 0L
     transit.curr_redl = 0L
     (*self.transits)[i]->set,transit
     transit = 0L
     tags = 0L
  endfor
          ;stop
          
  
  if keyword_set(throw) then begin
     self.message = 'throwing parameters!'
     self->message
     ;; dothrow =  ['Period','Inclination','a/R*','Rp/R*','Mid
     ;; Transit','Linear LD','Quadratic LD','Eccentricity','Omega']
     case self.parameterize_id of
        'basic': dothrow = ['Period','Inclination',$
                            'a/R*','Rp/R*','Mid Transit',$
                            'Eccentricity','Omega']
        'adv1': dothrow=['Period','b','tau_o','Rp/R*','Mid Transit']
        'adv2': dothrow=['Period','b','T','Rp/R*','Mid Transit']
     endcase
;;loop over parameters
     for j=0,n_elements(((*self.transits)[0]->get()).params)-1,1 do begin
        if max(strcmp(dothrow,((*self.transits)[0]->get()).params[j].param)) then begin
           sets = (*self.sets)[*,j]
           uniqsets =(sets[sort(sets)])[uniq(sets[sort(sets)])]
           uniqvals = uniqsets*0d0+sqrt(-1)

         ;  if j eq 4 then stop
           for i=0,n_elements(uniqsets)-1,1 do begin
              if abs(((*self.transits)[(where(sets eq uniqsets[0]))[0]]->get()).params[j].fixed-1) and $
                 ((*self.transits)[(where(sets eq uniqsets[0]))[0]]->get()).params[j].prior[0] ne 2 then begin
                 uniqvals[i] = ((*self.transits)[(where(sets eq uniqsets[i]))[0]]->get()).params[j].curr_link+$
                               (self.rand_obj->getrandomnumbers(1,/normal,/double)*$
                                (((*self.transits)[(where(sets eq uniqsets[i]))[0]]->get()).params[j].beta)*throw)
                 if ((*self.transits)[(where(sets eq uniqsets[i]))[0]]->get()).params[j].limited[0] then $
                    if uniqvals[i] lt ((*self.transits)[(where(sets eq uniqsets[i]))[0]]->get()).params[j].limits[0] then $
                       uniqvals[i] = ((*self.transits)[(where(sets eq uniqsets[i]))[0]]->get()).params[j].limits[0]
                 if ((*self.transits)[(where(sets eq uniqsets[i]))[0]]->get()).params[j].limited[1] then $
                    if uniqvals[i] gt ((*self.transits)[(where(sets eq uniqsets[i]))[0]]->get()).params[j].limits[1] then $
                       uniqvals[i] = ((*self.transits)[(where(sets eq uniqsets[i]))[0]]->get()).params[j].limits[1]
              endif
              
         ;     if j eq 4 then stop
           endfor
         ;  if j eq 4 then stop
           for i=0,n_elements(*self.transits)-1,1 do begin
              transit = (*self.transits)[i]->get()
              if abs(transit.params[j].fixed-1) then transit.params[j].curr_link =  uniqvals[where(uniqsets eq transit.params[j].set)]
              (*self.transits)[i]->set,transit
              transit = 0L
           endfor
        endif
     endfor  
  endif
  self.jumpcount*=0
  self.jumptot*=0
  self.iter = 0d0
  self->likelihood,/current,/setchain
  self->updatemod
  self->plot
  ;stop
  
end

pro tapmcmc::updatemod
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.params.value = transit.params.curr_link
     rednoise = dblarr(n_elements(transit.transit_tday))
    
     transit = self->disentangle_parameterization(transit,/val)
     
                                ;if self.deconLD then $
     dum = TAP_MPFit_function(transit.basic_params.value,time=*transit.model_tfine,flux=*transit.model_ffine,$
                              dfit=dfit,finefit=ffit,redfit=rednoise,tdat=transit.transit_tday,fdat=transit.transit_flux,$
                              rebin=transit.rebin,smooth=transit.t_int,deconLD=self.deconLD,adv_opt=transit.adv_opt) 
                                ;    else $
                                ;       dum = TAP_MPFit_function(transit.basic_params.value,time=*transit.model_tfine,flux=*transit.model_ffine,$
                                ;                                dfit=dfit,finefit=ffit,redfit=rednoise,tdat=transit.transit_tday,fdat=transit.transit_flux,$
                                ;                                rebin=transit.rebin,smooth=transit.t_int)
     
     transit.residuals=transit.transit_flux - dfit
     dum = 0L
     transit.rednoise=rednoise
     *transit.model_f = dfit
     transit.residuals = transit.transit_flux - dfit
     *transit.model_ffine = ffit
     redfit=0L
     ffit = 0L
     dfit = 0L
     
     (*self.transits)[i]->set,transit
     rednoise = 0L
     transit = 0L
  endfor
end

pro tapmcmc::guessbetas
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.params.beta = transit.params.value/1d3
     transit.params.beta = [6d-5,0.07,0.05,0.0004,0.0004,0.001,0.001,0.0001,0.0001,0.0001,1.8d-05,0.001,0.0001,7d-05,0d0]
     transit.params[4].beta = 1d-2
     if (where(transit.params.beta eq 0))[0] ne -1 then transit.params[where(transit.params.beta eq 0)].beta = 1d-4
     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
end

pro tapmcmc::plot
  xtitle= 'Parameter Value'
  histplot = -1
  bin = 25
  maxthick = 6

  if (self.iter mod 1000) eq 0 and self.iter gt 250 then begin
                                ;  print,'plotmode!
     if (self.iter mod 3000) eq 0 then self.plot+=1
                                ;if (self.iter mod 1000) eq 0 then begin
     if self.which_plot gt self.num_transits then self.which_plot = -1
     self.which_plot +=1
                                ;endif
     ;if self.plot ne 14 then begin
        k = 0
        while k lt 15 do begin
           if self.plot ge 13 then self.plot = 0
           if (self.jumpcount)[self.plot] ge 20 and self.plot ne 4 then k = 15 else begin
              k++
              self.plot +=1
           endelse
        endwhile
        k = 0L
     ;endif
     histplot = self.plot
  endif
  

  if histplot ge 0 then begin
     xr = [1d50,-1d50]
     yr = [0,0]
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        r = [.1*n_elements(*transit.params[self.plot].mcmc_chain),n_elements(*transit.params[self.plot].mcmc_chain)-1]
        xr = [min([xr[0],min((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]])]),max([xr[1],max((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]])])]
        transit = 0L
        r = 0L
     endfor
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        r = [.1*n_elements(*transit.params[self.plot].mcmc_chain),n_elements(*transit.params[self.plot].mcmc_chain)-1]
        if robust_sigma((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]]) ne 0 then begin
           plothist,(*transit.params[self.plot].mcmc_chain)[r[0]:r[1]],tx,ty,bin=(max(xr)-min(xr))/bin,/noplot
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
     !p.multi=[0,1,1]
     window = self.base_widget->window()
     wset,window.pix_window
     first = 0
     for i=0,n_elements(*self.transits)-1,1 do begin
        transit = (*self.transits)[i]->get()
        if transit.params[self.plot].prior[0] eq 1 then extra = ' (Gaussian Penalty)' else if transit.params[self.plot].prior[0] eq 2 then extra = ' (Gaussian Prior)' else extra=''
        
        r = [.1*n_elements(*transit.params[self.plot].mcmc_chain),$
             n_elements(*transit.params[self.plot].mcmc_chain)-1]
        if first eq 0 then begin
           if robust_sigma((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]]) ne 0 then begin
              plothist,(*transit.params[self.plot].mcmc_chain)[r[0]:r[1]], $
                       bin=(max(xr)-min(xr))/bin,color=((*self.transits)[i]->get()).color, $
                       background=(*self.colors).white,axiscolor=(*self.colors).black,xtitle=xtitle, $
                       title=((*self.transits)[i]->get()).params[self.plot].param+extra,xrange=xr,/xs, $
                       thick=maxthick,xticks=3,xthick=2,ythick=2,yrange=yr,/ys
              first = 1
           endif
        endif else begin
           if robust_sigma((*transit.params[self.plot].mcmc_chain)[r[0]:r[1]]) ne 0 then begin
              plothist,(*transit.params[self.plot].mcmc_chain)[r[0]:r[1]], $
                       bin=(max(xr)-min(xr))/bin,color=((*self.transits)[i]->get()).color, $
                       background=(*self.colors).white,axiscolor=(*self.colors).black,xtitle=xtitle, $
                       title=((*self.transits)[i]->get()).params[self.plot].param,/overplot,xrange=xr, $
                       /xs, xticks=3,thick=maxthick*(n_elements(*self.transits)-i)/n_elements(*self.transits),yrange=yr,/ys
           endif
        endelse
        if transit.params[self.plot].prior[0] gt 0 then begin
           vline,transit.params[self.plot].prior[1],thick=3,color=((*self.transits)[i]->get()).color
           vline,transit.params[self.plot].prior[1]+transit.params[self.plot].prior[2],thick=3,color=((*self.transits)[i]->get()).color,linestyle=2
           vline,transit.params[self.plot].prior[1]-transit.params[self.plot].prior[2],thick=3,color=((*self.transits)[i]->get()).color,linestyle=2
        endif
        
        transit = 0L
        r = 0L
     endfor
     xr = 0L
     yr = 0L
     wset,window.w_id
     device,copy=[0,0,window.x,window.y,0,0,window.pix_window]
     window = 0L
  endif
  bin = 0L
  histplot = 0L
  
  
  self->lcplot
                            
  
  
;stop
end

pro tapmcmc::likelihood,current=current,new=new,pass=pass,setchain=setchain
  likelihood = 0L
  if keyword_set(new) then *self.newlikes *= 0d0 else  *self.currlikes *= 0d0
  for i=0,n_elements(*self.transits)-1,1 do begin
     calc_likelihood = 1
     if keyword_set(pass) then if pass[i] eq 0 then calc_likelihood = 0
     if calc_likelihood then begin
        transit = (*self.transits)[i]->get()
        if keyword_set(current) then $
           transit = self->disentangle_parameterization(transit,/curr) else transit = self->disentangle_parameterization(transit,/new) ;; this has been moved to newlink. 
        if keyword_set(current) then $
           params = transit.basic_params.curr_link else $
              params=transit.basic_params.new_link
        dum = TAP_MPFit_function(params,time=*transit.model_tfine,flux=*transit.model_ffine,dfit=dfit,finefit=ffit,redfit=rednoise,tdat=transit.transit_tday,fdat=transit.transit_flux,rebin=transit.rebin,smooth=transit.t_int,deconLD=self.deconLD,adv_opt=transit.adv_opt) ;$
        
        likelihood = waveletlike((transit.transit_flux-dfit),params[12],params[13],/zeropad)
                                ;transit = 0L
                                ;transit = (*self.transits)[i]->get()
        
        if (where(transit.params.prior[0] eq 1))[0] ne -1 then begin
           these = where(transit.params.prior[0] eq 1)
           for j=0,n_elements(these)-1,1 do likelihood -= ((params[these[j]]-transit.params[these[j]].prior[1])^2d0/(2*(transit.params[these[j]].prior[2])^2d0))
           these = 0L
        endif
     endif else likelihood = 0d0 ;; this parameter for this transit did not pass, automatic fail, likelihood = 0
     
                                ; stop
     
     ;;  stop
     
     if self.reduce_like eq 0 then begin
        if keyword_set(new) then begin
                                ;transit.new_redl = likelihood 
           (*self.newlikes)[i] = likelihood
           if i eq 0 then self.new_redl = likelihood else self.new_redl += likelihood
        endif else begin
                                ;transit.curr_redl = likelihood
           (*self.currlikes)[i] = likelihood
                                ;if finite(transit.curr_redl) eq 0 then transit.curr_redl = -1d5
           if finite((*self.currlikes)[i]) eq 0 then  (*self.currlikes)[i] = -1d5
           if i eq 0 then self.curr_redl = likelihood else self.curr_redl += likelihood
        endelse
     endif else begin
        dd = 2
        if keyword_set(new) then begin
                                ;transit.new_redl = likelihood 
           (*self.newlikes)[i] = likelihood/(dd*(n_elements(transit.transit_tday)-n_elements(where(transit.params.fixed ne 1))))
           if i eq 0 then self.new_redl = likelihood/(dd*(n_elements(transit.transit_tday)-n_elements(where(transit.params.fixed ne 1)))) else self.new_redl += likelihood /(dd*(n_elements(transit.transit_tday)-n_elements(where(transit.params.fixed ne 1))))
        endif else begin
                                ;transit.curr_redl = likelihood
           (*self.currlikes)[i] = likelihood /(dd*(n_elements(transit.transit_tday)-n_elements(where(transit.params.fixed ne 1))))
                                ;if finite(transit.curr_redl) eq 0 then transit.curr_redl = -1d5
           if finite((*self.currlikes)[i]) eq 0 then  (*self.currlikes)[i] = -1d5
           
           if i eq 0 then self.curr_redl = likelihood /(dd*(n_elements(transit.transit_tday)-n_elements(where(transit.params.fixed ne 1)))) else self.curr_redl += likelihood /(dd*(n_elements(transit.transit_tday)-n_elements(where(transit.params.fixed ne 1))))
        endelse
     endelse
     
     if keyword_set(setchain) then begin
        ptr_free,transit.likelihood_chain
        transit.likelihood_chain = ptr_new((*self.currlikes)[i])
        (*self.transits)[i]->set,transit
     endif
     transit = 0L 
     likelihood = 0L
  endfor  
end

pro tapmcmc::testheap,title
  print,title
  title=0L
  help,/heap
 stop
  
end

pro tapmcmc::burnbetas,event,guessbetas=guessbetas
  if self.testheap then self->testheap,'burnbetas 0'
     
  lim1 = 1.15d0
  lim2 = .85d0

  coarse_adjust = 1
  if keyword_set(fine) then coarse_adjust = 0
  fine_adjust=2
  if keyword_set(guessbetas) then if coarse_adjust ne 0 or fine_adjust ne 0 then self->guessbetas
  self->setuprun
  initial_lock = 0d
  stabilized   = 0d 
  self.message = ' Stabilizing Betas'
  if 1-keyword_set(fine) then   self->message
 ; if self.testheap then self->testheap,'burnbetas 1'

  if coarse_adjust then begin
     tot_iter = 2d3
     while not initial_lock do begin
        tot_iter--
      ;  if self.testheap then self->testheap,'pre newlink'
        self->newlink
      ;  if self.testheap then self->testheap,'post newlink'
        if (self.iter mod 500d) eq 0 then begin
           for i=0,n_elements(*self.transits)-1,1 do begin
              ptr_free,self.adjust,self.dec,self.inc
              transit = (*self.transits)[i]->get()
              self.adjust = ptr_new(transit.params.beta*0d0)
              self.dec = ptr_new(where((transit.params.jumpct/transit.params.jumptot) gt transit.params.accept*1.25d0))
              self.inc = ptr_new(where((transit.params.jumpct/transit.params.jumptot) le transit.params.accept*0.75d0))
              fixit = ptr_new(where((transit.params.fixed) eq 1))     
              transit = 0L
              if (*self.dec)[0] ge 0 then (*self.adjust)[(*self.dec)] = 1.d0
              if (*self.inc)[0] ge 0 then (*self.adjust)[(*self.inc)] = -1.d0
              if (*fixit)[0] ge 0 then (*self.adjust)[(*fixit)] = 0.d0
              if (where(*self.adjust ne 0))[0] ge 0 then begin
                 ;;self.beta_adjusttype= 'crude' ; crude adjustment
                 
                 ;;if self.testheap then self->testheap,'burnbetas: Pre crude adhyst'
                 self->AdjustBeta,'crude',i
                 ;;if self.testheap then self->testheap,'burnbetas: Post crude adhyst'
                                  
                 
                                ;  print,'  adjusting:', (*self.adjust)
                                ;  print,'  NEW BETAS:', (*self.betas_store)
              endif else begin
                 initial_lock=1
              endelse
              ptr_free,fixit
           endfor
        endif
        if tot_iter le 0 then initial_lock = 1
      ;;; MORE
     endwhile

     tot_iter = 0L
     self.message = '  Coarse Stabilization Done'
     self->message
   ;  if self.testheap then self->testheap,'burnbetas 2: post coarse'


 ;;; MORE?
  endif else initial_lock = 1
 
  
  while fine_adjust gt 0 and self.iter le 1d5 do begin
     if self.testheap then self->testheap,'burnbetas 3: fineloop'
     
     while stabilized lt n_elements(*self.transits) and self.iter le 2d4 do begin 
       ; if self.testheap then self->testheap,'burnbetas 3.1: pre newlink'
     
        for z=0,100,1 do self->newlink       
        stabilized = 0
        for i=0,n_elements(*self.transits)-1,1 do begin
           cs = 0
           transit = (*self.transits)[i]->get()
           check = where(transit.params.fixed eq 0 and transit.params.prior[0] ne 2)
           if self.debug then begin
              if i eq 0 then begin
                 ;if (self.iter mod 250d0) eq 0 or self.iter lt 3 then begin
                    print,'-1? : ',(where((transit.params[check].accept*(lim1+(.1*(max([self.iter-6d3,0])/6d3))) lt (transit.params[check].jumpct/transit.params[check].jumptot)) or $
                                          (transit.params[check].accept*(lim2-(.1*(max([self.iter-6d3,0])/6d3))) gt (transit.params[check].jumpct/transit.params[check].jumptot))))
                    print,'1d4< ',self.iter       
                    print,transit.params[check].param
                    print,transit.params[check].jumptot
                    print,'-1?  ',(where(transit.params[check].jumptot le (200d0-min([150d0,(max([self.iter-6d3,0])/250d0)]))))
                 
                 ;endif
              endif
           endif
           cs =  (where((transit.params[check].accept*(lim1+(.1*(max([self.iter-6d4,0])/6d3))) lt (transit.params[check].jumpct/transit.params[check].jumptot)) or (transit.params[check].accept*(lim2-(.1*(max([self.iter-6d3,0])/6d3)))  gt (transit.params[check].jumpct/transit.params[check].jumptot))))
           cs2 = where(transit.params[check].jumptot le (200d0-min([150d0,(max([self.iter-6d3,0])/250)])))
           if self.debug and (self.iter mod 500) eq 0 then print,cs
           if cs[0] eq -1 and  self.iter gt 1d4 and cs2[0] eq -1 then begin
              if self.debug then if i eq 0 then begin
                 print,'GOOD!: fine_adj(1)',fine_adjust
                 print,i,stabilized
                                ;   stop
              endif
              stabilized += 1
              if fine_adjust eq 1 and i eq 0 then begin
                 self.message = '  Fine Stabilization Done'
                 self->message
              endif
           endif 
           if cs[0] ne -1 then begin
              ptr_free,self.adjust,self.need_adjust
              self.adjust = ptr_new(transit.params.beta*0d0)
              self.need_adjust = ptr_new(where(((transit.params.jumpct/transit.params.jumptot) - transit.params.accept)^2d0 gt $
                                               (transit.params.sval*((transit.params.accept*(1-transit.params.accept))/transit.params.jumptot)) and $
                                               transit.params.jumptot gt 1d2))
              if (*self.need_adjust)[0] ne -1 then (*self.adjust)[(*self.need_adjust)] = 1.d0
              
          
              if (where((*self.adjust) eq 1))[0] ne -1 then begin
                 if self.debug then print,stabilized
                 self->AdjustBeta,'fine',i
                 transit = 0L
                 transit = (*self.transits)[i]->get()
                 
                 if (where(transit.params.beta gt 1))[0] ne -1 then begin
                     
                    if self.debug then begin
                       print,'RESETTING a BETA!'
                       print,transit.params.beta
                       
                    endif
                    transit.params[where(transit.params.beta gt 1)].jumpct = 0d0
                    transit.params[where(transit.params.beta gt 1)].jumptot = 0d0
                    transit.params[where(transit.params.beta gt 1)].beta = .1d0
                    if self.debug then print,transit.params.beta
                    (*self.transits)[i]->set,transit
                ;    if self.debug then stop
                 endif
                 
                                ;self->message
              endif             ;else if  stabilized += 1
           endif   
           transit = 0L
           check= 0L
           
        endfor
     endwhile
     stabilized = 0L
     fine_adjust -= 1
  endwhile
  if self.testheap then self->testheap,'burnbetas 5: post fine'
  if  self.iter le 1d5 then self.message=' Betas have stabilized'  else $
     self.message = ' Issue with Beta stabilization--check outputs.'
  self->message
                                ;stop
end

pro tapmcmc::adjustbeta,type,active
  lim1 = 1.2d0
  lim2 = 0.8d0
  transit = (*self.transits)[active]->get()
  if self.debug then begin
     print,active
     print,'  rates:    ', (transit.params.jumpct/transit.params.jumptot)
     print,'   tots:    ', (transit.params.jumptot)
     print,'  betas:    ', transit.params.beta
     print,'  change:   ', *self.adjust
     print,'  iters:    ', self.iter, .44*(lim1+(.1*(max([self.iter-6d3,0])/6d3))),.44*(lim2-(.1*(max([self.iter-6d3,0])/1d4))),200d0-min([150d0,(max([self.iter-6d3,0])/250)])
     print,''
  endif
  case type of 
     'crude': begin
        if (where((*self.adjust) eq  1))[0] ge 0 then transit.params[where((*self.adjust) eq  1)].beta *= 1.5d0
        if (where((*self.adjust) eq -1))[0] ge 0 then transit.params[where((*self.adjust) eq -1)].beta /= 1.5d0
     end
     'fine':begin
        ptr_free,self.phi
        self.phi= ptr_new(dblarr(n_elements(transit.params.beta))*0d)
        area = where(transit.params.jumpct/transit.params.jumptot lt .1)
        if area[0] ne -1 then (*self.phi)[area] = .5d0
        area = where((0.5d*(transit.params.accept)) lt (transit.params.jumpct/transit.params.jumptot))
        if area[0] ne -1 then (*self.phi)[area] = 1d0
        area = where(((0.2d*transit.params.accept lt (transit.params.jumpct/transit.params.jumptot)) and $
                      ((0.5d*transit.params.accept) ge (transit.params.jumpct/transit.params.jumptot))))
        if area[0] ne -1 then  (*self.phi)[area] = 1.5d0
        area = where(((0.1d*transit.params.accept lt (transit.params.jumpct/transit.params.jumptot)) and $
                      ((0.2d*transit.params.accept) ge (transit.params.jumpct/transit.params.jumptot))))
        if area[0] ne -1 then   (*self.phi)[area] = 2d0
        
        transit.params[where((*self.adjust) eq 1)].beta = (transit.params.beta*((((transit.params.jumpct+1)/transit.params.jumptot)/transit.params.accept)^(*self.phi)))[where((*self.adjust) eq 1)]
        

        inc = 0
        dec = 0
        
        for i=0,n_elements(transit.params.beta)-1,1 do begin
           if (*self.adjust)[i] then begin
              inc = -1                                                                                                ; assume decrease
              if ((transit.params[i].jumpct/transit.params[i].jumptot)/(transit.params[i].accept)) le 1 then inc = 1 ; increase
              
              case transit.params[i].b_last of
                 0: transit.params[i].b_last = inc
                 -1: if inc eq 1 then transit.params[i].sval++
                 1: if inc eq -1 then transit.params[i].sval++
              endcase
              transit.params[i].b_last = inc     
           endif
        endfor
        transit.params[where(*self.adjust eq 1)].jumpct = 0d
        transit.params[where(*self.adjust eq 1)].jumptot = 0d
     end
  endcase
  (*self.transits)[active]->set,transit
  transit=0L
end

;; is a new parameter within physical bounds?
function tapmcmc::test_param,transit,pick
  ;; depends on what kind of param
  
  ;; is the parameter a limb darkening param?
  if max(strcmp(transit.params[pick].param,['Linear LD','Quadratic LD'])) eq 1 then begin ;; if u1, u2, make sure physical
     ;; are we using a deconvolved parameterization?  
     
     case self.deconLD of
        0: begin
           u1 = transit.params[where(strcmp(transit.params.param,'Linear LD'))].new_link
           u2 = transit.params[where(strcmp(transit.params.param,'Quadratic LD'))].new_link
        end
        1: begin
           u1 = (transit.params[where(strcmp(transit.params.param,'Quadratic LD'))].new_link $
                 +2d0*transit.params[where(strcmp(transit.params.param,'Linear LD'))].new_link)/5d0
           u2 = (transit.params[where(strcmp(transit.params.param,'Linear LD'))].new_link $
                 -2d0*transit.params[where(strcmp(transit.params.param,'Quadratic LD'))].new_link)/5d0
        end
        2: begin
           u1 = 2*sqrt(transit.params[where(strcmp(transit.params.param,'Linear LD'))].new_link)*$
                transit.params[where(strcmp(transit.params.param,'Quadratic LD'))].new_link
           u2 = sqrt(transit.params[where(strcmp(transit.params.param,'Linear LD'))].new_link)*$
                (1-2*transit.params[where(strcmp(transit.params.param,'Quadratic LD'))].new_link)
        end
     endcase
     if u1 lt 0 or u1+u2 gt 1d0 or u1+2*u2 lt 0d0 then return,0 ;; FAIL!
     ;; if u1 ge 1 or u2 ge 1 or u2 le -1 then pass[i] = 0
  endif
  
  ;; is the paramter one that must be positive, but is negative?
  if max(strcmp(transit.params[pick].param,['Sigma Red','Sigma White','Eccentricity'])) eq 1 then $
     if transit.params[pick].new_link lt 0 then return,0 ;; FAIL
  
  ;; is the paramter above or below a user set limit for the parameter?
  
  ;; below limit?
  if transit.params[pick].limited[0] eq 1 then if transit.params[pick].new_link lt transit.params[pick].limits[0] then return,0 ;; FAIL
  ;; above limit?
  if transit.params[pick].limited[1] eq 1 then if transit.params[pick].new_link ge transit.params[pick].limits[1] then return,0 ;; FAIL
  
  return,1 ;; if we get here, the parameter PASSES 
end


;; newlink
;; adds a link to the MCMC chain
;;   - two possible outcomes, either the current link is duplicated or
;;     a new_link is judged as "better" than the current link and
;;     added to the chain.
;;
;; in section A:
;;    Select the parameter to adjust.  TAP changes one parameter at a
;;    time (Gibbs sampler), a critical component in analyzing more
;;    than one transit at a time.  The loop there randomly picks until
;;    a non locked parameter is chosen.
;;
;;
;;
;;
;;
;;
;;
pro tapmcmc::newlink

  ;;  SECTION A:  parameter to adjust...
 
  ;; if all transits have locked the parameter selected, lock will be
  ;; switched to 1 and a new parameter will be chosen.
  
  ;; store parameter lock state for each transit
  lock = lonarr(n_elements(*self.transits))+1
  ;; store parameter pass state (automatic likelihood = 0 if, for
  ;; example, a parameter is beyond a user set LIMIT
  pass = lonarr(n_elements(*self.transits))
  
  while (where(lock eq 0))[0] eq -1 do begin
     pick = (sort(self.rand_obj2->getrandomnumbers(n_elements(((*self.transits)[0]->get()).params),/uniform)))[0]
     ;; setup test arrays for linked sets and locked values 
     linklock= [-1]
     linkval = [-1]
     linkpass= [-1]
     
     ;; reset lock and pass arrays to all zeros
     lock*=0
     pass*=0
     
     ;; parameter to permutate is in pick, now loop over all transits
     ;; and adjust the parameter.  There are a couple of options.  
     ;;
     for i=0,n_elements(*self.transits)-1,1 do begin
        ;; load the information of the transit and pull the current
        ;; link into the new one
        transit = (*self.transits)[i]->get()
        transit.params.new_link = transit.params.curr_link
        
        ;; is PICK locked for this specific transit?  If so, note this
        ;; and do not adjust new_link in any way.
        if transit.params[pick].fixed then lock[i] = 1 else begin
           ;; ok, the parameter value is not locked.  
           
           if (where(linklock eq transit.params[pick].set))[0] ne -1 then begin ;; this transit is part of a set, and a value has already been calculated
              ;; set this transit's new_link to the same value.
              transit.params[pick].new_link = linkval[where(linklock eq transit.params[pick].set)]
              pass[i] = linkpass[where(linklock eq transit.params[pick].set)]  ;; this is a set so must pass or fail in the same way.
           endif else begin  ;; either this transit and param is not in a set or the set has no good value yet
              ;; OK, how to adjust this parameter???
              
              ;; is a strict prior set?
              if transit.params[pick].prior[0] eq 2 then $ ;; if a PRIOR is set, pick from it
                 transit.params[pick].new_link = transit.params[pick].prior[1]+transit.params[pick].prior[2]* $
                                                 self.rand_obj->getrandomnumbers(1,/normal,/double) else $   ;; no prior:
                                                    transit.params[pick].new_link += (self.rand_obj->getrandomnumbers(1,/normal,/double)*(transit.params[pick].beta))
              
              
              ;; adjust the parameter value by a random gaussian
              ;; value of characteristic width beta
              
              ;; OK, is the value physical and can it pass?
              pass[i] = self->test_param(transit,pick)
              ;;endelse
              ;; now store the value and set in case other transits
              ;; are in the same linked set.
              linklock = [linklock,transit.params[pick].set]
              linkval  = [linkval,transit.params[pick].new_link]
              linkpass = [linkpass,pass[i]]
           endelse ;; end pick a new value for the param
        endelse    ;; end parameter not locked snippet
        (*self.transits)[i]->set,transit 
     endfor        ;; end loop over transits  
  endwhile         ;; repeat until a non locked parameter is picked.
  
  ;; now we have a new_link with a parameter which can evolve, and we
  ;; have an array of PASS values.  We now need to calculate
  ;; likelihoods ....
  self.jumptot[pick]++
  self->likelihood,/new,pass=pass
  
  tjump = self.rand_obj2->getrandomnumbers(1,/uniform)
  sets = (*self.sets)[*,pick]
  uniqsets = uniq(sets[sort(sets)])
  for i=0,n_elements(uniqsets)-1,1 do begin
     set = (sets[sort(sets)])[uniqsets[i]]
     if ((*self.prior)[where(sets eq (sets[sort(sets)])[uniqsets[i]]),pick])[0] eq 2 then prob = ((*self.accrate)[*,pick])[set] else $
        prob = exp(total((*self.newlikes)[where(sets eq set)])-total((*self.currlikes)[where(sets eq set)]))
     if tjump le min([prob,1d0]) then begin
        ;;(*self.jumps)[where(sets eq set)] = 1
        if total((*self.newlikes)[where(sets eq set)]) gt 0 then (*self.jumps)[where(sets eq set)] = 1 else $
           if total((*self.currlikes)[where(sets eq set)]) le 0 then (*self.jumps)[where(sets eq set)] = 1
     endif
     if (pass[where(sets eq set)])[0] eq 0 then (*self.jumps)[where(sets eq set)] = 0d0 ;; automatic fail regardless if param selection fails.
  endfor
  ;; if (self.iter mod 100) eq 0 then print,total(*self.currlikes)
  sets = 0L
  uniqsets = 0L
  acc = 0L
  
  if max(*self.jumps) then self.jumpcount[pick]++
  self.iter++
  prob = 0L
  tjump = 0L
  pass =0L
  lock =0L

  self->storelink,pick=pick
  pick=0L
end


pro tapmcmc::savechain
  for i=0,n_elements(*self.transits)-1,1 do begin
     transit = (*self.transits)[i]->get()
     transit.mcmc_complete = 1
     (*self.transits)[i]->set,transit
     transit = 0L
  endfor
  
  tap_state = (*self.transits)
  
;  stop
  save,TAP_state,filename='MCMC_chains/'+self.save_header+'_'+string(self.active_chain,format='(i2.2)')+'of'+string(self.num_chains,format='(i2.2)')+'.idlsav'
  TAP_state = 0L

  restore,'TAP_setup.idlsav'
  adjust= TAP_state[0]->get()
  if adjust.mcmc_complete eq 0 then *adjust.mcmc_files='MCMC_chains/'+self.save_header+'_'+string(self.active_chain,format='(i2.2)')+'of'+string(self.num_chains,format='(i2.2)')+'.idlsav' else $
     *adjust.mcmc_files=[*adjust.mcmc_files,'MCMC_chains/'+self.save_header+'_'+string(self.active_chain,format='(i2.2)')+'of'+string(self.num_chains,format='(i2.2)')+'.idlsav']
  adjust.mcmc_complete++
  TAP_State[0]->set,adjust
  adjust=0L
  save,TAP_state,filename='TAP_setup.idlsav'
  
end

function tapmcmc::info
  return,{savefile: self.save_header+'/TAP_setup.idlsav',$
          version: self.version}
end

pro tapmcmc::start,event
  centertlb,self.mcmc_base
  widget_control,self.mcmc_base,/realize
  
  widget_control, (*self.plot_windows)[0].window, Get_Value=wid
  (*self.plot_windows)[0].w_id=wid
   window,xsize=(*self.plot_windows)[0].x, $
         ysize=(*self.plot_windows)[0].y, $
         /pixmap,/free
  (*self.plot_windows)[0].pix_window = !d.window
  

  self.message = 'Executing '+self.save_header
  self.base_widget->message,message=self.message
  self->message
  
  ;; self.message = self.save_header
  ;; self->message
  
  self->ExecuteMCMC
end

pro tapmcmc::effective_length
  for i=0,self.num_transits-1,1 do begin
     transit = (*self.transits)[i]->get()
     range = self->burnrange(tnum=i)
    ; range = fillarr(1,.1*n_elements(*transit.params[0].mcmc_chain),n_elements(*transit.params[0].mcmc_chain)-1)
     tot_links = n_elements(range)*10d0
     
     for j=0,n_elements(transit.params)-1,1 do $
        if transit.params[j].fixed eq 0 and transit.params[j].prior[0] ne 2 then $ 
           if (*self.eff_len)[i,j] le self.min_eff[1] or finite((*self.eff_len)[i,j]) eq 0 then begin
        c_j  = [1d0]
        jact = [0d0]
        while c_j[n_elements(c_j)-1] gt 0.48d0 and jact[n_elements(jact)-1] le tot_links/2d0 do begin
           if n_elements(c_j) ge 2 then begin
              del =  c_j[n_elements(c_j)-2]-c_j[n_elements(c_j)-1] 
              step = max([jact[n_elements(c_j)-1]-jact[n_elements(c_j)-2],10])^(.1d0/del)
              step -= (step mod 10)
              if step lt 10 then step = 10
              if step gt 1000 then step = 1000
              if finite(step) eq 0 then step = 100
              jact = [jact,jact[n_elements(jact)-1]+step]
           endif else jact = [jact,jact[n_elements(jact)-1]+10d0]
           juse = jact[n_elements(jact)-1]/10d0
           avg0 = mean((*transit.params[j].mcmc_chain)[range],/double)
           avg1 = ((*transit.params[j].mcmc_chain)[range])[0:n_elements((*transit.params[j].mcmc_chain)[range])-(juse+1)]-avg0
           avg2 = ((*transit.params[j].mcmc_chain)[range])[(juse):n_elements((*transit.params[j].mcmc_chain)[range])-1]-avg0
           avg3 = mean(((*transit.params[j].mcmc_chain)[range]-avg0)^2d0,/double)
           c_j = [c_j,mean(avg1*avg2,/double)/avg3]
        endwhile
        (*self.corr_len)[i,j] =  interpol(jact,c_j,.5d0)
        (*self.eff_len)[i,j] =  tot_links/(*self.corr_len)[i,j]
     endif
     transit=0L
  endfor

  avg0 = 0L
  avg1 = 0L
  avg2 = 0L
  avg3 = 0L
  jact = 0L
  c_j = 0L
 
  self.message = '...min='+string(min((*self.eff_len)[where(finite(*self.eff_len))]),format='(i4)')+', req. '+string(self.min_eff[1],format='(i4)')
  if min((*self.eff_len)[where(finite(*self.eff_len))]) ge self.min_eff[1] then begin
     self.effpass = 1  
     self.message += ' good.'
  endif else self.message += ' extend 100k.'
  self->message 

  if 0 then begin
  window,23
  !p.multi=[0,1,2]
  plot,*((*self.transits)[1]->get()).params[1].mcmc_chain,thick=1,/ys,title=((*self.transits)[1]->get()).params[1].param
  plot,*((*self.transits)[1]->get()).params[0].likelihoods,thick=1,/ys,title='Likelihood'
  spawn,'say check this check this!'
  stop
  !p.multi=[0,1,1]
  colors=tap_colors()
  plot,*((*self.transits)[1]->get()).params[1].mcmc_chain,*((*self.transits)[1]->get()).params[2].mcmc_chain,/ys,xtitle=((*self.transits)[1]->get()).params[1].param,psym=8,symsize=.4,ytitle=((*self.transits)[1]->get()).params[2].param,color=colors.black,background=colors.white
  vline,88.89,color=colors.red
  hline,15.2,color=colors.red
  range = self->burnrange(tnum=1)
  vline,median((*((*self.transits)[1]->get()).params[1].mcmc_chain)[range]),color=colors.blue
  hline,median((*((*self.transits)[1]->get()).params[2].mcmc_chain)[range]),color=colors.blue
  stop


  plot,*((*self.transits)[1]->get()).params[5].mcmc_chain,*((*self.transits)[1]->get()).params[6].mcmc_chain,/ys,xtitle=((*self.transits)[1]->get()).params[5].param,psym=8,symsize=.4,ytitle=((*self.transits)[1]->get()).params[6].param,color=colors.black,xrange=[-.2,1.2],yrange=[-1.2,1.2],background=colors.white
  oplot,*((*self.transits)[0]->get()).params[5].mcmc_chain,*((*self.transits)[0]->get()).params[6].mcmc_chain,psym=8,symsize=.4,color=colors.blue
  oplot,*((*self.transits)[2]->get()).params[5].mcmc_chain,*((*self.transits)[2]->get()).params[6].mcmc_chain,psym=8,symsize=.4,color=colors.green

  vline,0d0,color=colors.red
  vline,1d0,color=colors.red
  hline,-1d0,color=colors.red
  hline,1d0,color=colors.red

stop

  plot,*((*self.transits)[1]->get()).params[5].mcmc_chain+*((*self.transits)[1]->get()).params[6].mcmc_chain,/ys,yrange=[-.2,1.2],psym=8,color=colors.black,background=colors.white
  oplot,*((*self.transits)[0]->get()).params[5].mcmc_chain+*((*self.transits)[0]->get()).params[6].mcmc_chain,psym=8,symsize=.4,color=colors.blue
  oplot,*((*self.transits)[2]->get()).params[5].mcmc_chain+*((*self.transits)[2]->get()).params[6].mcmc_chain,psym=8,symsize=.4,color=colors.green

  hline,0d0,color=colors.red
  hline,1d0,color=colors.red
  stop
endif
  ;stop
end


function tapmcmc::disentangle_parameterization,transit,val=val,curr=curr,new=new,pick=pick
  case self.parameterize_id of
     'basic': begin
        if keyword_set(val)  then transit.basic_params.value = transit.params.value
        if keyword_set(curr) then transit.basic_params.curr_link = transit.params.curr_link
        if keyword_set(new)  then transit.basic_params.new_link = transit.params.new_link
     end
     'adv1': begin
        stop
        print,'NOT USING THIS PARAMETERIZATION, REACHING HERE IS A CODE BUG'
     end
     'adv2': begin
                                ; stop
        if keyword_set(val) then gamma1 = 1d0+transit.params[where(strcmp(transit.params.param,'Eccentricity'))].value*sin(transit.params[where(strcmp(transit.params.param,'Omega'))].value)
        if keyword_set(val) then gamma2 = sqrt(1-transit.params[where(strcmp(transit.params.param,'Eccentricity'))].value^2d0)

        if keyword_set(curr) then gamma1c = 1d0+transit.params[where(strcmp(transit.params.param,'Eccentricity'))].curr_link*sin(transit.params[where(strcmp(transit.params.param,'Omega'))].curr_link)
        if keyword_set(curr) then gamma2c = sqrt(1-transit.params[where(strcmp(transit.params.param,'Eccentricity'))].curr_link^2d0)
        
        if keyword_set(new) then gamma1n = 1d0+transit.params[where(strcmp(transit.params.param,'Eccentricity'))].new_link*sin(transit.params[where(strcmp(transit.params.param,'Omega'))].new_link)
        if keyword_set(new) then gamma2n = sqrt(1-transit.params[where(strcmp(transit.params.param,'Eccentricity'))].new_link^2d0)
        ;; first just do the unchanged ones
        for i=0,n_elements(transit.basic_params)-1 do begin
           param = transit.basic_params[i].param
           if max(strcmp(param,['Period','Rp/R*','Mid Transit','Linear LD','Quadratic LD','Eccentricity',$
                                'Omega','OOT t^0','OOT t^1','OOT t^2','Sigma Red','Sigma White','Delta Light'])) then begin
              if keyword_set(val) then  transit.basic_params[i].value = transit.params[where(strcmp(transit.params.param,param))].value 
              if keyword_set(curr) then transit.basic_params[i].curr_link = transit.params[where(strcmp(transit.params.param,param))].curr_link
              if keyword_set(new) then transit.basic_params[i].new_link = transit.params[where(strcmp(transit.params.param,param))].new_link
           endif
           if strcmp(param,'Inclination') then begin
              ;;                 b         pi*T
              ;;  i = acos[ ----------- *  ---- ]
              ;;            sqrt(1-b^2)      P
              
              if keyword_set(val) then  transit.basic_params[i].value = (180/!dpi)*acos((transit.params[where(strcmp(transit.params.param,'b'))].value/$
                                                                                         sqrt(1-(transit.params[where(strcmp(transit.params.param,'b'))].value)^2)) * $
                                                                                        (!dpi*transit.params[where(strcmp(transit.params.param,'T'))].value)/$
                                                                                        transit.params[where(strcmp(transit.params.param,'Period'))].value)
              
              if keyword_set(curr) then  transit.basic_params[i].curr_link =  (180/!dpi)*acos((transit.params[where(strcmp(transit.params.param,'b'))].curr_link/$
                                                                                               sqrt(1-(transit.params[where(strcmp(transit.params.param,'b'))].curr_link)^2)) * $
                                                                                              (!dpi*transit.params[where(strcmp(transit.params.param,'T'))].curr_link)/$
                                                                                              transit.params[where(strcmp(transit.params.param,'Period'))].curr_link)
              
              if keyword_set(new) then  transit.basic_params[i].new_link =  (180/!dpi)*acos((transit.params[where(strcmp(transit.params.param,'b'))].new_link/$
                                                                                             sqrt(1-(transit.params[where(strcmp(transit.params.param,'b'))].new_link)^2)) * $
                                                                                            (!dpi*transit.params[where(strcmp(transit.params.param,'T'))].new_link)/$
                                                                                            transit.params[where(strcmp(transit.params.param,'Period'))].new_link)
              
           endif
           if strcmp(param,'a/R*') then begin
              ;;   a        P       gamma2^2 
              ;;  --- =  ------- * ---------- * sqrt(1-b^2)
              ;;   R*     pi* T      gamma1
              
              if keyword_set(val) then transit.basic_params[i].value = ( (gamma2^2/gamma1) * $
                                                                         (transit.params[where(strcmp(transit.params.param,'Period'))].value / $
                                                                          (!dpi*transit.params[where(strcmp(transit.params.param,'T'))].value)) * $
                                                                         sqrt(1-(transit.params[where(strcmp(transit.params.param,'b'))].value)^2))
              
              if keyword_set(curr) then transit.basic_params[i].curr_link = ( (gamma2c^2/gamma1c) * $
                                                                              (transit.params[where(strcmp(transit.params.param,'Period'))].curr_link / $
                                                                               (!dpi*transit.params[where(strcmp(transit.params.param,'T'))].curr_link)) * $
                                                                              sqrt(1-(transit.params[where(strcmp(transit.params.param,'b'))].curr_link)^2))
              
              
              if keyword_set(new) then transit.basic_params[i].new_link = ( (gamma2n^2/gamma1n) * $
                                                                            (transit.params[where(strcmp(transit.params.param,'Period'))].new_link / $
                                                                             (!dpi*transit.params[where(strcmp(transit.params.param,'T'))].new_link)) * $
                                                                            sqrt(1-(transit.params[where(strcmp(transit.params.param,'b'))].new_link)^2))
              
            endif
        endfor
     end
  endcase
  return,transit
end


pro tapmcmc::runchain,event
  if self.initial_run eq 1 then begin
     if self.active_chain eq 1 then self->BurnBetas,/guessbetas else self->burnbetas
     if self.testheap then self->testheap,'postbb'
     
     self->SetupRun,throw=5d0
     if self.testheap then self->testheap,'postbb post setup'
     
  endif else self->setuprun

  ;; 2.22 modifications... check effective length of chains after
  ;; reaching minimum paramter length, then extend if
  ;; necessary... extend by 1d4 at a time.

  
  self.effpass = 1
  if self.min_eff[0] then begin
     self.effpass = 0
     *self.corr_len *= 0
     *self.eff_len *= 0
     *self.eff_len += sqrt(-1)
  endif
  chain_length = self.chain_length[0]
  
  if self.testheap then self->testheap,'begin init run'
  
  for i=1d0,chain_length,1d0 do begin
     if (self.iter mod 1000d0) eq 0 then begin
        widget_control,self.progbar1,set_value=(self.iter/chain_length)
        widget_control,self.progbar2,set_value=((self.iter+(self.chain_length[0]*(self.active_chain-1)))/(self.num_chains*self.chain_length[0]))
     endif
     self->newlink
  endfor

  if self.effpass eq 0 then begin 
     self.message = 'Effective Length Loop...'
     self->message
     self->effective_length
     while self.effpass ne 1 do begin
        new_length = min([chain_length+1d5,self.chain_length[1]])
        widget_control,self.progbar1,set_value=(self.iter/new_length)
        widget_control,self.progbar2,set_value=((self.iter+(new_length*(self.active_chain-1)))/(self.num_chains*new_length))
        for i=chain_length,new_length,1d0 do begin
           if (self.iter mod 1000d0) eq 0 then begin
              widget_control,self.progbar1,set_value=(self.iter/new_length)
              widget_control,self.progbar2,set_value=((self.iter+(new_length*(self.active_chain-1)))/(self.num_chains*new_length))
           endif
           self->newlink
        endfor
        chain_length = new_length
        self->effective_length
        if self.effpass eq 0 then if chain_length ge self.chain_length[1] then begin
           self.effpass=1
           self.message='Maximum links reached.'
           self->message
        endif
    endwhile
  endif
  

  if self.testheap then self->testheap,'end init run'
  
  
  self->updatemod
  widget_control,self.progbar1,set_value=0d0
  widget_control,self.progbar2,set_value=((self.iter+(self.chain_length[0]*(self.active_chain-1)))/(self.num_chains*self.chain_length[0]))
end

pro tapmcmc::ExecuteMCMC,event
  self->lcplot
  ;; run the initial loops:
  
  self.initial_run = 1
  set = n_elements(*((*self.transits)[0]->get()).mcmc_files)
  if strcmp((*((*self.transits)[0]->get()).mcmc_files)[0],'-1') eq 0 then set++
  
  for i=set,self.num_chains,1 do begin
     self.active_chain = i
     self.message = 'Running Chain ' + string(self.active_chain,format='(i2.2)') + ' of ' + string(self.num_chains,format='(i2.2)')
     self->message
     self->RunChain

  ;   widget_control,self.progbar1,set_value=((self.iter/self.chain_length[0]))
  ;   widget_control,self.progbar2,set_value=((self.iter+(self.chain_length[0]*(self.active_chain-1)))/(self.num_chains*self.chain_length[0]))
     self->savechain
     self.base_widget->message,message=self.save_header+': Chain '+string(i,format='(i2.2)')+' of '+string(self.num_chains,format='(i2.2)')+' completed and stored.'
  ;   print,'ZAch-- try to get inference done after each chain...'
  ;   stop
  endfor
  
  self.initial_run = 0
  
  self.message = 'MCMC Run Complete'
  self->message
  ;if self.deconLD then begin
  ;    self.message = 'Reconstructing u1, u2'
  ;        
  ;    
  ;endif
  

  cd,'../'
end

function tapmcmc::burnrange,tnum=tnum
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
        if cross ne -1 then range[where(range lt cross),*] = 0d0
     endif 
  endfor
  
  ;;print,mm(range[*,0])
  ;;print,mm(range[where(range[*,1])])
  
  return,range[where(range[*,1])] 
end


pro tapmcmc::curr_analyze
  self.runval++
  if self.runval gt self.num_transits then self.runval = 1
  transit = (*self.transits)[self.runval-1]->get()
  for j=0,n_elements(transit.params)-1,1 do begin
     num = n_elements(*transit.params[j].mcmc_chain)
     range = mm(self->burnrange(tnum=self.runval-1))
     if range[1] eq n_elements(*transit.params[j].mcmc_chain) then range[1]--
     sorted = ((*transit.params[j].mcmc_chain)[range[0]:range[1]])[sort((*transit.params[j].mcmc_chain)[range[0]:range[1]])]
     range = n_elements(sorted)/100d0
     transit.params[j].runval = [sorted[50d0*range],$
                                 sorted[84.135d0*range]-sorted[50d0*range],$
                                 sorted[50d0*range]-sorted[15.865d0*range]]
     if transit.params[j].fixed then transit.params[j].runval[1:2] = -1d0
     
     num = 0L
     range = 0L
     sorted = 0L
  endfor
  (*self.transits)[self.runval-1]->set,transit
  transit = 0L
end


function tapmcmc::INIT,$
   base_parent=base_parent,$
   input_transits=input_transits,$
   restart=restart, version=version,$
   input_parameterization=input_parameterization,$
   _EXTRA=_extra
  

  ; SET TO 0 to deactivate c1, c2 parameterization of limb darkening. 
  ;self.deconLD = 1

  if keyword_set(restart) then begin
     self.save_header = restart
     cd,self.save_header
     self.save_header = (strsplit(self.save_header,'/',/extract))[n_elements(strsplit(self.save_header,'/'))-1]
;     print,self.save_header
  endif else begin
     self.save_header = 'TAPmcmc_'+curr_date(format='yyyymmdd_hhmm')
     if file_test(self.save_header) then begin
        t = self.save_header
        z = 2
        while(file_test(self.save_header)) do self.save_header = t+'_'+string(z,format='(i2.2)')   
        t = 0L
        z = 0L
     endif
     spawn,'mkdir '+self.save_header
     cd,self.save_header
     print,self.save_header
     spawn,'mv ../transit_setup.ascii .'
     spawn,'mkdir MCMC_chains'
  endelse
  
  ptr_free,self.colors
  self.colors = ptr_new(tap_colors())
  obj_destroy,self.base_widget
  self.base_widget = base_parent
  !EXCEPT = 0
  self.debug = 0
  self.testheap = 0
  self.version = version
  
  self.mcmc_base =  widget_base(frame=1,/column,title=strmid(self.save_header,8)+' '+self.version) 
  
  work_base = widget_base(self.mcmc_base,/column,frame=1)
  Xmanager, 'TAPmcmc',$
    self.mcmc_base,$
    /no_block  
  
  ptr_free,self.plot_windows
  self.plot_windows = ptr_new(replicate({x: 0d, y: 0d, w_id: 0L, pw_id: 0L, window:0L, pix_window: 0L, xrange: [0d0,0d0], yrange: [0d0,0d0]},1))
  (*self.plot_windows)[0].x = 300d0
  (*self.plot_windows)[0].y = 300d0
  (*self.plot_windows)[0].window = $
    widget_draw(work_base,$
                xsize=(*self.plot_windows)[0].x,$
                ysize=(*self.plot_windows)[0].y,$
                uvalue = 'mcmc plot 1',/align_center,frame=1)
  
  obj_destroy,self.progbarobj1,self.progbarobj2
  self.progbar1 = cw_progress(work_base,obj_ref=self.progbarobj1,/red,ysize=10d,xsize=(*self.plot_windows)[0].x)                 
  self.progbar2 = cw_progress(work_base,obj_ref=self.progbarobj2,/red,ysize=10d,xsize=(*self.plot_windows)[0].x)
  
  ptr_free,self.widget_bases
  self.widget_bases = ptr_new(lonarr(2))
  
  message = '('+ curr_date(format='hh:mm:ss yyyymmdd') +') TAP MCMC '+self.version
  (*self.widget_bases)[0] = widget_text(work_base, $
                                        font = textfont, $
                                        value = message, $
                                        /scroll, $
                                        ysize=5,scr_xsize=300)
  message = 0L

;  self.reduce_like = 1
 ; print,'REDUCE LIKE'

  ptr_free,self.transits
  self.transits = input_transits
  self.num_transits = n_elements(*self.transits)

  self.which_plot = self.num_transits  
  
  transit = (*self.transits)[0]->get()
  transit.mcmc_version = self.version
  self.deconLD = transit.ldtype
  (*self.transits)[0]->set,transit
  transit = 0L

  TAP_state = (*self.transits)
  save,TAP_state,filename='TAP_setup.idlsav'
  tap_state = 0L
  
  self.num_chains = ((*self.transits)[0]->get()).mcmc_params[0]
  self.chain_length =((*self.transits)[0]->get()).mcmc_params[1:2]
  
  self.min_eff = ((*self.transits)[0]->get()).min_effective_length
  ptr_free,self.corr_len,self.eff_len
;  self.min_eff = [1,100]
  self.corr_len = ptr_new(dblarr(self.num_transits,$
                                 n_elements(((*self.transits)[0]->get()).params)))
  self.eff_len = ptr_new(*self.corr_len+sqrt(-1))
  
  

  ;;
  ;; options:
  ;;  'basic' : "as is" for other TAP versions...
  ;;  'agol' : rp/r*, T, b, P, e, w, mu1, mu2, noise and corrections.
  self.parameterize_id = input_parameterization
  

 ; stop
  
  
  

  if 1-keyword_set(restart) then begin
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        ptr_free,transit.params.mcmc_chain
        for j=0,n_elements(transit.params)-1,1 do begin
           ptr_free,transit.params[j].mcmc_chain
           ptr_free,transit.params[j].likelihoods
           transit.params[j].mcmc_chain=ptr_new()         
           transit.params[j].likelihoods=ptr_new()      
        endfor        
        transit.mcmc_complete = 0
        (*self.transits)[i]->set,transit
        transit = 0L
     endfor
  endif else self.restart = restart
 
  ptr_free,self.sets,self.prior,self.accrate
  self.sets = ptr_new(lonarr(self.num_transits,n_elements(((*self.transits)[0]->get()).params)))
  self.prior = ptr_new(lonarr(self.num_transits,n_elements(((*self.transits)[0]->get()).params)))
  self.accrate = ptr_new(dblarr(self.num_transits,n_elements(((*self.transits)[0]->get()).params)))
  for i=0,self.num_transits-1,1 do (*self.sets)[i,*] =  ((*self.transits)[i]->get()).params.set
  for i=0,self.num_transits-1,1 do (*self.accrate)[i,*] =  ((*self.transits)[i]->get()).params.accept
  for i=0,self.num_transits-1,1 do (*self.prior)[i,*] = ((*self.transits)[i]->get()).params.prior[0]
  ptr_free,self.jumps,self.newlikes,self.currlikes
  self.jumps = ptr_new(lonarr(self.num_transits))
  self.newlikes = ptr_new(dblarr(self.num_transits))
  self.currlikes = ptr_new(dblarr(self.num_transits))

  self.rand_obj = obj_new('tap_RandomNumberGenerator')
  self.rand_obj2 = obj_new('tap_RandomNumberGenerator')
  
  if self.deconLD ne 0 then begin
     if self.deconLD eq 1 then print,"user requests u1+2u2,2u1-u2"
     if self.deconLD eq 2 then print,"user requests (u1+u2)^2, 0.5u1(u1+u2)^-1"
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        if transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].fixed eq 1 or $
           transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].fixed eq 1 or $ 
           transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].prior[0] eq 2 or $
           transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].prior[0] eq 2 or $
           transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].prior[0] eq 1 or $
           transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].prior[0] eq 1 then self.deconLD = 0
        transit=0L
     endfor
     if self.deconLD eq 0 then begin
        print,"switching to u1, u2"
        for i=0,self.num_transits-1,1 do begin
           transit = (*self.transits)[i]->get()
           transit.ldtype = 0
           (*self.transits)[i]->set,transit
        endfor
     endif
  endif
  if self.deconLD eq 1 then begin
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        
        u1 = transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].value
        u2 = transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].value
        stop_ld = 0
        if u1+u2 ge 1d0 $
           or u1+u2 le 0d0 $
           or u1 ge 1d0 $
           or u1 le 0d0 $
           or u2 le -1d0 $
           or u2 ge 1d0 then begin
           self.message = 'Transit '+transit.fname+' LD terms unphysical... adjusting'
           self->message
           if transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].fixed $
              and transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].fixed then stop_ld = 1
           if transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].fixed $ 
              and (u1 ge 1d0 $
                   or u1 le 0d0) then stop_ld = 1
           if transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].fixed $ 
              and (u2 le -1d0 $
                   or u2 ge 1d0) then stop_ld = 1
           
           if stop_ld then begin
              self.message = 'LD terms LOCKED at unphysical... stopping execution!'
              stop
           endif else begin
              u1 = .2d0
              u2 = .2d0
           endelse 
        endif 
        
        transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].value = 2*u1+u2
        transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].value   = u1-2*u2
        
        (*self.transits)[i]->set,transit
        
        transit = 0L
     endfor
     endif
  if self.deconLD eq 2 then begin
     for i=0,self.num_transits-1,1 do begin
        transit = (*self.transits)[i]->get()
        
        u1 = transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].value
        u2 = transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].value
        stop_ld = 0
        if u1+u2 ge 1d0 $
           or u1+u2 le 0d0 $
           or u1 ge 1d0 $
           or u1 le 0d0 $
           or u2 le -1d0 $
           or u2 ge 1d0 then begin
           self.message = 'Transit '+transit.fname+' LD terms unphysical... adjusting'
           self->message
           if transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].fixed $
              and transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].fixed then stop_ld = 1
           if transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].fixed $ 
              and (u1 ge 1d0 $
                   or u1 le 0d0) then stop_ld = 1
           if transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].fixed $ 
              and (u2 le -1d0 $
                   or u2 ge 1d0) then stop_ld = 1
           
           if stop_ld then begin
              self.message = 'LD terms LOCKED at unphysical... stopping execution!'
              stop
           endif else begin
              u1 = .2d0
              u2 = .2d0
           endelse 
        endif 
        
        transit.params[where(strcmp(transit.params.param,'Linear LD') eq 1)].value = (u1+u2)^2d0
        transit.params[where(strcmp(transit.params.param,'Quadratic LD') eq 1)].value   = 0.5*u1*(u1+u2)^(-1)
        
        (*self.transits)[i]->set,transit
        
        transit = 0L
     endfor
  endif
  
                                ;if self.deconLD then self.version += ' c1 c2 LD'
  

  return,1
end



pro tapmcmc__define
  struct = { tapmcmc, $
             colors: ptr_new(),$
             version: '',$
             debug: 0L,$
             testheap: 0L,$
             restart: '',$
             $ ;; widget stuff
             base_widget: obj_new(),$
             mcmc_base: 0L, $
             widget_bases: ptr_new(),$
             message_window: 0L, $
             message: '',$
             plot_windows: ptr_new(),$
             phased: 0L,$
             phasect: 0L,$
             progbar1: 0L,$
             progbar2: 0L,$
             progbarobj1: obj_new(),$
             progbarobj2: obj_new(),$
             $ ;; RaNDOM Numers!
             rand_obj: obj_new(),$
             rand_obj2: obj_new(),$
             $ 
             sets: ptr_new(),$
             prior: ptr_new(),$   
             accrate: ptr_new(),$
             jumps: ptr_new(),$
             newlikes: ptr_new(),$
             currlikes: ptr_new(),$
             $
             initial_run: 0,$
             active_chain: 0,$
             num_chains: 0,$
             chain_length: dblarr(2),$
             $
             parameterize_id: '',$
             $
             jump:0L ,$
             jumpcount: dblarr(15),$
             jumptot:   dblarr(15),$
             corr_len: ptr_new(),$
             eff_len: ptr_new(),$
             min_eff: dblarr(2),$
             effpass: 0L,$
             $
             adjust: ptr_new(),$
             need_adjust:  ptr_new(),$
             dec: ptr_new(),$
             inc: ptr_new(),$
             which_plot: 0,$
             $
             save_header: '',$
             $
             reduce_like: 0,$
             runval: 0L,$
             plot: 0d0,$
             iter: 0d0,$
             phi: ptr_new(),$
             num_transits: 0,$
             curr_redl: 0d0,$
             new_redl: 0d0,$
             deconLD: 0,$
             transits: ptr_new()  $
           }
  

end

function tapmcmc,base_parent=base_parent,input_transits=input_transits,only_version=only_version,_REF_EXTRA=_extra
  version = 'v2.52'
  ;; v2.40: -added advanced parameterization b^2, T
  ;;        -added additional LD parameterization
  ;; v2.50: added delta_light parameter (from johnson+ 2011ApJ...730...79J equation 9)
  ;; v2.51: -fixed advanced parameterization.... b^2->b
  ;;        -added likelihood tracking and outputting to mcmc chain files
  ;; v2.52: removed factor of 2 bug from conversion of likelihood into
  ;;        jump probability (in tapmcmc::newlink)

  if keyword_set(only_version) then return,{version: version} else $
     return,obj_new('tapmcmc',base_parent=base_parent,input_transits=input_transits,version=version,_EXTRA=_extra)
end

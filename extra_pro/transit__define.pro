;guess parameters object for widget

function transit::init
  self.lc = ptr_new(/allocate)
  return,1
end

function transit::cleanup
  ptr_free,self.lc
  return,1
end

pro transit::destroy
  ptr_free,self.lc
  obj_destroy,self
end

pro transit::set,value
; if data given, insert into pointer location
  if n_elements(value) ne 0 then *(self.lc)=value
  return
end

function transit::get,value
  if n_elements(*(self.lc)) ne 0 then value=*(self.lc)
  return,value
end

pro transit::load,value
  
end

pro transit::plot,model=model,final=final,emod=emod

end

pro transit__define
  void={transit,$
        lc: ptr_new()}
  return
end


Function generatearray,min,max,elements

  min *= 1d0
  max *= 1d0
  elements *= 1d0
  array = findgen(elements)
  array = (array*((max-min)/(elements-1d0))) + min
  
  return,array

END

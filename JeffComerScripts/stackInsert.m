function ret = stackInsert(c,stack,index)
if index <= 1
   ret = [{c} stack];
   return
end

if index > length(stack)
	ret = [stack {c}];   
end

ret = [{stack{1:index}} {c} {stack{index+1:end}}];




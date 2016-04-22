function [c,ret] = stackPop(stack)

if length(stack) == 0
   c = [];
   ret = {};
   return
end

if length(stack) == 1
   c = stack{1};
   ret = {};
   return
end

c = stack{1};
ret = {stack{2:end}};




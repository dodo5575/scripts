function result = stackEvaluate(stack)

buffer = {};
while length(stack) > 0
   [item stack] = stackPop(stack);
   
   if strcmp(item,'+')
      [b buffer] = stackPop(buffer);
      [a buffer] = stackPop(buffer);
      if (length(a)==0 | length(b)==0)
         result = [];
         return
      end
      
      buffer = stackPush(a+b, buffer);
   elseif strcmp(item,'-')
      [b buffer] = stackPop(buffer);
      [a buffer] = stackPop(buffer);
      if (length(a)==0 | length(b)==0)
         result = [];
         return
      end
      
      buffer = stackPush(a-b, buffer);
   elseif strcmp(item,'*')
      [b buffer] = stackPop(buffer);
      [a buffer] = stackPop(buffer);
      if (length(a)==0 | length(b)==0)
         result = [];
         return
      end

      buffer = stackPush(a*b, buffer);
   elseif strcmp(item,'/')
      [b buffer] = stackPop(buffer);
      [a buffer] = stackPop(buffer);
      if (length(a)==0 | length(b)==0 | b==0)
         result = [];
         return
      end

      buffer = stackPush(a/b, buffer);
   elseif strcmp(item,'^')   
      [b buffer] = stackPop(buffer);
      [a buffer] = stackPop(buffer);
      if (length(a)==0 | length(b)==0)
         result = [];
         return
      end
      
      buffer = stackPush(a^b, buffer);
elseif strcmp(item,'_')
	   [a buffer] = stackPop(buffer);
      if (length(a)==0)
         result = [];
         return
      end
      
      buffer = stackPush(-a, buffer);
else
   buffer = stackPush(item, buffer);
end
end

if length(buffer) ~= 1
   % Invalid stack.
   result = [];
else
   [result buffer] = stackPop(buffer);
end


   




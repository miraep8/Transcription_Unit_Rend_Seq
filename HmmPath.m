classdef HmmPath
   properties
      Cost = 0
      Path = []
   end
   methods
      function r = pathLen(obj)
          r = length(obj.Path);
      end
   end
end
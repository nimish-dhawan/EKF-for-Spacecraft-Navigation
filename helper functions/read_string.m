function val = read_string(str)
%Reads the numeric value val of the string str
%Returns zero if str is blank

if isempty(sscanf(str,'%f',1))
   val = 0;
else
   val = sscanf(str,'%f',1);
end
return

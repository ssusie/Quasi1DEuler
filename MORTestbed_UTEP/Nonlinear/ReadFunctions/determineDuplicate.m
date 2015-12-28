function  [] = determineDuplicate(text)

delim  = [1,findstr(text,'//')];
eqsign1 = regexp(text,'([^@].[=])','end');
eqsign2 = regexp(text,'([^,].[=])','end');
eqsign = intersect(eqsign1,eqsign2);

% eqsign = findstr(text,'=');

temp = cell(length(eqsign),1);
for i = 1:length(eqsign)
    temp{i} = text(delim(i):eqsign(i));
end

if length(temp) ~= length(unique(temp))
    error('There was a property duplicated in the FOM section and a CONFIG section.  This is not allowed, even if the duplication is consistent');
end

end
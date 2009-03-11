function d=finddim(v,div)
d=v;
while mod(d,div),
  d=d+1;
end;
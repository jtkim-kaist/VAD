function hz=meldm2hz(m)
hz=m;
sel=(m>10);
hz(sel)=1000*(10.^(0.0602*(m(sel)-10)));
hz(~sel)=100*m(~sel);

function m=hz2meldm(hz)
sel=(hz > 1000.0);
m=hz;
m(sel)=10+(log10(hz(sel))-3)/0.0602;
m(~sel)=hz(~sel)/100;

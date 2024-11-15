function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.mapping = data(2);
parms.alpha = data(3);
parms.revs_per_sec = data(4);
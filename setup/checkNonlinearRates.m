% check my rates. Do they agree with published curves?

r = 0:1e3:110e3;

nd = MSISatmosphere1(r'/1000);

rates = getNonlinearRates(r'/1000,nd,0);


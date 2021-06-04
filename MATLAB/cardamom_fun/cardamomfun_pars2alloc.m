function AF=cardamomfun_pars2alloc(CBR);

AF.fauto = CBR.PARS(:,2);
AF.ffol = CBR.PARS(:,3).*(1- AF.fauto);
AF.flab= CBR.PARS(:,13).*(1- AF.fauto - AF.ffol);
 AF.froo= CBR.PARS(:,4).*(1- AF.fauto - AF.ffol - AF.flab);
AF.fwoo= 1 - AF.fauto - AF.ffol - AF.flab - AF.froo;

AF.all=[AF.fauto,AF.ffol,AF.flab,AF.froo,AF.fwoo];


end
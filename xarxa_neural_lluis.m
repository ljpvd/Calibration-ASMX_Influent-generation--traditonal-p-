function net=xarxa_neural_lluis(in,out)

% in:  input vector
% out: output we want to generate
% net: nova xarxa neural entrenada amb in & out
% net: new neural network trained with in & out
% ara per ara s'han de passar les dades com a una fila i no com una columna
% important: you should pass the variables as single rows, not columns

% nova xarxa neural de 100 neurones
neurones=100;
net=newff(in,out,[neurones 1],{'tansig','purelin'}); 

% entrenar la xarxa amb les variables in & out
net.trainParam.epochs=10000;
net=train(net,in,out);

